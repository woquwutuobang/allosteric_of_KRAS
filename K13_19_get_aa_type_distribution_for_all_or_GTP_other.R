library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)
library(ggrepel)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"


classify_mutations_comprehensive <- function(input, anno, assay_sele, wt_aa) {
  #' 综合突变分类函数
  #'
  #' @param input 包含ddG数据的文件路径或数据框
  #' @param anno 包含注释信息的文件路径或数据框
  #' @param assay_sele 要分析的assay名称（如"K13", "RAF1"等）
  #' @param wt_aa 野生型氨基酸序列
  #' @return 包含突变分类的数据框
  
  # Read and process ddG data
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  # Data processing
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Read annotation data
  if (is.character(anno)) {
    anno_data <- fread(anno)
  } else {
    anno_data <- as.data.table(anno)
  }
  
  # Merge with annotation data
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  # Define binding site types
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  gxpmg_col <- paste0("GXPMG_scHAmin_ligand_", assay_sele)
  
  data_plot[get(scHAmin_col) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  # Calculate regulatory threshold
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  # Classify site types
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # Merge mutation data with site type information
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # Identify allosteric mutations
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # Classify mutation types
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "Orthosteric site huge differences"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "Orthosteric site small differences"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "GTP binding allosteric mutation"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "GTP binding other mutation"]
  
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == TRUE, mutation_type := "Allosteric mutation"]
  
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == FALSE, mutation_type := "Other mutation"]
  
  # Set factor levels for mutation types
  data_plot_mutation[, mutation_type := factor(mutation_type,
                                               levels = c("Orthosteric site huge differences",
                                                          "Orthosteric site small differences",
                                                          "GTP binding allosteric mutation",
                                                          "GTP binding other mutation",
                                                          "Allosteric mutation",
                                                          "Other mutation"))]
  
  # 添加有用的列
  data_plot_mutation[, `:=`(
    assay = assay_sele,
    reg_threshold = reg_threshold,
    is_significant = allosteric_mutation
  )]
  
  # 返回关键列
  result_cols <- c("Pos", "Pos_real", "wt_codon", "mt_codon", "mt", 
                   "mean_kcal/mol", "std_kcal/mol", "mean", "std",
                   "site_type", "allosteric_mutation", "mutation_type",
                   "assay", "reg_threshold", "is_significant")
  
  result_cols <- intersect(result_cols, names(data_plot_mutation))
  
  return(data_plot_mutation[, ..result_cols])
}



#### K13
mutation_K13 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "K13", 
  wt_aa )
#names(mutation_K13)
#mutation_K13 <- mutation_K13[, c(2:4, 23:26, 28)]

#### RAF1
mutation_RAF1 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "RAF1", 
  wt_aa )
#mutation_RAF1 <- mutation_RAF1[, c(2:4, 23:26, 28)]

#### K27
mutation_K27 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "K27", 
  wt_aa )
#mutation_K27 <- mutation_K27[, c(2:4, 23:26, 28)]

# 合并数据
mutation_merge1 <- merge(mutation_K27, mutation_RAF1, 
                         by = c("Pos","Pos_real", "wt_codon", "mt_codon", "mt"), 
                         suffixes = c("_K27", "_RAF1"))

mutation_merge2 <- merge(mutation_K13, mutation_RAF1, 
                         by = c("Pos","Pos_real", "wt_codon", "mt_codon", "mt"), 
                         suffixes = c("_K13", "_RAF1"))

# 识别共同的显著突变（同时是GTP binding allosteric mutation和Allosteric mutation）
common_mutations_K27_RAF1 <- mutation_merge1[
  (mutation_type_K27 == "GTP binding allosteric mutation" | 
     mutation_type_K27 == "Allosteric mutation") &
    (mutation_type_RAF1 == "GTP binding allosteric mutation" | 
       mutation_type_RAF1 == "Allosteric mutation")
]

common_mutations_K13_RAF1 <- mutation_merge2[
  (mutation_type_K13 == "GTP binding allosteric mutation" | 
     mutation_type_K13 == "Allosteric mutation") &
    (mutation_type_RAF1 == "GTP binding allosteric mutation" | 
       mutation_type_RAF1 == "Allosteric mutation")
]



# 定义化学类型函数
define_chemotype <- function(aa) {
  aromatic <- c("F","W","Y")
  aliphatic <- c("A","V","I","L","M")
  polar_uncharged <- c("S","T","N","Q","C")
  positive <- c("K","R","H")
  negative <- c("D","E")
  special <- c("G","P")
  
  if (aa %in% aromatic) return("Aromatic")
  if (aa %in% aliphatic) return("Aliphatic")
  if (aa %in% polar_uncharged) return("Polar uncharged")
  if (aa %in% positive) return("Positive")
  if (aa %in% negative) return("Negative")
  if (aa %in% special) return("Special")
  return(NA)
}



# 定义象限分类函数（修复版本）
classify_quadrant_vector2 <- function(x, y) {
  result <- character(length(x))
  for (i in seq_along(x)) {
    if (x[i] > 0 & y[i] > 0) {
      result[i] <- "Bp1>0, Bp2>0"
    } else if (x[i] < 0 & y[i] > 0) {
      result[i] <- "Bp1<0, Bp2>0"
    } else if (x[i] < 0 & y[i] < 0) {
      result[i] <- "Bp1<0, Bp2<0"
    } else if (x[i] > 0 & y[i] < 0) {
      result[i] <- "Bp1>0, Bp2<0"
    } else {
      result[i] <- "On axis"
    }
  }
  return(result)
}


# 提取突变信息
extract_mutation_info <- function(mutations_df) {
  result <- data.table()
  
  for (i in 1:nrow(mutations_df)) {
    mt <- mutations_df[i, mt]
    # 解析突变，格式如 "A15T"
    original_aa <- substr(mt, 1, 1)
    mutated_aa <- substr(mt, nchar(mt), nchar(mt))
    
    result <- rbind(result, data.table(
      mutation = mt,
      original_aa = original_aa,
      mutated_aa = mutated_aa,
      original_type = define_chemotype(original_aa),
      mutated_type = define_chemotype(mutated_aa)
    ))
  }
  
  return(result)
}

# 分析象限化学类型偏好
analyze_chemotype_by_quadrant <- function(common_mutations, data, x_col, y_col, comparison_name) {
  if (nrow(common_mutations) == 0) return(NULL)
  
  # 合并ddG值
  merged_data <- merge(common_mutations, 
                       data[, .(mt, x_ddG = get(x_col), y_ddG = get(y_col))], 
                       by = "mt")
  
  # 为每个突变添加象限分类
  merged_data[, quadrant := classify_quadrant_vector2(x_ddG, y_ddG)]
  
  # 提取突变信息
  mutation_info <- extract_mutation_info(merged_data)
  
  # 合并信息
  final_data <- merge(merged_data[, .(mt, quadrant)], mutation_info, by.x = "mt", by.y = "mutation")
  
  # 添加比较名称
  final_data[, comparison := comparison_name]
  
  return(final_data)
}

# 分析两组数据
k27_data <- analyze_chemotype_by_quadrant(common_mutations_K27_RAF1, mutation_merge1,
                                          "mean_kcal/mol_K27", "mean_kcal/mol_RAF1",
                                          "K27(Bp1) vs RAF1(Bp2)")

k13_data <- analyze_chemotype_by_quadrant(common_mutations_K13_RAF1, mutation_merge2,
                                          "mean_kcal/mol_K13", "mean_kcal/mol_RAF1",
                                          "K13(Bp1) vs RAF1(Bp2)")

# 合并数据
all_data <- rbind(k27_data, k13_data, fill = TRUE)

# 定义化学类型顺序
chemotype_order <- c("Aromatic", "Aliphatic", "Polar uncharged", "Positive", "Negative", "Special")
all_data[, original_type := factor(original_type, levels = chemotype_order)]
all_data[, mutated_type := factor(mutated_type, levels = chemotype_order)]

# 定义颜色映射
chemotype_colors <- c(
  "Aromatic" = "#FFB0A5",        
  "Aliphatic" = "#75C2F6",       
  "Polar uncharged" = "#09B636", 
  "Positive" = "#C68EFD",        
  "Negative" = "#F4AD0C",        
  "Special" = "#F1DD10"          
)


original_counts <- all_data[, .(count = .N), by = .(comparison, quadrant, type = original_type)]
mutated_counts <- all_data[, .(count = .N), by = .(comparison, quadrant, type = mutated_type)]

# 添加类型标签
original_counts[, time := "WT AA"]
mutated_counts[, time := "MT AA"]
combined_counts <- rbind(original_counts, mutated_counts)

# 方案1：堆叠百分比条形图
p1 <- ggplot(combined_counts, aes(x = quadrant, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  facet_grid(comparison ~ time) +
  scale_fill_manual(values = chemotype_colors, name = "Chemical Type") +
  labs(x = "Functional Pattern", y = "Proportion", 
       title = "Chemical Type Distribution by Functional Pattern") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 13, hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(labels = scales::percent_format())

print(p1)



ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251217/Chemical Type Distribution by Quadrant (Stacked) 2.png", 
       p1, width = 8, height = 6, dpi = 300)






#########分GTP和other进行绘图
# 给all_data加上 mutation_class 列
all_data[, mutation_class := ifelse(
  grepl("GTP binding allosteric mutation", comparison) | 
    (comparison == "K27(Bp1) vs RAF1(Bp2)" & mt %in% common_mutations_K27_RAF1[mutation_type_K27 == "GTP binding allosteric mutation", mt]) |
    (comparison == "K13(Bp1) vs RAF1(Bp2)" & mt %in% common_mutations_K13_RAF1[mutation_type_K13 == "GTP binding allosteric mutation", mt]),
  "GTP binding allosteric mutation",
  "Allosteric mutation"
)]

# 按 mutation_class 分组绘图
plot_by_mutation_class <- function(all_data_subset, mutation_class_name) {
  plot_data <- all_data_subset[mutation_class == mutation_class_name]
  
  original_counts <- plot_data[, .(count = .N), by = .(comparison, quadrant, type = original_type)]
  mutated_counts <- plot_data[, .(count = .N), by = .(comparison, quadrant, type = mutated_type)]
  
  original_counts[, time := "WT AA"]
  mutated_counts[, time := "MT AA"]
  
  combined_counts <- rbind(original_counts, mutated_counts)
  
  ggplot(combined_counts, aes(x = quadrant, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    facet_grid(comparison ~ time) +
    scale_fill_manual(values = chemotype_colors, name = "Chemical Type") +
    labs(x = "Functional Pattern", y = "Proportion", 
         title = paste0("Chemical Type Distribution: ", mutation_class_name)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, hjust = 0.5),
      strip.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    scale_y_continuous(labels = scales::percent_format())
}

# 绘制 GTP binding allosteric mutation 图
p_gtp <- plot_by_mutation_class(all_data, "GTP binding allosteric mutation")
print(p_gtp)

# 绘制 Allosteric mutation 图
p_allo <- plot_by_mutation_class(all_data, "Allosteric mutation")
print(p_allo)

# 保存图片
ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251217/Chemical_Type_Distribution_GTP_binding_allosteric_mutation2.png", 
       p_gtp, width = 8, height = 6, dpi = 300)

ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251217/Chemical_Type_Distribution_Allosteric_mutation2.png", 
       p_allo, width = 8, height = 6, dpi = 300)

