#########################
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

# 添加target位点标记
target_positions <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
common_mutations_K27_RAF1[, is_target_position := Pos_real %in% target_positions]
common_mutations_K13_RAF1[, is_target_position := Pos_real %in% target_positions]

# 打印共同的突变（按目标位点分组）
print_common_mutations <- function(common_mutations, comparison_name) {
  cat("\n", comparison_name, ":\n", sep = "")
  
  # 目标位点突变
  target_muts <- common_mutations[is_target_position == TRUE]
  if (nrow(target_muts) > 0) {
    cat("Target position mutations (Positions:", paste(unique(target_muts$Pos_real), collapse = ", "), "):\n")
    for (i in 1:nrow(target_muts)) {
      if (comparison_name == "K27 vs RAF1") {
        cat(sprintf("  %s (Position: %d, ddG_K27: %.3f, ddG_RAF1: %.3f)\n", 
                    target_muts[i, mt], target_muts[i, Pos_real],
                    target_muts[i, get("mean_kcal/mol_K27")], 
                    target_muts[i, get("mean_kcal/mol_RAF1")]))
      } else {
        cat(sprintf("  %s (Position: %d, ddG_K13: %.3f, ddG_RAF1: %.3f)\n", 
                    target_muts[i, mt], target_muts[i, Pos_real],
                    target_muts[i, get("mean_kcal/mol_K13")], 
                    target_muts[i, get("mean_kcal/mol_RAF1")]))
      }
    }
  }
  
  # 非目标位点突变
  non_target_muts <- common_mutations[is_target_position == FALSE]
  if (nrow(non_target_muts) > 0) {
    cat("\nNon-target position mutations (Positions:", paste(unique(non_target_muts$Pos_real), collapse = ", "), "):\n")
    for (i in 1:nrow(non_target_muts)) {
      if (comparison_name == "K27 vs RAF1") {
        cat(sprintf("  %s (Position: %d, ddG_K27: %.3f, ddG_RAF1: %.3f)\n", 
                    non_target_muts[i, mt], non_target_muts[i, Pos_real],
                    non_target_muts[i, get("mean_kcal/mol_K27")], 
                    non_target_muts[i, get("mean_kcal/mol_RAF1")]))
      } else {
        cat(sprintf("  %s (Position: %d, ddG_K13: %.3f, ddG_RAF1: %.3f)\n", 
                    non_target_muts[i, mt], non_target_muts[i, Pos_real],
                    non_target_muts[i, get("mean_kcal/mol_K13")], 
                    non_target_muts[i, get("mean_kcal/mol_RAF1")]))
      }
    }
  }
  
  # 统计信息
  cat(sprintf("\nTotal: %d mutations (Target: %d, Non-target: %d)\n", 
              nrow(common_mutations), nrow(target_muts), nrow(non_target_muts)))
}

#print_common_mutations(common_mutations_K27_RAF1, "K27 vs RAF1")
#print_common_mutations(common_mutations_K13_RAF1, "K13 vs RAF1")


create_scatter_plot <- function(
    data, x_var, y_var, x_lab, y_lab, common_mutations, title
) {
  
  # ===============================
  # 1. 颜色映射（名字要全程一致）
  # ===============================
  color_palette <- c(
    "Target position" = "#F4AD0C",      # nucleotide pocket allosteric mutation
    "Non-target position" = "#C68EFD",  # allosteric mutation
    "Other" = "grey80"
  )
  
  # ===============================
  # 2. 构造绘图数据
  # ===============================
  plot_data <- copy(data)
  plot_data[, highlight_type := "Other"]
  
  if (nrow(common_mutations) > 0) {
    common_mt_list <- common_mutations$mt
    
    plot_data[
      mt %in% common_mt_list,
      highlight_type := ifelse(
        Pos_real %in% target_positions,
        "Target position",
        "Non-target position"
      )
    ]
  }
  
  # ===============================
  # 3. Legend 位置
  # ===============================
  legend_x <- 2.8
  legend_y <- 2.8
  
  # ===============================
  # 4. 主图
  # ===============================
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    
    geom_point(
      data = plot_data[highlight_type == "Other"],
      color = color_palette["Other"],
      alpha = 0.5, size = 1
    ) +
    
    geom_point(
      data = plot_data[highlight_type == "Non-target position"],
      color = color_palette["Non-target position"],
      size = 3, alpha = 0.8
    ) +
    
    geom_point(
      data = plot_data[highlight_type == "Target position"],
      color = color_palette["Target position"],
      size = 3.5, alpha = 0.9
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "black", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "black", alpha = 0.5) +
    
    labs(x = x_lab, y = y_lab, title = title) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    coord_cartesian(xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3))
  
  # ===============================
  # 5. 标签（只给共同突变）
  # ===============================
  if (nrow(common_mutations) > 0) {
    
    label_data <- plot_data[
      highlight_type %in% c("Target position", "Non-target position")
    ]
    
    p <- p +
      ggrepel::geom_text_repel(
        data = label_data,
        aes(label = mt, color = highlight_type),
        size = 3.2,
        max.overlaps = 50,
        box.padding = 0.5,
        min.segment.length = 0.1,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = color_palette,
        guide = "none"
      )
  }
  
  # ===============================
  # 6. 手动画 Color Legend（✅ 已修复）
  # ===============================
  p <- p +
    
    annotate(
      "rect",
      xmin = legend_x - 0.6, xmax = legend_x + 0.6,
      ymin = legend_y - 0.8, ymax = legend_y + 0.1,
      fill = "white", color = "black",
      alpha = 0.85, linewidth = 0.3
    ) +
    
    annotate(
      "text",
      x = legend_x, y = legend_y,
      label = "Color Legend",size = 3.5
    ) +
    
    annotate(
      "point",
      x = legend_x - 0.4, y = legend_y - 0.25,
      color = color_palette["Target position"], size = 3
    ) +
    annotate(
      "text",
      x = legend_x - 0.3, y = legend_y - 0.25,
      label = "Nucleotide pocket",
      hjust = 0, size = 3.2
    ) +
    
    annotate(
      "point",
      x = legend_x - 0.4, y = legend_y - 0.45,
      color = color_palette["Non-target position"], size = 3
    ) +
    annotate(
      "text",
      x = legend_x - 0.3, y = legend_y - 0.45,
      label = "Other region",
      hjust = 0, size = 3.2
    )
  
  return(p)
}

# 创建散点图
p1 <- create_scatter_plot(
  data = mutation_merge1,
  x_var = "mean_kcal/mol_K27",
  y_var = "mean_kcal/mol_RAF1", 
  x_lab = "Binding ΔΔG (K27) (kcal/mol)",
  y_lab = "Binding ΔΔG (RAF1) (kcal/mol)",
  common_mutations = common_mutations_K27_RAF1,
  title = "K27 vs RAF1 - Common Allosteric Mutations"
)

p2 <- create_scatter_plot(
  data = mutation_merge2,
  x_var = "mean_kcal/mol_K13",
  y_var = "mean_kcal/mol_RAF1",
  x_lab = "Binding ΔΔG (K13) (kcal/mol)", 
  y_lab = "Binding ΔΔG (RAF1) (kcal/mol)",
  common_mutations = common_mutations_K13_RAF1,
  title = "K13 vs RAF1 - Common Allosteric Mutations"
)

p1
p2

# 组合图形
combined_plot <- p1 + p2 + 
  plot_layout(ncol = 2)

# 显示图形
print(combined_plot)

# 保存图形
#ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251220/Common_Allosteric_Mutations_Scatter_with_Legend.png", 
#       p1, width = 13, height = 6, dpi = 300)

ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251220/K27 VS RAF1 highlight allosteric mutations.pdf", 
       p1,device = cairo_pdf, width = 6, height = 6, dpi = 300)



# 保存图形
#ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251220/Common_Allosteric_Mutations_Scatter_with_Legend.png", 
#       p2, width = 13, height = 6, dpi = 300)

ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251220/K13 VS RAF1 highlight allosteric mutations.pdf", 
       p2,device = cairo_pdf, width = 6, height = 6, dpi = 300)
