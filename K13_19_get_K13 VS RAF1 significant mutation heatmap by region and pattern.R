library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(krasddpcams)

# 读取数据和注释
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 定义函数，进行K13和RAF1突变分析
classify_mutations_comprehensive <- function(input, anno, assay_sele, wt_aa) {
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  if (is.character(anno)) {
    anno_data <- fread(anno)
  } else {
    anno_data <- as.data.table(anno)
  }
  
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  gxpmg_col <- paste0("GXPMG_scHAmin_ligand_", assay_sele)
  
  data_plot[get(scHAmin_col) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  reg_threshold <- data_plot[binding_type == "binding site", 
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE), 
                             .SDcols = c("mean", "sigma")]
  
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
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
  
  data_plot_mutation[, mutation_type := factor(mutation_type,
                                               levels = c("Orthosteric site huge differences",
                                                          "Orthosteric site small differences",
                                                          "GTP binding allosteric mutation",
                                                          "GTP binding other mutation",
                                                          "Allosteric mutation",
                                                          "Other mutation"))]
  
  result_cols <- c("Pos", "Pos_real", "wt_codon", "mt_codon", "mt", 
                   "mean_kcal/mol", "std_kcal/mol", "mean", "std",
                   "site_type", "allosteric_mutation", "mutation_type",
                   "assay", "reg_threshold", "is_significant")
  
  result_cols <- intersect(result_cols, names(data_plot_mutation))
  
  return(data_plot_mutation[, ..result_cols])
}


# 分析K13和RAF1的突变
mutation_K13 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "K13", 
  wt_aa
)

mutation_RAF1 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "RAF1", 
  wt_aa
)

# 合并K13和RAF1的显著突变信息
mutation_merge <- merge(mutation_K13, mutation_RAF1, 
                        by = c("Pos", "Pos_real", "wt_codon", "mt_codon", "mt"), 
                        suffixes = c("_K13", "_RAF1"))

# 过滤显著突变（Allosteric mutation和GTP binding allosteric mutation）
common_mutations <- mutation_merge[
  (mutation_type_K13 %in% c("Allosteric mutation", "GTP binding allosteric mutation")) &
    (mutation_type_RAF1 %in% c("Allosteric mutation", "GTP binding allosteric mutation"))
]


library(data.table)
library(ggplot2)

# ============================================================
# Step 1. 给 K13 vs RAF1 的显著突变定义 pattern
# ============================================================

common_mutations[, pattern := fifelse(
  sign(mean_K13) == sign(mean_RAF1),
  "Same-direction",
  fifelse(
    sign(mean_K13) != sign(mean_RAF1),
    "Opposite-direction",
    "Other"
  )
)]

# ============================================================
# Step 2. 定义区域：GTP-binding vs Other regions
# （用 site_type，而不是硬编码位置）
# ============================================================

common_mutations[, region := fifelse(
  site_type_K13 == "GTP binding interface site" |
    site_type_RAF1 == "GTP binding interface site",
  "GTP-binding region",
  "Other regions"
)]

# ============================================================
# Step 3. 构建热图用的 long-format 数据
# ============================================================

heatmap_df <- melt(
  common_mutations,
  measure.vars = c("mean_K13", "mean_RAF1"),
  variable.name = "binder",
  value.name = "ddG"
)

heatmap_df[, binder := fifelse(binder == "mean_K13", "K13", "RAF1")]
heatmap_df[, region := NULL]

# 合并注释信息（pattern / region / Pos_real）
heatmap_df <- merge(
  heatmap_df,
  unique(common_mutations[, .(mt, Pos_real, pattern, region)]),
  by = "mt"
)

# ============================================================
# Step 4. Y 轴排序规则
# 在同一 region × pattern 内，按 Pos_real 从大到小
# ============================================================

heatmap_df[, mt := factor(
  mt,
  levels = heatmap_df[
    order(region, pattern, -Pos_real),
    unique(mt)
  ]
)]

# ============================================================
# Step 5. 画图函数（facet = pattern）
# ============================================================

plot_heatmap_by_region <- function(df, region_name, title_text) {
  ggplot(
    df[region == region_name],
    aes(x = binder, y = mt, fill = ddG)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(
      pattern ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_fill_gradient2(
      low = "#4575b4",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      name = expression(Delta*Delta*G)
    ) +
    labs(
      x = NULL,
      y = "Mutation",
      title = title_text
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    )
}

# ============================================================
# Step 6. 分别画两张图
# ============================================================

p_gtp <- plot_heatmap_by_region(
  heatmap_df,
  region_name = "GTP-binding region",
  title_text = "GTP-binding region: distinct K13 vs RAF1 mutation patterns"
)

p_other <- plot_heatmap_by_region(
  heatmap_df,
  region_name = "Other regions",
  title_text = "Non-GTP regions: distinct K13 vs RAF1 mutation patterns"
)

# 显示
p_gtp
p_other

ggsave(
  "Figure4_K13_vs_RAF1_GTP.pdf",
  p_gtp,
  width = 4,
  height = 6
)

ggsave(
  "Figure4_K13_vs_RAF1_Other.pdf",
  p_other,
  width = 4,
  height = 8
)






################################# 这个对
library(data.table)
library(ggplot2)

# ============================================================
# Step 0. 确保 common_mutations 是 data.table
# ============================================================
setDT(common_mutations)

# ============================================================
# Step 1. 定义 4-class functional pattern (K13 vs RAF1)
# ============================================================
common_mutations[, pattern := fifelse(
  mean_K13 > 0 & mean_RAF1 > 0, "Both disrupted",
  fifelse(
    mean_K13 < 0 & mean_RAF1 < 0, "Both enhanced",
    fifelse(
      mean_K13 < 0 & mean_RAF1 > 0, "K13 enhanced / RAF1 disrupted",
      fifelse(
        mean_K13 > 0 & mean_RAF1 < 0, "RAF1 enhanced / K13 disrupted",
        NA_character_
      )
    )
  )
)]

# 固定 pattern 顺序（非常重要）
common_mutations[, pattern := factor(
  pattern,
  levels = c(
    "Both disrupted",
    "K13 enhanced / RAF1 disrupted",
    "Both enhanced",
    "RAF1 enhanced / K13 disrupted"
  )
)]

# ============================================================
# Step 2. 定义区域（用 site_type，而不是硬编码位置）
# ============================================================
common_mutations[, region := fifelse(
  site_type_K13 == "GTP binding interface site" |
    site_type_RAF1 == "GTP binding interface site",
  "GTP-binding region",
  "Other regions"
)]

# ============================================================
# Step 3. 构建热图用 long-format 数据
# ============================================================
heatmap_df <- melt(
  common_mutations,
  measure.vars = c("mean_K13", "mean_RAF1"),
  variable.name = "binder",
  value.name = "ddG"
)

heatmap_df[, binder := fifelse(binder == "mean_K13", "K13", "RAF1")]

# 防止 merge 后出现 .x / .y
heatmap_df[, region := NULL]
heatmap_df[, Pos_real := NULL]

# 删除 melt 后可能重复的列
heatmap_df[, c("pattern", "region", "Pos_real") := NULL]

# 再 merge 注释
heatmap_df <- merge(
  heatmap_df,
  unique(common_mutations[, .(mt, Pos_real, pattern, region)]),
  by = "mt"
)

# ============================================================
# Step 4. 排序 & 设置 mt factor
# ============================================================
# 先排序：region -> pattern -> Pos_real (降序)
setorder(heatmap_df, region, pattern, -Pos_real)

# 再设置 mt factor 顺序
heatmap_df[, mt := factor(mt, levels = unique(mt))]

# ============================================================
# Step 5. 画图函数（pattern 纵向分面）
# ============================================================
plot_heatmap_by_region <- function(df, region_name, title_text) {
  ggplot(
    df[region == region_name],
    aes(x = binder, y = mt, fill = ddG)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(
      pattern ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_fill_gradient2(
      low = "#1B38A6",
      mid = "white",
      high = "#F4270C",
      midpoint = 0,
      name = expression(Delta*Delta*G)
    ) +
    labs(
      x = NULL,
      y = "Mutation",
      title = title_text
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      #strip.text.y = element_text(size = 6),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    )
}

# ============================================================
# Step 6. 分别画两张图
# ============================================================
p_gtp <- plot_heatmap_by_region(
  heatmap_df,
  region_name = "GTP-binding region",
  title_text = "GTP-binding region"
)

p_other <- plot_heatmap_by_region(
  heatmap_df,
  region_name = "Other regions",
  title_text = "Non-GTP regions"
)

# ============================================================
# Step 7. 显示
# ============================================================
p_gtp
p_other



ggsave(
  "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251225/Figure4_K13_vs_RAF1_GTP.pdf",
  device=cairo_pdf,
  p_gtp,
  width = 4,
  height = 6
)

ggsave(
  "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251225/Figure4_K13_vs_RAF1_Other.pdf",
  p_other,
  device=cairo_pdf,
  width = 4,
  height = 8
)
