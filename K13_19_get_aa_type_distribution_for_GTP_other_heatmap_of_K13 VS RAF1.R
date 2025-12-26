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

common_mutations[, pattern := factor(
  pattern,
  levels = c(
    "Both disrupted",
    "K13 enhanced / RAF1 disrupted",
    "Both enhanced",
    "RAF1 enhanced / K13 disrupted"
  )
)]


sig_levels <- c("Allosteric mutation", "GTP binding allosteric mutation")

sig_mut <- common_mutations[
  mutation_type_K13 %in% sig_levels &
    mutation_type_RAF1 %in% sig_levels
]


long_sig <- rbindlist(list(
  sig_mut[, .(
    Pos,
    chemotype = sapply(mt_codon, define_chemotype),
    mean_ddG = mean_K13,
    assay = "K13",
    pattern
  )],
  sig_mut[, .(
    Pos,
    chemotype = sapply(mt_codon, define_chemotype),
    mean_ddG = mean_RAF1,
    assay = "RAF1",
    pattern
  )]
))


heat_sig <- long_sig[
  !is.na(chemotype),
  .(mean_ddG = mean(mean_ddG, na.rm = TRUE)),
  by = .(Pos, chemotype, assay, pattern)
]



chemotype_levels <- c(
  "Aromatic",
  "Aliphatic",
  "Polar uncharged",
  "Positive",
  "Negative",
  "Special"
)

assay_levels <- c("K13", "RAF1")

heat_sig[, chemotype := factor(chemotype, levels = chemotype_levels)]
heat_sig[, assay := factor(assay, levels = assay_levels)]



heat_sig[, x_var := paste(chemotype, assay, sep = " | ")]

x_levels <- as.vector(
  unlist(
    lapply(
      chemotype_levels,
      function(ct) paste(ct, assay_levels, sep = " | ")
    )
  )
)

heat_sig[, x_var := factor(x_var, levels = x_levels)]




GTP_pocket <- c(
  12, 13, 14, 15, 16, 17, 18,
  28, 29, 30, 32, 34, 35,
  57, 60, 61,
  116, 117, 119, 120,
  145, 146, 147
)



# 准备星号数据
star_df <- unique(heat_sig[Pos %in% GTP_pocket, .(Pos, pattern)])

# 使用geom_point的版本（星号形状不如字符好看）
p_sig <- ggplot(
  heat_sig,
  aes(x = x_var, y = factor(Pos), fill = mean_ddG)
) +
  geom_tile(color = "grey80") +
  
  geom_point(
    data = star_df,
    aes(x = 0.3, y = factor(Pos)),  # 放在绘图区域的左侧
    shape = 8,  # 星号形状
    color = "#FFD700",
    size = 1,  # 增加大小
    stroke = 1.5,  # 增加边框宽度（相当于加粗）
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ pattern, nrow = 1) +
  scale_fill_gradient2(
    low = "#1B38A6",
    mid = "white",
    high = "#F4270C",
    midpoint = 0,
    name = expression(Delta*Delta*G)
  ) +
  scale_x_discrete(expand = expansion(add = c(1, 0.5))) +  # 左侧扩展空间
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5.5, 5.5, 5.5, 40)
  ) +
  labs(
    x = "Mutation chemotype | Assay",
    y = "KRAS position",
    title = "Significant allosteric mutations: chemotype-averaged effects"
  )

p_sig


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251226/sig_muts_heatmap by patterns and region.pdf", 
                device = cairo_pdf, height = 10, width = 16, dpi = 300)
