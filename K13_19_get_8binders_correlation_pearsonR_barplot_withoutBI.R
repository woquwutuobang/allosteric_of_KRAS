##########without BI
library(data.table)
library(ggplot2)
library(dplyr)

# =========================================================
# 0. 定义 binding interface map（每个 assay 的 interface）
# =========================================================

binding_sites_map <- list(
  RAF1 = c(21, 25, 29,31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67),
  PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73),
  SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153),
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)


# =========================================================
# 1. 定义 p 值计算函数（之前版本）
# =========================================================
krasddpcams__pvalue <- function(x, std) {
  2 * pnorm(-abs(x / std))
}

# =========================================================
# 2. 定义 classify_mutation_type 函数
# =========================================================
classify_mutation_type <- function(input, assay_sele, anno, wt_aa) {
  ddG <- fread(input)
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]]))
  )
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
  
  anno_data <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all.x = TRUE)
  data_plot[, binding_type := "allosteric site"]
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type == "binding site", site_type := "Binding interface site"]
  
  data_plot_mutation <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation[Pos > 1 & !is.na(mt)]
  data_plot_mutation[, mutation_type := "Other mutation"]
  
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(`mean_kcal/mol`) - reg_threshold, `std_kcal/mol`),
    method = "BH") < 0.05 & (abs(`mean_kcal/mol`) - reg_threshold) > 0]
  
  data_plot_mutation[site_type == "Binding interface site" & allosteric_mutation == TRUE,
                     mutation_type := "Orthosteric site huge differences"]
  data_plot_mutation[site_type != "Binding interface site" & allosteric_mutation == TRUE,
                     mutation_type := "Allosteric mutation"]
  
  data_plot_mutation[, assay := assay_sele]
  return(data_plot_mutation)
}



calculate_correlation <- function(input_x, input_y, assay_x, assay_y, anno, wt_aa) {
  
  # 先分类、计算ΔΔG
  data_x <- classify_mutation_type(input_x, assay_x, anno, wt_aa)
  data_y <- classify_mutation_type(input_y, assay_y, anno, wt_aa)
  
  # -----------------------------
  # 第一步：过滤掉 assay_x 的 binding interface 位点突变
  # -----------------------------
  if (assay_x %in% names(binding_sites_map)) {
    data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
  }
  
  # -----------------------------
  # 第二步：过滤掉 assay_y 的 binding interface 位点突变
  # -----------------------------
  if (assay_y %in% names(binding_sites_map)) {
    data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
  }
  
  # -----------------------------
  # 第三步：只保留共享突变用于 Pearson 计算
  # -----------------------------
  data_x_clean <- data_x[, .(mt, Pos_real, `mean_kcal/mol`)]
  setnames(data_x_clean, "mean_kcal/mol", paste0("mean_", assay_x))
  
  data_y_clean <- data_y[, .(mt, Pos_real, `mean_kcal/mol`)]
  setnames(data_y_clean, "mean_kcal/mol", paste0("mean_", assay_y))
  
  merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
  
  # -----------------------------
  # 第四步：计算 Pearson R
  # -----------------------------
  cor_test <- cor.test(
    merged_data[[paste0("mean_", assay_x)]],
    merged_data[[paste0("mean_", assay_y)]]
  )
  
  return(list(
    r = round(cor_test$estimate, 3),
    p = cor_test$p.value
  ))
}




# =========================================================
# 4. 计算所有 assay 对的 R/p
# =========================================================
calculate_all_correlations <- function(input_files, anno, wt_aa) {
  assays <- names(input_files)
  res_list <- list()
  
  for (i in 1:(length(assays)-1)) {
    for (j in (i+1):length(assays)) {
      cor_res <- calculate_correlation(
        input_x = input_files[[assays[i]]],
        input_y = input_files[[assays[j]]],
        assay_x = assays[i],
        assay_y = assays[j],
        anno = anno,
        wt_aa = wt_aa
      )
      
      res_list[[paste0(assays[i], "_vs_", assays[j])]] <- data.table(
        assay1 = assays[i],
        assay2 = assays[j],
        R = cor_res$r,
        p = cor_res$p
      )
    }
  }
  
  do.call(rbind, res_list)
}

# =========================================================
# 5. 构建排序和显著性标签
# =========================================================
prepare_correlation_data <- function(cor_dt) {
  BI1 <- c("RAF1", "RALGDS", "PI3KCG", "SOS1","K55", "K27")
  BI2 <- c("K13", "K19")
  darpin <- c("K55", "K27", "K13", "K19")
  active <- c("RAF1", "RALGDS", "PI3KCG", "K55")
  
  cor_dt[, Comparison := paste0(assay1, " vs ", assay2)]
  
  cor_dt[, Type1 := ifelse(assay1 %in% darpin, "DARPin", "non-DARPin")]
  cor_dt[, Type2 := ifelse(assay2 %in% darpin, "DARPin", "non-DARPin")]
  cor_dt[, Act1 := ifelse(assay1 %in% active, "Active", "Inactive")]
  cor_dt[, Act2 := ifelse(assay2 %in% active, "Active", "Inactive")]
  
  cor_dt[, DARPinClass := fifelse(Type1=="non-DARPin" & Type2=="non-DARPin", "nonDARPin_vs_nonDARPin",
                                  fifelse(Type1=="DARPin" & Type2=="DARPin","DARPin_vs_DARPin",
                                          "nonDARPin_vs_DARPin"))]
  cor_dt[, ActivityClass := fifelse(Act1=="Active" & Act2=="Active", "Active_vs_Active",
                                    fifelse(Act1=="Inactive" & Act2=="Inactive","Inactive_vs_Inactive",
                                            "Active_vs_Inactive"))]
  
  cor_dt[, comparison_type := fifelse(assay1 %in% BI1 & assay2 %in% BI1, "BI1_vs_BI1",
                                      fifelse(assay1 %in% BI2 & assay2 %in% BI2, "BI2_vs_BI2",
                                              "BI1_vs_BI2"))]
  
  cor_dt[, significance := fcase(
    p < 0.001, "***",
    p < 0.01,  "**",
    p < 0.05,  "*",
    default = "NS"
  )]
  
  
  cor_dt <- cor_dt[order(
    factor(comparison_type, levels=c("BI1_vs_BI1","BI1_vs_BI2","BI2_vs_BI2")),
    factor(DARPinClass, levels=c("nonDARPin_vs_nonDARPin","nonDARPin_vs_DARPin","DARPin_vs_DARPin")),
    factor(ActivityClass, levels=c("Active_vs_Active","Active_vs_Inactive","Inactive_vs_Inactive")),
    factor(Comparison, levels=unique(Comparison))
  )]
  
  cor_dt[, Comparison := factor(Comparison, levels=Comparison)]
  return(cor_dt)
}

# =========================================================
# 6. 绘制柱状图
# =========================================================
plot_correlation_bar <- function(cor_dt) {
  ggplot(cor_dt, aes(x=Comparison, y=R, fill=comparison_type)) +
    geom_bar(stat="identity", width=0.7, alpha=0.85) +
    geom_text(aes(label=significance), vjust=-0.5, size=3, color="#F1DD10") +
    scale_fill_manual(values=c(
      "BI1_vs_BI1"="#F4270C",
      "BI1_vs_BI2"="#F4AD0C",
      "BI2_vs_BI2"="#1B38A6"
    )) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=10),
          legend.position = "none") +
    ylab("Pearson R")
}

# =========================================================
# 7. 示例运行
# =========================================================
input_files <- list(
  RAF1   = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  RALGDS = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAL.txt",
  PI3KCG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_PI3.txt",
  SOS1   = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_SOS.txt",
  K55    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K55.txt",
  K27    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",
  K13    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  K19    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"
)

anno <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
wt_aa <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQYRMKKLNSSDDGTQGCMGLPCVVM"

cor_dt <- calculate_all_correlations(input_files, anno, wt_aa)
cor_dt <- prepare_correlation_data(cor_dt)
p_bar <- plot_correlation_bar(cor_dt)
print(p_bar)



ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figure 3/28Rvalue_barplot_withoutBI_8binder.pdf",
  plot = p_bar,
  device = cairo_pdf,
  width = 12,
  height = 8
)

