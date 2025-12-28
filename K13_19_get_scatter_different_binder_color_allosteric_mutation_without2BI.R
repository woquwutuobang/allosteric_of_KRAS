library(data.table)
library(ggplot2)
library(dplyr)
library(krasddpcams)

# p值计算函数
krasddpcams__pvalue <- function(x, std) {
  2 * pnorm(-abs(x/std))
}

# 判断mutation_type的函数
classify_mutation_type <- function(input, assay_sele, anno, wt_aa) {
  # Read and process ddG data
  ddG <- fread(input)
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
  
  # Merge with annotation data
  anno_data <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  # Define binding site types
  data_plot[, binding_type := "allosteric site"]
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
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
    krasddpcams__pvalue(abs(`mean_kcal/mol`) - reg_threshold, `std_kcal/mol`),
    method = "BH") < 0.05 & (abs(`mean_kcal/mol`) - reg_threshold) > 0]
  
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
  data_plot_mutation <- within(data_plot_mutation,
                               mutation_type <- factor(mutation_type,
                                                       levels = c("Orthosteric site huge differences",
                                                                  "Orthosteric site small differences",
                                                                  "GTP binding allosteric mutation",
                                                                  "GTP binding other mutation",
                                                                  "Allosteric mutation",
                                                                  "Other mutation")))
  
  # 添加assay标识
  data_plot_mutation[, assay := assay_sele]
  
  return(data_plot_mutation)
}

plot_assay_correlation_no_interface <- function(input1, input2, assay1, assay2, anno, wt_aa, 
                                                point_size = 2, alpha = 0.7, base_size = 10,
                                                xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3),
                                                include_GTP_binding = TRUE) {
  
  # 获取两个assay的mutation type数据
  data_assay1 <- classify_mutation_type(input1, assay1, anno, wt_aa)
  data_assay2 <- classify_mutation_type(input2, assay2, anno, wt_aa)
  
  # 选择需要的列并重命名
  data1_clean <- data_assay1[, .(mt, Pos_real, `mean_kcal/mol`, mutation_type, site_type)]
  setnames(data1_clean, "mean_kcal/mol", paste0("mean_", assay1))
  setnames(data1_clean, "mutation_type", paste0("mutation_type_", assay1))
  setnames(data1_clean, "site_type", paste0("site_type_", assay1))
  
  data2_clean <- data_assay2[, .(mt, Pos_real, `mean_kcal/mol`, mutation_type, site_type)]
  setnames(data2_clean, "mean_kcal/mol", paste0("mean_", assay2))
  setnames(data2_clean, "mutation_type", paste0("mutation_type_", assay2))
  setnames(data2_clean, "site_type", paste0("site_type_", assay2))
  
  # 合并两个assay的数据
  merged_data <- merge(data1_clean, data2_clean, by = c("mt", "Pos_real"), all = TRUE)
  
  # ===============================
  # 筛选：去掉结合界面（Binding interface site），保留核苷酸结合口袋
  # ===============================
  filtered_data <- merged_data[
    # 保留所有非Binding interface site的位点
    !(get(paste0("site_type_", assay1)) == "Binding interface site") &
      !(get(paste0("site_type_", assay2)) == "Binding interface site")
  ]
  
  print(paste("原始数据点数:", nrow(merged_data)))
  print(paste("去结合界面后点数:", nrow(filtered_data)))
  
  # 统计不同类型的位点数量
  site_type_counts <- filtered_data[, .(
    N_total = .N,
    N_GTP1 = sum(get(paste0("site_type_", assay1)) == "GTP binding interface site", na.rm = TRUE),
    N_GTP2 = sum(get(paste0("site_type_", assay2)) == "GTP binding interface site", na.rm = TRUE),
    N_Reminder1 = sum(get(paste0("site_type_", assay1)) == "Reminder", na.rm = TRUE),
    N_Reminder2 = sum(get(paste0("site_type_", assay2)) == "Reminder", na.rm = TRUE)
  )]
  
  print("过滤后位点类型统计:")
  print(site_type_counts)
  
  # ===============================
  # 计算相关性（去除NA）
  # ===============================
  x_col <- paste0("mean_", assay1)
  y_col <- paste0("mean_", assay2)
  
  cor_data <- filtered_data[!is.na(get(x_col)) & !is.na(get(y_col))]
  
  print(paste("用于相关性计算的点数:", nrow(cor_data)))
  
  cor_test <- cor.test(cor_data[[x_col]], cor_data[[y_col]], method = "pearson")
  R_value <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # 格式化显示
  R_label <- sprintf("R = %.2f", R_value)
  p_label <- ifelse(p_value < 2.2e-16, "p < 2.2e-16", 
                    ifelse(p_value < 0.001, "p < 0.001", 
                           sprintf("p = %.3f", p_value)))
  cor_label <- paste(R_label, p_label, sep = ", ")
  
  # ===============================
  # 定义颜色分组（包括核苷酸结合口袋的变构突变）
  # ===============================
  filtered_data[, color_group := "Other"]
  
  # 1. assay1的变构突变（包括普通变构和GTP结合变构）
  filtered_data[
    (get(paste0("mutation_type_", assay1)) %in% c("Allosteric mutation", "GTP binding allosteric mutation")) &
      !(get(paste0("mutation_type_", assay2)) %in% c("Allosteric mutation", "GTP binding allosteric mutation")), 
    color_group := paste0(assay1, " allosteric")
  ]
  
  # 2. assay2的变构突变（包括普通变构和GTP结合变构）
  filtered_data[
    (get(paste0("mutation_type_", assay2)) %in% c("Allosteric mutation", "GTP binding allosteric mutation")) &
      !(get(paste0("mutation_type_", assay1)) %in% c("Allosteric mutation", "GTP binding allosteric mutation")), 
    color_group := paste0(assay2, " allosteric")
  ]
  
  # 3. 重叠的变构突变
  filtered_data[
    (get(paste0("mutation_type_", assay1)) %in% c("Allosteric mutation", "GTP binding allosteric mutation")) &
      (get(paste0("mutation_type_", assay2)) %in% c("Allosteric mutation", "GTP binding allosteric mutation")), 
    color_group := "Overlap allosteric"
  ]
  
  # 设置颜色分组和颜色映射
  color_levels <- c("Other", 
                    paste0(assay1, " allosteric"), 
                    paste0(assay2, " allosteric"), 
                    "Overlap allosteric")
  
  color_map <- setNames(
    c("grey", "#FF0066", "#A31300", "#C68EFD"),
    color_levels
  )
  
  filtered_data[, color_group := factor(color_group, levels = color_levels)]
  filtered_data <- filtered_data[order(as.numeric(color_group))]
  
  # 统计各组的数量
  print("颜色分组统计:")
  print(table(filtered_data$color_group))
  
  # 详细统计不同类型的突变
  mutation_detail <- filtered_data[, .(
    N_assay1_allosteric = sum(get(paste0("mutation_type_", assay1)) == "Allosteric mutation", na.rm = TRUE),
    N_assay1_GTP_allosteric = sum(get(paste0("mutation_type_", assay1)) == "GTP binding allosteric mutation", na.rm = TRUE),
    N_assay2_allosteric = sum(get(paste0("mutation_type_", assay2)) == "Allosteric mutation", na.rm = TRUE),
    N_assay2_GTP_allosteric = sum(get(paste0("mutation_type_", assay2)) == "GTP binding allosteric mutation", na.rm = TRUE)
  )]
  
  print("变构突变详细统计:")
  print(mutation_detail)
  
  # ===============================
  # 绘图 + 相关性标注
  # ===============================
  p <- ggplot(filtered_data, aes(x = .data[[x_col]], 
                                 y = .data[[y_col]], 
                                 color = color_group)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = color_map, name = "Allosteric Mutations",
                       breaks = color_levels[color_levels != "Other"],
                       labels = c(
                         paste0(assay1, " allosteric\n(inc. GTP binding)"),
                         paste0(assay2, " allosteric\n(inc. GTP binding)"),
                         "Overlap allosteric\n(inc. GTP binding)"
                       )) +
    labs(
      x = bquote("Binding"~Delta*Delta*"G ("*.(assay1)*") (kcal/mol)"), 
      y = bquote("Binding"~Delta*Delta*"G ("*.(assay2)*") (kcal/mol)"),
      subtitle = paste0("Excluding binding interface sites (n = ", nrow(cor_data), ")")
    ) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = base_size - 2),
      axis.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 3),
      legend.title = element_text(size = base_size - 1),
      legend.key.size = unit(0.5, "cm"),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, face = "italic")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    
    # === 添加相关性文字 ===
    annotate("text", 
             x = xlim[1] + 0.05 * diff(xlim), 
             y = ylim[2] - 0.05 * diff(ylim), 
             label = cor_label,
             hjust = 0, vjust = 1, size = 3.5)
  
  # 调试信息
  print(paste("去结合界面后Pearson R:", R_value))
  print(paste("去结合界面后P value:", p_value))
  
  return(list(
    plot = p,
    data = filtered_data,
    R = R_value,
    P = p_value,
    n_points = nrow(cor_data),
    color_distribution = table(filtered_data$color_group)
  ))
}


# 使用示例
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

########=====
# 比较K13和RAF1（去结合界面，保留核苷酸口袋）
result_k13_raf1 <- plot_assay_correlation_no_interface(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K13",
  assay2 = "RAF1",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  wt_aa = wt_aa
)

# 显示图形
print(result_k13_raf1$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f3/20251127/K13_vs_RAF1_non_binding_interface_correlation.pdf", 
       plot = result_k13_raf1$plot, 
       width = 4, 
       height = 4.5,  # 稍微增加高度以适应更长的图例
       device = cairo_pdf)

########=====
# 比较K13和K19（去结合界面，保留核苷酸口袋）
result_k13_k19 <- plot_assay_correlation_no_interface(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay1 = "K13",
  assay2 = "K19",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa
)

# 显示图形
print(result_k13_k19$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f3/20251127/K13_vs_K19_non_binding_interface_correlation.pdf", 
       plot = result_k13_k19$plot, 
       width = 4, 
       height = 4.5,
       device = cairo_pdf)

########=====
# 比较K27和RAF1（去结合界面，保留核苷酸口袋）
result_k27_raf1 <- plot_assay_correlation_no_interface(
  input1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  input2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay1 = "K27",
  assay2 = "RAF1",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  wt_aa = wt_aa
)

# 显示图形
print(result_k27_raf1$plot)

# 保存图形
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f3/20251127/K27_vs_RAF1_non_binding_interface_correlation.pdf", 
       plot = result_k27_raf1$plot, 
       width = 4, 
       height = 4.5,
       device = cairo_pdf)

