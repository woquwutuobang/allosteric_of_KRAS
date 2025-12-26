#### 20251119   without Binding interface

### version 5 - filtered by binding interface sites
# =========================================================
# KRAS: 下三角散点 + 右上角 R+p 面板（去掉 binding interface 位点后比较）
# =========================================================
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtext)

# =========================================================
# binding interface 位点映射（按你提供的列表）
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
# 1. p 值计算函数
# =========================================================
krasddpcams__pvalue <- function(x, std) {
  2 * pnorm(-abs(x / std))
}

# =========================================================
# 2. 分类函数 (保持原样)
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

# =========================================================
# 3. 计算相关系数和p值（在合并前去除 binding interface 位点）
# =========================================================
calculate_correlation <- function(input_x, input_y, assay_x, assay_y, anno, wt_aa) {
  data_x <- classify_mutation_type(input_x, assay_x, anno, wt_aa)
  data_y <- classify_mutation_type(input_y, assay_y, anno, wt_aa)
  
  # 去掉各自 assay 的 binding interface 位点
  if (!is.null(binding_sites_map[[assay_x]])) {
    data_x <- data_x[ !(Pos_real %in% binding_sites_map[[assay_x]]) ]
  }
  if (!is.null(binding_sites_map[[assay_y]])) {
    data_y <- data_y[ !(Pos_real %in% binding_sites_map[[assay_y]]) ]
  }
  
  data_x_clean <- data_x[, .(mt, Pos_real, `mean_kcal/mol`)]
  setnames(data_x_clean, "mean_kcal/mol", paste0("mean_", assay_x))
  data_y_clean <- data_y[, .(mt, Pos_real, `mean_kcal/mol`)]
  setnames(data_y_clean, "mean_kcal/mol", paste0("mean_", assay_y))
  
  merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
  
  # 若合并后数据点不足，返回 NA
  if (nrow(merged_data) < 3) {
    return(list(r = NA, p = NA))
  }
  
  cor_test <- cor.test(merged_data[[paste0("mean_", assay_x)]], 
                       merged_data[[paste0("mean_", assay_y)]])
  
  return(list(
    r = round(cor_test$estimate, 3),
    p = cor_test$p.value
  ))
}

# =========================================================
# 4. 散点图函数（左下三角）- 修复坐标轴标题并去除 binding interface 位点
# =========================================================
plot_assay_correlation_fixed <- function(input_x, input_y, assay_x, assay_y, anno, wt_aa,
                                         point_size = 1.5, alpha = 0.6, base_size = 9,
                                         xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3),
                                         show_x_title = FALSE, show_y_title = FALSE) {
  data_x <- classify_mutation_type(input_x, assay_x, anno, wt_aa)
  data_y <- classify_mutation_type(input_y, assay_y, anno, wt_aa)
  
  # 去掉 binding interface 位点
  if (!is.null(binding_sites_map[[assay_x]])) {
    data_x <- data_x[ !(Pos_real %in% binding_sites_map[[assay_x]]) ]
  }
  if (!is.null(binding_sites_map[[assay_y]])) {
    data_y <- data_y[ !(Pos_real %in% binding_sites_map[[assay_y]]) ]
  }
  
  data_x_clean <- data_x[, .(mt, Pos_real, `mean_kcal/mol`, mutation_type)]
  setnames(data_x_clean, "mean_kcal/mol", paste0("mean_", assay_x))
  setnames(data_x_clean, "mutation_type", paste0("mutation_type_", assay_x))
  data_y_clean <- data_y[, .(mt, Pos_real, `mean_kcal/mol`, mutation_type)]
  setnames(data_y_clean, "mean_kcal/mol", paste0("mean_", assay_y))
  setnames(data_y_clean, "mutation_type", paste0("mutation_type_", assay_y))
  
  merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
  
  merged_data[, color_group := "Other"]
  merged_data[get(paste0("mutation_type_", assay_x)) == "Allosteric mutation" &
                !(get(paste0("mutation_type_", assay_y)) == "Allosteric mutation"), color_group := "X_allosteric"]
  merged_data[!(get(paste0("mutation_type_", assay_x)) == "Allosteric mutation") &
                get(paste0("mutation_type_", assay_y)) == "Allosteric mutation", color_group := "Y_allosteric"]
  merged_data[get(paste0("mutation_type_", assay_x)) == "Allosteric mutation" &
                get(paste0("mutation_type_", assay_y)) == "Allosteric mutation", color_group := "Overlap_allosteric"]
  
  color_levels <- c("Other","X_allosteric","Y_allosteric","Overlap_allosteric")
  color_map <- c("grey","#FF0066","#A31300","#C68EFD")
  names(color_map) <- color_levels
  merged_data[, color_group := factor(color_group, levels=color_levels)]
  
  p <- ggplot(merged_data, aes(
    x = .data[[paste0("mean_", assay_x)]],
    y = .data[[paste0("mean_", assay_y)]],
    color = color_group
  )) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = color_map) +
    labs(x = if(show_x_title) paste0("ΔΔG (", assay_x, ")") else NULL,
         y = if(show_y_title) paste0("ΔΔG (", assay_y, ")") else NULL) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = base_size-2, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = base_size-2),
      axis.title.x = if(show_x_title) element_text(size = base_size-1, margin = margin(t = 5)) else element_blank(),
      axis.title.y = if(show_y_title) element_text(size = base_size-1, margin = margin(r = 5)) else element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    coord_cartesian(xlim = xlim, ylim = ylim)
  
  return(p)
}

# =========================================================
# 5. R值和p值文本图函数（右上三角） - 保持使用 calculate_correlation（已过滤）
# =========================================================
plot_correlation_text <- function(input_x, input_y, assay_x, assay_y, anno, wt_aa, base_size = 9) {
  cor_result <- calculate_correlation(input_x, input_y, assay_x, assay_y, anno, wt_aa)
  
  p_value_text <- ifelse(is.na(cor_result$p), "p = NA", ifelse(cor_result$p < 0.001, "p < 0.001", 
                                                               paste0("p = ", round(cor_result$p, 3))))
  r_text <- ifelse(is.na(cor_result$r), "R = NA", paste0("R = ", cor_result$r))
  
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.6, 
             label = r_text, 
             size = base_size / 2.5, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.4, 
             label = p_value_text, 
             size = base_size / 3) +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 5),
          plot.background = element_rect(fill = "white", color = NA))
  
  return(p)
}

# =========================================================
# 6. 创建assay名称标签函数 - 修改为对齐版本
# =========================================================
create_assay_label <- function(assay_name, orientation = "horizontal", base_size = 9) {
  if (orientation == "horizontal") {
    p <- ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, 
               fill = "#87CEEB", alpha = 0.8, color = "grey80") +
      annotate("text", x = 0.5, y = 0.5, label = assay_name, 
               size = base_size / 2.5) +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 2))
  } else { # vertical
    p <- ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, 
               fill = "#87CEEB", alpha = 0.8, color = "grey80") +
      annotate("text", x = 0.5, y = 0.5, label = assay_name, 
               size = base_size / 2.5, angle = 90) +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 2))
  }
  return(p)
}

# =========================================================
# 7. 创建完整矩阵图 - 修复对齐问题（并使用过滤版绘图/计算）
# =========================================================
create_all_assay_comparisons_full <- function() {
  input_files <- list(
    RAF1="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
    RALGDS="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt",
    PI3KCG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt",
    SOS1="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt",
    K55="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt",
    K27="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
    K13="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
    K19="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt"
  )
  anno <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
  wt_aa <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPARTVETRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQYRMKKLNSSDDGTQGCMGLPCVVM"
  
  assays <- names(input_files)
  n <- length(assays)
  
  # 创建空图
  blank_plot <- ggplot() + theme_void() + theme(plot.margin = margin(1, 1, 1, 1))
  
  # 创建所有图形组件
  all_plots <- list()
  
  # 1. 创建主矩阵图形
  for (i in 1:n) {
    for (j in 1:n) {
      if (i > j) { # 左下三角：散点图 (filtered)
        # 显示坐标轴标题的条件：最后一行显示X轴标题，第一列显示Y轴标题
        show_x_title <- (i == n)  # 最后一行显示X轴标题
        show_y_title <- (j == 1)  # 第一列显示Y轴标题
        
        all_plots[[paste0("scatter_", i, "_", j)]] <- plot_assay_correlation_fixed(
          input_x = input_files[[assays[j]]], 
          input_y = input_files[[assays[i]]],
          assay_x = assays[j], 
          assay_y = assays[i],
          anno = anno, 
          wt_aa = wt_aa,
          show_x_title = show_x_title,
          show_y_title = show_y_title
        )
      } else if (i < j) { # 右上三角：R值文本 (filtered via calculate_correlation)
        all_plots[[paste0("text_", i, "_", j)]] <- plot_correlation_text(
          input_x = input_files[[assays[j]]], 
          input_y = input_files[[assays[i]]],
          assay_x = assays[j], 
          assay_y = assays[i],
          anno = anno, 
          wt_aa = wt_aa
        )
      } else { # 对角：空白
        all_plots[[paste0("diag_", i, "_", j)]] <- blank_plot
      }
    }
  }
  
  # 2. 创建标签图形
  top_labels <- lapply(assays, function(x) create_assay_label(x, "horizontal"))
  bottom_labels <- lapply(assays, function(x) create_assay_label(x, "horizontal"))
  left_labels <- lapply(assays, function(x) create_assay_label(x, "vertical"))
  right_labels <- lapply(assays, function(x) create_assay_label(x, "vertical"))
  
  # 3. 构建主矩阵
  main_matrix_rows <- list()
  for (i in 1:n) {
    row_plots <- list()
    for (j in 1:n) {
      if (i > j) {
        row_plots[[j]] <- all_plots[[paste0("scatter_", i, "_", j)]]
      } else if (i < j) {
        row_plots[[j]] <- all_plots[[paste0("text_", i, "_", j)]]
      } else {
        row_plots[[j]] <- all_plots[[paste0("diag_", i, "_", j)]]
      }
    }
    main_matrix_rows[[i]] <- wrap_plots(row_plots, ncol = n)
  }
  
  # 4. 验证散点图数量
  scatter_count <- sum(sapply(names(all_plots), function(x) grepl("^scatter", x)))
  cat("散点图数量:", scatter_count, "/", choose(n, 2), "\n")
  
  # 5. 组合所有部分 - 完全修正版本
  # 顶部标签行
  top_row <- wrap_plots(
    c(list(blank_plot), top_labels, list(blank_plot)), 
    ncol = n + 2, 
    widths = c(0.4, rep(0.9, n), 0.4)
  )
  
  # 底部标签行
  bottom_row <- wrap_plots(
    c(list(blank_plot), bottom_labels, list(blank_plot)), 
    ncol = n + 2, 
    widths = c(0.4, rep(0.9, n), 0.4)
  )
  
  # 中间行（左侧标签 + 主网格 + 右侧标签）
  middle_rows <- list()
  for (i in 1:n) {
    middle_rows[[i]] <- wrap_plots(
      c(list(left_labels[[i]]), list(main_matrix_rows[[i]]), list(right_labels[[i]])), 
      ncol = 3, 
      widths = c(0.05, 0.9, 0.05)
    )
  }
  middle_grid <- wrap_plots(middle_rows, ncol = 1)
  
  # 最终组合
  final_plot <- wrap_plots(
    list(top_row, middle_grid, bottom_row), 
    ncol = 1, 
    heights = c(0.05, 0.9, 0.05)
  )
  
  return(final_plot)
}

# =========================================================
# 8. 运行
# =========================================================
final_plot <- create_all_assay_comparisons_full()
print(final_plot)


cat("Saving PNG...\n")
png("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/sf4/20251119/sf4b_correlation_matrix_complete_aligned_manual_binding.png",
    width = 20, height = 20, units = "in", res = 300)
print(final_plot)
dev.off()
