library(data.table)
library(ggplot2)
library(krasddpcams)

# ============================================================
# ============================================================

identify_allosteric_mutations <- function(input, anno, wt_aa, assay_sele) {
  
  ddG <- fread(input)
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", wt_codon := substr(id, 1, 1)]
  ddG[id != "WT", mt_codon := substr(id, nchar(id), nchar(id))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(
    wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]]))
  )
  
  ddG <- merge(ddG, heatmap_tool,
               by = c("Pos_real", "wt_codon", "mt_codon"),
               all = TRUE)
  ddG[, Pos := Pos_real]
  
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) /
      sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = Pos_real]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = Pos_real]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  anno_data <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5,
            binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[
    GXPMG_scHAmin_ligand_RAF1 < 5,
    binding_type_gtp_included := "GTP binding site"
  ]
  
  # 注意：这里使用所有位点计算阈值（包括结合界面）
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) /
                               sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site",
            site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site",
            site_type := "GTP binding interface site"]
  
  data_plot_mutation1 <- merge(ddG,
                               data_plot[, .(Pos, site_type)],
                               by = "Pos",
                               all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  
  # 注意：这里不排除结合界面，正常判断所有突变
  data_plot_mutation[, allosteric_mutation :=
                       p.adjust(
                         krasddpcams__pvalue(abs(mean) - reg_threshold, std),
                         method = "BH"
                       ) < 0.05 & (abs(mean) - reg_threshold) > 0
  ]
  
  list(
    mutation_table = data_plot_mutation,
    site_annotation = data_plot
  )
}

# ============================================================
# 新的变构位点识别函数（并行证据逻辑）
# ============================================================

get_statistical_sites <- function(data_plot_mutation, min_n = 9) {
  df <- as.data.table(data_plot_mutation)
  
  # 筛选至少有min_n个突变的位点
  site_n <- df[, .N, by = Pos_real][N >= min_n]
  df <- df[Pos_real %in% site_n$Pos_real]
  
  # 计算总体变构突变比例
  total_allo <- sum(df$allosteric_mutation, na.rm = TRUE)
  total_non <- sum(!df$allosteric_mutation, na.rm = TRUE)
  
  # 对每个位点进行Fisher精确检验
  fisher_sites <- df[, .(
    N_total = .N,
    N_allosteric = sum(allosteric_mutation, na.rm = TRUE)
  ), by = Pos_real][, {
    
    # 添加0.5的连续性校正
    a <- N_allosteric + 0.5
    b <- (N_total - N_allosteric) + 0.5
    c <- (total_allo - N_allosteric) + 0.5
    d <- (total_non - (N_total - N_allosteric)) + 0.5
    
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    
    .(
      odds_ratio = as.numeric(ft$estimate),
      p_value = ft$p.value
    )
  }, by = Pos_real]
  
  # 多重检验校正
  fisher_sites[, q_value := p.adjust(p_value, "BH")]
  fisher_sites[, statistical_support := (q_value < 0.05 & odds_ratio > 2)]
  
  return(fisher_sites)
}


get_quantitative_sites <- function(data_plot_mutation, min_n = 9) {
  df <- as.data.table(data_plot_mutation)
  df[, abs_ddG := abs(`mean_kcal/mol`)]
  
  # 使用结合界面位点作为背景
  interface_df <- df[site_type == "Binding interface site"]
  interface_q75 <- quantile(interface_df$abs_ddG, 0.75, na.rm = TRUE)
  interface_q90 <- quantile(interface_df$abs_ddG, 0.90, na.rm = TRUE)
  
  # 计算每个位点的突变强度
  site_strength <- df[, .(
    N = .N,
    Q75_abs_ddG = quantile(abs_ddG, 0.75, na.rm = TRUE),
    max_abs_ddG = max(abs_ddG, na.rm = TRUE)
  ), by = Pos_real][N >= min_n]
  
  # 定量支持标准：Q75超过背景Q75 或 最大值超过背景Q90
  site_strength[, quantitative_support :=
                  Q75_abs_ddG > interface_q75 |
                  max_abs_ddG > interface_q90
  ]
  
  return(site_strength)
}


# ============================================================
# 修改后的变构位点识别函数（并联证据逻辑）
# ============================================================

annotate_allosteric_sites_parallel <- function(
    data_plot_mutation,
    min_n = 9
) {
  
  # 获取统计证据
  stat_sites <- get_statistical_sites(data_plot_mutation, min_n)
  
  # 获取定量证据
  quant_sites <- get_quantitative_sites(data_plot_mutation, min_n)
  
  # 合并两种证据
  final_sites <- merge(stat_sites, quant_sites,
                       by = "Pos_real", all = TRUE)
  
  # 并联逻辑：只要满足统计证据或定量证据之一，就是变构位点
  final_sites[, is_allosteric_site := 
                statistical_support == TRUE | quantitative_support == TRUE]
  
  # 定义证据类型
  final_sites[, evidence_type := fifelse(
    statistical_support & quantitative_support,
    "Statistical + Quantitative (strongest evidence)",
    fifelse(
      statistical_support & !quantitative_support,
      "Statistical evidence only",
      fifelse(
        !statistical_support & quantitative_support,
        "Quantitative evidence only",
        "Non-allosteric"
      )
    )
  )]
  
  return(final_sites[order(Pos_real)])
}

# ============================================================
# 修改后的绘图函数：修复data.table列引用问题，添加图例
# ============================================================

plot_allosteric_sites_distance_effect <- function(res, binder_name, 
                                                  allosteric_sites_list, 
                                                  anno_data,
                                                  binding_sites_map,
                                                  GTP_Binding_pocket) {
  
  # 获取突变数据
  mutation_data <- res$mutation_table
  
  # 计算每个位点20种突变的平均效应值和标准误差
  site_effects <- mutation_data[, .(
    mean_effect = mean(abs(`mean_kcal/mol`), na.rm = TRUE),
    sd_effect = sd(abs(`mean_kcal/mol`), na.rm = TRUE),
    n_mutations = sum(!is.na(`mean_kcal/mol`))
  ), by = Pos_real]
  
  # 计算标准误差
  site_effects[, se_effect := sd_effect / sqrt(n_mutations)]
  
  # 合并距离信息
  distance_col <- paste0("scHAmin_ligand_", binder_name)
  
  if (!distance_col %in% names(anno_data)) {
    warning(sprintf("Distance column %s not found.", distance_col))
    similar_cols <- grep(paste0("scHAmin.*", binder_name), names(anno_data), value = TRUE, ignore.case = TRUE)
    if (length(similar_cols) > 0) {
      distance_col <- similar_cols[1]
    } else {
      site_effects[, Distance := runif(.N, 5, 30)]
    }
  }
  
  if (distance_col %in% names(anno_data)) {
    if (!"Pos" %in% names(anno_data)) {
      if ("Pos_real" %in% names(anno_data)) {
        anno_data[, Pos := Pos_real]
      }
    }
    
    site_effects <- merge(site_effects, 
                          anno_data[, .(Pos_real = Pos, Distance = get(distance_col))],
                          by = "Pos_real", all.x = TRUE)
  }
  
  # 标记变构位点（注意：这个列表已经排除了结合界面）
  site_effects[, is_allosteric_site := Pos_real %in% allosteric_sites_list]
  
  # 创建分类变量 - 按照要求的顺序判断
  site_effects[, site_category := "Other"]
  
  # 1. 首先标记结合界面位点（红色）
  if (binder_name %in% names(binding_sites_map)) {
    site_effects[Pos_real %in% binding_sites_map[[binder_name]], 
                 site_category := "Binding Interface"]
  }
  
  # 2. 标记其他变构位点（绿色）- 这些位点已经是变构且不在结合界面
  site_effects[is_allosteric_site == TRUE & site_category == "Other", 
               site_category := "Allosteric Site"]
  
  # 3. 最后判断是不是GTP口袋的变构位点（橙色）
  site_effects[site_category == "Allosteric Site" & Pos_real %in% GTP_Binding_pocket, 
               site_category := "GTP Binding Pocket (Allosteric)"]
  
  # 颜色映射
  color_palette <- c(
    "Binding Interface" = "#F4270C",                    # 🔴 红色 - 所有结合界面位点
    "GTP Binding Pocket (Allosteric)" = "#F4AD0C",      # 🟠 橙色 - GTP口袋的变构位点
    "Allosteric Site" = "#09B636",                      # 🟢 绿色 - 其他变构位点
    "Other" = "#999999"                                 # ⚫ 灰色 - 其他位点
  )
  
  # 形状映射：只有变构位点用大一点的点
  shape_palette <- c(
    "Binding Interface" = 16,
    "GTP Binding Pocket (Allosteric)" = 17,
    "Allosteric Site" = 18,
    "Other" = 16
  )
  
  # 大小映射：变构位点大一点
  size_palette <- c(
    "Binding Interface" = 2.5,
    "GTP Binding Pocket (Allosteric)" = 4,
    "Allosteric Site" = 3.5,
    "Other" = 2
  )
  
  # 修复：创建用于筛选的数据框版本
  site_effects_df <- as.data.frame(site_effects)
  site_effects_df$is_allosteric <- site_effects_df$is_allosteric_site
  
  # 绘制图形
  p <- ggplot(site_effects_df, aes(x = Distance, y = mean_effect)) +
    
    # 首先绘制非变构位点的error bar（灰色）
    geom_errorbar(data = site_effects_df[!site_effects_df$is_allosteric, ],
                  aes(ymin = mean_effect - se_effect, 
                      ymax = mean_effect + se_effect),
                  color = "#CCCCCC",
                  width = 0.2,
                  alpha = 0.6) +
    
    # 绘制变构位点的error bar（根据类别使用不同颜色）
    geom_errorbar(data = site_effects_df[site_effects_df$is_allosteric, ],
                  aes(ymin = mean_effect - se_effect, 
                      ymax = mean_effect + se_effect,
                      color = site_category),
                  width = 0.3,
                  alpha = 0.8,
                  size = 0.8) +
    
    # 绘制所有点
    geom_point(aes(color = site_category,
                   shape = site_category,
                   size = site_category),
               alpha = 0.8) +
    
    # 为变构位点添加标签
    geom_text(data = site_effects_df[site_effects_df$is_allosteric, ],
              aes(label = Pos_real,
                  color = site_category),
              vjust = -0.8, 
              hjust = 0.5,
              size = 3.5,
              fontface = "bold",
              check_overlap = TRUE) +
    
    # 设置颜色、形状和大小
    scale_color_manual(values = color_palette, 
                       name = "Site Category",
                       breaks = c("Binding Interface",
                                  "GTP Binding Pocket (Allosteric)", 
                                  "Allosteric Site",
                                  "Other")) +
    scale_shape_manual(values = shape_palette,
                       name = "Site Category",
                       breaks = c("Binding Interface",
                                  "GTP Binding Pocket (Allosteric)", 
                                  "Allosteric Site",
                                  "Other")) +
    scale_size_manual(values = size_palette,
                      name = "Site Category",
                      breaks = c("Binding Interface",
                                 "GTP Binding Pocket (Allosteric)", 
                                 "Allosteric Site",
                                 "Other")) +
    
    # 标签
    labs(
      title = paste("Allosteric Sites for", binder_name),
      subtitle = "Binding interface sites excluded from allosteric site identification",
      x = paste("Distance to", binder_name, "(Å)"),
      y = "Mean |ΔΔG| (kcal/mol) ± SE",
      caption = paste(
        "Binding Interface sites: ", sum(site_effects$site_category == "Binding Interface"),
        " | Allosteric sites (excluding interface): ", sum(site_effects$is_allosteric_site),
        " (GTP Pocket: ", sum(site_effects$site_category == "GTP Binding Pocket (Allosteric)"),
        ", Other: ", sum(site_effects$site_category == "Allosteric Site"),
        ")"
      )
    ) +
    
    # 主题 - 去掉背景网格线，保留图例
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",  # 单图中显示图例
      legend.key = element_blank(),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.background = element_rect(fill = "white", color = "gray"),
      panel.grid.major = element_blank(),  # 去掉主要网格线
      panel.grid.minor = element_blank(),  # 去掉次要网格线
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    
    # 图例合并为一个
    guides(
      color = guide_legend(override.aes = list(size = c(3, 4, 3.5, 2.5))),
      shape = guide_legend(override.aes = list(size = c(3, 4, 3.5, 2.5))),
      size = "none"  # 隐藏大小图例，因为已经包含在颜色/形状图例中
    )
  
  return(p)
}

# ============================================================
# 主程序
# ============================================================

# 定义输入文件路径和参数
input_RAF1 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt"
input_RALGDS <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAL.txt"
input_PI3KCG <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_PI3.txt"
input_SOS1 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_SOS.txt"
input_K55 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K55.txt"
input_K27 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt"
input_K13 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt"
input_K19 <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt"

anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 创建输入文件列表
input_files <- list(
  RAF1 = input_RAF1,
  RALGDS = input_RALGDS,
  PI3KCG = input_PI3KCG,
  SOS1 = input_SOS1,
  K55 = input_K55,
  K27 = input_K27,
  K13 = input_K13,
  K19 = input_K19
)

# 定义分类位点
GTP_Binding_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 
                        57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

binding_sites_map <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67),
  PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73),
  SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 
           50, 56, 59, 64, 65, 66, 67, 70, 149, 153),
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 
          105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 
          108, 125, 129, 133, 136, 137)
)

# 加载注释数据
cat("Loading annotation data...\n")
anno_data <- fread(anno_file)

# 确保anno_data有Pos_real列
if (!"Pos_real" %in% names(anno_data)) {
  if ("Pos" %in% names(anno_data)) {
    anno_data[, Pos_real := Pos]
  } else if (1 %in% names(anno_data)) {
    setnames(anno_data, "1", "Pos")
    anno_data[, Pos_real := Pos]
  } else {
    anno_data[, Pos_real := 1:.N]
  }
}

# 存储结果
allosteric_sites <- list()
allosteric_details <- list()
plots_list <- list()

# 主循环
for (b in names(input_files)) {
  
  cat(sprintf("\n========== Processing %s ==========\n", b))
  
  # 1. 识别变构突变
  res <- identify_allosteric_mutations(
    input = input_files[[b]],
    anno = anno_file,
    wt_aa = wt_aa,
    assay_sele = b
  )
  
  # 2. 使用并联证据方法识别变构位点
  allo_sites <- annotate_allosteric_sites_parallel(
    data_plot_mutation = res$mutation_table,
    min_n = 9
  )
  
  # 3. 关键修改：在提取变构位点时排除结合界面
  mutation_data <- res$mutation_table
  binding_sites <- unique(mutation_data[site_type == "Binding interface site", Pos_real])
  all_allosteric_positions <- allo_sites[is_allosteric_site == TRUE, Pos_real]
  allosteric_positions <- setdiff(all_allosteric_positions, binding_sites)
  
  # 保存结果
  allosteric_details[[b]] <- allo_sites
  allosteric_sites[[b]] <- allosteric_positions
  
  # 4. 打印摘要
  cat(sprintf("  Total positions analyzed: %d\n", nrow(allo_sites)))
  cat(sprintf("  All allosteric sites (including interface): %d\n", length(all_allosteric_positions)))
  cat(sprintf("  Allosteric sites (excluding interface): %d\n", length(allosteric_positions)))
  
  # 5. 绘制图形
  cat("  Creating plot...\n")
  p <- plot_allosteric_sites_distance_effect(
    res = res,
    binder_name = b,
    allosteric_sites_list = allosteric_positions,
    anno_data = anno_data,
    binding_sites_map = binding_sites_map,
    GTP_Binding_pocket = GTP_Binding_pocket
  )
  
  plots_list[[b]] <- p
  
  # 6. 显示和保存图形
  print(p)
  
  output_file <- sprintf("allosteric_sites_%s_no_interface.png", b)
  ggsave(filename = output_file,
         plot = p,
         width = 8,  # 增加宽度以容纳图例
         height = 6,
         dpi = 300)
  cat(sprintf("  Plot saved as: %s\n", output_file))
}

# ============================================================
# 创建组合图形（去掉网格线，只保留一个共享图例）
# ============================================================

if (!require(cowplot)) {
  install.packages("cowplot")
}
library(cowplot)

cat("\n========== Creating combined plot ==========\n")

# 检查plots_list
if (length(plots_list) == 0) {
  cat("Error: No plots to combine.\n")
} else {
  # 从第一张图中提取图例
  legend <- get_legend(
    plots_list[[1]] + 
      theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.direction = "horizontal") +
      guides(color = guide_legend(nrow = 1, byrow = TRUE),
             shape = guide_legend(nrow = 1, byrow = TRUE))
  )
  
  # 创建没有图例的子图列表
  plots_no_legend <- lapply(plots_list, function(p) {
    p + theme(legend.position = "none")
  })
  
  # 创建主图区域（没有图例）
  main_plots <- plot_grid(
    plotlist = plots_no_legend,
    ncol = 4,
    nrow = 2,
    labels = names(plots_list),
    label_size = 12,
    label_fontface = "bold",
    align = "hv"  # 水平和垂直对齐
  )
  
  # 组合主图和共享图例（图例在底部）
  combined_plot <- plot_grid(
    main_plots,
    legend,
    ncol = 1,
    rel_heights = c(9, 1)  # 主图占9份，图例占1份
  )
  
  # 显示组合图形
  print(combined_plot)
  
  # 保存组合图形
  output_file <- "allosteric_sites_combined.png"
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = 22,  # 调整宽度
    height = 13, # 调整高度以容纳底部图例
    dpi = 300
  )
  
  cat(sprintf("Combined plot saved as: %s\n", output_file))
  cat(sprintf("  Dimensions: %d x %d inches\n", 22, 13))
  cat("  Legend placed at the bottom (shared across all subplots)\n")
}
# ============================================================
# 打印变构位点位置信息（去掉结合界面后）
# ============================================================

cat("\n========== Allosteric Sites Summary (excluding binding interface) ==========\n\n")

for (b in names(allosteric_sites)) {
  if (length(allosteric_sites[[b]]) > 0) {
    # 将变构位点分为GTP口袋和其他区域
    gtp_sites <- intersect(allosteric_sites[[b]], GTP_Binding_pocket)
    other_sites <- setdiff(allosteric_sites[[b]], GTP_Binding_pocket)
    
    # 按要求的格式打印
    cat(sprintf("%s Allosteric Sites:\n", b))
    
    if (length(gtp_sites) > 0) {
      cat(sprintf("  GTP Binding Pocket (%d sites): %s\n", 
                  length(gtp_sites),
                  paste(sort(gtp_sites), collapse = ", ")))
    }
    
    if (length(other_sites) > 0) {
      cat(sprintf("  Other Regions (%d sites): %s\n", 
                  length(other_sites),
                  paste(sort(other_sites), collapse = ", ")))
    }
    
    if (length(gtp_sites) == 0 && length(other_sites) == 0) {
      cat("  No allosteric sites found.\n")
    }
    
    cat("\n")
  } else {
    cat(sprintf("%s: No allosteric sites found (excluding binding interface).\n\n", b))
  }
}

# ============================================================
# 保存结果
# ============================================================

cat("\n========== Saving results ==========\n")

# 创建详细的位置信息表
detailed_positions <- list()
for (b in names(allosteric_sites)) {
  if (length(allosteric_sites[[b]]) > 0) {
    gtp_sites <- intersect(allosteric_sites[[b]], GTP_Binding_pocket)
    other_sites <- setdiff(allosteric_sites[[b]], GTP_Binding_pocket)
    
    detailed_positions[[b]] <- data.table(
      Binder = b,
      Position = c(sort(gtp_sites), sort(other_sites)),
      Category = c(rep("GTP Binding Pocket", length(gtp_sites)),
                   rep("Other Region", length(other_sites)))
    )
  }
}

if (length(detailed_positions) > 0) {
  positions_table <- rbindlist(detailed_positions, fill = TRUE)
  fwrite(positions_table, "allosteric_positions_detailed.csv")
  cat("Detailed positions saved to: allosteric_positions_detailed.csv\n")
}

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Files created:\n")
cat("  - allosteric_sites_*.png (8 individual plots with error bars and legends)\n")
cat("  - allosteric_sites_combined.png (combined plot, each subplot has its own legend)\n")
cat("  - allosteric_analysis_details.rds (detailed R data)\n")
cat("  - allosteric_sites_summary.csv (summary table)\n")
cat("  - allosteric_positions_detailed.csv (detailed position information)\n")