# ============================================================
# Allosteric site OR analysis with parallel evidence logic
# ============================================================

library(data.table)
library(ggplot2)

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
  fisher_sites[, statistical_support := (q_value < 0.05 & odds_ratio > 1)]
  
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
  
  # 修改1: 从"and"逻辑改为"or"逻辑 - 符合任一条件即可
  final_sites[, evidence_type := fifelse(
    statistical_support | quantitative_support,  # 改为 OR 逻辑
    "Significant (either statistical or quantitative)",
    "Non-significant"
  )]
  
  return(final_sites[order(Pos_real)])
}

# ============================================================
# ============================================================

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

# ============================================================
# ============================================================

allosteric_sites <- list()
allosteric_details <- list()

for (b in names(input_files)) {
  
  res <- identify_allosteric_mutations(
    input = input_files[[b]],
    anno = anno_file,
    wt_aa = wt_aa,
    assay_sele = b
  )
  
  # 使用新的并行证据方法识别变构位点
  allo_sites <- annotate_allosteric_sites_parallel(
    data_plot_mutation = res$mutation_table,
    min_n = 9
  )
  
  # 保存详细结果
  allosteric_details[[b]] <- allo_sites
  
  # 修改2: 提取显著位点（统计或定量证据）
  significant_sites <- allo_sites[
    evidence_type == "Significant (either statistical or quantitative)",
    Pos_real
  ]
  
  # 修改3: 排除结合界面位点
  # 获取结合界面位点信息
  binding_sites <- res$mutation_table[site_type == "Binding interface site", unique(Pos_real)]
  gtp_binding_sites <- res$mutation_table[site_type == "GTP binding interface site", unique(Pos_real)]
  all_binding_sites <- unique(c(binding_sites, gtp_binding_sites))
  
  # 只保留非结合界面的显著位点
  allosteric_only_sites <- significant_sites[!significant_sites %in% all_binding_sites]
  
  allosteric_sites[[b]] <- allosteric_only_sites
  
  # 可选：打印每个配体的变构位点信息
  cat(sprintf("\n%s - Allosteric sites:\n", b))
  cat(sprintf("  Total significant sites: %d\n", length(significant_sites)))
  cat(sprintf("  After removing binding sites: %d\n", length(allosteric_only_sites)))
  cat(sprintf("  Binding interface sites removed: %s\n", 
              paste(sort(all_binding_sites), collapse = ", ")))
  cat(sprintf("  Final allosteric sites: %s\n", 
              paste(sort(allosteric_only_sites), collapse = ", ")))
}

universe_sites <- sort(unique(unlist(allosteric_sites)))
cat(sprintf("\nTotal unique allosteric sites across all binders: %d\n", 
            length(universe_sites)))

# ===========================================================
# 构建列联表并计算OR
# ===========================================================

build_contingency <- function(A, B, universe) {
  a <- universe %in% A
  b <- universe %in% B
  matrix(
    c(sum(a & b),
      sum(a & !b),
      sum(!a & b),
      sum(!a & !b)),
    nrow = 2,
    byrow = TRUE
  )
}

compute_or <- function(mat, cc = 0.5) {
  m <- mat + cc
  (m[1,1] * m[2,2]) / (m[1,2] * m[2,1])
}

binders <- names(allosteric_sites)
res_or <- list()
k <- 1

for (i in 1:(length(binders) - 1)) {
  for (j in (i + 1):length(binders)) {
    
    mat <- build_contingency(
      allosteric_sites[[binders[i]]],
      allosteric_sites[[binders[j]]],
      universe_sites
    )
    
    res_or[[k]] <- data.frame(
      assay1 = binders[i],
      assay2 = binders[j],
      OR = compute_or(mat),
      p_value = fisher.test(mat)$p.value
    )
    k <- k + 1
  }
}

or_results <- rbindlist(res_or)

or_results$significance <- cut(
  or_results$p_value,
  c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

# ============================================================
# 分类、排序、绘图
# ============================================================

binders_all <- c("RAF1","RALGDS","PI3KCG","SOS1","K55","K27","K13","K19")
non_darpin <- c("RAF1","RALGDS","PI3KCG","SOS1")
darpin     <- c("K55","K27","K13","K19")
active     <- c("RAF1","RALGDS","PI3KCG","K55")
inactive   <- c("SOS1","K27","K13","K19")
BI1 <- c("RAF1", "RALGDS", "PI3KCG", "SOS1","K55", "K27")
BI2 <- c("K13", "K19")

# 添加额外的分类信息
or_results$DARPinClass <- apply(or_results, 1, function(x) {
  if (x["assay1"] %in% non_darpin & x["assay2"] %in% non_darpin) {
    return("nonDARPin_vs_nonDARPin")
  } else if (x["assay1"] %in% darpin & x["assay2"] %in% darpin) {
    return("DARPin_vs_DARPin")
  } else {
    return("nonDARPin_vs_DARPin")
  }
})

or_results$ActivityClass <- apply(or_results, 1, function(x) {
  if (x["assay1"] %in% active & x["assay2"] %in% active) {
    return("Active_vs_Active")
  } else if (x["assay1"] %in% inactive & x["assay2"] %in% inactive) {
    return("Inactive_vs_Inactive")
  } else {
    return("Active_vs_Inactive")
  }
})

or_results$Comparison <- paste0(or_results$assay1, " vs ", or_results$assay2)

or_results$comparison_type <- apply(or_results, 1, function(x) {
  assay1_in_BI1 <- x["assay1"] %in% BI1
  assay1_in_BI2 <- x["assay1"] %in% BI2
  assay2_in_BI1 <- x["assay2"] %in% BI1
  assay2_in_BI2 <- x["assay2"] %in% BI2
  
  if (assay1_in_BI1 && assay2_in_BI1) {
    return("BI1_vs_BI1")
  } else if (assay1_in_BI2 && assay2_in_BI2) {
    return("BI2_vs_BI2")
  } else {
    return("BI1_vs_BI2")
  }
})

# 排序
or_results <- or_results[order(
  # 第一层：comparison_type
  factor(or_results$comparison_type, levels = c(
    "BI1_vs_BI1",
    "BI1_vs_BI2", 
    "BI2_vs_BI2"
  )),
  # 第二层：DARPin分类
  factor(or_results$DARPinClass, levels = c(
    "nonDARPin_vs_nonDARPin",
    "nonDARPin_vs_DARPin",
    "DARPin_vs_DARPin"
  )),
  # 第三层：活性分类
  factor(or_results$ActivityClass, levels = c(
    "Active_vs_Active",
    "Active_vs_Inactive", 
    "Inactive_vs_Inactive"
  )),
  # 最后按比较名称排序以确保一致性
  factor(or_results$Comparison, levels = unique(or_results$Comparison))
), ]

# 设置因子顺序
or_results$Comparison <- factor(or_results$Comparison, levels = or_results$Comparison)
or_results$comparison_type <- factor(
  or_results$comparison_type, 
  levels = c("BI1_vs_BI1", "BI1_vs_BI2", "BI2_vs_BI2")
)

# 计算log10(OR)
or_results$OR_log10 <- log10(or_results$OR)

# ============================================================
# 可视化
# ============================================================

p_or <- ggplot(or_results,
               aes(x = Comparison,
                   y = OR_log10,
                   fill = comparison_type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("log10(OR)=%.2f", OR_log10)),
            vjust = -0.5, size = 2) +
  geom_text(aes(label = significance),
            vjust = -1.4, size = 3, color = "#F1DD10") +
  scale_fill_manual(values = c(
    "BI1_vs_BI1" = "#F4270C",
    "BI1_vs_BI2" = "#F4AD0C",
    "BI2_vs_BI2" = "#1B38A6"
  )) +
  labs(
    title = "Allosteric Site Overlap Analysis (OR Logic: Statistical OR Quantitative)",
    subtitle = "Significant sites: Statistical OR Quantitative evidence, excluding binding interfaces",
    x = "Comparison",
    y = "log10(Odds Ratio)",
    fill = "Binding Interface Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_or)
