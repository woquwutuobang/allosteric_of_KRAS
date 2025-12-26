library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)

# =====================================================
# assay 命名映射（核心修正点）
# =====================================================
assay_map <- data.table(
  assay_dis  = c("RAF1","RALGDS","PI3KCG","SOS1","K55","K27","K13","K19"),
  assay_file = c("RAF1","RAL","PI3","SOS","K55","K27","K13","K19")
)

# =====================================================
# Single panel: raw + median exponential fit
# =====================================================
make_energy_distance_panel <- function(df, x_col, y_col,
                                       x_label = NULL,   # 新增
                                       y_label = NULL,   # 新增
                                       interface_cutoff = 5,
                                       xlim = c(0, 35),
                                       ylim = c(0, 2.5),
                                       show_x = FALSE,
                                       show_y = FALSE) {
  
  if (!(x_col %in% colnames(df)) || !(y_col %in% colnames(df))) {
    return(ggplot() + theme_void())
  }
  
  df <- df[, .(x = get(x_col), y = abs(get(y_col)))]
  df <- df[complete.cases(df)]
  
  df_bi    <- df[x < interface_cutoff]
  df_nonbi <- df[x >= interface_cutoff]
  
  if (nrow(df_nonbi) < 5) return(ggplot() + theme_void())
  
  # =========================
  # median (non-interface)
  # =========================
  df_median <- df_nonbi %>%
    group_by(x) %>%
    summarise(y = median(y), .groups = "drop")
  
  # =========================
  # exponential fits
  # =========================
  fit_raw <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_nonbi,
        start = list(a = 1, b = -0.1),
        control = nls.control(warnOnly = TRUE)),
    error = function(e) NULL
  )
  
  fit_med <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_median,
        start = list(a = 1, b = -0.1),
        control = nls.control(warnOnly = TRUE)),
    error = function(e) NULL
  )
  
  # =========================
  # ⭐ 提取 a / b / p-value
  # =========================
  stat_label <- NULL
  
  if (!is.null(fit_raw)) {
    coef_raw <- coef(summary(fit_raw))
    
    a_val <- coef_raw["a", "Estimate"]
    b_val <- coef_raw["b", "Estimate"]
    
    # 线性化后算相关显著性（稳健）
    cor_test <- cor.test(df_nonbi$x, log(df_nonbi$y))
    
    p_val <- cor_test$p.value
    
    stat_label <- sprintf(
      "a = %.2f\nb = %.3f\np = %.1e",
      a_val, b_val, p_val
    )
  }
  
  # =========================
  # fit lines
  # =========================
  fit_df <- rbind(
    if (!is.null(fit_raw)) {
      data.frame(
        x = seq(min(df_nonbi$x), max(df_nonbi$x), length.out = 200),
        y = predict(fit_raw,
                    newdata = data.frame(x = seq(min(df_nonbi$x),
                                                 max(df_nonbi$x),
                                                 length.out = 200))),
        type = "Raw"
      )
    },
    if (!is.null(fit_med)) {
      data.frame(
        x = seq(min(df_median$x), max(df_median$x), length.out = 200),
        y = predict(fit_med,
                    newdata = data.frame(x = seq(min(df_median$x),
                                                 max(df_median$x),
                                                 length.out = 200))),
        type = "Median"
      )
    }
  )
  
  # =========================
  # base plot
  # =========================
  p <- ggplot() +
    geom_point(data = df_nonbi, aes(x, y),
               color = "#75C2F6", alpha = 0.15, size = 1.2) +
    geom_point(data = df_bi, aes(x, y),
               color = "#FFB6C1", alpha = 0.15, size = 1.2) +
    geom_point(data = df_median, aes(x, y),
               color = "#1B38A6", size = 1.4) +
    geom_line(data = fit_df,
              aes(x, y, linetype = type, color = type),
              linewidth = 0.6) +
    geom_vline(xintercept = interface_cutoff,
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_manual(values = c(Raw = "#1B38A6", Median = "#FF6A56")) +
    scale_linetype_manual(values = c(Raw = "solid", Median = "dashed")) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = 7) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_line(color = "black", linewidth = 0.3)
    )
  
  # =========================
  # ⭐ 右上角统计标注
  # =========================
  if (!is.null(stat_label)) {
    p <- p +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = stat_label,
        hjust = 1.05,
        vjust = 1.1,
        size = 3.8,      # 👈 panel 内已经算“明显偏大”
        lineheight = 1.1
      )
  }
  
  # =========================
  # axes control（你原来的）
  # =========================
  if (show_x) {
    p <- p +
      scale_x_continuous(breaks = seq(0, 35, 10)) +
      theme(
        axis.text.x  = element_text(size = 6),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.x = element_text(size = 7)
      ) +
      labs(x = x_label)
  }
  
  if (show_y) {
    p <- p +
      scale_y_continuous(breaks = seq(0, 2.5, 0.5)) +
      theme(
        axis.text.y  = element_text(size = 6),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_text(size = 7)
      ) +
      labs(y = y_label)
  }
  
  return(p)
}


# =====================================================
# 8x8 matrix plot（修正版）
# =====================================================
plot_energy_distance_8x8 <- function(input_template, anno_file,
                                     interface_cutoff = 5) {
    
  assays_dis  <- assay_map$assay_dis
  assays_file <- assay_map$assay_file
  
  anno <- fread(file = anno_file)
  anno[, Pos_real := Pos]
  
  plot_mat <- list()
  
  for (i in seq_along(assays_dis)) {
    
    assay_energy_dis  <- assays_dis[i]
    assay_energy_file <- assays_file[i]
    
    input <- gsub("ASSAY", assay_energy_file, input_template)
    dat   <- fread(file = input)
    
    dat[, Pos_real := Pos + 1]
    dat <- dat[, c(20:23)]
    colnames(dat)[1:3] <- paste0(colnames(dat)[1:3], "_", assay_energy_dis)
    
    df <- merge(anno, dat, by = "Pos_real")
    
    for (j in seq_along(assays_dis)) {
      assay_dist <- assays_dis[j]
      
      x_col <- paste0("scHAmin_ligand_", assay_dist)
      y_col <- paste0("mean_kcal/mol_", assay_energy_dis)
      
      p <- make_energy_distance_panel(
        df,
        x_col, y_col,
        x_label = paste0("Distance to binder ", assay_dist, " (Å)"),
        y_label = paste0("|ΔΔG| (", assay_energy_dis, ") kcal/mol"),
        interface_cutoff = interface_cutoff,
        show_x = (i == length(assays_dis)),  # 最底一行
        show_y = (j == 1)                    # 最左一列
      )
      
      plot_mat[[length(plot_mat) + 1]] <- p
    }
  }
  
  wrap_plots(plot_mat, ncol = 8)
}


# =====================================================
# Run
# =====================================================
p <- plot_energy_distance_8x8(
  input_template = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_ASSAY.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  interface_cutoff = 5
)

p

ggsave(
  filename = "./some results data/20251225/energy_distance_decay_8x8_raw_median2.png",
  plot = p,
  width = 13,
  height = 13,
  units = "in"
)





















