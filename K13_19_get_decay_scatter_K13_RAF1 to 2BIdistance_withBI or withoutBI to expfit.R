########## with BI to expfit

# Energy-Distance Decay Analysis with Exponential Fitting (Including Binding Interface) using BOTH raw and median values
# This script analyzes the relationship between binding energy changes and distance from binding partners using both raw and median values

library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)        
library(gridExtra)  

# ===============================
# Function: Distance Effect Fitting (Exponential Model) - Including Binding Interface using BOTH raw and median
# ===============================
plot_energy_distance_decay_expfit_both_with_interface <- function(input, assay_sele, anno_file,
                                                                  x_cols = c("scHAmin_ligand_K13", "scHAmin_ligand_RAF1"),
                                                                  titles = c("Distance to K13(Å)", "Distance to RAF1(Å)"),
                                                                  x_range = c(0, 35), y_range = c(0, 3)) {
  
  # 1. Read data
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  plot_list <- list()
  
  # 2. Iterate through all x_cols
  for (i in seq_along(x_cols)) {
    xvector <- anno_final[[x_cols[i]]]
    yvector <- abs(anno_final[[y_col]])  # Take absolute value of energy
    
    df <- data.frame(x = xvector, y = yvector)
    df <- df[complete.cases(df), ]
    
    # --- Separate binding interface (BI) and non-interface points ---
    df_bi <- df[df$x < 5, ]      # Binding interface: distance < 5Å
    df_non_bi <- df[df$x >= 5, ] # Non-interface: distance >= 5Å
    
    # --- Calculate median values for binding interface points ---
    df_median_bi <- df_bi %>%
      group_by(x) %>%
      summarise(
        y_median = median(y, na.rm = TRUE),
        n = n()
      ) %>%
      filter(n >= 1)
    
    # --- Calculate median values for non-interface points ---
    df_median_non_bi <- df_non_bi %>%
      group_by(x) %>%
      summarise(
        y_median = median(y, na.rm = TRUE),
        n = n()
      ) %>%
      filter(n >= 1)
    
    # --- Combine all data for fitting ---
    df_all <- df  # All data points
    
    # --- Calculate median values for all distances (for completeness) ---
    df_median_all <- df_all %>%
      group_by(x) %>%
      summarise(
        y_median = median(y, na.rm = TRUE),
        n = n()
      ) %>%
      filter(n >= 1)
    
    # --- Fit exponential model using RAW data (all data) ---
    residual_sum_of_squares_raw <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess_raw <- c(a = 1, b = -0.1)
    opt_params_raw <- tryCatch(
      optim(initial_guess_raw, residual_sum_of_squares_raw, 
            x = df_all$x, y = df_all$y)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    fit_model_raw <- tryCatch(
      nls(y ~ a * exp(b * x), data = df_all, 
          start = list(a = opt_params_raw[1], b = opt_params_raw[2])),
      error = function(e) NULL
    )
    
    # --- Fit exponential model using MEDIAN data (all data) ---
    residual_sum_of_squares_median <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess_median <- c(a = 1, b = -0.1)
    opt_params_median <- tryCatch(
      optim(initial_guess_median, residual_sum_of_squares_median, 
            x = df_median_all$x, y = df_median_all$y_median)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    fit_model_median <- tryCatch(
      nls(y_median ~ a * exp(b * x), data = df_median_all, 
          start = list(a = opt_params_median[1], b = opt_params_median[2])),
      error = function(e) NULL
    )
    
    # --- Fit results for both models ---
    fit_df_raw <- data.frame()
    fit_df_median <- data.frame()
    annotation_text_raw <- NULL
    annotation_text_median <- NULL
    
    # Raw data fitting results
    if (!is.null(fit_model_raw)) {
      x_seq <- seq(min(df_all$x, na.rm = TRUE), max(df_all$x, na.rm = TRUE), length.out = 200)
      y_fit_raw <- predict(fit_model_raw, newdata = data.frame(x = x_seq))
      fit_df_raw <- data.frame(x = x_seq, y = y_fit_raw, type = "Raw")
      
      fit_summary_raw <- summary(fit_model_raw)
      coefs_raw <- fit_summary_raw$coefficients
      a_val_raw <- round(coefs_raw["a", "Estimate"], 3)
      b_val_raw <- round(coefs_raw["b", "Estimate"], 3)
      p_val_b_raw <- coefs_raw["b", "Pr(>|t|)"]
      
      p_text_raw <- if (is.na(p_val_b_raw)) {
        "p = NA"
      } else if (p_val_b_raw < 0.001) {
        "p < 0.001"
      } else if (p_val_b_raw < 0.05) {
        "p < 0.05"
      } else {
        paste0("p = ", round(p_val_b_raw, 3))
      }
      
      annotation_text_raw <- paste0("Raw: a = ", a_val_raw, ", b = ", b_val_raw, "\n", p_text_raw)
    }
    
    # Median data fitting results
    if (!is.null(fit_model_median)) {
      x_seq <- seq(min(df_median_all$x, na.rm = TRUE), max(df_median_all$x, na.rm = TRUE), length.out = 200)
      y_fit_median <- predict(fit_model_median, newdata = data.frame(x = x_seq))
      fit_df_median <- data.frame(x = x_seq, y = y_fit_median, type = "Median")
      
      fit_summary_median <- summary(fit_model_median)
      coefs_median <- fit_summary_median$coefficients
      a_val_median <- round(coefs_median["a", "Estimate"], 3)
      b_val_median <- round(coefs_median["b", "Estimate"], 3)
      p_val_b_median <- coefs_median["b", "Pr(>|t|)"]
      
      p_text_median <- if (is.na(p_val_b_median)) {
        "p = NA"
      } else if (p_val_b_median < 0.001) {
        "p < 0.001"
      } else if (p_val_b_median < 0.05) {
        "p < 0.05"
      } else {
        paste0("p = ", round(p_val_b_median, 3))
      }
      
      annotation_text_median <- paste0("Median: a = ", a_val_median, ", b = ", b_val_median, "\n", p_text_median)
    }
    
    # --- Combine fit data ---
    fit_df_combined <- rbind(fit_df_raw, fit_df_median)
    
    # --- Shared x breaks ---
    x_breaks <- pretty(x_range, n = 5)
    
    # ========== Main plot ==========
    p_main <- ggplot() +
      # Plot non-interface individual points (light blue)
      geom_point(data = df_non_bi, aes(x = x, y = y),
                 alpha = 0.1, size = 3, color = "#75C2F6") +
      # Plot binding interface individual points (light red)
      geom_point(data = df_bi, aes(x = x, y = y),
                 alpha = 0.1, size = 3, color = "#FFB6C1") +
      # Plot non-interface median points (dark blue)
      geom_point(data = df_median_non_bi, aes(x = x, y = y_median),
                 color = "#1B38A6", size = 2, shape = 19) +
      # Plot binding interface median points (dark red)
      geom_point(data = df_median_bi, aes(x = x, y = y_median),
                 color = "#8B0000", size = 2, shape = 19) +
      # Plot fitted curves
      geom_line(data = fit_df_combined, aes(x = x, y = y, color = type, linetype = type),
                linewidth = 0.8) +
      # Add vertical line to indicate interface boundary (5Å)
      geom_vline(xintercept = 5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      # Add text annotation for interface boundary
      annotate("text", x = 5, y = max(y_range) * 0.95, 
               label = "Interface boundary", hjust = -0.1, size = 3, color = "gray50") +
      scale_x_continuous(
        limits = x_range,
        breaks = x_breaks,
        expand = c(0.02, 0)
      ) +
      scale_y_continuous(
        limits = y_range,
        expand = c(0, 0)
      ) +
      scale_color_manual(
        name = "Fitting",
        values = c("Raw" = "#1B38A6", "Median" = "#FF6A56"),
        labels = c("Raw data fitting", "Median data fitting")
      ) +
      scale_linetype_manual(
        name = "Fitting",
        values = c("Raw" = "solid", "Median" = "dashed"),
        labels = c("Raw data fitting", "Median data fitting")
      ) +
      theme_classic(base_size = 14) +
      theme(
        text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        
        # X-axis ticks - Force 90 degree rotation
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          margin = margin(t = 2, b = 2),
          lineheight = 0.9
        ),
        
        # Y-axis ticks
        axis.text.y = element_text(
          margin = margin(l = 2, r = 2),
          hjust = 1,
          lineheight = 0.9
        ),
        axis.title.y = if(i > 1) element_blank() else element_text(size = 14),
        
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(linewidth = 0.6, colour = "black"),
        plot.margin = margin(5, 5, 15, 5),
        
        # Legend position
        legend.position = c(0.7, 0.85),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(0.4, "cm")
      ) +
      labs(
        x = titles[i],
        y = if(i == 1) paste0("Binding ΔΔG (", assay_sele,")(kcal/mol)") else ""
      ) +
      guides(
        color = guide_legend(title = "Fitting Type"),
        linetype = guide_legend(title = "Fitting Type")
      )
    
    # --- Annotation text with matching colors ---
    if (!is.null(annotation_text_raw)) {
      p_main <- p_main +
        annotation_custom(
          grob = textGrob(
            label = annotation_text_raw,
            x = unit(0.97, "npc"),
            y = unit(0.95, "npc"),
            just = c("right", "top"),
            gp = gpar(fontsize = 11, col = "#FF6A56")  # Same color as raw fitting line
          )
        )
    }
    
    if (!is.null(annotation_text_median)) {
      p_main <- p_main +
        annotation_custom(
          grob = textGrob(
            label = annotation_text_median,
            x = unit(0.97, "npc"),
            y = unit(0.85, "npc"),
            just = c("right", "top"),
            gp = gpar(fontsize = 11, col = "#1B38A6")  # Same color as median fitting line
          )
        )
    }
    
    # --- Add legend for point types ---
    p_main <- p_main +
      # Add point color explanation
      annotate("text", x = Inf, y = Inf, 
               label = "Points: Light=raw, Dark=median\nBlue=Non-interface, Red=Interface", 
               hjust = 1.05, vjust = 1.5, size = 2.8, color = "gray30")
    
    # --- Add box around axis titles ---
    p_main <- p_main +
      theme(
        axis.title.x = element_text(
          margin = margin(t = 10),
          color = "black"
        ),
        axis.title.y = if(i == 1) element_text(
          margin = margin(r = 10),
          color = "black"
        ) else element_blank()
      )
    
    plot_list[[i]] <- p_main
  }
  
  # 3. Combine all plots with adjusted layout
  final_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      theme = theme(
        plot.margin = margin(10, 10, 10, 10)
      )
    )
  
  return(final_plot)
}

# ===============================
# Analysis for RAF1 (Including Binding Interface) using BOTH raw and median
# ===============================
plot_result_raf1_both_with_interface <- plot_energy_distance_decay_expfit_both_with_interface(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
  assay_sele = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_result_raf1_both_with_interface)

# Save plot using ggsave (outside the function)
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251119/RAF1_expfit_both_withBI.pdf",
  plot = plot_result_raf1_both_with_interface,
  device = cairo_pdf,
  height = 4,
  width = 8
)

# ===============================
# Analysis for K13 (Including Binding Interface) using BOTH raw and median
# ===============================
plot_result_k13_both_with_interface <- plot_energy_distance_decay_expfit_both_with_interface(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_result_k13_both_with_interface)

# Save plot using ggsave (outside the function)
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251119/K13_expfit_both_withBI.pdf",
  plot = plot_result_k13_both_with_interface,
  device = cairo_pdf,
  height = 4,
  width = 8
)













########## without BI to expfit

# Energy-Distance Decay Analysis with Exponential Fitting (Including Binding Interface) using BOTH raw and median values
# This script analyzes the relationship between binding energy changes and distance from binding partners using both raw and median values

library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)        
library(gridExtra)  

# ===============================
# Function: Distance Effect Fitting (Exponential Model) - Fitting without binding interface but plotting all points
# ===============================
plot_energy_distance_decay_expfit_both_without_interface_plot_with_interface <- function(input, assay_sele, anno_file,
                                                                                         x_cols = c("scHAmin_ligand_K13", "scHAmin_ligand_RAF1"),
                                                                                         titles = c("Distance to K13(Å)", "Distance to RAF1(Å)"),
                                                                                         x_range = c(0, 35), y_range = c(0, 3)) {
  
  # 1. Read data
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  plot_list <- list()
  
  # 2. Iterate through all x_cols
  for (i in seq_along(x_cols)) {
    xvector <- anno_final[[x_cols[i]]]
    yvector <- abs(anno_final[[y_col]])  # Take absolute value of energy
    
    df <- data.frame(x = xvector, y = yvector)
    df <- df[complete.cases(df), ]
    
    # --- Separate binding interface (BI) and non-interface points ---
    df_bi <- df[df$x < 5, ]      # Binding interface: distance < 5Å
    df_non_bi <- df[df$x >= 5, ] # Non-interface: distance >= 5Å
    
    # --- Calculate median values for each distance (for non-interface only for fitting) ---
    df_median_non_bi <- df_non_bi %>%
      group_by(x) %>%
      summarise(
        y_median = median(y, na.rm = TRUE),
        n = n()
      ) %>%
      filter(n >= 1)  # Ensure we have at least one data point per distance
    
    # Also calculate median for BI points (for plotting only)
    df_median_bi <- df_bi %>%
      group_by(x) %>%
      summarise(
        y_median = median(y, na.rm = TRUE),
        n = n()
      ) %>%
      filter(n >= 1)
    
    # --- Fit exponential model using RAW data (non-interface only) ---
    residual_sum_of_squares_raw <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess_raw <- c(a = 1, b = -0.1)
    opt_params_raw <- tryCatch(
      optim(initial_guess_raw, residual_sum_of_squares_raw, 
            x = df_non_bi$x, y = df_non_bi$y)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    fit_model_raw <- tryCatch(
      nls(y ~ a * exp(b * x), data = df_non_bi, 
          start = list(a = opt_params_raw[1], b = opt_params_raw[2])),
      error = function(e) NULL
    )
    
    # --- Fit exponential model using MEDIAN data (non-interface only) ---
    residual_sum_of_squares_median <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess_median <- c(a = 1, b = -0.1)
    opt_params_median <- tryCatch(
      optim(initial_guess_median, residual_sum_of_squares_median, 
            x = df_median_non_bi$x, y = df_median_non_bi$y_median)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    fit_model_median <- tryCatch(
      nls(y_median ~ a * exp(b * x), data = df_median_non_bi, 
          start = list(a = opt_params_median[1], b = opt_params_median[2])),
      error = function(e) NULL
    )
    
    # --- Fit results for both models ---
    fit_df_raw <- data.frame()
    fit_df_median <- data.frame()
    annotation_text_raw <- NULL
    annotation_text_median <- NULL
    
    # Raw data fitting results (non-interface)
    if (!is.null(fit_model_raw)) {
      x_seq <- seq(min(df_non_bi$x, na.rm = TRUE), max(df_non_bi$x, na.rm = TRUE), length.out = 200)
      y_fit_raw <- predict(fit_model_raw, newdata = data.frame(x = x_seq))
      fit_df_raw <- data.frame(x = x_seq, y = y_fit_raw, type = "Raw")
      
      fit_summary_raw <- summary(fit_model_raw)
      coefs_raw <- fit_summary_raw$coefficients
      a_val_raw <- round(coefs_raw["a", "Estimate"], 3)
      b_val_raw <- round(coefs_raw["b", "Estimate"], 3)
      p_val_b_raw <- coefs_raw["b", "Pr(>|t|)"]
      
      p_text_raw <- if (is.na(p_val_b_raw)) {
        "p = NA"
      } else if (p_val_b_raw < 0.001) {
        "p < 0.001"
      } else if (p_val_b_raw < 0.05) {
        "p < 0.05"
      } else {
        paste0("p = ", round(p_val_b_raw, 3))
      }
      
      annotation_text_raw <- paste0("Raw: a = ", a_val_raw, ", b = ", b_val_raw, "\n", p_text_raw)
    }
    
    # Median data fitting results (non-interface)
    if (!is.null(fit_model_median)) {
      x_seq <- seq(min(df_median_non_bi$x, na.rm = TRUE), max(df_median_non_bi$x, na.rm = TRUE), length.out = 200)
      y_fit_median <- predict(fit_model_median, newdata = data.frame(x = x_seq))
      fit_df_median <- data.frame(x = x_seq, y = y_fit_median, type = "Median")
      
      fit_summary_median <- summary(fit_model_median)
      coefs_median <- fit_summary_median$coefficients
      a_val_median <- round(coefs_median["a", "Estimate"], 3)
      b_val_median <- round(coefs_median["b", "Estimate"], 3)
      p_val_b_median <- coefs_median["b", "Pr(>|t|)"]
      
      p_text_median <- if (is.na(p_val_b_median)) {
        "p = NA"
      } else if (p_val_b_median < 0.001) {
        "p < 0.001"
      } else if (p_val_b_median < 0.05) {
        "p < 0.05"
      } else {
        paste0("p = ", round(p_val_b_median, 3))
      }
      
      annotation_text_median <- paste0("Median: a = ", a_val_median, ", b = ", b_val_median, "\n", p_text_median)
    }
    
    # --- Combine fit data ---
    fit_df_combined <- rbind(fit_df_raw, fit_df_median)
    
    # --- Shared x breaks ---
    x_breaks <- pretty(x_range, n = 5)
    
    # ========== Main plot ==========
    p_main <- ggplot() +
      # Plot non-interface individual points (blue, transparent)
      geom_point(data = df_non_bi, aes(x = x, y = y),
                 alpha = 0.1, size = 3, color = "#75C2F6") +
      # Plot binding interface individual points (light red, transparent)
      geom_point(data = df_bi, aes(x = x, y = y),
                 alpha = 0.1, size = 3, color = "#FFB6C1") +
      # Plot non-interface median points (dark blue)
      geom_point(data = df_median_non_bi, aes(x = x, y = y_median),
                 color = "#1B38A6", size = 2, shape = 19) +
      # Plot binding interface median points (dark red)
      geom_point(data = df_median_bi, aes(x = x, y = y_median),
                 color = "#8B0000", size = 2, shape = 19) +
      # Plot fitted curves (only fitted on non-interface data)
      geom_line(data = fit_df_combined, aes(x = x, y = y, color = type, linetype = type),
                linewidth = 0.8) +
      # Add vertical line to indicate interface boundary (5Å)
      geom_vline(xintercept = 5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      # Add text annotation for interface boundary
      annotate("text", x = 5, y = max(y_range) * 0.95, 
               label = "Interface boundary", hjust = -0.1, size = 3, color = "gray50") +
      scale_x_continuous(
        limits = x_range,
        breaks = x_breaks,
        expand = c(0.02, 0)
      ) +
      scale_y_continuous(
        limits = y_range,
        expand = c(0, 0)
      ) +
      scale_color_manual(
        name = "Fitting (Non-interface)",
        values = c("Raw" = "#1B38A6", "Median" = "#FF6A56"),
        labels = c("Raw data fitting", "Median data fitting")
      ) +
      scale_linetype_manual(
        name = "Fitting (Non-interface)",
        values = c("Raw" = "solid", "Median" = "dashed"),
        labels = c("Raw data fitting", "Median data fitting")
      ) +
      theme_classic(base_size = 14) +
      theme(
        text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour = "black"),
        
        # X-axis ticks - Force 90 degree rotation
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          margin = margin(t = 2, b = 2),
          lineheight = 0.9
        ),
        
        # Y-axis ticks
        axis.text.y = element_text(
          margin = margin(l = 2, r = 2),
          hjust = 1,
          lineheight = 0.9
        ),
        axis.title.y = if(i > 1) element_blank() else element_text(size = 14),
        
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(linewidth = 0.6, colour = "black"),
        plot.margin = margin(5, 5, 15, 5),
        
        # Legend position
        legend.position = c(0.7, 0.85),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(0.4, "cm")
      ) +
      labs(
        x = titles[i],
        y = if(i == 1) paste0("Binding ΔΔG (", assay_sele,")(kcal/mol)") else ""
      ) +
      guides(
        color = guide_legend(title = "Fitting Type"),
        linetype = guide_legend(title = "Fitting Type")
      )
    
    # --- Annotation text with matching colors ---
    if (!is.null(annotation_text_raw)) {
      p_main <- p_main +
        annotation_custom(
          grob = textGrob(
            label = annotation_text_raw,
            x = unit(0.97, "npc"),
            y = unit(0.95, "npc"),
            just = c("right", "top"),
            gp = gpar(fontsize = 11, col = "#FF6A56")  # Same color as raw fitting line
          )
        )
    }
    
    if (!is.null(annotation_text_median)) {
      p_main <- p_main +
        annotation_custom(
          grob = textGrob(
            label = annotation_text_median,
            x = unit(0.97, "npc"),
            y = unit(0.85, "npc"),
            just = c("right", "top"),
            gp = gpar(fontsize = 11, col = "#1B38A6")  # Same color as median fitting line
          )
        )
    }
    
    # --- Add color legend for point types ---
    # Create custom legend for point types
    point_legend <- data.frame(
      label = c("Non-interface (fitted)", "Interface (excluded)"),
      color = c("#75C2F6", "#FFB6C1"),
      shape = c(19, 19)
    )
    
    # Add legend using annotate with dummy points
    p_main <- p_main +
      annotate("point", x = -Inf, y = -Inf, color = "#75C2F6", size = 3, alpha = 0.5) +
      annotate("point", x = -Inf, y = -Inf, color = "#FFB6C1", size = 3, alpha = 0.5) +
      annotate("text", x = Inf, y = Inf, 
               label = "Points: Light=raw, Dark=median\nBlue=Non-interface, Red=Interface", 
               hjust = 1.05, vjust = 1.5, size = 2.8, color = "gray30")
    
    # --- Add box around axis titles ---
    p_main <- p_main +
      theme(
        axis.title.x = element_text(
          margin = margin(t = 10),
          color = "black"
        ),
        axis.title.y = if(i == 1) element_text(
          margin = margin(r = 10),
          color = "black"
        ) else element_blank()
      )
    
    plot_list[[i]] <- p_main
  }
  
  # 3. Combine all plots with adjusted layout
  final_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      theme = theme(
        plot.margin = margin(10, 10, 10, 10)
      )
    )
  
  return(final_plot)
}

# ===============================
# Analysis for RAF1 (Fitting without binding interface but plotting with interface)
# ===============================
plot_result_raf1_without_interface_plot_with_interface <- plot_energy_distance_decay_expfit_both_without_interface_plot_with_interface(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
  assay_sele = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
)

print(plot_result_raf1_without_interface_plot_with_interface)

# Save plot using ggsave (outside the function)
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251119/RAF1_expfit_withoutBI_plot_withBI.pdf",
  plot = plot_result_raf1_without_interface_plot_with_interface,
  device = cairo_pdf,
  height = 4,
  width = 8
)

# ===============================
# Analysis for K13 (Fitting without binding interface but plotting with interface)
# ===============================
plot_result_k13_without_interface_plot_with_interface <- plot_energy_distance_decay_expfit_both_without_interface_plot_with_interface(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_result_k13_without_interface_plot_with_interface)

# Save plot using ggsave (outside the function)
ggsave(
  filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251119/K13_expfit_withoutBI_plot_withBI.pdf",
  plot = plot_result_k13_without_interface_plot_with_interface,
  device = cairo_pdf,
  height = 4,
  width = 8
)