# Energy-Distance Decay Analysis with Exponential Fitting (Excluding Binding Interface)
# Using BOTH raw data and median values per position with BI1/BI2 color scheme
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)

# ===============================
# Function: Calculate Exponential Fit Parameters for BOTH raw and median (with interface exclusion)
# ===============================
calculate_exp_fit_parameters_both_exclude_interface <- function(input, assay_sele, anno_file, exclude_interface = TRUE, interface_cutoff = 5) {
  
  # 1. Read data
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  
  # Define secondary structures
  rects_sheet <- data.frame(xstart = c(3, 38, 51, 77, 109, 139),
                            xend = c(9, 44, 57, 84, 115, 143),
                            col = c("β1", "β2", "β3", "β4", "β5", "β6"))
  
  rects_helix <- data.frame(xstart = c(15, 67, 87, 127, 148),
                            xend = c(24, 73, 104, 136, 166),
                            col = c("α1", "α2", "α3", "α4", "α5"))
  
  # Prepare results data frame
  results <- data.frame()
  
  # Function to calculate exponential fit for both raw and median data
  calculate_exp_fit_both <- function(xvector, yvector, pos_vector, region_name, data_type) {
    # Remove NA values
    valid_idx <- complete.cases(data.frame(x = xvector, y = yvector, pos = pos_vector))
    xvector <- xvector[valid_idx]
    yvector <- abs(yvector[valid_idx])  # Take absolute value of energy
    pos_vector <- pos_vector[valid_idx]
    
    if (length(xvector) < 5) {  # Need enough points for fitting
      return(data.frame(
        region = region_name,
        data_type = data_type,
        a_value = NA,
        b_value = NA,
        p_value = NA,
        n_points = length(xvector),
        status = "insufficient_data"
      ))
    }
    
    # For median data: calculate median y for each position
    if (data_type == "median") {
      pos_median <- data.frame(
        pos = pos_vector,
        x = xvector,
        y = yvector
      ) %>%
        group_by(pos) %>%
        summarise(
          x_median = median(x, na.rm = TRUE),
          y_median = median(y, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (nrow(pos_median) < 5) {
        return(data.frame(
          region = region_name,
          data_type = data_type,
          a_value = NA,
          b_value = NA,
          p_value = NA,
          n_points = nrow(pos_median),
          status = "insufficient_data"
        ))
      }
      
      xvector <- pos_median$x_median
      yvector <- pos_median$y_median
    }
    
    # Initial parameter optimization
    residual_sum_of_squares <- function(params, x, y) {
      a <- params[1]; b <- params[2]
      predicted <- a * exp(b * x)
      sum((y - predicted)^2, na.rm = TRUE)
    }
    
    initial_guess <- c(a = 1, b = -0.1)
    opt_params <- tryCatch(
      optim(initial_guess, residual_sum_of_squares, x = xvector, y = yvector)$par,
      error = function(e) c(a = 1, b = -0.1)
    )
    
    # Fit exponential model
    fit_model <- tryCatch(
      nls(y ~ a * exp(b * x), 
          data = data.frame(x = xvector, y = yvector), 
          start = list(a = opt_params[1], b = opt_params[2]),
          control = nls.control(warnOnly = TRUE, maxiter = 100)),
      error = function(e) NULL
    )
    
    if (!is.null(fit_model) && !is.na(coef(fit_model)["b"])) {
      fit_summary <- summary(fit_model)
      coefs <- fit_summary$coefficients
      
      return(data.frame(
        region = region_name,
        data_type = data_type,
        a_value = coefs["a", "Estimate"],
        b_value = coefs["b", "Estimate"],
        p_value = coefs["b", "Pr(>|t|)"],
        n_points = length(xvector),
        status = "success"
      ))
    } else {
      return(data.frame(
        region = region_name,
        data_type = data_type,
        a_value = NA,
        b_value = NA,
        p_value = NA,
        n_points = length(xvector),
        status = "fit_failed"
      ))
    }
  }
  
  # Get all distance columns for different assays
  distance_assays <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
  distance_cols <- paste0("scHAmin_ligand_", distance_assays)
  
  # Calculate for each distance assay
  for (dist_assay in distance_assays) {
    x_col <- paste0("scHAmin_ligand_", dist_assay)
    
    if (!x_col %in% colnames(anno_final)) next
    
    # Prepare main data
    xvector <- anno_final[[x_col]]
    yvector <- anno_final[[y_col]]
    pos_real <- anno_final$Pos_real
    
    main_df <- data.frame(
      x = xvector,
      y = yvector,
      pos = pos_real
    )
    main_df <- main_df[complete.cases(main_df), ]
    
    # Exclude binding interface residues if requested
    if (exclude_interface) {
      main_df <- main_df[main_df$x > interface_cutoff, ]
    }
    
    if (nrow(main_df) < 5) next
    
    # Calculate for overall data - BOTH raw and median
    overall_result_raw <- calculate_exp_fit_both(main_df$x, main_df$y, main_df$pos, "overall", "raw")
    overall_result_median <- calculate_exp_fit_both(main_df$x, main_df$y, main_df$pos, "overall", "median")
    
    overall_results <- rbind(overall_result_raw, overall_result_median)
    overall_results$assay <- assay_sele
    overall_results$distance_assay <- dist_assay
    overall_results$exclude_interface <- exclude_interface
    overall_results$interface_cutoff <- ifelse(exclude_interface, interface_cutoff, NA)
    results <- rbind(results, overall_results)
    
    # Calculate for sheet regions - BOTH raw and median
    sheet_positions <- unique(unlist(lapply(1:nrow(rects_sheet), function(i) {
      rects_sheet$xstart[i]:rects_sheet$xend[i]
    })))
    
    sheet_df <- main_df[main_df$pos %in% sheet_positions, ]
    if (nrow(sheet_df) >= 5) {
      sheet_result_raw <- calculate_exp_fit_both(sheet_df$x, sheet_df$y, sheet_df$pos, "sheet", "raw")
      sheet_result_median <- calculate_exp_fit_both(sheet_df$x, sheet_df$y, sheet_df$pos, "sheet", "median")
      
      sheet_results <- rbind(sheet_result_raw, sheet_result_median)
      sheet_results$assay <- assay_sele
      sheet_results$distance_assay <- dist_assay
      sheet_results$exclude_interface <- exclude_interface
      sheet_results$interface_cutoff <- ifelse(exclude_interface, interface_cutoff, NA)
      results <- rbind(results, sheet_results)
    }
    
    # Calculate for helix regions - BOTH raw and median
    helix_positions <- unique(unlist(lapply(1:nrow(rects_helix), function(i) {
      rects_helix$xstart[i]:rects_helix$xend[i]
    })))
    
    helix_df <- main_df[main_df$pos %in% helix_positions, ]
    if (nrow(helix_df) >= 5) {
      helix_result_raw <- calculate_exp_fit_both(helix_df$x, helix_df$y, helix_df$pos, "helix", "raw")
      helix_result_median <- calculate_exp_fit_both(helix_df$x, helix_df$y, helix_df$pos, "helix", "median")
      
      helix_results <- rbind(helix_result_raw, helix_result_median)
      helix_results$assay <- assay_sele
      helix_results$distance_assay <- dist_assay
      helix_results$exclude_interface <- exclude_interface
      helix_results$interface_cutoff <- ifelse(exclude_interface, interface_cutoff, NA)
      results <- rbind(results, helix_results)
    }
    
    # Calculate for other regions - BOTH raw and median
    other_positions <- setdiff(main_df$pos, c(sheet_positions, helix_positions))
    other_df <- main_df[main_df$pos %in% other_positions, ]
    if (nrow(other_df) >= 5) {
      other_result_raw <- calculate_exp_fit_both(other_df$x, other_df$y, other_df$pos, "other", "raw")
      other_result_median <- calculate_exp_fit_both(other_df$x, other_df$y, other_df$pos, "other", "median")
      
      other_results <- rbind(other_result_raw, other_result_median)
      other_results$assay <- assay_sele
      other_results$distance_assay <- dist_assay
      other_results$exclude_interface <- exclude_interface
      other_results$interface_cutoff <- ifelse(exclude_interface, interface_cutoff, NA)
      results <- rbind(results, other_results)
    }
  }
  
  return(results)
}

# ===============================
# Main Analysis for All Assays (with interface exclusion option) - BOTH raw and median
# ===============================
analyze_all_assays_both_exclude_interface <- function(input_template, anno_file, output_file = NULL, exclude_interface = TRUE, interface_cutoff = 5) {
  
  assays <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
  
  all_results <- data.frame()
  
  for (assay in assays) {
    input_file <- gsub("ASSAY", assay, input_template)
    cat("Processing", assay, "...\n")
    
    result <- calculate_exp_fit_parameters_both_exclude_interface(
      input = input_file,
      assay_sele = assay,
      anno_file = anno_file,
      exclude_interface = exclude_interface,
      interface_cutoff = interface_cutoff
    )
    
    # Add group information
    result$group <- ifelse(assay %in% c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27"), "BI1", "BI2")
    
    all_results <- rbind(all_results, result)
  }
  
  # Save results if output file specified
  if (!is.null(output_file)) {
    fwrite(all_results, output_file)
    cat("Results saved to:", output_file, "\n")
  }
  
  return(all_results)
}

# ===============================
# Helper Functions for BOTH data types
# ===============================

# Helper function to convert p-value to significance stars
get_significance_stars <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# ===============================
# CODE 1: Overall Structure Only (One Plot) - Excluding Interface - BOTH raw and median with BI1/BI2 colors
# ===============================

# Function to create overall structure plots (excluding interface) with BI1/BI2 color scheme
create_overall_plot_both_exclude_interface <- function(results_data) {
  
  # Filter data - ONLY KEEP DIAGONAL (assay == distance_assay) and overall region
  plot_data <- results_data %>%
    filter(region == "overall" & !is.na(a_value) & !is.na(b_value) & assay == distance_assay) %>%
    mutate(
      assay = factor(assay, levels = c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")),
      group = factor(group, levels = c("BI1", "BI2")),
      data_type = factor(data_type, levels = c("raw", "median")),
      # Create combined group_data_type for coloring
      group_data_type = paste0(group, "_", data_type),
      # Add significance stars
      significance_a = sapply(p_value, get_significance_stars),
      significance_b = sapply(p_value, get_significance_stars)
    )
  
  if (nrow(plot_data) == 0) {
    cat("No data available for overall structure\n")
    return(NULL)
  }
  
  # Define BI1/BI2 color scheme
  group_data_type_colors <- c(
    "BI1_raw" = "#1B38A6",      # Deep Blue for BI1 raw
    "BI1_median" = "#75C2F6",   # Light Blue for BI1 median
    "BI2_raw" = "#F4270C",      # Deep Red for BI2 raw
    "BI2_median" = "#FFB0A5"    # Light Red for BI2 median
  )
  
  group_data_type_labels <- c(
    "BI1_raw" = "BI1 - Raw Data",
    "BI1_median" = "BI1 - Median Data", 
    "BI2_raw" = "BI2 - Raw Data",
    "BI2_median" = "BI2 - Median Data"
  )
  
  # Calculate y positions for significance annotations (above bars)
  plot_data <- plot_data %>%
    group_by(assay, group_data_type) %>%
    mutate(
      y_position_a = a_value * 1.15,
      y_position_b = ifelse(b_value >= 0, b_value * 1.15, b_value * 0.85)
    ) %>%
    ungroup()
  
  # Adjust y positions to be within reasonable limits
  max_a <- max(plot_data$a_value, na.rm = TRUE) * 1.2
  plot_data$y_position_a <- pmin(plot_data$y_position_a, max_a)
  
  # Create a-value plot with BI1/BI2 colored stars
  p_a <- ggplot(plot_data, aes(x = assay, y = a_value, fill = group_data_type)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_a, label = significance_a, color = group_data_type), 
              position = position_dodge(0.8), 
              vjust = -0.5, size = 3, fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = group_data_type_colors, labels = group_data_type_labels) +
    scale_color_manual(values = group_data_type_colors) +
    scale_y_continuous(limits = c(0, max_a)) +
    labs(
      title = "Initial energy (a) Values - Overall Structure\n(Excluding Binding Interface < 5Å)",
      x = "Assay",
      y = "Initial energy (a) value",
      fill = "Group & Data Type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    )
  
  # Create b-value plot with BI1/BI2 colored stars
  p_b <- ggplot(plot_data, aes(x = assay, y = b_value, fill = group_data_type)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_b, label = significance_b, color = group_data_type), 
              position = position_dodge(0.8), 
              vjust = ifelse(plot_data$b_value >= 0, -0.5, 1.5), 
              size = 3, fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = group_data_type_colors, labels = group_data_type_labels) +
    scale_color_manual(values = group_data_type_colors) +
    scale_y_continuous(limits = c(-0.15, 0.05)) +
    labs(
      title = "Decay Rate (b) Values - Overall Structure\n(Excluding Binding Interface < 5Å)",
      x = "Assay",
      y = "Decay Rate (b) value",
      fill = "Group & Data Type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    )
  
  # Combine plots
  final_plot <- p_a + p_b + plot_layout(ncol = 2)
  
  return(final_plot)
}

# ===============================
# CODE 2: Secondary Structure Regions (One Plot) - Excluding Interface - BOTH raw and median with BI1/BI2 colors
# ===============================

# Function to create secondary structure plots (excluding interface) with BI1/BI2 color scheme
create_secondary_structure_plot_both_exclude_interface <- function(results_data) {
  
  # Define regions and their titles
  regions <- list(
    sheet = "β-Sheet Regions", 
    helix = "α-Helix Regions",
    other = "Other Regions"
  )
  
  # Prepare data for all secondary structure regions
  plot_data_list <- list()
  
  for (region_name in names(regions)) {
    region_data <- results_data %>%
      filter(region == region_name & !is.na(a_value) & !is.na(b_value) & assay == distance_assay) %>%
      mutate(
        assay = factor(assay, levels = c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")),
        group = factor(group, levels = c("BI1", "BI2")),
        data_type = factor(data_type, levels = c("raw", "median")),
        # Create combined group_data_type for coloring
        group_data_type = paste0(group, "_", data_type),
        region_label = regions[[region_name]],
        # Add significance stars
        significance_a = sapply(p_value, get_significance_stars),
        significance_b = sapply(p_value, get_significance_stars)
      )
    
    if (nrow(region_data) > 0) {
      plot_data_list[[region_name]] <- region_data
    }
  }
  
  if (length(plot_data_list) == 0) {
    cat("No data available for secondary structure regions\n")
    return(NULL)
  }
  
  # Combine all region data
  all_plot_data <- bind_rows(plot_data_list)
  
  # Define BI1/BI2 color scheme
  group_data_type_colors <- c(
    "BI1_raw" = "#1B38A6",      # Deep Blue for BI1 raw
    "BI1_median" = "#75C2F6",   # Light Blue for BI1 median
    "BI2_raw" = "#F4270C",      # Deep Red for BI2 raw
    "BI2_median" = "#FFB0A5"    # Light Red for BI2 median
  )
  
  group_data_type_labels <- c(
    "BI1_raw" = "BI1 - Raw Data",
    "BI1_median" = "BI1 - Median Data", 
    "BI2_raw" = "BI2 - Raw Data",
    "BI2_median" = "BI2 - Median Data"
  )
  
  # Calculate y positions for significance annotations
  all_plot_data <- all_plot_data %>%
    group_by(region, assay, group_data_type) %>%
    mutate(
      y_position_a = a_value * 1.15,
      y_position_b = ifelse(b_value >= 0, b_value * 1.15, b_value * 0.85)
    ) %>%
    ungroup()
  
  # Create a-value plot for all secondary structures with BI1/BI2 colored stars
  p_a <- ggplot(all_plot_data, aes(x = assay, y = a_value, fill = group_data_type)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_a, label = significance_a, color = group_data_type), 
              position = position_dodge(0.8), 
              vjust = -0.5, size = 3, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ region_label, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = group_data_type_colors, labels = group_data_type_labels) +
    scale_color_manual(values = group_data_type_colors) +
    labs(
      title = "Initial energy (a) Values by Secondary Structure\n(Excluding Binding Interface < 5Å)",
      x = "Assay",
      y = "Initial energy (a) value",
      fill = "Group & Data Type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 10)
    )
  
  # Create b-value plot for all secondary structures with BI1/BI2 colored stars
  p_b <- ggplot(all_plot_data, aes(x = assay, y = b_value, fill = group_data_type)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_text(aes(y = y_position_b, label = significance_b, color = group_data_type), 
              position = position_dodge(0.8), 
              vjust = ifelse(all_plot_data$b_value >= 0, -0.5, 1.5), 
              size = 3, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ region_label, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = group_data_type_colors, labels = group_data_type_labels) +
    scale_color_manual(values = group_data_type_colors) +
    labs(
      title = "Decay Rate (b) Values by Secondary Structure\n(Excluding Binding Interface < 5Å)",
      x = "Assay",
      y = "Decay Rate (b) value",
      fill = "Group & Data Type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 10)
    )
  
  # Combine plots
  final_plot <- p_a / p_b + plot_layout(heights = c(1, 1))
  
  return(final_plot)
}


# ===============================
# Execute Analysis - Excluding Binding Interface - BOTH raw and median with BI1/BI2 colors
# ===============================

# Run analysis excluding binding interface (<5Å) with BOTH data types
results_both_exclude_interface <- analyze_all_assays_both_exclude_interface(
  input_template = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901_2/task_901/weights/weights_Binding_ASSAY.txt",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",  
  #output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251117/exponential_fit_results_both_diagonal_exclude_interface2.csv",
  exclude_interface = TRUE,
  interface_cutoff = 5
)

# CODE 1: Create overall structure plot (excluding interface) with BI1/BI2 color scheme
overall_plot_both_exclude_interface <- create_overall_plot_both_exclude_interface(results_both_exclude_interface)

overall_plot_both_exclude_interface

# Save overall plot using ggsave (outside the function)
if (!is.null(overall_plot_both_exclude_interface)) {
  ggsave(
    filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251117/overall_structure_both_BI_colors_exclude_interface2.pdf",
    plot = overall_plot_both_exclude_interface,
    width = 16,
    height = 6,
    device = cairo_pdf
  )
  cat("Overall structure plot (BI1/BI2 color scheme, excluding interface) saved.\n")
}

# CODE 2: Create secondary structure plot (excluding interface) with BI1/BI2 color scheme
secondary_plot_both_exclude_interface <- create_secondary_structure_plot_both_exclude_interface(results_both_exclude_interface)

secondary_plot_both_exclude_interface

# Save secondary structure plot using ggsave (outside the function)
if (!is.null(secondary_plot_both_exclude_interface)) {
  ggsave(
    filename = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f5/20251117/secondary_structure_both_BI_colors_exclude_interface2.pdf",
    plot = secondary_plot_both_exclude_interface,
    width = 16,
    height = 10,
    device = cairo_pdf
  )
  cat("Secondary structure plot (BI1/BI2 color scheme, excluding interface) saved.\n")
}

# Print summary statistics for both data types
print_summary_statistics_both(results_both_exclude_interface)

cat("=== Analysis Complete ===\n")
cat("1. Results for both raw and median data (excluding binding interface <5Å) saved\n")
cat("2. Bar plots with BI1/BI2 color scheme created\n")
cat("3. Color Scheme:\n")
cat("   - BI1 Raw: Deep Blue (#1B38A6)\n")
cat("   - BI1 Median: Light Blue (#75C2F6)\n")
cat("   - BI2 Raw: Deep Red (#F4270C)\n")
cat("   - BI2 Median: Light Red (#FF6A56)\n")








