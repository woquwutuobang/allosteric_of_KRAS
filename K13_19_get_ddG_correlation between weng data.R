library(data.table)
library(ggplot2)
library(ggpubr)

# ============================================================
# Function: plot_ddG_correlation
# Description:
#   Perform correlation analysis between two ddG datasets and create
#   a scatter plot with correlation coefficient and p-value.
#   Datasets are merged by position and ID.
# 
# Parameters:
#   data1: First ddG dataset (data.frame with "Pos", "id", and "mean_kcal/mol")
#   data2: Second ddG dataset (data.frame with "Pos", "id", and "mean_kcal/mol")  
#   data1_name: Name for first dataset (for X-axis label)
#   data2_name: Name for second dataset (for Y-axis label)
#   output_file: Optional path to save the plot
#   limits: Axis limits (default: c(-1.6, 2.8))
#
# Returns:
#   ggplot object with correlation plot
# ============================================================

plot_ddG_correlation <- function(data1, data2, data1_name = "Dataset 1", data2_name = "Dataset 2", output_file = NULL, limits = c(-1.6, 2.8)) {
  
  # Convert data to data.table format
  data1_dt <- as.data.table(data1)
  data2_dt <- as.data.table(data2)
  
  # Check for duplicate positions
  cat("Checking for duplicate positions...\n")
  cat("Data1 duplicate positions:", sum(duplicated(data1_dt$Pos)), "\n")
  cat("Data2 duplicate positions:", sum(duplicated(data2_dt$Pos)), "\n")
  
  # Merge datasets by position and ID
  merged_data <- merge(data1_dt, data2_dt, by = c("Pos", "id"), suffixes = c("_1", "_2"))
  
  # Extract mean ddG values
  x_values <- merged_data[["mean_kcal/mol_1"]]
  y_values <- merged_data[["mean_kcal/mol_2"]]
  
  # Calculate correlation statistics
  cor_test <- cor.test(x_values, y_values, method = "pearson")
  r_value <- round(cor_test$estimate, 3)
  p_value <- cor_test$p.value
  
  # Format p-value for display
  p_text <- ifelse(p_value < 0.001, "p < 0.001", 
                   ifelse(p_value < 0.01, "p < 0.01",
                          ifelse(p_value < 0.05, "p < 0.05",
                                 paste0("p = ", round(p_value, 3)))))
  
  # Create correlation plot
  p <- ggplot(merged_data, aes(x = x_values, y = y_values)) +
    geom_point(alpha = 0.35, size = 0.6, color = "#75C2F6") +
    geom_smooth(method = "lm", se = TRUE, color = "#FF6A56", linewidth = 0.6) +
    # Add correlation statistics annotation
    annotate("text", 
             x = -Inf, y = Inf, 
             label = paste0("R = ", r_value, "\n", p_text),
             hjust = -0.1, vjust = 1.5, 
             size = 8/.pt, color = "black") +
    # Set fixed axis limits
    scale_x_continuous(
      limits = limits,
      breaks = scales::pretty_breaks(n = 6)
    ) +
    scale_y_continuous(
      limits = limits,
      breaks = scales::pretty_breaks(n = 6)
    ) +
    labs(
      x = paste(data1_name, "ΔΔG (kcal/mol)"),
      y = paste(data2_name, "ΔΔG (kcal/mol)"),
      title = "ΔΔG Correlation Analysis"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8, color = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4)
    ) +
    coord_fixed(ratio = 1)  # Ensure equal aspect ratio
  
  # Print correlation statistics
  cat("Correlation Analysis:\n")
  cat("Number of positions:", nrow(merged_data), "\n")
  cat("Pearson R:", r_value, "\n")
  cat("P-value:", p_value, "\n")
  cat("P-value display:", p_text, "\n")
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, device = cairo_pdf, width = 4, height = 4, units = "in", dpi = 300)
    cat("Plot saved to:", output_file, "\n")
  }
  
  return(p)
}



# ============================================================
# Example usage:
# ============================================================

### folding
ddG_folding_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Folding.txt")

#names(ddG_folding_new)
ddG_folding_new<-ddG_folding_new[,c(1,3,20:22)]

ddG_folding_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Folding.txt")
#names(ddG_folding_weng)

ddG_folding_weng<-ddG_folding_weng[,c(1,3,20:22)]

data1 = ddG_folding_new
data2 = ddG_folding_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_folding_new,
  data2 = ddG_folding_weng,
  data1_name = "Folding_this study",
  data2_name = "Folding_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between folding_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)





### RAF1
ddG_RAF1_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt")

#names(ddG_folding_new)
ddG_RAF1_new<-ddG_RAF1_new[,c(1,3,20:22)]

ddG_RAF1_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_RAF.txt")
#names(ddG_folding_weng)

ddG_RAF1_weng<-ddG_RAF1_weng[,c(1,3,20:22)]

data1 = ddG_RAF1_new
data2 = ddG_RAF1_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_RAF1_new,
  data2 = ddG_RAF1_weng,
  data1_name = "RAF1_this study",
  data2_name = "RAF1_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between RAF1_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)







### SOS1
ddG_SOS1_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_SOS.txt")

#names(ddG_folding_new)
ddG_SOS1_new<-ddG_SOS1_new[,c(1,3,20:22)]

ddG_SOS1_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_SOS.txt")
#names(ddG_folding_weng)

ddG_SOS1_weng<-ddG_SOS1_weng[,c(1,3,20:22)]

data1 = ddG_SOS1_new
data2 = ddG_SOS1_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_SOS1_new,
  data2 = ddG_SOS1_weng,
  data1_name = "SOS1_this study",
  data2_name = "SOS1_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between SOS1_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)








### K55
ddG_K55_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K55.txt")

#names(ddG_folding_new)
ddG_K55_new<-ddG_K55_new[,c(1,3,20:22)]

ddG_K55_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_K55.txt")
#names(ddG_folding_weng)

ddG_K55_weng<-ddG_K55_weng[,c(1,3,20:22)]

data1 = ddG_K55_new
data2 = ddG_K55_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_K55_new,
  data2 = ddG_K55_weng,
  data1_name = "K55_this study",
  data2_name = "K55_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between K55_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)



### K27
ddG_K27_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt")

#names(ddG_folding_new)
ddG_K27_new<-ddG_K27_new[,c(1,3,20:22)]

ddG_K27_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_K27.txt")
#names(ddG_folding_weng)

ddG_K27_weng<-ddG_K27_weng[,c(1,3,20:22)]

data1 = ddG_K27_new
data2 = ddG_K27_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_K27_new,
  data2 = ddG_K27_weng,
  data1_name = "K27_this study",
  data2_name = "K27_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between K27_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)



### RALGDS
ddG_RALGDS_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAL.txt")

#names(ddG_folding_new)
ddG_RALGDS_new<-ddG_RALGDS_new[,c(1,3,20:22)]

ddG_RALGDS_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_RAL.txt")
#names(ddG_folding_weng)

ddG_RALGDS_weng<-ddG_RALGDS_weng[,c(1,3,20:22)]

data1 = ddG_RALGDS_new
data2 = ddG_RALGDS_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_RALGDS_new,
  data2 = ddG_RALGDS_weng,
  data1_name = "RALGDS_this study",
  data2_name = "RALGDS_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between RALGDS_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)



### PIK3CG
ddG_PIK3CG_new<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_PI3.txt")

#names(ddG_folding_new)
ddG_PIK3CG_new<-ddG_PIK3CG_new[,c(1,3,20:22)]

ddG_PIK3CG_weng<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/Data/energy_article_data/weights_Binding_PI3.txt")
#names(ddG_folding_weng)

ddG_PIK3CG_weng<-ddG_PIK3CG_weng[,c(1,3,20:22)]

data1 = ddG_PIK3CG_new
data2 = ddG_PIK3CG_weng

# Perform correlation analysis and create plot
correlation_plot <- plot_ddG_correlation(
  data1 = ddG_PIK3CG_new,
  data2 = ddG_PIK3CG_weng,
  data1_name = "PIK3CG_this study",
  data2_name = "PIK3CG_Weng",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s2/20251010/ddG correlation between PIK3CG_new and folding_weng test.pdf",
  limits = c(-1.6, 2.8)
)

# Display the plot
print(correlation_plot)
