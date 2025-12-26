library(ggplot2)
library(dplyr)
library(patchwork)
library(wlab.block)

# ============================================================
# Function: plot_assay_fitness_density
# Description:
#   Generate density plots of normalized fitness distributions for 
#   synonymous, missense, and stop mutations across experimental blocks
#   for any given assay type (abundance, binding, etc.).
#   The top panel shows the global distribution across all blocks,
#   while the bottom panel shows block-specific distributions.
# 
# Parameters:
#   assay_type: Name of the assay (e.g., "abundance", "RAF1", "SOS1", etc.)
#   block1: Path to block1 fitness data file
#   block2: Path to block2 fitness data file  
#   block3: Path to block3 fitness data file

# Returns:
#   Combined ggplot object with two panels
# ============================================================


plot_fitness_density_overall <- function(
    assay_type,
    block1,
    block2,
    block3
) {
  
  # Load normalized fitness data from three experimental blocks
  nor_fit <- nor_fitness(
    block1 = block1,
    block2 = block2,
    block3 = block3
  )
  
  # Classify mutation types
  nor_fit_classified <- nor_fit %>%
    mutate(
      mut_type = case_when(
        # Synonymous
        Nham_aa == 0 & Nham_nt > 0 ~ "Synonymous",
        # Stop
        STOP == TRUE | STOP_readthrough == TRUE ~ "Stop",
        # Missense
        Nham_aa > 0 & indel == FALSE &
          STOP == FALSE & STOP_readthrough == FALSE ~ "Missense"
      )
    ) %>%
    filter(!is.na(mut_type))
  
  # Print mutation type distribution (useful QC)
  cat("Mutation type distribution for", assay_type, ":\n")
  print(table(nor_fit_classified$mut_type))
  
  # Prepare plotting data
  nor_fit_plot <- nor_fit_classified %>%
    mutate(
      mut_type = factor(
        mut_type,
        levels = c("Synonymous", "Missense", "Stop")
      )
    )
  
  # Overall density plot (all blocks combined)
  p_overall <- ggplot(
    nor_fit_plot,
    aes(x = nor_fitness, color = mut_type)
  ) +
    geom_density(size = 1) +
    scale_color_manual(
      values = c(
        "Synonymous" = "#09B636",
        "Missense"   = "#F4AD0C",
        "Stop"       = "#FF6A56"
      ),
      name = "Mutation Type"
    ) +
    labs(
      title = paste(toupper(assay_type), "Fitness Distribution"),
      x = "Normalized Fitness",
      y = "Density"
    ) +
    xlim(-1.5, 0.5) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10),
      text = element_text(size = 10),
      axis.text = element_text(size = 10)
    )
  
  return(p_overall)
}



p_K13 <- plot_fitness_density_overall(
  assay_type = "K13",
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData"
)

print(p_K13)

ggsave(
  filename = "RAF1_fitness_density_overall.pdf",
  plot = p_K13,
  width = 4,
  height = 3,
  units = "in"
)



