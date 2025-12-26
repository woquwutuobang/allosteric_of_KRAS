# Fitness Correlation Analysis Across Replicates
# This script creates correlation plots to assess reproducibility between experimental replicates

library(data.table)
library(krasddpcams)
library(ggplot2)
library(GGally)
library(ggpubr)

# Define color scheme
colour_scheme<-list(
  "blue"="#1B38A6",#rgb(27, 56, 166)
  "red"="#F4270C",#rgb(244, 39, 12)
  "orange"="#F4AD0C",#rgb(244, 173, 12)
  "green"="#09B636",#rgb(9, 182, 54)
  "yellow"="#F1DD10",#rgb(241, 221, 16)
  "purple"="#C68EFD",#rgb(198, 142, 253)
  "hot pink"="#FF0066",#rgb(255, 0, 102)
  "light blue"="#75C2F6",#rgb(117, 194, 246)
  "light red"="#FF6A56",#rgb(255, 106, 86)       # The red ones are unified with this, but the heatmap ones remain unchanged.
  "dark red"="#A31300",#rgb(163, 19, 0)
  "dark green"="#007A20",#rgb(0, 122, 32)
  "pink"="#FFB0A5" #rgb(255, 176, 165)
)


#' Create fitness correlation plot across experimental blocks
#'
#' @param input Normalized fitness data
#' @param assay_name Name of the assay for plot title
#' @param colour_scheme Color scheme for different blocks
#' @return GGally plot object showing correlations between replicates

krasddpcams__plot_fitness_correlation_blocks <- function(
    input,
    assay_name,
    colour_scheme
){
  d <- GGally::ggpairs(
    input,
    columns = c(18,19,20), ##K13/K19
    #columns = c(18, 20, 22),##abundance/RAF1
    columnLabels = c("replicate 1", "replicate 2", "replicate 3"),
    mapping = ggplot2::aes(color = block),
    lower = list(continuous = function(data, mapping, ...){
      ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_bin2d(bins = 100, alpha = 0.2) +
        ggplot2::scale_fill_gradient(low = "white", high = "black")
    }),
    upper = list(continuous = GGally::wrap("cor", 
                                           mapping = ggplot2::aes(color = block), 
                                           size = 8 * 0.35)),
    diag = list(continuous = "blankDiag")
  ) +
    ggplot2::scale_color_manual(values = c(
      colour_scheme[["red"]], 
      colour_scheme[["green"]], 
      colour_scheme[["blue"]]
    )) +
    ggplot2::ggtitle(assay_name) +
    ggpubr::theme_classic2() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 8),
      strip.text.x = ggplot2::element_text(size = 8),
      strip.text.y = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8)
    )
  
  return(d)
}




K19_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K13_block1_0710/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData",
                                                        block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                                                        block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")

K13_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K19_block1/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates_cleaned.RData",
                                                        block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K19_block2/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates_cleaned.RData",
                                                        block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")


stability_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                              block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                                                              block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")

RAF1_nor_df<- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                                        block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
                                                        block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")



###
d_K13 <- krasddpcams__plot_fitness_correlation_blocks(K13_nor_df,"K13", colour_scheme)
print(d_K13)

ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/figureS1a_fitness_correlation_blocks_K13.pdf", d_K13, device = cairo_pdf,height = 4, width=4)


###
d_K19 <- krasddpcams__plot_fitness_correlation_blocks(K19_nor_df,"K19", colour_scheme)
print(d_K19)

ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/figureS1a_fitness_correlation_blocks_K19.pdf", d_K19, device = cairo_pdf,height = 4, width=4)


###
d_stability <- krasddpcams__plot_fitness_correlation_blocks(stability_nor_df,"stability", colour_scheme)
print(d_stability)

ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/figureS1a_fitness_correlation_blocks_stability.pdf", d_stability, device = cairo_pdf,height = 4, width=4)


###
d_RAF1 <- krasddpcams__plot_fitness_correlation_blocks(RAF1_nor_df,"RAF1", colour_scheme)
print(d_RAF1)

ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure_s1/20251009/figureS1a_fitness_correlation_blocks_RAF1.pdf", d_RAF1, device = cairo_pdf,height = 4, width=4)

