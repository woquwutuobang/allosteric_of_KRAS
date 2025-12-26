# KRAS DARPins Binding Fitness Analysis
# This script analyzes and visualizes the fitness effects of mutations on KRAS binding to DARPins K13 and K19

library(ggplot2)
library(data.table)
library(dplyr)
library(wlab.block)
library(krasddpcams)

# Wild-type KRAS sequence (residues 2-189)
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"



# Function to merge two DiMSum data frames
merge_dimsum_data <- function(merge_1, merge_2) {
  a1 <- as.character(substitute(merge_1))
  a2 <- as.character(substitute(merge_2))
  merge_1[, assay := a1]
  merge_2[, assay := a2]
  output <- rbind(merge_1, merge_2)
  return(output)
}

# Function to identify mutation positions
identify_mutation_positions <- function(input, wt_aa) {
  output <- input
  output[, AA_Pos1 := which(unlist(strsplit(aa_seq, "")) != 
                              unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[1], aa_seq]
  output[, AA_Pos2 := which(unlist(strsplit(aa_seq, "")) != 
                              unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[2], aa_seq]
  
  for (i in 1:188) {
    output[AA_Pos1 == i, mt1 := substr(aa_seq, i, i)]
  }
  for (i in 1:188) {
    output[AA_Pos2 == i, mt2 := substr(aa_seq, i, i)]
  }
  for (i in 1:188) {
    output[AA_Pos1 == i, wtcodon1 := substr(wt_aa, i, i)]
  }
  for (i in 1:188) {
    output[AA_Pos2 == i, wtcodon2 := substr(wt_aa, i, i)]
  }
  
  output[, codon1 := substr(aa_seq, AA_Pos1, AA_Pos1)]
  output[, codon2 := substr(aa_seq, AA_Pos2, AA_Pos2)]
  output[, AA_Pos1 := AA_Pos1 + 1]  # Convert to 1-based indexing
  output[, AA_Pos2 := AA_Pos2 + 1]
  output[, mt1 := paste0(wtcodon1, AA_Pos1, codon1)]
  output[, mt2 := paste0(wtcodon2, AA_Pos2, codon2)]
  
  return(output)
}

# ==============================
# General Function: Plot binder binding fitness
# ==============================
plot_binding_fitness <- function(input, assay_sele, anno) {
  
  # 1. Separate stability data (stab) and specified binding data
  input_abundance <- input[assay == "stab", ]
  input_abundance_single <- krasddpcams__nor_overlap_single_mt_fitness(input_abundance)
  
  input_binding <- input[assay == assay_sele, ]
  input_binding_single <- krasddpcams__nor_overlap_single_mt_fitness(input_binding)
  
  # 2. Merge both conditions
  input_long <- rbind(input_abundance_single, input_binding_single)
  input_dc <- dcast(input_long,
                    nt_seq + aa_seq + Nham_aa + AA_Pos1 + wtcodon1 ~ assay,
                    value.var = c("nor_fitness_nooverlap",
                                  "nor_fitness_nooverlap_sigma",
                                  "nor_gr_nooverlap",
                                  "nor_gr_nooverlap_sigma"),
                    drop = TRUE)
  
  input_single_pos <- input_dc
  input_single_pos[, position := AA_Pos1]
  input_single_pos[, WT_AA := wtcodon1]
  
  # 3. Merge with annotation file
  anno_single <- merge(input_single_pos, anno,
                       by.x = c("position", "WT_AA"),
                       by.y = c("Pos_real", "codon"),
                       all = TRUE)
  
  # 4. Label binding interface residues
  interface_col <- paste0("scHAmin_ligand_", assay_sele)
  anno_single[, type_bs := "others"]
  anno_single[get(interface_col) < 5, type_bs := "binding_interface"]
  
  # 5. Create plot
  p <- ggplot() +
    # Non-binding interface mutations
    geom_point(data = anno_single[position > 1 & type_bs == "others", ],
               aes(x = nor_fitness_nooverlap_stab,
                   y = get(paste0("nor_fitness_nooverlap_", assay_sele))),
               color = "black", alpha = 0.6, size = 1.5) +
    # Binding interface mutations
    geom_point(data = anno_single[position > 1 & type_bs == "binding_interface", ],
               aes(x = nor_fitness_nooverlap_stab,
                   y = get(paste0("nor_fitness_nooverlap_", assay_sele))),
               color = "#F4270C", alpha = 0.6, size = 2) +
    
    # Fixed coordinate ranges
    coord_cartesian(xlim = c(-1.3, 1.1), ylim = c(-1.3, 0.3)) +
    
    # Classic theme
    theme_classic() +
    
    # Rotate X-axis ticks + control font
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "right",
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(3.1, "mm"),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    
    labs(
      x = "Abundance Fitness",
      y = paste0(assay_sele, " Binding Fitness"),  # Automatically named based on binder
      color = NULL
    )
  
  return(p)
}



# Load and normalize experimental data
stability_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData"
)

K13_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block1_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData"
)

K19_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K19_block1_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData",
  block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData"
)


RAF1_nor_df <- krasddpcams__normalize_growthrate_fitness(
  block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
  block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData"
)


# Simplify variable names
stab <- stability_nor_df
K13 <- K13_nor_df  
K19 <- K19_nor_df
RAF1 <- RAF1_nor_df


# Merge abundance and binding data
all_data_K13 <- merge_dimsum_data(stab, K13)
all_data_K19 <- merge_dimsum_data(stab, K19)

all_data_RAF1 <- merge_dimsum_data(stab, RAF1)


# Identify mutation positions
all_data_pos_K13 <- identify_mutation_positions(all_data_K13, wt_aa)
all_data_pos_K19 <- identify_mutation_positions(all_data_K19, wt_aa)
all_data_pos_RAF1 <- identify_mutation_positions(all_data_RAF1, wt_aa)


# Load annotation data
anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv")
anno[, Pos_real := Pos] 
#names(anno)
summery<-ddG_data_assay(input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
                        wt_aa = wt_aa)
#names(summery)
anno<-merge(anno,summery,by="Pos_real",all=T)

# Generate binding fitness plots
plot_K13 <- plot_binding_fitness(
  input = all_data_pos_K13,
  assay_sele = "K13",
  anno = anno
)

plot_K19 <- plot_binding_fitness(
  input = all_data_pos_K19, 
  assay_sele = "K19",
  anno = anno
)


plot_RAF1 <- plot_binding_fitness(
  input = all_data_pos_RAF1, 
  assay_sele = "RAF1",
  anno = anno
)

# Display plots
print(plot_K13)
print(plot_K19)
print(plot_RAF1)
# Save plots
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/K13_binding_vs_abundance_fitness 3.pdf", plot = plot_K13, device = cairo_pdf, 
       height = 4, width = 4)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/K19_binding_vs_abundance_fitness 3.pdf", plot = plot_K19, device = cairo_pdf,
       height = 4, width = 4)
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251031_version_all_figure/f1/20251103/RAF1_binding_vs_abundance_fitness.pdf", plot = plot_RAF1, device = cairo_pdf,
       height = 4, width = 4)
