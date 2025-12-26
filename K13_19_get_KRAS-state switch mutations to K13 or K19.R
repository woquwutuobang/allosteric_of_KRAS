#############################=================================
#### KRAS-state switch mutations to K13 


library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Function to classify mutations based on new criteria
classify_mutations <- function(ddG_path, anno_path, assay_sele, wt_aa) {
  # Read input files as data.tables
  ddG <- fread(ddG_path)  # Make sure ddG is a data.table
  anno <- fread(anno_path)  # Read annotation file as data.table
  
  # Ensure ddG is treated as a data.table if not already
  setDT(ddG)
  
  # Data processing
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Merge with annotation data
  data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
  
  # Define binding site types
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  # Calculate regulatory threshold
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  # Classify site types
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # Merge mutation data with site type information
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # Identify allosteric mutations
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # Classify mutation types based on new criteria
  data_plot_mutation[, mutation_type := ifelse(
    site_type == "GTP binding interface site" & allosteric_mutation == TRUE, "GTP binding allosteric mutation",
    ifelse(site_type == "GTP binding interface site" & allosteric_mutation == FALSE, "GTP binding other mutation", 
           "Other mutation"))]
  
  return(data_plot_mutation)
}

####K13
mutation_K13 <- classify_mutations("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",
                                   "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
                                   assay_sele = "K13", wt_aa = wt_aa)

#names(mutation_K13)
mutation_K13<-mutation_K13[,c(2:4,23:26,28)]




###RAF1
mutation_RAF1 <- classify_mutations("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
                                    "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
                                    assay_sele = "RAF1", wt_aa = wt_aa)

#names(mutation_RAF1)
mutation_RAF1<-mutation_RAF1[,c(2:4,23:26,28)]



###K27
mutation_K27 <- classify_mutations("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt",
                                   "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
                                   assay_sele = "K27", wt_aa = wt_aa)

#names(mutation_K27)
mutation_K27<-mutation_K27[,c(2:4,23:26,28)]



merged_K27_RAF1 <- merge(mutation_K27, mutation_RAF1, by = c("Pos_real","wt_codon","mt_codon","mt"), suffixes = c("_K27", "_RAF1"))
names(merged_K27_RAF1)
# Step 1: 筛选所有 GTP binding allosteric mutation 的突变
GTP_allosteric_mutations <- merged_K27_RAF1 %>%
  dplyr::filter(
    mutation_type_K27 == "GTP binding allosteric mutation" &
      mutation_type_RAF1 == "GTP binding allosteric mutation"
  )

# Step 2: 进一步筛选位于第二象限的突变
KRAS_switch_mutations <- GTP_allosteric_mutations %>%
  dplyr::filter(`mean_kcal/mol_K27` < 0 & `mean_kcal/mol_RAF1` > 0)

# Step 3: 将这些突变定义为 KRAS_switch_mutation
KRAS_switch_mutations$mutation_type <- "KRAS_switch_mutation"

# Step 4: 打印出 KRAS_switch_mutation 对应的 `mt` 列（突变类型）
print(KRAS_switch_mutations$mt)

# Step 5: 合并 K13 和 RAF1 的数据框，准备绘制 K13 vs RAF1 的散点图
merged_K13_RAF1 <- merge(mutation_K13, mutation_RAF1, by = c("Pos_real","wt_codon","mt_codon","mt"), suffixes = c("_K13", "_RAF1"))

# Step 6: 为 K13 vs RAF1 散点图创建特殊的高亮显示列
# 如果突变在 KRAS_switch_mutations 中，标记为 "highlight"，否则为 "other"
merged_K13_RAF1 <- merged_K13_RAF1 %>%
  mutate(highlight = ifelse(mt %in% KRAS_switch_mutations$mt, "highlight", "other"))



ggplot(merged_K13_RAF1, aes(x = `mean_kcal/mol_K13`, y = `mean_kcal/mol_RAF1`, color = highlight)) +
  geom_point(data = subset(merged_K13_RAF1, highlight == "other"), 
             aes(x = `mean_kcal/mol_K13`, y = `mean_kcal/mol_RAF1`), 
             color = "gray", alpha = 0.6, size = 2) +
  
  geom_point(data = subset(merged_K13_RAF1, highlight == "highlight"), 
             aes(x = `mean_kcal/mol_K13`, y = `mean_kcal/mol_RAF1`), 
             color = "#007A20", alpha = 0.8, size = 2.5) +
  
  ggrepel::geom_text_repel(
    data = subset(merged_K13_RAF1, highlight == "highlight"),
    aes(label = mt),
    color = "#007A20",
    size = 3,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2
  )+
  
  labs(
    x = bquote("Binding" ~ Delta*Delta*"G ("*"K13"*") (kcal/mol)"),
    y = bquote("Binding" ~ Delta*Delta*"G ("*"RAF1"*") (kcal/mol)"),
    title = "K13 vs RAF1 with KRAS-state Switch Mutations in highlight"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  
  # ✅ 固定坐标范围
  coord_cartesian(
    xlim = c(-1.5, 3.3),
    ylim = c(-1.5, 3.3)
  )+

  annotate(
    "point",
    x = 0.9, y = 3,
    color = "#007A20",
    size = 3
  ) +
  
  # 文字说明
  annotate(
    "text",
    x = 1, y = 3,
    label = "KRAS-state switch mutations",
    color = "#007A20",
    size = 4,
    hjust = 0
  )




ggsave(
  "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251220/KRAS-state swithch mutations to K13.png",
  width = 6,height = 6, dpi = 300,
  limitsize = FALSE
)








#############################################======================================================
#### KRAS-state switch mutations to K13

library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# Function to classify mutations based on new criteria
classify_mutations <- function(ddG_path, anno_path, assay_sele, wt_aa) {
  # Read input files as data.tables
  ddG <- fread(ddG_path)  # Make sure ddG is a data.table
  anno <- fread(anno_path)  # Read annotation file as data.table
  
  # Ensure ddG is treated as a data.table if not already
  setDT(ddG)
  
  # Data processing
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Merge with annotation data
  data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
  
  # Define binding site types
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  # Calculate regulatory threshold
  reg_threshold <- data_plot[binding_type == "binding site",
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE),
                             .SDcols = c("mean", "sigma")]
  
  # Classify site types
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  # Merge mutation data with site type information
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # Identify allosteric mutations
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # Classify mutation types based on new criteria
  data_plot_mutation[, mutation_type := ifelse(
    site_type == "GTP binding interface site" & allosteric_mutation == TRUE, "GTP binding allosteric mutation",
    ifelse(site_type == "GTP binding interface site" & allosteric_mutation == FALSE, "GTP binding other mutation", 
           "Other mutation"))]
  
  return(data_plot_mutation)
}

####K19
mutation_K19 <- classify_mutations("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt",
                                   "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
                                   assay_sele = "K19", wt_aa = wt_aa)

#names(mutation_K13)
mutation_K19<-mutation_K19[,c(2:4,23:26,28)]




###RAF1
mutation_RAF1 <- classify_mutations("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
                                    "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
                                    assay_sele = "RAF1",wt_aa = wt_aa)

#names(mutation_RAF1)
mutation_RAF1<-mutation_RAF1[,c(2:4,23:26,28)]



###K27
mutation_K27 <- classify_mutations("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt",
                                   "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
                                   assay_sele = "K27", wt_aa = wt_aa)

#names(mutation_K27)
mutation_K27<-mutation_K27[,c(2:4,23:26,28)]




merged_K27_RAF1 <- merge(mutation_K27, mutation_RAF1, by = c("Pos_real","wt_codon","mt_codon","mt"), suffixes = c("_K27", "_RAF1"))
names(merged_K27_RAF1)
# Step 1: 筛选所有 GTP binding allosteric mutation 的突变
GTP_allosteric_mutations <- merged_K27_RAF1 %>%
  dplyr::filter(
    mutation_type_K27 == "GTP binding allosteric mutation" &
      mutation_type_RAF1 == "GTP binding allosteric mutation"
  )

# Step 2: 进一步筛选位于第二象限的突变
KRAS_switch_mutations <- GTP_allosteric_mutations %>%
  dplyr::filter(`mean_kcal/mol_K27` < 0 & `mean_kcal/mol_RAF1` > 0)

# Step 3: 将这些突变定义为 KRAS_switch_mutation
KRAS_switch_mutations$mutation_type <- "KRAS_switch_mutation"

# Step 4: 打印出 KRAS_switch_mutation 对应的 `mt` 列（突变类型）
print(KRAS_switch_mutations$mt)

# Step 5: 合并 K19 和 RAF1 的数据框，准备绘制 K13 vs RAF1 的散点图
merged_K19_RAF1 <- merge(mutation_K19, mutation_RAF1, by = c("Pos_real","wt_codon","mt_codon","mt"), suffixes = c("_K19", "_RAF1"))

# Step 6: 为 K19 vs RAF1 散点图创建特殊的高亮显示列
# 如果突变在 KRAS_switch_mutations 中，标记为 "highlight"，否则为 "other"
merged_K19_RAF1 <- merged_K19_RAF1 %>%
  mutate(highlight = ifelse(mt %in% KRAS_switch_mutations$mt, "highlight", "other"))



ggplot(merged_K19_RAF1, aes(x = `mean_kcal/mol_K19`, y = `mean_kcal/mol_RAF1`, color = highlight)) +
  geom_point(data = subset(merged_K19_RAF1, highlight == "other"), 
             aes(x = `mean_kcal/mol_K19`, y = `mean_kcal/mol_RAF1`), 
             color = "gray", alpha = 0.6, size = 2) +
  
  geom_point(data = subset(merged_K19_RAF1, highlight == "highlight"), 
             aes(x = `mean_kcal/mol_K19`, y = `mean_kcal/mol_RAF1`), 
             color = "#007A20", alpha = 0.8, size = 2.5) +
  
  ggrepel::geom_text_repel(
    data = subset(merged_K19_RAF1, highlight == "highlight"),
    aes(label = mt),
    color = "#007A20",
    size = 3,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2
  )+
  
  labs(
    x = bquote("Binding" ~ Delta*Delta*"G ("*"K19"*") (kcal/mol)"),
    y = bquote("Binding" ~ Delta*Delta*"G ("*"RAF1"*") (kcal/mol)"),
    title = "K19 vs RAF1 with KRAS-state Switch Mutations in highlight"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  
  # ✅ 固定坐标范围
  coord_cartesian(
    xlim = c(-1.5, 3.3),
    ylim = c(-1.5, 3.3)
  )+
  
  annotate(
    "point",
    x = 0.9, y = 3,
    color = "#007A20",
    size = 3
  ) +
  
  # 文字说明
  annotate(
    "text",
    x = 1, y = 3,
    label = "KRAS-state switch mutations",
    color = "#007A20",
    size = 4,
    hjust = 0
  )




ggsave(
  "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251220/KRAS-state swithch mutations to K19.png",
  width = 6,height = 6, dpi = 300,
  limitsize = FALSE
)
