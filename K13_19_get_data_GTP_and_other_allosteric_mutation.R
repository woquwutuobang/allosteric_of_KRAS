library(openxlsx)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(krasddpcams)
library(ggrepel)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"


classify_mutations_comprehensive <- function(input, anno, assay_sele, wt_aa) {
  #' 综合突变分类函数
  #'
  #' @param input 包含ddG数据的文件路径或数据框
  #' @param anno 包含注释信息的文件路径或数据框
  #' @param assay_sele 要分析的assay名称（如"K13", "RAF1"等）
  #' @param wt_aa 野生型氨基酸序列
  #' @return 包含突变分类的数据框
  
  # Read and process ddG data
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
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
  
  # Read annotation data
  if (is.character(anno)) {
    anno_data <- fread(anno)
  } else {
    anno_data <- as.data.table(anno)
  }
  
  # Merge with annotation data
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  # Define binding site types
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  gxpmg_col <- paste0("GXPMG_scHAmin_ligand_", assay_sele)
  
  data_plot[get(scHAmin_col) < 5, binding_type := "binding site"]
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
  
  # Classify mutation types
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "Orthosteric site huge differences"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "Orthosteric site small differences"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == TRUE, mutation_type := "GTP binding allosteric mutation"]
  
  data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] &
                       allosteric_mutation == FALSE, mutation_type := "GTP binding other mutation"]
  
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == TRUE, mutation_type := "Allosteric mutation"]
  
  data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") &
                       allosteric_mutation == FALSE, mutation_type := "Other mutation"]
  
  # Set factor levels for mutation types
  data_plot_mutation[, mutation_type := factor(mutation_type,
                                               levels = c("Orthosteric site huge differences",
                                                          "Orthosteric site small differences",
                                                          "GTP binding allosteric mutation",
                                                          "GTP binding other mutation",
                                                          "Allosteric mutation",
                                                          "Other mutation"))]
  
  # 添加有用的列
  data_plot_mutation[, `:=`(
    assay = assay_sele,
    reg_threshold = reg_threshold,
    is_significant = allosteric_mutation
  )]
  
  # 返回关键列
  result_cols <- c("Pos", "Pos_real", "wt_codon", "mt_codon", "mt", 
                   "mean_kcal/mol", "std_kcal/mol", "mean", "std",
                   "site_type", "allosteric_mutation", "mutation_type",
                   "assay", "reg_threshold", "is_significant")
  
  result_cols <- intersect(result_cols, names(data_plot_mutation))
  
  return(data_plot_mutation[, ..result_cols])
}



#### K13
mutation_K13 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "K13", 
  wt_aa )
#names(mutation_K13)
#mutation_K13 <- mutation_K13[, c(2:4, 23:26, 28)]

#### RAF1
mutation_RAF1 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "RAF1", 
  wt_aa )
#mutation_RAF1 <- mutation_RAF1[, c(2:4, 23:26, 28)]

#### K27
mutation_K27 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K27.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "K27", 
  wt_aa )
#mutation_K27 <- mutation_K27[, c(2:4, 23:26, 28)]

# 合并数据
mutation_merge1 <- merge(mutation_K27, mutation_RAF1, 
                         by = c("Pos","Pos_real", "wt_codon", "mt_codon", "mt"), 
                         suffixes = c("_K27", "_RAF1"))

mutation_merge2 <- merge(mutation_K13, mutation_RAF1, 
                         by = c("Pos","Pos_real", "wt_codon", "mt_codon", "mt"), 
                         suffixes = c("_K13", "_RAF1"))

# 识别共同的显著突变（同时是GTP binding allosteric mutation和Allosteric mutation）
common_mutations_K27_RAF1 <- mutation_merge1[
  (mutation_type_K27 == "GTP binding allosteric mutation" | 
     mutation_type_K27 == "Allosteric mutation") &
    (mutation_type_RAF1 == "GTP binding allosteric mutation" | 
       mutation_type_RAF1 == "Allosteric mutation")
]

common_mutations_K13_RAF1 <- mutation_merge2[
  (mutation_type_K13 == "GTP binding allosteric mutation" | 
     mutation_type_K13 == "Allosteric mutation") &
    (mutation_type_RAF1 == "GTP binding allosteric mutation" | 
       mutation_type_RAF1 == "Allosteric mutation")
]




# 设置输出文件夹
out_dir <- "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/"

# 定义保留列
#keep_cols <- c("Pos_real", "mt")

# ---------------- K27 vs RAF1 ----------------
# GTP binding allosteric mutation
K27_RAF1_GTP <- common_mutations_K27_RAF1[
  mutation_type_K27 == "GTP binding allosteric mutation" | 
    mutation_type_RAF1 == "GTP binding allosteric mutation",
  .(Pos_real, mt, `mean_kcal/mol_K27`, `std_kcal/mol_K27`,
    `mean_kcal/mol_RAF1`, `std_kcal/mol_RAF1`)
]

# Allosteric mutation
K27_RAF1_Allo <- common_mutations_K27_RAF1[
  mutation_type_K27 == "Allosteric mutation" | 
    mutation_type_RAF1 == "Allosteric mutation",
  .(Pos_real, mt, `mean_kcal/mol_K27`, `std_kcal/mol_K27`,
    `mean_kcal/mol_RAF1`, `std_kcal/mol_RAF1`)
]

# ---------------- K13 vs RAF1 ----------------
# GTP binding allosteric mutation
K13_RAF1_GTP <- common_mutations_K13_RAF1[
  mutation_type_K13 == "GTP binding allosteric mutation" | 
    mutation_type_RAF1 == "GTP binding allosteric mutation",
  .(Pos_real, mt, `mean_kcal/mol_K13`, `std_kcal/mol_K13`,
    `mean_kcal/mol_RAF1`, `std_kcal/mol_RAF1`)
]

# Allosteric mutation
K13_RAF1_Allo <- common_mutations_K13_RAF1[
  mutation_type_K13 == "Allosteric mutation" | 
    mutation_type_RAF1 == "Allosteric mutation",
  .(Pos_real, mt, `mean_kcal/mol_K13`, `std_kcal/mol_K13`,
    `mean_kcal/mol_RAF1`, `std_kcal/mol_RAF1`)
]

# ---------------- 保存 CSV ----------------
fwrite(K27_RAF1_GTP, file.path(out_dir, "K27_RAF1_GTP_binding_allosteric.csv"))
fwrite(K27_RAF1_Allo, file.path(out_dir, "K27_RAF1_allosteric.csv"))
fwrite(K13_RAF1_GTP, file.path(out_dir, "K13_RAF1_GTP_binding_allosteric.csv"))
fwrite(K13_RAF1_Allo, file.path(out_dir, "K13_RAF1_allosteric.csv"))

# ---------------- 保存 Excel ----------------
wb <- createWorkbook()

addWorksheet(wb, "K27_RAF1_GTP")
writeData(wb, "K27_RAF1_GTP", K27_RAF1_GTP)

addWorksheet(wb, "K27_RAF1_Allo")
writeData(wb, "K27_RAF1_Allo", K27_RAF1_Allo)

addWorksheet(wb, "K13_RAF1_GTP")
writeData(wb, "K13_RAF1_GTP", K13_RAF1_GTP)

addWorksheet(wb, "K13_RAF1_Allo")
writeData(wb, "K13_RAF1_Allo", K13_RAF1_Allo)

saveWorkbook(wb, file.path(out_dir, "common_mutations_split.xlsx"), overwrite = TRUE)

cat("四个 CSV 文件和一个 Excel 文件已保存完成！\n")
