library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(krasddpcams)

# 读取数据和注释
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# 定义函数，进行K13和RAF1突变分析
classify_mutations_comprehensive <- function(input, anno, assay_sele, wt_aa) {
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  ddG[, `:=`(Pos_real, Pos_ref + 1)]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
  
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  if (is.character(anno)) {
    anno_data <- fread(anno)
  } else {
    anno_data <- as.data.table(anno)
  }
  
  data_plot <- merge(weighted_mean_ddG, anno_data, by = "Pos", all = TRUE)
  
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  gxpmg_col <- paste0("GXPMG_scHAmin_ligand_", assay_sele)
  
  data_plot[get(scHAmin_col) < 5, binding_type := "binding site"]
  data_plot[, binding_type_gtp_included := binding_type]
  data_plot[get("GXPMG_scHAmin_ligand_RAF1") < 5, binding_type_gtp_included := "GTP binding site"]
  
  reg_threshold <- data_plot[binding_type == "binding site", 
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE), 
                             .SDcols = c("mean", "sigma")]
  
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", site_type := "GTP binding interface site"]
  
  data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
  data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
  data_plot_mutation[, mutation_type := "Reminder"]
  
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH") < 0.05 & (abs(mean) - reg_threshold) > 0]
  
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
  
  data_plot_mutation[, mutation_type := factor(mutation_type,
                                               levels = c("Orthosteric site huge differences",
                                                          "Orthosteric site small differences",
                                                          "GTP binding allosteric mutation",
                                                          "GTP binding other mutation",
                                                          "Allosteric mutation",
                                                          "Other mutation"))]
  
  result_cols <- c("Pos", "Pos_real", "wt_codon", "mt_codon", "mt", 
                   "mean_kcal/mol", "std_kcal/mol", "mean", "std",
                   "site_type", "allosteric_mutation", "mutation_type",
                   "assay", "reg_threshold", "is_significant")
  
  result_cols <- intersect(result_cols, names(data_plot_mutation))
  
  return(data_plot_mutation[, ..result_cols])
}


# 分析K13和RAF1的突变
mutation_K13 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "K13", 
  wt_aa
)

mutation_RAF1 <- classify_mutations_comprehensive(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv",
  assay_sele = "RAF1", 
  wt_aa
)



pvalue <- function(
    av,
    se,
    degreesFreedom = 5,
    mu = 0,
    testType = "ztest"
){
  
  #Perform z test
  if(testType=="ztest"){
    zscore <- (av - mu)/se
    pval <- 2*pnorm(abs(zscore), lower.tail = FALSE)
    return(pval)
  }
  
  #Perform t test
  if(testType=="ttest"){
    tstat <- (av - mu)/se
    pval <- 2*pt(abs(tstat), degreesFreedom, lower = FALSE)
    return(pval)
  }
  
}

mutation_K13[, pval := pvalue(mean, std)]
mutation_K13[, fdr := p.adjust(pval, method = "BH")]

mutation_RAF1[, pval := pvalue(mean, std)]
mutation_RAF1[, fdr := p.adjust(pval, method = "BH")]

define_chemotype <- function(aa) {
  if (aa %in% c("F","W","Y")) return("Aromatic")
  if (aa %in% c("A","V","I","L","M")) return("Aliphatic")
  if (aa %in% c("S","T","N","Q","C")) return("Polar uncharged")
  if (aa %in% c("K","R","H")) return("Positive")
  if (aa %in% c("D","E")) return("Negative")
  if (aa %in% c("G","P")) return("Special")
  return(NA)
}

mutation_K13[, chemotype := sapply(mt_codon, define_chemotype)]
mutation_RAF1[, chemotype := sapply(mt_codon, define_chemotype)]


heat_K13 <- mutation_K13[fdr < 0.05, .(
  mean_ddG = mean(`mean_kcal/mol`, na.rm = TRUE)
), by = .(Pos, chemotype)]
heat_K13[, assay := "K13"]

heat_RAF1 <- mutation_RAF1[fdr < 0.05, .(
  mean_ddG = mean(`mean_kcal/mol`, na.rm = TRUE)
), by = .(Pos, chemotype)]
heat_RAF1[, assay := "RAF1"]


heat_sig <- rbind(heat_K13, heat_RAF1)

chemotype_levels <- c("Aromatic","Aliphatic","Polar uncharged","Positive","Negative","Special")
assay_levels <- c("K13", "RAF1")
heat_sig[, chemotype := factor(chemotype, levels = chemotype_levels)]
heat_sig[, assay := factor(assay, levels = assay_levels)]
heat_sig[, x_var := paste(chemotype, assay, sep = " | ")]

x_levels <- as.vector(
  unlist(
    lapply(chemotype_levels, function(ct) paste(ct, assay_levels, sep = " | "))
  )
)
heat_sig[, x_var := factor(x_var, levels = x_levels)]



# GTP pocket 位点
GTP_pocket <- c(
  12, 13, 14, 15, 16, 17, 18,
  28, 29, 30, 32, 34, 35,
  57, 60, 61,
  116, 117, 119, 120,
  145, 146, 147
)

# 星号数据
star_df <- data.table(Pos = GTP_pocket)

p_sig <- ggplot(
  heat_sig,
  aes(y = x_var, x = factor(Pos), fill = mean_ddG)
) +
  geom_tile(color = "grey80") +
  
  # ⭐ 星号放在 X 轴旁边
  geom_text(
    data = star_df,
    aes(x = factor(Pos), y = 0),   # 0 表示在最下方或最左侧
    label = "★",
    color = "#FFD700",
    size = 4,
    inherit.aes = FALSE
  ) +
  
  scale_fill_gradient2(
    low = "#1B38A6",
    mid = "white",
    high = "#F4270C",
    midpoint = 0,
    name = expression(Delta*Delta*G)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  coord_fixed() +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(5.5, 5.5, 5.5, 40)
  ) +
  labs(
    y = "Mutation chemotype | Assay",
    x = "KRAS position",
    title = "All mutations (FDR-corrected)"
  )

p_sig



ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251226/all_muts_heatmap aa type.pdf", 
                device = cairo_pdf, height = 6, width = 20, dpi = 300)


