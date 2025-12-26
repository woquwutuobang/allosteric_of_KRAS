library(ggplot2)
library(data.table)
library(dplyr)
library(wlab.block)
library(krasddpcams)
library(ggpubr)
library(rlang) 
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"



# ================================
# Function: Scatter plot (ΔΔGf vs ΔΔGb)
# ================================
krasddpcams__plot_scatter_ddGb_ddGf <- function(ddG1, assay1, ddG2, assay2,
                                                anno, binder = c("K13", "K19","RAF1"),
                                                colour_scheme,
                                                xlim = c(-0.7, 2.7),
                                                ylim = c(-1.7, 2.7)) {
  binder <- match.arg(binder)
  
  # ---- 1. 读取并合并 ddG 数据 ----
  ddG1 <- krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1)
  ddG2 <- krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2)
  all_ddG <- rbind(ddG1, ddG2)
  all_ddG_dc <- dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ assay, value.var = "mean_kcal/mol")
  all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", by.y = "Pos")
  
  # ---- 2. 定义 binding interface ----
  dist_col <- paste0("scHAmin_ligand_", binder)
  all_ddG_dc_anno[, `:=`(binding_type, "others")]
  all_ddG_dc_anno[get(dist_col) <= 5, `:=`(binding_type, "binding interface")]
  
  # ---- 3. 绘图 ----
  p <- ggplot() +
    geom_point(
      data = all_ddG_dc_anno[Pos_real > 1 & binding_type == "others", ],
      aes(x = !!sym(assay1), y = !!sym(binder)),
      color = "black", alpha = 0.6, size = 1.5
    ) +
    geom_point(
      data = all_ddG_dc_anno[Pos_real > 1 & binding_type == "binding interface"],
      aes(x = !!sym(assay1), y = !!sym(binder)),
      color = "#F4270C", alpha = 0.6, size = 2) +
    
    # ---- 固定坐标范围 ----
  scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    
    theme_classic() +
    labs(
      x = "Folding ΔΔG (kcal/mol)",
      y = paste0(binder, " Binding ΔΔG (kcal/mol)")
    ) +
    theme(
      text = element_text(size = 8),
      legend.position = "right",
      legend.text = element_text(size = 8),
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1, size = 8, colour = "black"
      ),
      axis.text.y = element_text(
        size = 8, colour = "black",
        vjust = 0.5, hjust = 0.5,
        margin = margin(0, 0, 0, 0, "mm")
      ),
      axis.title = element_text(size = 8),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(3.1, "mm"),
      legend.key.size = unit(1, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    coord_fixed(ratio = 1, expand = TRUE, clip = "on")
  
  return(p)
}



anno <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv")

# K13
p_k13 <- krasddpcams__plot_scatter_ddGb_ddGf(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Folding.txt",
  assay1 = "folding",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K13.txt",
  assay2 = "K13",
  anno = anno,
  binder = "K13"
)
p_k13

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/figure2a_scatter_ddGb_ddGf_K13.pdf", p_k13, device = cairo_pdf, height = 4, width = 4)

# K19
p_k19 <- krasddpcams__plot_scatter_ddGb_ddGf(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Folding.txt",
  assay1 = "folding",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_K19.txt",
  assay2 = "K19",
  anno = anno,
  binder = "K19"
)
p_k19
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/figure2a_scatter_ddGb_ddGf_K19.pdf", p_k19, device = cairo_pdf, height = 4, width = 4)





# RAF1
p_RAF1 <- krasddpcams__plot_scatter_ddGb_ddGf(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Folding.txt",
  assay1 = "folding",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Binding_RAF1.txt",
  assay2 = "RAF1",
  anno = anno,
  binder = "RAF1"
)
p_RAF1
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/figure1/20251029/figure2a_scatter_ddGb_ddGf_RAF1.pdf", p_RAF1, device = cairo_pdf, height = 4, width = 4)




