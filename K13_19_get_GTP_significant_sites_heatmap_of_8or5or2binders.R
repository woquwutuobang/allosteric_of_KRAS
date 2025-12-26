library(krasddpcams)
library(data.table)
library(ggplot2)

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



krasddpcams__plot_ddG_all_allosteric_heatmap <- function (ddG1, assay1, ddG2, assay2, ddG3, assay3, ddG4, assay4, 
                                                          ddG5, assay5, ddG6, assay6, ddG7, assay7, ddG8, assay8, anno, wt_aa, colour_scheme) 
{
  # 指定关心的位点
  Pos_list <- c(15, 16, 18, 29, 32, 57, 116, 117, 145)  # 你关注的突变位点
  
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  
  # 合并突变信息和热图工具数据
  heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, 
                                                            "")), each = 20), Pos_real = rep(2:188, each = 20), mt_codon = unlist(aa_list))
  
  # 读取 ddG 数据并合并
  ddG1 <- krasddpcams__read_ddG(ddG1, assay1)
  input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG2 <- krasddpcams__read_ddG(ddG2, assay2)
  input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG3 <- krasddpcams__read_ddG(ddG3, assay3)
  input3_heatmap <- merge(ddG3, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG4 <- krasddpcams__read_ddG(ddG4, assay4)
  input4_heatmap <- merge(ddG4, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG5 <- krasddpcams__read_ddG(ddG5, assay5)
  input5_heatmap <- merge(ddG5, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG6 <- krasddpcams__read_ddG(ddG6, assay6)
  input6_heatmap <- merge(ddG6, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG7 <- krasddpcams__read_ddG(ddG7, assay7)
  input7_heatmap <- merge(ddG7, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG8 <- krasddpcams__read_ddG(ddG8, assay8)
  input8_heatmap <- merge(ddG8, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  # 合并所有输入数据
  ddG <- rbind(input1_heatmap, input2_heatmap, input3_heatmap, 
               input4_heatmap, input5_heatmap, input6_heatmap, input7_heatmap, input8_heatmap)
  
  # 对于wt和mt相同的突变，设置mean_kcal/mol 为 0
  ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
  
  # 合并注释数据
  output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", 
                  all = T)
  
  # 筛选关心的位点
  output <- output[Pos_real %in% Pos_list, ]
  
  # 设置突变氨基酸的顺序
  output <- within(output, mt_codon <- factor(mt_codon, levels = c("D", 
                                                                   "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", "P", 
                                                                   "A", "V", "I", "L", "M", "F", "W", "Y")))
  
  # 设置assay的顺序
  output <- within(output, assay <- factor(assay, levels = c("RAF1", 
                                                             "RALGDS", "PI3KCG", "SOS1", "K27", "K55", "K13", "K19")))
  
  # 计算wt_codon和Pos_real的组合，作为位置的唯一标识
  output[, `:=`(wtcodon_pos, paste0(wt_codon, Pos_real))]
  
  # 排序wtcodon_pos
  output <- within(output, wtcodon_pos <- factor(wtcodon_pos, 
                                                 levels = unique(output[order(Pos_real), wtcodon_pos])))
  
  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = output, mapping = ggplot2::aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)) + 
    ggplot2::geom_text(data = output[Pos_real > 1 & wt_codon == mt_codon, ], 
                       mapping = ggplot2::aes(x = assay, y = mt_codon, label = "-"), size = 5 * 5/14) + 
    ggplot2::scale_fill_gradient2(low = colour_scheme[["blue"]], 
                                  mid = "gray", high = colour_scheme[["red"]], na.value = "white") + 
    ggplot2::facet_wrap(~wtcodon_pos, nrow = 2) + 
    ggplot2::ylab("Mutant AA") + 
    ggplot2::xlab("Binding partners") + 
    ggplot2::labs(fill = NULL) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 10, family = "Arial"), 
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(), 
      legend.position = "bottom", 
      legend.text = ggplot2::element_text(size = 10), 
      strip.text.x = ggplot2::element_text(size = 10), 
      strip.background = ggplot2::element_rect(colour = "white", fill = "white"), 
      axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
      axis.text.y = ggplot2::element_text(family = "Arial", size = 10, vjust = 0.5, hjust = 1, 
                                          margin = ggplot2::margin(0, 5, 0, 0)),  # Ensure proper spacing for Y axis
      panel.spacing.y = ggplot2::unit(3, "mm"), 
      panel.spacing.x = ggplot2::unit(1, "mm"), 
      legend.key.height = ggplot2::unit(3.1, "mm"), 
      legend.key.width = ggplot2::unit(2, "mm"), 
      legend.key.size = ggplot2::unit(1, "mm")) + 
    ggplot2::coord_fixed()
  
}

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

krasddpcams__plot_ddG_all_allosteric_heatmap(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt", assay1 = "RAF1",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAL.txt", assay2 = "RALGDS",
  ddG3 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_PI3.txt", assay3 = "PI3KCG",
  ddG4 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_SOS.txt", assay4 = "SOS1",
  ddG5 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K55.txt", assay5 = "K55",
  ddG6 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt", assay6 = "K27",
  ddG7 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt", assay7 = "K13",
  ddG8 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt", assay8 = "K19",
  anno = anno, wt_aa = wt_aa,
  colour_scheme)





ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f4/20251216/GTP_significantsites_heatmap_of_8binders.pdf", device = cairo_pdf,height = 8, width=6)











########################## 5binders RAF1/K55/K27/K13/K19

library(krasddpcams)
library(data.table)
library(ggplot2)

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


krasddpcams__plot_ddG_all_allosteric_heatmap <- function (ddG1, assay1, ddG2, assay2, ddG3, assay3, ddG4, assay4, 
                                                          ddG5, assay5,anno, wt_aa, colour_scheme,legend_limits=c(-1.6,2.8)) 
{
  
  Pos_list <- c(15, 16, 18, 29, 32, 57, 116, 117, 145)  
  
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  
  # 合并突变信息和热图工具数据
  heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, 
                                                            "")), each = 20), Pos_real = rep(2:188, each = 20), mt_codon = unlist(aa_list))
  
  # 读取 ddG 数据并合并
  ddG1 <- krasddpcams__read_ddG(ddG1, assay1)
  input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG2 <- krasddpcams__read_ddG(ddG2, assay2)
  input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG3 <- krasddpcams__read_ddG(ddG3, assay3)
  input3_heatmap <- merge(ddG3, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG4 <- krasddpcams__read_ddG(ddG4, assay4)
  input4_heatmap <- merge(ddG4, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG5 <- krasddpcams__read_ddG(ddG5, assay5)
  input5_heatmap <- merge(ddG5, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  # 合并所有输入数据
  ddG <- rbind(input1_heatmap, input2_heatmap, input3_heatmap, 
               input4_heatmap, input5_heatmap)
  
  # 对于wt和mt相同的突变，设置mean_kcal/mol 为 0
  ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
  
  # 合并注释数据
  output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", 
                  all = T)
  
  # 筛选关心的位点
  output <- output[Pos_real %in% Pos_list, ]
  
  # 设置突变氨基酸的顺序
  output <- within(output, mt_codon <- factor(mt_codon, levels = c("D", 
                                                                   "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", "P", 
                                                                   "A", "V", "I", "L", "M", "F", "W", "Y")))
  
  # 设置assay的顺序
  output <- within(output, assay <- factor(assay, levels = c("RAF1", "K55","K27" ,"K13", "K19")))
  
  # 计算wt_codon和Pos_real的组合，作为位置的唯一标识
  output[, `:=`(wtcodon_pos, paste0(wt_codon, Pos_real))]
  
  # 排序wtcodon_pos
  output <- within(output, wtcodon_pos <- factor(wtcodon_pos, 
                                                 levels = unique(output[order(Pos_real), wtcodon_pos])))
  
  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = output, mapping = ggplot2::aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)) + 
    ggplot2::geom_text(data = output[Pos_real > 1 & wt_codon == mt_codon, ], 
                       mapping = ggplot2::aes(x = assay, y = mt_codon, label = "-"), size = 5 * 5/14) + 
    ggplot2::scale_fill_gradient2(limits = legend_limits,low = colour_scheme[["blue"]], 
                                  mid = "gray", high = colour_scheme[["red"]], na.value = "white") + 
    ggplot2::facet_wrap(~wtcodon_pos, nrow = 2) + 
    ggplot2::ylab("Mutant AA") + 
    ggplot2::xlab("Binding partners") + 
    ggplot2::labs(fill = NULL) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 10, family = "Arial"), 
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(), 
      legend.position = "bottom", 
      legend.text = ggplot2::element_text(size = 10), 
      strip.text.x = ggplot2::element_text(size = 10), 
      strip.background = ggplot2::element_rect(colour = "white", fill = "white"), 
      axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
      axis.text.y = ggplot2::element_text(family = "Arial", size = 10, vjust = 0.5, hjust = 1, 
                                          margin = ggplot2::margin(0, 5, 0, 0)),  # Ensure proper spacing for Y axis
      panel.spacing.y = ggplot2::unit(3, "mm"), 
      panel.spacing.x = ggplot2::unit(1, "mm"), 
      legend.key.height = ggplot2::unit(3.1, "mm"), 
      legend.key.width = ggplot2::unit(2, "mm"), 
      legend.key.size = ggplot2::unit(1, "mm")) + 
    ggplot2::coord_fixed()
}

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

krasddpcams__plot_ddG_all_allosteric_heatmap(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt", assay1 = "RAF1",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K55.txt", assay2 = "K55",
  ddG3 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt", assay3 = "K27",
  ddG4 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt", assay4 = "K13",
  ddG5 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt", assay5 = "K19",
  anno = anno, wt_aa = wt_aa,
  colour_scheme,legend_limits = c(-1.6,2.8))




ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f4/20251216/GTP_sites_significantsites_heatmap_of_5binders.pdf", device = cairo_pdf,height = 8, width=6)





########################## 3binders RAF1/K27/K13

library(krasddpcams)
library(data.table)
library(ggplot2)

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


krasddpcams__plot_ddG_all_allosteric_heatmap <- function (ddG1, assay1, ddG2, assay2, ddG3, assay3,anno, wt_aa, colour_scheme,legend_limits=c(-1.6,2.8)) 
{
  
  Pos_list <- c(15, 16, 18, 29, 32, 57, 116, 117, 145)  
  
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  
  # 合并突变信息和热图工具数据
  heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, 
                                                            "")), each = 20), Pos_real = rep(2:188, each = 20), mt_codon = unlist(aa_list))
  
  # 读取 ddG 数据并合并
  ddG1 <- krasddpcams__read_ddG(ddG1, assay1)
  input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG2 <- krasddpcams__read_ddG(ddG2, assay2)
  input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG3 <- krasddpcams__read_ddG(ddG3, assay3)
  input3_heatmap <- merge(ddG3, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  
  # 合并所有输入数据
  ddG <- rbind(input1_heatmap, input2_heatmap, input3_heatmap)
  
  # 对于wt和mt相同的突变，设置mean_kcal/mol 为 0
  ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
  
  # 合并注释数据
  output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", 
                  all = T)
  
  # 筛选关心的位点
  output <- output[Pos_real %in% Pos_list, ]
  
  # 设置突变氨基酸的顺序
  output <- within(output, mt_codon <- factor(mt_codon, levels = c("D", 
                                                                   "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", "P", 
                                                                   "A", "V", "I", "L", "M", "F", "W", "Y")))
  
  # 设置assay的顺序
  output <- within(output, assay <- factor(assay, levels = c("RAF1", "K27" ,"K13")))
  
  # 计算wt_codon和Pos_real的组合，作为位置的唯一标识
  output[, `:=`(wtcodon_pos, paste0(wt_codon, Pos_real))]
  
  # 排序wtcodon_pos
  output <- within(output, wtcodon_pos <- factor(wtcodon_pos, 
                                                 levels = unique(output[order(Pos_real), wtcodon_pos])))
  
  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = output, mapping = ggplot2::aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)) + 
    ggplot2::geom_text(data = output[Pos_real > 1 & wt_codon == mt_codon, ], 
                       mapping = ggplot2::aes(x = assay, y = mt_codon, label = "-"), size = 5 * 5/14) + 
    ggplot2::scale_fill_gradient2(limits = legend_limits,low = colour_scheme[["blue"]], 
                                  mid = "gray", high = colour_scheme[["red"]], na.value = "white") + 
    ggplot2::facet_wrap(~wtcodon_pos, nrow = 2) + 
    ggplot2::ylab("Mutant AA") + 
    ggplot2::xlab("Binding partners") + 
    ggplot2::labs(fill = NULL) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 10, family = "Arial"), 
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(), 
      legend.position = "bottom", 
      legend.text = ggplot2::element_text(size = 10), 
      strip.text.x = ggplot2::element_text(size = 10), 
      strip.background = ggplot2::element_rect(colour = "white", fill = "white"), 
      axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
      axis.text.y = ggplot2::element_text(family = "Arial", size = 10, vjust = 0.5, hjust = 1, 
                                          margin = ggplot2::margin(0, 5, 0, 0)),  # Ensure proper spacing for Y axis
      panel.spacing.y = ggplot2::unit(3, "mm"), 
      panel.spacing.x = ggplot2::unit(1, "mm"), 
      legend.key.height = ggplot2::unit(3.1, "mm"), 
      legend.key.width = ggplot2::unit(2, "mm"), 
      legend.key.size = ggplot2::unit(1, "mm")) + 
    ggplot2::coord_fixed()
}

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

krasddpcams__plot_ddG_all_allosteric_heatmap(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt", assay1 = "RAF1",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt", assay2 = "K27",
  ddG3 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt", assay3 = "K13",
  anno = anno, wt_aa = wt_aa,
  colour_scheme,legend_limits = c(-1.6,2.8))




ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/探究中/mutation invoved site check on structure/results_data/20251216/GTP_sites_significantsites_heatmap_of_3binders_K27_RAF1_K13.pdf", device = cairo_pdf,height = 8, width=5)





########################## 2binders RAF1/K13

library(krasddpcams)
library(data.table)
library(ggplot2)

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


krasddpcams__plot_ddG_all_allosteric_heatmap <- function (ddG1, assay1, ddG2, assay2,anno, wt_aa, colour_scheme,legend_limits) 
{
  
  Pos_list <- c(15, 16, 18, 29, 32, 57, 116, 117, 145)  
  
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  
  # 合并突变信息和热图工具数据
  heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, 
                                                            "")), each = 20), Pos_real = rep(2:188, each = 20), mt_codon = unlist(aa_list))
  
  # 读取 ddG 数据并合并
  ddG1 <- krasddpcams__read_ddG(ddG1, assay1)
  input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG2 <- krasddpcams__read_ddG(ddG2, assay2)
  input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  
  # 合并所有输入数据
  ddG <- rbind(input1_heatmap, input2_heatmap)
  
  # 对于wt和mt相同的突变，设置mean_kcal/mol 为 0
  ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
  
  # 合并注释数据
  output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", 
                  all = T)
  
  # 筛选关心的位点
  output <- output[Pos_real %in% Pos_list, ]
  
  # 设置突变氨基酸的顺序
  output <- within(output, mt_codon <- factor(mt_codon, levels = c("D", 
                                                                   "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", "P", 
                                                                   "A", "V", "I", "L", "M", "F", "W", "Y")))
  
  # 设置assay的顺序
  output <- within(output, assay <- factor(assay, levels = c("RAF1","K13")))
  
  # 计算wt_codon和Pos_real的组合，作为位置的唯一标识
  output[, `:=`(wtcodon_pos, paste0(wt_codon, Pos_real))]
  
  # 排序wtcodon_pos
  output <- within(output, wtcodon_pos <- factor(wtcodon_pos, 
                                                 levels = unique(output[order(Pos_real), wtcodon_pos])))
  
  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = output, mapping = ggplot2::aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)) + 
    ggplot2::geom_text(data = output[Pos_real > 1 & wt_codon == mt_codon, ], 
                       mapping = ggplot2::aes(x = assay, y = mt_codon, label = "-"), size = 5 * 5/14) + 
    ggplot2::scale_fill_gradient2(limits = legend_limits,low = colour_scheme[["blue"]], 
                                  mid = "gray", high = colour_scheme[["red"]], na.value = "white") + 
    ggplot2::facet_wrap(~wtcodon_pos, nrow = 2) + 
    ggplot2::ylab("Mutant AA") + 
    ggplot2::xlab("Binding partners") + 
    ggplot2::labs(fill = NULL) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 10, family = "Arial"), 
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(), 
      legend.position = "bottom", 
      legend.text = ggplot2::element_text(size = 10), 
      strip.text.x = ggplot2::element_text(size = 10), 
      strip.background = ggplot2::element_rect(colour = "white", fill = "white"), 
      axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
      axis.text.y = ggplot2::element_text(family = "Arial", size = 10, vjust = 0.5, hjust = 1, 
                                          margin = ggplot2::margin(0, 5, 0, 0)),  # Ensure proper spacing for Y axis
      panel.spacing.y = ggplot2::unit(3, "mm"), 
      panel.spacing.x = ggplot2::unit(1, "mm"), 
      legend.key.height = ggplot2::unit(3.1, "mm"), 
      legend.key.width = ggplot2::unit(2, "mm"), 
      legend.key.size = ggplot2::unit(1, "mm")) + 
    ggplot2::coord_fixed()
}

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

krasddpcams__plot_ddG_all_allosteric_heatmap(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt", assay1 = "RAF1",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt", assay2 = "K13",
  anno = anno, wt_aa = wt_aa,
  colour_scheme,legend_limits = c(-1.6,2.8))




ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f4/20251216/GTP_sites_significantsites_heatmap_of_2binders_RAF1_K13.pdf", device = cairo_pdf,height = 15, width=20)







########################## 2binders RAF1/K19

library(krasddpcams)
library(data.table)
library(ggplot2)

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


krasddpcams__plot_ddG_all_allosteric_heatmap <- function (ddG1, assay1, ddG2, assay2,anno, wt_aa, colour_scheme,legend_limits) 
{
  
  Pos_list <- c(15, 16, 18, 29, 32, 57, 116, 117, 145)  
  
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", 
                                     "")))
  
  # 合并突变信息和热图工具数据
  heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, 
                                                            "")), each = 20), Pos_real = rep(2:188, each = 20), mt_codon = unlist(aa_list))
  
  # 读取 ddG 数据并合并
  ddG1 <- krasddpcams__read_ddG(ddG1, assay1)
  input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  ddG2 <- krasddpcams__read_ddG(ddG2, assay2)
  input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", 
                                                     "Pos_real", "mt_codon"), all = T)
  
  
  # 合并所有输入数据
  ddG <- rbind(input1_heatmap, input2_heatmap)
  
  # 对于wt和mt相同的突变，设置mean_kcal/mol 为 0
  ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
  
  # 合并注释数据
  output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", 
                  all = T)
  
  # 筛选关心的位点
  output <- output[Pos_real %in% Pos_list, ]
  
  # 设置突变氨基酸的顺序
  output <- within(output, mt_codon <- factor(mt_codon, levels = c("D", 
                                                                   "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", "P", 
                                                                   "A", "V", "I", "L", "M", "F", "W", "Y")))
  
  # 设置assay的顺序
  output <- within(output, assay <- factor(assay, levels = c("RAF1","K19")))
  
  # 计算wt_codon和Pos_real的组合，作为位置的唯一标识
  output[, `:=`(wtcodon_pos, paste0(wt_codon, Pos_real))]
  
  # 排序wtcodon_pos
  output <- within(output, wtcodon_pos <- factor(wtcodon_pos, 
                                                 levels = unique(output[order(Pos_real), wtcodon_pos])))
  
  ggplot2::ggplot() + 
    ggplot2::geom_tile(data = output, mapping = ggplot2::aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)) + 
    ggplot2::geom_text(data = output[Pos_real > 1 & wt_codon == mt_codon, ], 
                       mapping = ggplot2::aes(x = assay, y = mt_codon, label = "-"), size = 5 * 5/14) + 
    ggplot2::scale_fill_gradient2(limits = legend_limits,low = colour_scheme[["blue"]], 
                                  mid = "gray", high = colour_scheme[["red"]], na.value = "white") + 
    ggplot2::facet_wrap(~wtcodon_pos, nrow = 2) + 
    ggplot2::ylab("Mutant AA") + 
    ggplot2::xlab("Binding partners") + 
    ggplot2::labs(fill = NULL) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 10, family = "Arial"), 
      axis.ticks.x = ggplot2::element_blank(), 
      axis.ticks.y = ggplot2::element_blank(), 
      legend.position = "bottom", 
      legend.text = ggplot2::element_text(size = 10), 
      strip.text.x = ggplot2::element_text(size = 10), 
      strip.background = ggplot2::element_rect(colour = "white", fill = "white"), 
      axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
      axis.text.y = ggplot2::element_text(family = "Arial", size = 10, vjust = 0.5, hjust = 1, 
                                          margin = ggplot2::margin(0, 5, 0, 0)),  # Ensure proper spacing for Y axis
      panel.spacing.y = ggplot2::unit(3, "mm"), 
      panel.spacing.x = ggplot2::unit(1, "mm"), 
      legend.key.height = ggplot2::unit(3.1, "mm"), 
      legend.key.width = ggplot2::unit(2, "mm"), 
      legend.key.size = ggplot2::unit(1, "mm")) + 
    ggplot2::coord_fixed()
}

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
anno<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

krasddpcams__plot_ddG_all_allosteric_heatmap(
  ddG1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt", assay1 = "RAF1",
  ddG2 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt", assay2 = "K19",
  anno = anno, wt_aa = wt_aa,
  colour_scheme,legend_limits = c(-1.6,2.8))




ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f4/20251216/GTP_sites_significantsites_heatmap_of_2binders_RAF1_K19.pdf", device = cairo_pdf,height = 15, width=20)
