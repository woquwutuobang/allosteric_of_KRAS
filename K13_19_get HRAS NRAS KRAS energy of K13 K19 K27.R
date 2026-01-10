library(data.table)
library(krasddpcams)
library(dplyr)
library(purrr)
library(tidyr)

ddG_K13<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",assay_sele = "K13")
names(ddG_K13)
ddG_K13<-ddG_K13[,c(1:3,23:27)]

ddG_K19<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",assay_sele = "K19")
ddG_K19<-ddG_K19[,c(1:3,23:27)]
names(ddG_K19)

ddG_K27<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",assay_sele = "K27")
ddG_K27<-ddG_K27[,c(1:3,23:27)]
names(ddG_K27)

ddG_all <- bind_rows(ddG_K13, ddG_K19, ddG_K27)


mut_list <- list(
  Group1_HRAS_166 = c("H95Q","E107D","P121A","S122A","D126E","T127S","K128R","F141Y","D151E","K165Q"),
  Group2_NRAS_166 = c("T87S","E91A","H94N","H95L","E107D","S122T","Q131H","D132E","R135K","D151E","K165Q"),
  
  Group3_HRAS_188 = c("H95Q","E107D","P121A","S122A","D126E","T127S","K128R","F141Y","D151E","K165Q",
                      "E168L","K169R","M170K","S171L","K172N","D173P","G174P","K175D","K176E","K177S",
                      "K178G","K179P","K180G","S181C","K182M","T183S","K184C","C185K","V186C"),
  
  Group4_NRAS_188 = c("T87S","E91A","H94N","H95L","E107D","S122T","Q131H","D132E","R135K","D151E","K165Q",
                      "H166Y","K167R","E168M","M170K","S171L","K172N","D173S","G174S","K175D","K176D",
                      "K177G","K178T","K179Q","K180G","S181C","K182M","T183G","K184L","C185P","V186C"),
  
  Group5_HRAS_G12V_166 = c("G12V","H95Q","E107D","P121A","S122A","D126E","T127S","K128R","F141Y","D151E","K165Q"),
  Group6_NRAS_Q61H_166 = c("Q61H","T87S","E91A","H94N","H95L","E107D","S122T","Q131H","D132E","R135K","D151E","K165Q","H166Y"),
  
  Group7_HRAS_G12V_188 = c("G12V","H95Q","E107D","P121A","S122A","D126E","T127S","K128R","F141Y","D151E","K165Q",
                           "E168L","K169R","M170K","S171L","K172N","D173P","G174P","K175D","K176E","K177S",
                           "K178G","K179P","K180G","S181C","K182M","T183S","K184C","C185K","V186C"),
  
  Group8_NRAS_Q61H_188 = c("Q61H","T87S","E91A","H94N","H95L","E107D","S122T","Q131H","D132E","R135K","D151E","K165Q",
                           "H166Y","K167R","E168M","M170K","S171L","K172N","D173S","G174S","K175D","K176D",
                           "K177G","K178T","K179Q","K180G","S181C","K182M","T183G","K184L","C185P","V186C"),
  
  Group9_H95Q  = c("H95Q"),
  Group10_H95L = c("H95L"),
  Group11_E107D = c("E107D"),
  Group12_G12D = c("G12D"),
  Group13_S17N = c("S17N")
)




group_ddG_sum <- map_dfr(names(mut_list), function(g) {
  
  muts <- mut_list[[g]]
  
  ddG_all %>%
    filter(mt %in% muts) %>%
    group_by(assay) %>%
    summarise(
      group = g,
      sum_mean_kcal = sum(`mean_kcal/mol`, na.rm = TRUE),
      .groups = "drop"
    )
})


final_table <- group_ddG_sum %>%
  mutate(assay = paste0(assay, "_energy")) %>%   # K13 → K13_energy
  pivot_wider(
    names_from  = assay,
    values_from = sum_mean_kcal
  ) %>%
  rename(Assay = group)


write.csv(
  final_table,
  file = "KRAS_group_energy_summary.csv",
  row.names = FALSE
)
