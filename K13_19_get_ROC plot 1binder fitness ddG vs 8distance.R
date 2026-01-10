library(krasddpcams)
library(data.table)
library(ggplot2)
library(ROCR)


wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

##########################============================================
krasddpcams__normalize_growthrate_fitness_5blocks <- function(
    block1_dimsum_df,
    block2a_dimsum_df, block2b_dimsum_df,
    block3a_dimsum_df, block3b_dimsum_df
) {
  
  # 1. 读取数据
  load(block1_dimsum_df); block1 <- as.data.table(all_variants)
  load(block2a_dimsum_df); block2a <- as.data.table(all_variants)
  load(block2b_dimsum_df); block2b <- as.data.table(all_variants)
  load(block3a_dimsum_df); block3a <- as.data.table(all_variants)
  load(block3b_dimsum_df); block3b <- as.data.table(all_variants)
  
  # 2. 合并 block2 和 block3
  block2 <- rbind(block2a, block2b, fill = TRUE)
  block3 <- rbind(block3a, block3b, fill = TRUE)
  
  # 3. 合并所有 block
  data_before_nor <- rbind(
    block1 = block1,
    block2 = block2,
    block3 = block3,
    idcol = "block",
    fill = TRUE
  )
  
  # 4. 计算权重
  data_before_nor$gr_over_sigmasquared <- data_before_nor$growthrate / (data_before_nor$growthrate_sigma)^2
  data_before_nor$one_over_sigmasquared <- 1 / (data_before_nor$growthrate_sigma)^2
  data_before_nor$fitness_over_sigmasquared <- data_before_nor$fitness / (data_before_nor$sigma)^2
  data_before_nor$one_over_fitness_sigmasquared <- 1 / (data_before_nor$sigma)^2
  
  # 5. 初始化存储
  block_names <- c("block1", "block2", "block3")
  stop_means_gr <- numeric(length(block_names))
  wt_means_gr <- numeric(length(block_names))
  stop_means_fitness <- numeric(length(block_names))
  wt_means_fitness <- numeric(length(block_names))
  
  # 6. 计算加权平均值
  for(i in seq_along(block_names)) {
    blk <- block_names[i]
    dead_gr_blk <- data_before_nor[STOP == TRUE & block == blk, ]
    wt_gr_blk <- data_before_nor[WT == TRUE & block == blk, ]
    stop_means_gr[i] <- sum(dead_gr_blk$gr_over_sigmasquared, na.rm = TRUE) / sum(dead_gr_blk$one_over_sigmasquared, na.rm = TRUE)
    wt_means_gr[i] <- sum(wt_gr_blk$gr_over_sigmasquared, na.rm = TRUE) / sum(wt_gr_blk$one_over_sigmasquared, na.rm = TRUE)
    
    dead_fit_blk <- data_before_nor[STOP == TRUE & block == blk, ]
    wt_fit_blk <- data_before_nor[WT == TRUE & block == blk, ]
    stop_means_fitness[i] <- sum(dead_fit_blk$fitness_over_sigmasquared, na.rm = TRUE) / sum(dead_fit_blk$one_over_fitness_sigmasquared, na.rm = TRUE)
    wt_means_fitness[i] <- sum(wt_fit_blk$fitness_over_sigmasquared, na.rm = TRUE) / sum(wt_fit_blk$one_over_fitness_sigmasquared, na.rm = TRUE)
  }
  
  # 7. 缩放矩阵
  scaling_data_gr <- data.frame(
    block1 = c(stop_means_gr[1], wt_means_gr[1]),
    block2 = c(stop_means_gr[2], wt_means_gr[2]),
    block3 = c(stop_means_gr[3], wt_means_gr[3])
  )
  
  scaling_data_fitness <- data.frame(
    block1 = c(stop_means_fitness[1], wt_means_fitness[1]),
    block2 = c(stop_means_fitness[2], wt_means_fitness[2]),
    block3 = c(stop_means_fitness[3], wt_means_fitness[3])
  )
  
  # 8. 计算缩放系数
  coef_gr <- list(); coef_fitness <- list()
  for(blk in block_names[-1]) {
    lm_gr <- lm(as.formula(paste("block1 ~", blk)), data = scaling_data_gr)
    coef_gr[[blk]] <- list(slope = coef(lm_gr)[2], intercept = coef(lm_gr)[1])
    
    lm_fit <- lm(as.formula(paste("block1 ~", blk)), data = scaling_data_fitness)
    coef_fitness[[blk]] <- list(slope = coef(lm_fit)[2], intercept = coef(lm_fit)[1])
  }
  
  # 9. 应用归一化
  data_after_nor <- copy(data_before_nor)
  
  # 生长率
  data_after_nor[block == "block1", `:=`(nor_gr = growthrate, nor_gr_sigma = growthrate_sigma)]
  for(blk in block_names[-1]) {
    slope <- coef_gr[[blk]]$slope
    intercept <- coef_gr[[blk]]$intercept
    data_after_nor[block == blk, `:=`(nor_gr = growthrate * slope + intercept, nor_gr_sigma = growthrate_sigma * slope)]
  }
  
  # fitness
  data_after_nor[block == "block1", `:=`(nor_fitness = fitness, nor_fitness_sigma = sigma)]
  for(blk in block_names[-1]) {
    slope <- coef_fitness[[blk]]$slope
    intercept <- coef_fitness[[blk]]$intercept
    data_after_nor[block == blk, `:=`(nor_fitness = fitness * slope + intercept, nor_fitness_sigma = sigma * slope)]
  }
  
  return(data_after_nor)
}


#########################################==========================================================
krasddpcams__merge_dimsum_df9 <- function(
    merge_1, merge_2, merge_3, merge_4, merge_5, merge_6, merge_7, merge_8, merge_9
) {
  # 获取每个输入的变量名
  a1 <- as.character(substitute(merge_1))
  a2 <- as.character(substitute(merge_2))
  a3 <- as.character(substitute(merge_3))
  a4 <- as.character(substitute(merge_4))
  a5 <- as.character(substitute(merge_5))
  a6 <- as.character(substitute(merge_6))
  a7 <- as.character(substitute(merge_7))
  a8 <- as.character(substitute(merge_8))
  a9 <- as.character(substitute(merge_9))
  
  # 给每个数据表添加 assay 列
  merge_1[, `:=`(assay, a1)]
  merge_2[, `:=`(assay, a2)]
  merge_3[, `:=`(assay, a3)]
  merge_4[, `:=`(assay, a4)]
  merge_5[, `:=`(assay, a5)]
  merge_6[, `:=`(assay, a6)]
  merge_7[, `:=`(assay, a7)]
  merge_8[, `:=`(assay, a8)]
  merge_9[, `:=`(assay, a9)]
  
  # 合并
  output <- rbind(
    merge_1, merge_2, merge_3, merge_4, merge_5,
    merge_6, merge_7, merge_8, merge_9
  )
  
  return(output)
}


###################################################=======================================================
krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface <- function(
    fitness_all,
    weighted_meab_abs_ddGf,
    weighted_meab_abs_ddGb,
    anno_input,
    bind,
    bind_ligands
){
  
  
  ## -----------------------------
  ## 1. fitness
  ## -----------------------------
  abundance_fitness_mean <- fitness_all[
    Nham_aa == 1 & assay == "stab",
    .(mean_abundance_fitness = mean(nor_fitness, na.rm = TRUE)),
    by = "AA_Pos1"
  ]
  setnames(abundance_fitness_mean, "AA_Pos1", "Pos_real")
  
  binding_fitness_mean <- fitness_all[
    Nham_aa == 1 & assay == bind,
    .(mean_binding_fitness = mean(nor_fitness, na.rm = TRUE)),
    by = "AA_Pos1"
  ]
  setnames(binding_fitness_mean, "AA_Pos1", "Pos_real")
  
  ## -----------------------------
  ## 2. ddG
  ## -----------------------------
  ddGf <- weighted_meab_abs_ddGf[, 1:2]
  setnames(ddGf, c("Pos_real", "weighted_mean_ddGf"))
  
  ddGb <- weighted_meab_abs_ddGb[, 1:2]
  setnames(ddGb, c("Pos_real", "weighted_mean_ddGb"))
  
  ## -----------------------------
  ## 3. anno
  ## -----------------------------
  anno <- copy(anno_input)
  
  ## -----------------------------
  ## 4. merge once（重要）
  ## -----------------------------
  base_dt <- merge(abundance_fitness_mean, binding_fitness_mean, by = "Pos_real")
  base_dt <- merge(base_dt, ddGf, by = "Pos_real")
  base_dt <- merge(base_dt, ddGb, by = "Pos_real")
  
  metric_names_plot <- c("mean_binding_fitness", "weighted_mean_ddGb")
  
  perf_list <- list()
  
  ## =============================
  ## 5. ligand 循环
  ## =============================
  for(bind_ligand in bind_ligands){
    
    anno_tmp <- copy(anno)
    anno_tmp[, boostDM := 0L]
    anno_tmp[get(bind_ligand) < 5, boostDM := 1L]
    anno_tmp <- anno_tmp[, .(Pos, boostDM)]
    
    dt <- merge(
      base_dt,
      anno_tmp,
      by.x = "Pos_real",
      by.y = "Pos"
    )
    
    dt <- dt[order(Pos_real)][!duplicated(Pos_real)]
    
    for(metric in metric_names_plot){
      
      dt[, plot_metric := get(metric)]
      
      roc_df <- data.frame(
        predictions = dt$plot_metric,
        labels = dt$boostDM
      )
      
      ## 防止全 0 / 全 1
      if(length(unique(roc_df$labels)) < 2){
        message("Skipping ", bind_ligand, " - ", metric, " (single class)")
        next
      }
      
      roc_df$predictions_lm <-
        lm(labels ~ predictions, data = roc_df)$fitted.values
      
      pred <- ROCR::prediction(roc_df$predictions_lm, roc_df$labels)
      perf <- ROCR::performance(pred, "tpr", "fpr")
      auc  <- round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
      
      perf_list[[paste(bind_ligand, metric, sep = "_")]] <-
        data.table(
          FPR = perf@x.values[[1]],
          TPR = perf@y.values[[1]],
          ligand = bind_ligand,
          measure = metric,
          auc = auc
        )
    }
  }
  
  plot_dt <- rbindlist(perf_list)
  
  ## -----------------------------
  ## 6. AUC label
  ## -----------------------------
  plot_dt[, measure := factor(
    measure,
    levels = c("weighted_mean_ddGb", "mean_binding_fitness")
  )]
  
  auc_dt <- plot_dt[
    , .(auc = unique(auc)),
    by = .(ligand, measure)
  ]
  
  auc_dt[, `:=`(
    FPR = 0.05,
    TPR = ifelse(
      measure == "weighted_mean_ddGb", 0.30, 0.60
    )
  )]
  
  plot_cols <- c(
    weighted_mean_ddGb = "black",
    mean_binding_fitness = "#F4270C"
  )
  
  ## -----------------------------
  ## 7. plot
  ## -----------------------------
  p <- p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = FPR, y = TPR, color = measure)
  ) +
    ggplot2::geom_line(linewidth = 0.35) +
    ggplot2::geom_abline(
      linetype = 2,
      linewidth = 0.35,
      color = "grey50"
    ) +
    
    ## ⭐ 关键：facet
    ggplot2::facet_wrap(~ ligand, ncol = 4) +
    
    ## AUC text（每个 facet 一份）
    ggplot2::geom_text(
      data = auc_dt,
      ggplot2::aes(
        x = FPR,
        y = TPR,
        label = paste0("AUC = ", auc),
        color = measure
      ),
      size = 3,
      hjust = 0,
      vjust = 0,
      inherit.aes = FALSE,
      show.legend = FALSE
    )+
    
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    
    ggplot2::scale_colour_manual(
      values = plot_cols,
      name = NULL,
      labels = c(
        "weighted mean of absolute ddGb",
        "mean binding fitness"
      )
    ) +
    
    ggplot2::theme_classic(base_size = 7) +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "vertical",
      legend.text = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(size = 7)
    ) +
    
    ggplot2::coord_fixed()
  
  return(p)
}







stability_nor_df <- krasddpcams__normalize_growthrate_fitness_5blocks(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                                      block2a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                                      block2b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                                                                      block3a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                                                      block3b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")


RAF_nor_df <- krasddpcams__normalize_growthrate_fitness_5blocks(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                                                block2a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                                                block2b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
                                                                block3a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
                                                                block3b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")


RALGDS_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                                                           block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                                                           block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData")


PI3KCG_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                                           block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                                           block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData")


SOS1_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                                                         block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                                                         block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData")


K55_nor_df <- krasddpcams__normalize_growthrate_fitness_5blocks(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                                                                block2a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                                                                block2b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K55_block2_Q20_rbg2_filter1_fitness_replicates.RData",
                                                                block3a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
                                                                block3b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K55_block3_Q20_rbg_1_filter1_fitness_replicates.RData")



K27_nor_df <- krasddpcams__normalize_growthrate_fitness_5blocks(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                                                                block2a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                                                                block2b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K27_block2_Q20_rbg3_filter1_fitness_replicates.RData",
                                                                block3a_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
                                                                block3b_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K27_block3_Q20_rbg3_filter1_fitness_replicates.RData")


K13_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData",
                                                        block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                                                        block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")



K19_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates_cleaned.RData",
                                                        block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates_cleaned.RData",
                                                        block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")




stab <- stability_nor_df
RAF <- RAF_nor_df
RALGDS <- RALGDS_nor_df
PI3KCG <- PI3KCG_nor_df
SOS1 <- SOS1_nor_df
K55 <- K55_nor_df
K27 <- K27_nor_df
K13 <- K13_nor_df
K19 <- K19_nor_df


all_data <- krasddpcams__merge_dimsum_df9(stab,RAF,RALGDS,PI3KCG,SOS1,K55,K27,K13,K19)
all_data_pos<-krasddpcams__pos_id(all_data,wt_aa)




bind_ligands <- c(
  "scHAmin_ligand_K13",
  "scHAmin_ligand_K19",
  "scHAmin_ligand_K27",
  "scHAmin_ligand_K55",
  "scHAmin_ligand_SOS1",
  "scHAmin_ligand_PI3KCG",
  "scHAmin_ligand_RALGDS",
  "scHAmin_ligand_RAF1"
)





############## K13
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_K13<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",
                                                                  assay_sele = "K13")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_K13,
  anno_input = anno_input,
  bind = "K13",
  bind_ligands = bind_ligands
)

p + ggtitle("K13 binding interface")




############## K19
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_K19<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt",
                                                                  assay_sele = "K19")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
#names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_K19,
  anno_input = anno_input,
  bind = "K19",
  bind_ligands = bind_ligands
)

p + ggtitle("K19 binding interface")





############## RAF1
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_RAF1<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",
                                                                  assay_sele = "RAF1")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_RAF1,
  anno_input = anno_input,
  bind = "RAF",
  bind_ligands = bind_ligands
)

p + ggtitle("RAF1 binding interface")





############## RALGDS
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_RALGDS<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAL.txt",
                                                                   assay_sele = "RALGDS")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
#names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_RALGDS,
  anno_input = anno_input,
  bind = "RALGDS",
  bind_ligands = bind_ligands
)

p + ggtitle("RALGDS binding interface")






############## PI3KCG
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_PI3KCG<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_PI3.txt",
                                                                     assay_sele = "PI3KCG")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
#names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_PI3KCG,
  anno_input = anno_input,
  bind = "PI3KCG",
  bind_ligands = bind_ligands
)

p + ggtitle("PI3KCG binding interface")




############## K55
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_K55<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K55.txt",
                                                                     assay_sele = "K55")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
#names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_K55,
  anno_input = anno_input,
  bind = "K55",
  bind_ligands = bind_ligands
)

p + ggtitle("K55 binding interface")







############## K27
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_K27<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt",
                                                                  assay_sele = "K27")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
#names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_K27,
  anno_input = anno_input,
  bind = "K27",
  bind_ligands = bind_ligands
)

p + ggtitle("K27 binding interface")







############## SOS1
weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Folding.txt",
                                                                   assay_sele = "folding")
weighted_mean_abs_ddG_SOS1<-krasddpcams__get_weighted_mean_abs_ddG(ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_SOS.txt",
                                                                  assay_sele = "SOS1")


anno_input<-fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")
#names(anno_input)

p <- krasddpcams__plot_multiROCAUC_ddG_fitness_binding_interface(
  fitness_all = all_data_pos,
  weighted_meab_abs_ddGf = weighted_mean_abs_ddG_fold,
  weighted_meab_abs_ddGb = weighted_mean_abs_ddG_SOS1,
  anno_input = anno_input,
  bind = "SOS1",
  bind_ligands = bind_ligands
)

p + ggtitle("SOS1 binding interface")

