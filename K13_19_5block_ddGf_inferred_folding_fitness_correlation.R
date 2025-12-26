#### 20250917
library(data.table)
library(krasddpcams)
## for 0901 mochi model
## Gf


## S1f 2D density plots showing non-linear relationships (global epistasis) between observed AbundancePCA fitness and changes in free energy of folding. 


krasddpcams__merge_ddGf_fitness_5blocks_for_0901_mochi_model <- function(prediction = prediction, folding_ddG = folding_ddG, 
                                                                         block1_dimsum_df = block1_dimsum_df, 
                                                                         block2_dimsum_df = block2_dimsum_df, 
                                                                         block3_dimsum_df = block3_dimsum_df,
                                                                         block4_dimsum_df = block4_dimsum_df,
                                                                         block5_dimsum_df = block5_dimsum_df,
                                                                         wt_aa_input = wt_aa_input) {
  
  # 读取预测和折叠自由能数据
  pre <- fread(prediction)
  folding_ddG <- fread(folding_ddG)
  pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
  
  # 加载所有5个block的数据
  load(block1_dimsum_df)
  block1 <- as.data.table(all_variants)
  load(block2_dimsum_df)
  block2 <- as.data.table(all_variants)
  load(block3_dimsum_df)
  block3 <- as.data.table(all_variants)
  load(block4_dimsum_df)
  block4 <- as.data.table(all_variants)
  load(block5_dimsum_df)
  block5 <- as.data.table(all_variants)
  
  # 合并所有block的数据
  data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3,
                           block4 = block4, block5 = block5, 
                           idcol = "block", fill = TRUE)
  
  # 计算加权fitness所需的中间量
  data_before_nor$fitness_over_sigmasquared <- data_before_nor$fitness/(data_before_nor$sigma)^2
  data_before_nor$one_over_fitness_sigmasquared <- 1/(data_before_nor$sigma)^2
  
  # 计算每个block的终止密码子fitness（加权平均）
  calculate_stop_fitness <- function(block_name) {
    dead_fitness <- data_before_nor[STOP == TRUE & block == block_name, ]
    stop_fitness <- sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(dead_fitness$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(stop_fitness)
  }
  
  stop1_fitness <- calculate_stop_fitness("block1")
  stop2_fitness <- calculate_stop_fitness("block2")
  stop3_fitness <- calculate_stop_fitness("block3")
  stop4_fitness <- calculate_stop_fitness("block4")
  stop5_fitness <- calculate_stop_fitness("block5")
  
  # 计算每个block的野生型fitness（加权平均）
  calculate_wt_fitness <- function(block_name) {
    wt_fitness_block <- data_before_nor[WT == TRUE & block == block_name, ]
    wt_fitness <- sum(wt_fitness_block$fitness_over_sigmasquared, na.rm = TRUE) /
      sum(wt_fitness_block$one_over_fitness_sigmasquared, na.rm = TRUE)
    return(wt_fitness)
  }
  
  wt1_fitness <- calculate_wt_fitness("block1")
  wt2_fitness <- calculate_wt_fitness("block2")
  wt3_fitness <- calculate_wt_fitness("block3")
  wt4_fitness <- calculate_wt_fitness("block4")
  wt5_fitness <- calculate_wt_fitness("block5")
  
  # 创建标准化数据
  scaling_data_fitness <- data.frame(
    block1 = c(stop1_fitness, wt1_fitness),
    block2 = c(stop2_fitness, wt2_fitness),
    block3 = c(stop3_fitness, wt3_fitness),
    block4 = c(stop4_fitness, wt4_fitness),
    block5 = c(stop5_fitness, wt5_fitness)
  )
  
  # 以block1为基准，计算其他block的转换参数
  calculate_scaling_params <- function(target_block) {
    lm_model <- lm(formula = block1 ~ get(target_block), data = scaling_data_fitness)
    return(list(
      slope = lm_model$coefficients[[2]],
      intercept = lm_model$coefficients[[1]]
    ))
  }
  
  params_block2 <- calculate_scaling_params("block2")
  params_block3 <- calculate_scaling_params("block3")
  params_block4 <- calculate_scaling_params("block4")
  params_block5 <- calculate_scaling_params("block5")
  
  d2 <- params_block2$slope
  e2 <- params_block2$intercept
  d3 <- params_block3$slope
  e3 <- params_block3$intercept
  d4 <- params_block4$slope
  e4 <- params_block4$intercept
  d5 <- params_block5$slope
  e5 <- params_block5$intercept
  
  # 处理预测数据
  pre_nor <- pre_pos
  
  # 提取预测fitness（需要根据实际数据结构调整索引）
  extract_prediction <- function(row) {
    return(row[78 + as.numeric(row[92])])
  }
  pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = extract_prediction)
  pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
  
  # 提取加性性状（需要根据实际数据结构调整索引）
  extract_additive_trait0 <- function(row) {
    return(row[92 + as.numeric(row[92]) * 2 - 1])
  }
  extract_additive_trait1 <- function(row) {
    return(row[92 + as.numeric(row[92]) * 2])
  }
  
  pre_nor$additive_trait0 <- apply(pre_nor, MARGIN = 1, FUN = extract_additive_trait0)
  pre_nor$additive_trait0 <- as.numeric(pre_nor$additive_trait0)
  pre_nor$additive_trait1 <- apply(pre_nor, MARGIN = 1, FUN = extract_additive_trait1)
  pre_nor$additive_trait1 <- as.numeric(pre_nor$additive_trait1)
  pre_nor[, `:=`(additive_trait, additive_trait0 + additive_trait1)]
  
  # 应用标准化转换（5个block）
  pre_nor[phenotype == 1, `:=`(pre_nor_mean_fitness, mean)]
  pre_nor[phenotype == 2, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
  pre_nor[phenotype == 3, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
  pre_nor[phenotype == 4, `:=`(pre_nor_mean_fitness, mean * d4 + e4)]
  pre_nor[phenotype == 5, `:=`(pre_nor_mean_fitness, mean * d5 + e5)]
  
  pre_nor[phenotype == 1, `:=`(pre_nor_fitness_sigma, std)]
  pre_nor[phenotype == 2, `:=`(pre_nor_fitness_sigma, std * d2)]
  pre_nor[phenotype == 3, `:=`(pre_nor_fitness_sigma, std * d3)]
  pre_nor[phenotype == 4, `:=`(pre_nor_fitness_sigma, std * d4)]
  pre_nor[phenotype == 5, `:=`(pre_nor_fitness_sigma, std * d5)]
  
  pre_nor[phenotype == 1, `:=`(ob_nor_fitness, fitness)]
  pre_nor[phenotype == 2, `:=`(ob_nor_fitness, fitness * d2 + e2)]
  pre_nor[phenotype == 3, `:=`(ob_nor_fitness, fitness * d3 + e3)]
  pre_nor[phenotype == 4, `:=`(ob_nor_fitness, fitness * d4 + e4)]
  pre_nor[phenotype == 5, `:=`(ob_nor_fitness, fitness * d5 + e5)]
  
  pre_nor[phenotype == 1, `:=`(ob_nor_fitness_sigma, sigma)]
  pre_nor[phenotype == 2, `:=`(ob_nor_fitness_sigma, sigma * d2)]
  pre_nor[phenotype == 3, `:=`(ob_nor_fitness_sigma, sigma * d3)]
  pre_nor[phenotype == 4, `:=`(ob_nor_fitness_sigma, sigma * d4)]
  pre_nor[phenotype == 5, `:=`(ob_nor_fitness_sigma, sigma * d5)]
  
  pre_nor[phenotype == 1, `:=`(pre_nor_fitness, predicted_fitness)]
  pre_nor[phenotype == 2, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
  pre_nor[phenotype == 3, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
  pre_nor[phenotype == 4, `:=`(pre_nor_fitness, predicted_fitness * d4 + e4)]
  pre_nor[phenotype == 5, `:=`(pre_nor_fitness, predicted_fitness * d5 + e5)]
  
  return(pre_nor)
}






wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

#### S1f 2D density plots showing non-linear relationships (global epistasis) between observed AbundancePCA fitness and changes in free energy of folding

Folding_pre_ob_fitness_for_0901_mochi_model<-krasddpcams__merge_ddGf_fitness_5blocks_for_0901_mochi_model(prediction="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/predictions/predicted_phenotypes_all.txt",
                                                                                                          folding_ddG="C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/weights_Folding.txt",
                                                                                                          block1_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                                                                          block2_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                                                                          block3_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                                                                                                          block4_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                                                                                          block5_dimsum_df="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                                                                                                          wt_aa_input=wt_aa)





## block1
krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness_for_0901_mochi_model,
                                 fold_n=1,
                                 mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Abundance1.txt",
                                 phenotypen = 1,RT=0.001987*(273+30),bin_input = 50)
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250917/0901/model_evalution/figureS1f1_block1_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

## block2-1
krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness_for_0901_mochi_model,
                                 fold_n=1,
                                 mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Abundance2_1.txt",
                                 phenotypen = 2,RT=0.001987*(273+30),bin_input = 50)
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250917/0901/model_evalution/figureS1f2_block2_1_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

## block2-2
krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness_for_0901_mochi_model,
                                 fold_n=1,
                                 mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Abundance2_2.txt",
                                 phenotypen = 3,RT=0.001987*(273+30),bin_input = 50)
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250917/0901/model_evalution/figureS1f3_block2_2_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")


## block3-1
krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness_for_0901_mochi_model,
                                 fold_n=1,
                                 mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Abundance3_1.txt",
                                 phenotypen = 4,RT=0.001987*(273+30),bin_input = 50)
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250917/0901/model_evalution/figureS1f3_block3_1_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")



## block3-2
krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness_for_0901_mochi_model,
                                 fold_n=1,
                                 mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Abundance3_2.txt",
                                 phenotypen = 5,RT=0.001987*(273+30),bin_input = 50)
ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250917/0901/model_evalution/figureS1f3_block3_2_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")







