library(krasddpcams)
library(data.table)

### 20250920

## Gb 
## for 0901 energy

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

## S1i  Non-linear relationships (global epistasis) between observed BindingPCA fitness and both free energies of binding and folding. 

krasddpcams__merge_ddGb_ob_pre_fitness_3blocks_for_0901_mochi_model <- function(prediction = prediction, 
                                                                                block1_dimsum_df = block1_dimsum_df, 
                                                                                block2_dimsum_df = block2_dimsum_df, 
                                                                                block3_dimsum_df = block3_dimsum_df,
                                                                                assay_sele = assay_sele, 
                                                                                wt_aa_input = wt_aa_input) {
  
  # 读取预测数据
  pre <- fread(prediction)
  pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
  
  # 加载所有3个block的数据
  load(block1_dimsum_df)
  block1 <- as.data.table(all_variants)
  load(block2_dimsum_df)
  block2 <- as.data.table(all_variants)
  load(block3_dimsum_df)
  block3 <- as.data.table(all_variants)
  
  # 合并所有block的数据
  data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3,
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
  
  # 创建标准化数据
  scaling_data_fitness <- data.frame(
    block1 = c(stop1_fitness, wt1_fitness),
    block2 = c(stop2_fitness, wt2_fitness),
    block3 = c(stop3_fitness, wt3_fitness)
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
  
  d2 <- params_block2$slope
  e2 <- params_block2$intercept
  d3 <- params_block3$slope
  e3 <- params_block3$intercept
  
  # 处理预测数据
  pre_nor <- pre_pos
  
  # 提取预测fitness
  extract_prediction <- function(row) {
    return(row[78 + as.numeric(row[92])])
  }
  pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = extract_prediction)
  pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
  
  # 简单的assay到phenotype映射
  assay_sele_df <- data.table(assay = c("K13", "K19", "K27", "K55", "PI3", "RAF1", "RAL", "SOS"), 
                              phenotype_base = c(6, 11, 16, 21, 26, 31, 36, 41))
  
  # 检查assay是否存在
  if (!assay_sele %in% assay_sele_df$assay) {
    stop("assay '", assay_sele, "' 不在可选的assay列表中。可选的assay有: ", 
         paste(assay_sele_df$assay, collapse = ", "))
  }
  
  # 获取基础phenotype编号
  base_pheno <- assay_sele_df[assay == assay_sele, phenotype_base]
  
  # 应用标准化转换（3个block对应3个phenotype）
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_mean_fitness, mean)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
  
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness_sigma, std)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness_sigma, std * d2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness_sigma, std * d3)]
  
  pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness, fitness)]
  pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness, fitness * d2 + e2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness, fitness * d3 + e3)]
  
  pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness_sigma, sigma)]
  pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness_sigma, sigma * d2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness_sigma, sigma * d3)]
  
  pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness, predicted_fitness)]
  pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
  pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
  
  return(pre_nor)
}

wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

## RAF1 - 使用3块版本的数据处理函数
#### S1i
RAF1_pre_ob_fitness_for_0901_mochi_model_3blocks <- krasddpcams__merge_ddGb_ob_pre_fitness_3blocks_for_0901_mochi_model(
  prediction = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/predictions/predicted_phenotypes_all.txt",
  block1_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
  block2_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
  block3_dimsum_df = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_sele = "RAF1",
  wt_aa_input = wt_aa)

## block1
Cairo::CairoPDF(file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250920/0901/model_evalution/RAF1_figureS1i1_block1_ddGf_ddGb_ob_pre_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input = RAF1_pre_ob_fitness_for_0901_mochi_model_3blocks,
                                             folding_assay = Abundance1,
                                             binding_assay = Binding1_RAF1,
                                             RT = 0.001987*(273+30),
                                             mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Binding1_RAF1.txt", colour_scheme)
dev.off()

## block2-1
Cairo::CairoPDF(file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250920/0901/model_evalution/RAF1_figureS1i2_block2_1_ddGf_ddGb_ob_pre_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input = RAF1_pre_ob_fitness_for_0901_mochi_model_3blocks,
                                             folding_assay = Abundance2_1,
                                             binding_assay = Binding2_1_RAF1,
                                             RT = 0.001987*(273+30),
                                             mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Binding2_1_RAF1.txt", colour_scheme)
dev.off()

## block3-1 (原block2-2)
Cairo::CairoPDF(file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI/Results/20250920/0901/model_evalution/RAF1_figureS1i3_block3_1_ddGf_ddGb_ob_pre_fitness.pdf")
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input = RAF1_pre_ob_fitness_for_0901_mochi_model_3blocks,
                                             folding_assay = Abundance2_2,
                                             binding_assay = Binding2_2_RAF1,
                                             RT = 0.001987*(273+30),
                                             mochi_parameters = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/weights/linears_weights_Binding2_2_RAF1.txt", colour_scheme)
dev.off()

# 注意：如果您需要为其他assay（如K13、K19）也创建3块版本，需要相应调整函数调用
# 下面是一个示例：

## K13 - 3块版本示例
# K13_pre_ob_fitness_for_0901_mochi_model_3blocks <- krasddpcams__merge_ddGb_ob_pre_fitness_3blocks_for_0901_mochi_model(
#   prediction = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20250901/task_901/predictions/predicted_phenotypes_all.txt",
#   block1_dimsum_df = "path/to/K13_block1.RData",
#   block2_dimsum_df = "path/to/K13_block2.RData",
#   block3_dimsum_df = "path/to/K13_block3.RData",
#   assay_sele = "K13",
#   wt_aa_input = wt_aa)
# 
# ## block1
# Cairo::CairoPDF(file = "path/to/K13_block1_ddGf_ddGb_ob_pre_fitness.pdf")
# krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input = K13_pre_ob_fitness_for_0901_mochi_model_3blocks,
#                                                     folding_assay = Abundance1,
#                                                     binding_assay = Binding1_K13,
#                                                     RT = 0.001987*(273+30),
#                                                     mochi_parameters = "path/to/linears_weights_Binding1_K13.txt", colour_scheme)
# dev.off()