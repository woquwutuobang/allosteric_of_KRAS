library(wlab.block)
library(data.table)
library(krasddpcams)

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

wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
####Abundance
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)

names(nor_fit_single)


# 使用 z-test 或 t-test 计算 'mean_kcal/mol' 的 p 值
nor_fit_single[, p_value := pvalue(`nor_fitness`, `nor_fitness_sigma`, testType = "ztest")]

# 筛选 p 值小于 0.05 的结果
significant_results <- nor_fit_single[p_value < 0.05]

# 使用 Benjamini-Hochberg 方法进行 FDR 校正
nor_fit_single[, p_adjusted := p.adjust(p_value, method = "BH")]

# 筛选 FDR < 0.05 的结果
significant_results_FDR <- nor_fit_single[p_adjusted < 0.05]

# 统计 'mean_kcal/mol' 大于 0 和小于 0 的数量
count_positive <- sum(significant_results_FDR$`nor_fitness` > 0)
count_negative <- sum(significant_results_FDR$`nor_fitness` < 0)

# 打印结果
cat("数量统计：\n")
cat("mean_kcal/mol > 0 的数量：", count_positive, "\n")
cat("mean_kcal/mol < 0 的数量：", count_negative, "\n")








####RAF1
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/final_fitness_for_plot/20250525_RDatas/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)









####K13
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block1_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)





####K19
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_cleaned_RData_20250829_2/K19_block1_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
