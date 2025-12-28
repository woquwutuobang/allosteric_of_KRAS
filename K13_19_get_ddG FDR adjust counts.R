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



ddG<-fread("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt")


names(ddG)


# 使用 z-test 或 t-test 计算 'mean_kcal/mol' 的 p 值
ddG[, p_value := pvalue(`mean_kcal/mol`, `std_kcal/mol`, testType = "ztest")]

# 筛选 p 值小于 0.05 的结果
significant_results <- ddG[p_value < 0.05]

# 使用 Benjamini-Hochberg 方法进行 FDR 校正
ddG[, p_adjusted := p.adjust(p_value, method = "BH")]

# 筛选 FDR < 0.05 的结果
significant_results_FDR <- ddG[p_adjusted < 0.05]

# 统计 'mean_kcal/mol' 大于 0 和小于 0 的数量
count_positive <- sum(significant_results_FDR$`mean_kcal/mol` > 0)
count_negative <- sum(significant_results_FDR$`mean_kcal/mol` < 0)

# 打印结果
cat("数量统计：\n")
cat("mean_kcal/mol > 0 的数量：", count_positive, "\n")
cat("mean_kcal/mol < 0 的数量：", count_negative, "\n")


