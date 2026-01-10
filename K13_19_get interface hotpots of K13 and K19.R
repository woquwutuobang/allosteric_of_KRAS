library(data.table)
library(krasddpcams)
library(dplyr)


#####K13
ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt"
assay_sele = "K13"

ddG <- krasddpcams__read_ddG(ddG = ddG, assay_sele = assay_sele)
names(significant_ddG)


# 1. Z检验计算p-value
pvalue <- function(av, se, degreesFreedom = 5, mu = 0, testType = "ztest") {
  # Perform z test
  if (testType == "ztest") {
    zscore <- (av - mu) / se
    pval <- 2 * pnorm(abs(zscore), lower.tail = FALSE)
    return(pval)
  }
  
  # Perform t test
  if (testType == "ttest") {
    tstat <- (av - mu) / se
    pval <- 2 * pt(abs(tstat), degreesFreedom, lower = FALSE)
    return(pval)
  }
}


ddG$pvalue <- mapply(pvalue, av = ddG$`mean_kcal/mol`, se = ddG$`std_kcal/mol`)

# 2. 对p-value进行FDR校正
ddG$pvalue_fdr <- p.adjust(ddG$pvalue, method = "fdr")

# 只保留FDR < 0.05的突变
significant_ddG <- ddG %>%
  filter(pvalue_fdr < 0.05)

# 假设interface_positions已经给定
interface_positions <- c(88, 91, 87, 129, 90, 133, 94, 137, 95, 68, 136, 99, 102, 101, 107, 98)

# 筛选出界面上的突变数据
interface_ddG <- significant_ddG %>%
  filter(Pos_real %in% interface_positions)

# 计算界面所有位点的 mean_kcal/mol 和 std_kcal/mol 的平均值
mean_ddG_avg <- mean(interface_ddG$`mean_kcal/mol`, na.rm = TRUE)  # 所有界面位点的mean_kcal/mol的平均值
std_ddG_avg <- mean(interface_ddG$`std_kcal/mol`, na.rm = TRUE)  # 所有界面位点的std_kcal/mol的平均值

# 设置阈值：mean + 2 * std
threshold <- mean_ddG_avg + 2 * std_ddG_avg

# 计算每个界面位点的 mean_kcal/mol 中位数
median_by_pos <- interface_ddG %>%
  group_by(Pos_real) %>%
  summarise(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE))  # 计算每个位点的mean_kcal/mol的中位数

# 判断每个界面位点的中位数是否大于阈值
interface_ddG <- median_by_pos %>%
  mutate(hotspot = ifelse(median_ddG >= threshold, TRUE, FALSE))

# 查看热点位置
hotspot_positions <- interface_ddG %>%
  filter(hotspot == TRUE) %>%
  select(Pos_real)

# 打印热点残基位置
print(hotspot_positions)







#####K19
ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt"
assay_sele = "K19"

ddG <- krasddpcams__read_ddG(ddG = ddG, assay_sele = assay_sele)
names(significant_ddG)


# 1. Z检验计算p-value
pvalue <- function(av, se, degreesFreedom = 5, mu = 0, testType = "ztest") {
  # Perform z test
  if (testType == "ztest") {
    zscore <- (av - mu) / se
    pval <- 2 * pnorm(abs(zscore), lower.tail = FALSE)
    return(pval)
  }
  
  # Perform t test
  if (testType == "ttest") {
    tstat <- (av - mu) / se
    pval <- 2 * pt(abs(tstat), degreesFreedom, lower = FALSE)
    return(pval)
  }
}


ddG$pvalue <- mapply(pvalue, av = ddG$`mean_kcal/mol`, se = ddG$`std_kcal/mol`)

# 2. 对p-value进行FDR校正
ddG$pvalue_fdr <- p.adjust(ddG$pvalue, method = "fdr")

# 只保留FDR < 0.05的突变
significant_ddG <- ddG %>%
  filter(pvalue_fdr < 0.05)

# 假设interface_positions已经给定
interface_positions <- c(88, 91, 87, 129, 90, 133, 94, 137, 95, 68, 136, 99, 102, 101, 107, 98)

# 筛选出界面上的突变数据
interface_ddG <- significant_ddG %>%
  filter(Pos_real %in% interface_positions)

# 计算界面所有位点的 mean_kcal/mol 和 std_kcal/mol 的平均值
mean_ddG_avg <- mean(interface_ddG$`mean_kcal/mol`, na.rm = TRUE)  # 所有界面位点的mean_kcal/mol的平均值
std_ddG_avg <- mean(interface_ddG$`std_kcal/mol`, na.rm = TRUE)  # 所有界面位点的std_kcal/mol的平均值

# 设置阈值：mean + 2 * std
threshold <- mean_ddG_avg + 2 * std_ddG_avg

# 计算每个界面位点的 mean_kcal/mol 中位数
median_by_pos <- interface_ddG %>%
  group_by(Pos_real) %>%
  summarise(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE))  # 计算每个位点的mean_kcal/mol的中位数

# 判断每个界面位点的中位数是否大于阈值
interface_ddG <- median_by_pos %>%
  mutate(hotspot = ifelse(median_ddG >= threshold, TRUE, FALSE))

# 查看热点位置
hotspot_positions <- interface_ddG %>%
  filter(hotspot == TRUE) %>%
  select(Pos_real)

# 打印热点残基位置
print(hotspot_positions)





#####K27
ddG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt"
assay_sele = "K27"

ddG <- krasddpcams__read_ddG(ddG = ddG, assay_sele = assay_sele)
names(significant_ddG)


# 1. Z检验计算p-value
pvalue <- function(av, se, degreesFreedom = 5, mu = 0, testType = "ztest") {
  # Perform z test
  if (testType == "ztest") {
    zscore <- (av - mu) / se
    pval <- 2 * pnorm(abs(zscore), lower.tail = FALSE)
    return(pval)
  }
  
  # Perform t test
  if (testType == "ttest") {
    tstat <- (av - mu) / se
    pval <- 2 * pt(abs(tstat), degreesFreedom, lower = FALSE)
    return(pval)
  }
}


ddG$pvalue <- mapply(pvalue, av = ddG$`mean_kcal/mol`, se = ddG$`std_kcal/mol`)

# 2. 对p-value进行FDR校正
ddG$pvalue_fdr <- p.adjust(ddG$pvalue, method = "fdr")

# 只保留FDR < 0.05的突变
significant_ddG <- ddG %>%
  filter(pvalue_fdr < 0.05)

# 假设interface_positions已经给定
interface_positions <- c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71)

# 筛选出界面上的突变数据
interface_ddG <- significant_ddG %>%
  filter(Pos_real %in% interface_positions)

# 计算界面所有位点的 mean_kcal/mol 和 std_kcal/mol 的平均值
mean_ddG_avg <- mean(interface_ddG$`mean_kcal/mol`, na.rm = TRUE)  # 所有界面位点的mean_kcal/mol的平均值
std_ddG_avg <- mean(interface_ddG$`std_kcal/mol`, na.rm = TRUE)  # 所有界面位点的std_kcal/mol的平均值

# 设置阈值：mean + 2 * std
threshold <- mean_ddG_avg + 2 * std_ddG_avg

# 计算每个界面位点的 mean_kcal/mol 中位数
median_by_pos <- interface_ddG %>%
  group_by(Pos_real) %>%
  summarise(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE))  # 计算每个位点的mean_kcal/mol的中位数

# 判断每个界面位点的中位数是否大于阈值
interface_ddG <- median_by_pos %>%
  mutate(hotspot = ifelse(median_ddG >= threshold, TRUE, FALSE))

# 查看热点位置
hotspot_positions <- interface_ddG %>%
  filter(hotspot == TRUE) %>%
  select(Pos_real)

# 打印热点残基位置
print(hotspot_positions)




