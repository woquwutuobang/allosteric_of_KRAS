
## 针对0901数据
# 筛选phenotype=1的数据
phenotype1_data <- K19_pre_ob_fitness_for_0901_mochi_model[phenotype == 11, ]

# 检查Abundance和Binding列中哪些值为1
abundance_binding_cols <- c("Abundance1","Abundance2_1","Abundance2_2","Abundance3_1",
                            "Abundance3_2","Binding1_K13","Binding2_1_K13","Binding2_2_K13",
                            "Binding3_1_K13","Binding3_2_K13","Binding1_K19","Binding2_1_K19",
                            "Binding2_2_K19","Binding3_1_K19","Binding3_2_K19","Binding1_K27",
                            "Binding2_1_K27","Binding2_2_K27","Binding3_1_K27","Binding3_2_K27",
                            "Binding1_K55","Binding2_1_K55","Binding2_2_K55","Binding3_1_K55",
                            "Binding3_2_K55","Binding1_PI3","Binding2_1_PI3","Binding2_2_PI3",
                            "Binding3_1_PI3","Binding3_2_PI3","Binding1_RAF1","Binding2_1_RAF1",
                            "Binding2_2_RAF1","Binding3_1_RAF1","Binding3_2_RAF1","Binding1_RAL",
                            "Binding2_1_RAL","Binding2_2_RAL","Binding3_1_RAL","Binding3_2_RAL",
                            "Binding1_SOS","Binding2_1_SOS","Binding2_2_SOS","Binding3_1_SOS","Binding3_2_SOS")

# 找出哪些列在phenotype=1时有值为1的观测
cols_with_value_1 <- abundance_binding_cols[sapply(abundance_binding_cols, function(col) {
  any(phenotype1_data[[col]] == 1, na.rm = TRUE)
})]

# 查看结果
print(cols_with_value_1)



##### 针对0830数据

phenotype1_data <- K13_pre_ob_fitness_for_0830_mochi_model[phenotype == 6, ]

# 检查Abundance和Binding列中哪些值为1
abundance_binding_cols <- c("Abundance1","Abundance2_1","Abundance2_2","Abundance3_1","Abundance3_2",
                            "Binding1_K13","Binding2_K13","Binding3_K13",
                            "Binding1_K19","Binding2_K19", "Binding3_K19",
                            "Binding1_K27", "Binding2_K27","Binding3_K27",
                            "Binding1_K55","Binding2_K55","Binding3_K55",
                            "Binding1_PI3","Binding2_PI3","Binding3_PI3",
                            "Binding1_RAF1","Binding2_1_RAF1","Binding2_2_RAF1","Binding3_1_RAF1","Binding3_2_RAF1",
                            "Binding1_RAL","Binding2_RAL","Binding3_RAL",
                            "Binding1_SOS","Binding2_SOS","Binding3_SOS")

# 找出哪些列在phenotype=1时有值为1的观测
cols_with_value_1 <- abundance_binding_cols[sapply(abundance_binding_cols, function(col) {
  any(phenotype1_data[[col]] == 1, na.rm = TRUE)
})]

# 查看结果
print(cols_with_value_1)





# 筛选Binding3_2_RAF1=1的数据
binding_1_data <- Folding_pre_ob_fitness[Binding3_2_RAF1 == 1, ]

# 查看这些观测的phenotype值分布
if(nrow(binding_1_data) > 0) {
  phenotype_values <- unique(binding_1_data$phenotype)
  print(paste("当Binding3_2_RAF1=1时，phenotype的值为:", paste(phenotype_values, collapse = ", ")))
  
  # 统计每个phenotype的出现次数
  phenotype_counts <- table(binding_1_data$phenotype)
  print("每个phenotype的出现次数:")
  print(phenotype_counts)
} else {
  print("没有找到Binding3_2_RAF1=1的观测")
}
