library(wlab.block)
library(dplyr)
library(ggplot2)
library(patchwork)

wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

####Abundance

## block2
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData")



## block3
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")




### RAF1

## block2
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData")




## block3
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")



### K13

## block1
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K13_block1_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block1_Q20_rbg_filter2_20250829_fitness_replicates.RData")


## block2
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData")

## block3

load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")



### K19
## block1

load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block1_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block1_Q20_rbg_filter2_20250829_fitness_replicates.RData")

## block2

load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData")

## block3

load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/rbg_filter2_Q20_RData_report_20250829/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

#names(all_variants)
#names(synonymous)

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")






#### 20251121

## K13-block1 0710

load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K13_block1_0710/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData")  

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData")





## K19-block1 1027+before

load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K19_block1/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates.RData")  

all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates.RData")






## K19-block2 1027
#aa_complement("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K19_block2/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates.RData",wt_aa)
load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K19_block2/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates.RData")  
#load("C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_20251027/K19_block1/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates_cleaned.RData")  



all_variants <- rbind(all_variants, synonymous)
all_variants <- all_variants %>% filter(sigma < 0.5)


save(all_variants, 
     file = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates.RData")

