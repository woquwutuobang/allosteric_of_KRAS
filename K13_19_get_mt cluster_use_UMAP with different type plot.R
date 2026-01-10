library(umap)     # UMAP
library(ggplot2)  # 可视化
library(data.table)
library(krasddpcams)
library(tidyr)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

#################=========================================================
## version 1 with folding to cluster
################==========================================================

ddG_RAF1<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",assay_sele = "RAF1")
ddG_RALGDS<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAL.txt",assay_sele = "RALGDS")
ddG_PI3KCG<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_PI3.txt",assay_sele = "PI3KCG")
ddG_SOS1<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_SOS.txt",assay_sele = "SOS1")
ddG_K55<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K55.txt",assay_sele = "K55")
ddG_K27<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",assay_sele = "K27")
ddG_K13<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",assay_sele = "K13")
ddG_K19<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",assay_sele = "K19")
ddG_folding<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Folding.txt",assay_sele = "folding")

ddG_list <- list(
  RAF1 = ddG_RAF1,
  RALGDS = ddG_RALGDS,
  PI3KCG = ddG_PI3KCG,
  SOS1 = ddG_SOS1,
  K55 = ddG_K55,
  K27 = ddG_K27,
  K13 = ddG_K13,
  K19 = ddG_K19,
  folding = ddG_folding
)

# 指定 assay 顺序
assay_order <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19","folding")

################################################################################
## 1. 数据预处理
################################################################################

# 合并所有数据，保留mt列
all_ddG_data <- rbindlist(lapply(names(ddG_list), function(assay_name) {
  dt <- ddG_list[[assay_name]]
  if ("mt" %in% names(dt)) {
    result <- dt[, .(assay = assay_name, Pos_real, mt, ddG = `mean_kcal/mol`)]
    return(result)
  } else {
    cat(paste("警告:", assay_name, "数据没有mt列\n"))
    return(NULL)
  }
}), fill = TRUE)

# 检查数据
cat("总数据行数:", nrow(all_ddG_data), "\n")
cat("唯一突变:", length(unique(all_ddG_data$mt)), "\n")

# 删除mt列为NA的行
all_ddG_data <- all_ddG_data[!is.na(mt)]

# 删除格式不正确的mt
all_ddG_data <- all_ddG_data[grepl("^[A-Z][0-9]+[A-Z]$", mt)]

# 将数据转换为宽格式
wide_ddG_data <- dcast(all_ddG_data, mt + Pos_real ~ assay, value.var = "ddG", fun.aggregate = mean)

cat("\n转换后数据行数:", nrow(wide_ddG_data), "\n")

# 检查缺失值
assay_columns <- setdiff(names(wide_ddG_data), c("mt", "Pos_real"))
missing_counts <- colSums(is.na(wide_ddG_data[, ..assay_columns]))
cat("\n各assay缺失值数量:\n")
print(missing_counts)

# 删除有缺失值的行
wide_ddG_data_complete <- wide_ddG_data[complete.cases(wide_ddG_data[, ..assay_columns])]
cat("\n完全数据行数:", nrow(wide_ddG_data_complete), "\n")

################################################################################
## 2. 解析突变信息
################################################################################

# 解析突变信息的函数
parse_mutation <- function(mt_string) {
  if (is.na(mt_string) || mt_string == "" || nchar(mt_string) < 3) {
    return(list(wt = NA, pos = NA, mut = NA, valid = FALSE))
  }
  
  if (!grepl("^[A-Z][0-9]+[A-Z]$", mt_string)) {
    return(list(wt = NA, pos = NA, mut = NA, valid = FALSE))
  }
  
  wt_aa <- substr(mt_string, 1, 1)
  mut_aa <- substr(mt_string, nchar(mt_string), nchar(mt_string))
  position_str <- substr(mt_string, 2, nchar(mt_string) - 1)
  pos <- as.numeric(position_str)
  
  return(list(wt = wt_aa, pos = pos, mut = mut_aa, valid = TRUE))
}

# 解析所有突变
parsed_mutations <- lapply(wide_ddG_data_complete$mt, parse_mutation)

# 只保留有效突变
valid_indices <- which(sapply(parsed_mutations, function(x) x$valid))
wide_ddG_data_complete <- wide_ddG_data_complete[valid_indices]
parsed_mutations <- parsed_mutations[valid_indices]

# 添加解析后的信息
wide_ddG_data_complete[, wt_aa := sapply(parsed_mutations, function(x) x$wt)]
wide_ddG_data_complete[, mut_aa := sapply(parsed_mutations, function(x) x$mut)]
wide_ddG_data_complete[, parsed_pos := sapply(parsed_mutations, function(x) x$pos)]

cat("\n有效突变数量:", nrow(wide_ddG_data_complete), "\n")

################################################################################
## 3. 氨基酸类型分类
################################################################################

# 氨基酸类型分类
amino_acid_categories <- list(
  aromatic = c("F", "W", "Y"),
  aliphatic = c("A", "V", "I", "L", "M"),
  polar_uncharged = c("S", "T", "N", "Q", "C"),
  positive = c("K", "R", "H"),
  negative = c("D", "E"),
  special = c("G", "P")
)

# 创建氨基酸到类别的映射
get_category <- function(aa) {
  if (is.na(aa)) return(NA_character_)
  
  for (category_name in names(amino_acid_categories)) {
    if (aa %in% amino_acid_categories[[category_name]]) {
      return(category_name)
    }
  }
  return("other")
}

# 添加氨基酸类别
wide_ddG_data_complete[, wt_category := sapply(wt_aa, get_category)]
wide_ddG_data_complete[, mut_category := sapply(mut_aa, get_category)]

# 删除类别为NA的行
wide_ddG_data_complete <- wide_ddG_data_complete[!is.na(wt_category) & !is.na(mut_category)]

# 创建突变类型标签
wide_ddG_data_complete[, mutation_type := paste(wt_category, "to", mut_category)]
wide_ddG_data_complete[, mutation_simple := paste(wt_aa, "→", mut_aa)]

cat("\n最终数据行数:", nrow(wide_ddG_data_complete), "\n")

################################################################################
## 4. UMAP 分析 (替换t-SNE部分)
################################################################################

# 准备UMAP输入数据
umap_input_data <- as.matrix(wide_ddG_data_complete[, ..assay_columns])

cat("\nUMAP输入数据维度:", dim(umap_input_data), "\n")

# 进行 UMAP 降维
set.seed(123)
umap_result <- umap(umap_input_data)

# 将 UMAP 结果合并到数据中
wide_ddG_data_complete[, umap_1 := umap_result$layout[, 1]]
wide_ddG_data_complete[, umap_2 := umap_result$layout[, 2]]

################################################################################
## 5. 功能区域和二级结构定义
################################################################################

# 功能区域定义
K13_binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 
                                101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_binding_interface_site <- c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71) 
GTP_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

# 二级结构区域
rects_sheet <- data.frame(
  xstart = c(3, 38, 51, 77, 109, 139),
  xend = c(9, 44, 57, 84, 115, 143),
  col = c("β1", "β2", "β3", "β4", "β5", "β6")
)

rects_helix <- data.frame(
  xstart = c(15, 67, 87, 127, 148),
  xend = c(24, 73, 104, 136, 166),
  col = c("α1", "α2", "α3", "α4", "α5")
)

# 创建二级结构标签的函数
create_secondary_structure_labels <- function(pos) {
  for (i in 1:nrow(rects_sheet)) {
    if (pos >= rects_sheet$xstart[i] && pos <= rects_sheet$xend[i]) {
      return(rects_sheet$col[i])
    }
  }
  
  for (i in 1:nrow(rects_helix)) {
    if (pos >= rects_helix$xstart[i] && pos <= rects_helix$xend[i]) {
      return(rects_helix$col[i])
    }
  }
  
  return("Other")
}

# 添加功能区域和二级结构标签
wide_ddG_data_complete[, functional_group := "Other"]
wide_ddG_data_complete[Pos_real %in% K13_binding_interface_site, functional_group := "K13 Interface"]
wide_ddG_data_complete[Pos_real %in% RAF1_binding_interface_site, functional_group := "RAF1 Interface"]
wide_ddG_data_complete[Pos_real %in% GTP_pocket, functional_group := "GTP Pocket"]
wide_ddG_data_complete[Pos_real %in% c(10:17, 25:40, 58:76), functional_group := "Functional Loop"]

wide_ddG_data_complete[, secondary_structure := sapply(Pos_real, create_secondary_structure_labels)]
secondary_structure_order <- c(paste0("β", 1:6), paste0("α", 1:5), "Other")
wide_ddG_data_complete[, secondary_structure := factor(secondary_structure, levels = secondary_structure_order)]

################################################################################
## 6. 绘制6个UMAP图 (替换t-SNE为UMAP)
################################################################################

########## 图1：按位置编号着色
p1 <- ggplot(wide_ddG_data_complete, aes(x = umap_1, y = umap_2, color = Pos_real)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "UMAP Plot (colored by Position)", 
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
  scale_color_gradient(low = "grey80", high = "#C68EFD",
                       name = "Position")

print(p1)

########## 图2：按ddG值着色，分面展示各assay
# 转换为长格式
long_ddG_data <- melt(wide_ddG_data_complete,
                      id.vars = c("mt", "Pos_real", "umap_1", "umap_2", 
                                  "wt_aa", "mut_aa", "wt_category", "mut_category",
                                  "mutation_type", "mutation_simple",
                                  "functional_group", "secondary_structure"),
                      measure.vars = assay_columns,
                      variable.name = "assay",
                      value.name = "ddG")

long_ddG_data[, assay := factor(assay, levels = assay_order)]

p2 <- ggplot(long_ddG_data, aes(x = umap_1, y = umap_2, color = ddG)) +
  geom_point(size = 1, alpha = 0.7) +
  labs(title = "UMAP Plot (colored by ddG value)", 
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.position = "right") +
  scale_color_gradient2(low = "#1B38A6", mid = "grey", high = "#F4270C", 
                        midpoint = 0,
                        name = "ddG (kcal/mol)",
                        limits = c(-1, 2)) +
  facet_wrap(~assay, nrow = 3)

print(p2)

########## 图3：按功能区域着色
functional_colors <- c(
  "K13 Interface" = "#1B38A6",    
  "RAF1 Interface" = "#F4270C",   
  "GTP Pocket" = "#F4AD0C",       
  "Functional Loop" = "#F1DD10",  
  "Other" = "grey80"             
)

p3 <- ggplot(wide_ddG_data_complete, aes(x = umap_1, y = umap_2, color = functional_group)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "UMAP Plot (colored by Functional Region)", 
       x = "UMAP 1", y = "UMAP 2",
       color = "Functional Region") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = functional_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p3)

########## 图4：按二级结构着色
secondary_colors <- c(
  "β1" = "#8DD3C7", "β2" = "#FFFFB3", "β3" = "#BEBADA", 
  "β4" = "#FB8072", "β5" = "#80B1D3", "β6" = "#FDB462",
  "α1" = "#B3DE69", "α2" = "#FCCDE5", "α3" = "#D9D9D9",
  "α4" = "#BC80BD", "α5" = "#CCEBC5",
  "Other" = "#999999"
)

p4 <- ggplot(wide_ddG_data_complete, aes(x = umap_1, y = umap_2, color = secondary_structure)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "UMAP Plot (colored by Secondary Structure)", 
       x = "UMAP 1", y = "UMAP 2",
       color = "Secondary Structure") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = secondary_colors) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))

print(p4)

########## 图5：按突变型氨基酸类别着色
mut_category_colors <- c(
  "aromatic" = "#E41A1C",        # 红色
  "aliphatic" = "#377EB8",       # 蓝色
  "polar_uncharged" = "#4DAF4A", # 绿色
  "positive" = "#984EA3",        # 紫色
  "negative" = "#FF7F00",        # 橙色
  "special" = "#FFFF33",         # 黄色
  "other" = "#A65628"            # 棕色
)

p5 <- ggplot(wide_ddG_data_complete, aes(x = umap_1, y = umap_2, color = mut_category)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "UMAP Plot (colored by Mutant Amino Acid Category)", 
       x = "UMAP 1", y = "UMAP 2",
       color = "Mutant AA Category") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = mut_category_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p5)

########## 图6：按突变类型（野生型→突变型）着色
mutation_type_counts <- sort(table(wide_ddG_data_complete$mutation_type), decreasing = TRUE)
top_mutation_types <- names(head(mutation_type_counts, 10))
wide_ddG_data_complete[, mutation_type_plot := ifelse(mutation_type %in% top_mutation_types, 
                                                      mutation_type, "Other")]

# 为常见突变类型分配颜色
mutation_type_colors <- c(
  "aliphatic to aliphatic" = "#66C2A5",
  "aliphatic to polar_uncharged" = "#FC8D62",
  "aliphatic to special" = "#8DA0CB",
  "polar_uncharged to aliphatic" = "#E78AC3",
  "polar_uncharged to polar_uncharged" = "#A6D854",
  "aliphatic to aromatic" = "#FFD92F",
  "polar_uncharged to special" = "#E5C494",
  "special to aliphatic" = "#B3B3B3",
  "aliphatic to positive" = "#BC80BD",
  "polar_uncharged to aromatic" = "#CCEBC5",
  "Other" = "#999999"
)

p6 <- ggplot(wide_ddG_data_complete, aes(x = umap_1, y = umap_2, color = mutation_type_plot)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "UMAP Plot (colored by Mutation Type)", 
       x = "UMAP 1", y = "UMAP 2",
       color = "Mutation Type\n(WT → Mutant)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = mutation_type_colors) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1)))

print(p6)

