library(Rtsne)    # t-SNE
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

ddG_list <- list(
  RAF1 = ddG_RAF1,
  RALGDS = ddG_RALGDS,
  PI3KCG = ddG_PI3KCG,
  SOS1 = ddG_SOS1,
  K55 = ddG_K55,
  K27 = ddG_K27,
  K13 = ddG_K13,
  K19 = ddG_K19
)

# 指定 assay 顺序
assay_order <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

# 计算每个 Pos_real 在每个 assay 的 mean_kcal/mol 中位数
median_dt <- rbindlist(lapply(names(ddG_list), function(assay_name) {
  dt <- ddG_list[[assay_name]]
  dt[, .(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE)), by = Pos_real][, assay := assay_name]
}))

# 将 assay 转为 factor，指定顺序
median_dt[, assay := factor(assay, levels = assay_order)]

pivot_dt <- dcast(median_dt, Pos_real ~ assay, value.var = "median_ddG")

# 检查缺失值
sum(is.na(pivot_dt))
# 删除含有缺失值的行
pivot_dt <- pivot_dt[complete.cases(pivot_dt), ]

################################################################################
## t-SNE 分析
################################################################################

# 进行 t-SNE 降维
set.seed(123)
tsne_result <- Rtsne(as.matrix(pivot_dt[, -1, with = FALSE]), dims = 2, pca = TRUE)

# 将 t-SNE 结果合并到数据中
pivot_dt[, tsne_1 := tsne_result$Y[, 1]]
pivot_dt[, tsne_2 := tsne_result$Y[, 2]]

########## 1. t-SNE 图（按位置编号着色）
p_tsne_pos <- ggplot(pivot_dt, aes(x = tsne_1, y = tsne_2, color = Pos_real)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Plot (colored by Position)", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12)) +
  scale_color_gradient(low = "grey80", high = "#C68EFD",
                       name = "Position")

print(p_tsne_pos)

########## 2. t-SNE 图（按ddG值着色，分面展示各assay）

# 转换为长格式用于分面展示
pivot_long <- melt(pivot_dt,
                   id.vars = c("Pos_real", "tsne_1", "tsne_2"),
                   measure.vars = assay_order,
                   variable.name = "assay",
                   value.name = "median_ddG")

# 确保assay顺序
pivot_long[, assay := factor(assay, levels = assay_order)]

p_tsne_ddG <- ggplot(pivot_long, aes(x = tsne_1, y = tsne_2, color = median_ddG)) +
  geom_point(size = 1.5, alpha = 0.8) +
  labs(title = "t-SNE Plot (colored by ddG value)", 
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10),
        legend.position = "right") +
  scale_color_gradient2(low = "#1B38A6", mid = "grey", high = "#F4270C", 
                        midpoint = 0,
                        name = "ddG (kcal/mol)",
                        breaks = seq(-1, 2, 1),
                        limits = c(-1, 2)) +
  facet_wrap(~assay, nrow = 2)

print(p_tsne_ddG)




########## 3. t-SNE 图（按功能区域着色）

################################################################################
## 定义功能区域
################################################################################
K13_binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 
                                101, 102, 105, 106, 107, 129, 133, 136, 137, 138)
RAF1_binding_interface_site <- c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71) 
GTP_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)



# 1). 结合界面和GTP口袋的着色图

# 定义颜色方案
binding_site_colors <- c(
  "K13 Interface" = "#1B38A6",    # 深蓝色
  "RAF1 Interface" = "#F4270C",   # 红色
  "GTP Pocket" = "#F4AD0C",       # 橙色
  "Other" = "grey80"             # 灰色
)

# 创建结合界面和GTP口袋的分类标签
pivot_dt[, binding_group := "Other"]
pivot_dt[Pos_real %in% K13_binding_interface_site, binding_group := "K13 Interface"]
pivot_dt[Pos_real %in% RAF1_binding_interface_site, binding_group := "RAF1 Interface"]
pivot_dt[Pos_real %in% GTP_pocket, binding_group := "GTP Pocket"]

# 检查分布
cat("结合界面和GTP口袋位置分布统计:\n")
print(table(pivot_dt$binding_group))

# 绘制图1：结合界面和GTP口袋
p_binding_sites <- ggplot(pivot_dt, aes(x = tsne_1, y = tsne_2, color = binding_group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE Plot (Binding Interfaces & GTP Pocket)", 
       x = "t-SNE 1", y = "t-SNE 2",
       color = "Binding Site / GTP Pocket") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = binding_site_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_binding_sites)

################################################################################

# 2). 功能环区的着色图

# 定义功能环区的颜色方案
functional_loop_colors <- c(
  "P-loop (10-17)" = "#2E8B57",      # 海绿色
  "switch I (25-40)" = "#8A2BE2",    # 蓝紫色
  "switch II (58-76)" = "#FF4500",   # 橙红色
  "Other" = "grey80"                 # 灰色
)

# 创建功能环区的分类标签
pivot_dt[, loop_group := "Other"]

# 分别标记每个功能环区
pivot_dt[Pos_real %in% 10:17, loop_group := "P-loop (10-17)"]
pivot_dt[Pos_real %in% 25:40, loop_group := "switch I (25-40)"]
pivot_dt[Pos_real %in% 58:76, loop_group := "switch II (58-76)"]

# 检查分布
cat("\n功能环区位置分布统计:\n")
print(table(pivot_dt$loop_group))

# 绘制图2：功能环区
p_functional_loops <- ggplot(pivot_dt, aes(x = tsne_1, y = tsne_2, color = loop_group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE Plot (Functional Loops)", 
       x = "t-SNE 1", y = "t-SNE 2",
       color = "Functional Loop") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = functional_loop_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_functional_loops)





########## 4. t-SNE 图（按二级结构着色）

################################################################################
## 定义二级结构区域
################################################################################

# β-sheet 区域
rects_sheet <- data.frame(
  xstart = c(3, 38, 51, 77, 109, 139),
  xend = c(9, 44, 57, 84, 115, 143),
  col = c("β1", "β2", "β3", "β4", "β5", "β6")
)

# α-helix 区域
rects_helix <- data.frame(
  xstart = c(15, 67, 87, 127, 148),
  xend = c(24, 73, 104, 136, 166),
  col = c("α1", "α2", "α3", "α4", "α5")
)

# 创建二级结构标签的函数
create_secondary_structure_labels <- function(pos) {
  # 检查是否在 β-sheet 中
  for (i in 1:nrow(rects_sheet)) {
    if (pos >= rects_sheet$xstart[i] && pos <= rects_sheet$xend[i]) {
      return(rects_sheet$col[i])
    }
  }
  
  # 检查是否在 α-helix 中
  for (i in 1:nrow(rects_helix)) {
    if (pos >= rects_helix$xstart[i] && pos <= rects_helix$xend[i]) {
      return(rects_helix$col[i])
    }
  }
  
  # 如果都不在，返回"Other"
  return("Other")
}

# 添加二级结构标签到数据
pivot_dt[, secondary_structure := sapply(Pos_real, create_secondary_structure_labels)]

# 将二级结构转为因子，指定顺序
secondary_structure_order <- c(paste0("β", 1:6), paste0("α", 1:5), "Other")
pivot_dt[, secondary_structure := factor(secondary_structure, levels = secondary_structure_order)]

# 检查分布
cat("二级结构分布统计:\n")
print(table(pivot_dt$secondary_structure))
# 定义二级结构的配色方案
# 使用Set3调色板，确保有足够多的颜色
secondary_colors <- c(
  "β1" = "#41A67E", "β2" = "#FFFFB3", "β3" = "#BEBADA", 
  "β4" = "#FB8072", "β5" = "#1581BF", "β6" = "#FDB462",
  "α1" = "#B3DE69", "α2" = "#FCCDE5", "α3" = "#FF0087",
  "α4" = "#9112BC", "α5" = "#CCEBC5",
  "Other" = "#CBCBCB"
)

p_tsne_secondary <- ggplot(pivot_dt, aes(x = tsne_1, y = tsne_2, color = secondary_structure)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE Plot (colored by Secondary Structure)", 
       x = "t-SNE 1", y = "t-SNE 2",
       color = "Secondary Structure") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10)) +
  scale_color_manual(values = secondary_colors) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))

print(p_tsne_secondary)
