library(data.table)
library(krasddpcams)
library(ggplot2)



wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

ddG_RAF1<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",assay_sele = "RAF1")
ddG_RALGDS<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAL.txt",assay_sele = "RALGDS")
ddG_PI3KCG<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_PI3.txt",assay_sele = "PI3KCG")
ddG_SOS1<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_SOS.txt",assay_sele = "SOS1")
ddG_K55<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K55.txt",assay_sele = "K55")
ddG_K27<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",assay_sele = "K27")
ddG_K13<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",assay_sele = "K13")
ddG_K19<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",assay_sele = "K19")

names(ddG_RAF1)



# 将所有 ddG 数据合并到一个列表
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



# 绘制热图
ggplot(median_dt, aes(x = assay, y = Pos_real, fill = median_ddG)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#1B38A6", mid = "grey", high = "#F4270C", midpoint = 0, name = "ΔΔG (kcal/mol)") +
  scale_y_reverse() +  # Pos_real 从上到下
  theme_minimal() +
  labs(x = "Assay", y = "Position", title = "Median ΔΔG Heatmap") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
