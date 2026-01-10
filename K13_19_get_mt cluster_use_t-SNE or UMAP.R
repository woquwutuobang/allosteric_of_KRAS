library(data.table)
library(krasddpcams)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(umap)


wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

ddG_RAF1<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAF1.txt",assay_sele = "RAF1")
ddG_RALGDS<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_RAL.txt",assay_sele = "RALGDS")
ddG_PI3KCG<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_PI3.txt",assay_sele = "PI3KCG")
ddG_SOS1<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_SOS.txt",assay_sele = "SOS1")
ddG_K55<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K55.txt",assay_sele = "K55")
ddG_K27<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K27.txt",assay_sele = "K27")
ddG_K13<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K13.txt",assay_sele = "K13")
ddG_K19<-krasddpcams__read_ddG("C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20251121/task_901/weights/weights_Binding_K19.txt",assay_sele = "K19")

names(ddG_RAF1)




make_mut_id <- function(dt) {
  dt[, mut_id := paste0(Pos_real, "_", mt)]
  dt
}

ddG_list <- list(
  K13    = make_mut_id(ddG_K13),
  K19    = make_mut_id(ddG_K19),
  K27    = make_mut_id(ddG_K27),
  K55    = make_mut_id(ddG_K55),
  RAF1   = make_mut_id(ddG_RAF1),
  RALGDS = make_mut_id(ddG_RALGDS),
  PI3KCG = make_mut_id(ddG_PI3KCG),
  SOS1   = make_mut_id(ddG_SOS1)
)


ddG_long <- rbindlist(
  lapply(names(ddG_list), function(a) {
    ddG_list[[a]][, .(
      mut_id,
      assay = a,
      ddG = `mean_kcal/mol`
    )]
  })
)


ddG_wide <- dcast(
  ddG_long,
  mut_id ~ assay,
  value.var = "ddG",
  fun.aggregate = mean,   # 或者 median
  na.rm = TRUE
)


ddG_mat <- data.table::as.data.table(na.omit(ddG_wide))


############################ t-SNE plot
set.seed(123)

X <- as.matrix(ddG_mat[, -1, with = FALSE])

tsne_res <- Rtsne(
  X,
  perplexity = 30,
  scale = TRUE,
  pca = TRUE
)

tsne_dt <- data.table(
  mut_id = ddG_mat$mut_id,
  TSNE1 = tsne_res$Y[,1],
  TSNE2 = tsne_res$Y[,2]
)



# melt 数据
tsne_long <- melt(
  cbind(tsne_dt, ddG_mat[, -1, with = FALSE]),
  id.vars = c("mut_id", "TSNE1", "TSNE2"),
  variable.name = "assay",
  value.name = "ddG"
)

# 指定 assay 顺序，确保是 factor
assay_order <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
tsne_long$assay <- factor(tsne_long$assay, levels = assay_order)

# 绘图
ggplot(tsne_long, aes(TSNE1, TSNE2, color = ddG)) +
  geom_point(size = 0.8) +
  scale_color_gradient2(low = "#1B38A6", mid = "grey", high = "#F4270C", midpoint = 0) +
  facet_wrap(~ assay, nrow = 2, ncol = 4, drop = FALSE) +  # drop=FALSE 保留所有 factor 水平
  theme_classic()





###### UMAP
umap_res <- umap(X)

umap_dt <- data.table(
  mut_id = ddG_mat$mut_id,
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2]
)

# melt 数据
umap_long <- melt(
  cbind(umap_dt, ddG_mat[, -1, with = FALSE]),
  id.vars = c("mut_id", "UMAP1", "UMAP2"),
  variable.name = "assay",
  value.name = "ddG"
)

# 指定 assay 顺序
assay_order <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
umap_long$assay <- factor(umap_long$assay, levels = assay_order)

# 绘图
ggplot(umap_long, aes(UMAP1, UMAP2, color = ddG)) +
  geom_point(size = 0.8) +
  scale_color_gradient2(low = "#1B38A6", mid = "grey", high = "#F4270C", midpoint = 0) +
  facet_wrap(~ assay, nrow = 2, ncol = 4, drop = FALSE) +  # 两行四列排列
  theme_classic()

