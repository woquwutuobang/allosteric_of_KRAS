library(dplyr)
library(ggplot2)

# ==============================
# 1. WT序列和位置
# ==============================
wt_aa <- unlist(strsplit("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM", ""))
positions <- 1:length(wt_aa)

# ==============================
# 2. 配体结合位点
# ==============================
K13_binding_interface <- c(63,68,87,88,90,91,92,94,95,96,97,98,99,101,102,105,106,107,129,133,136,137,138) 
K19_binding_interface <- c(68,87,88,90,91,92,94,95,97,98,99,101,102,105,107,108,125,129,133,136,137) 
RAF1_binding_interface <- c(21,25,29,31,33,36,37,38,39,40,41,67,71) 
K27_binding_interface <- c(21,24,25,27,29,30,31,32,33,34,35,36,38,39,40,41,43,52,54,67,70,71)

# ==============================
# 3. 构建数据框
# ==============================
df <- data.frame(
  Position = rep(positions, 4),
  Residue = rep(wt_aa, 4),
  Partner = rep(c("DARPin K13", "DARPin K19", "RAF1-RBD","K27"), each = length(wt_aa)),
  Contact = FALSE
)

# 设置接触位点
df$Contact[df$Partner == "DARPin K13" & df$Position %in% K13_binding_interface] <- TRUE
df$Contact[df$Partner == "DARPin K19" & df$Position %in% K19_binding_interface] <- TRUE
df$Contact[df$Partner == "RAF1-RBD" & df$Position %in% RAF1_binding_interface] <- TRUE
df$Contact[df$Partner == "K27" & df$Position %in% K27_binding_interface] <- TRUE

# 添加 ResLabel
df <- df %>%
  mutate(ResLabel = paste0(Residue, Position))

# ==============================
# 4. 只保留BI1/BI2配体接触的残基
# ==============================
plot_df <- df %>%
  filter(Partner %in% c("RAF1-RBD", "K27", "DARPin K13", "DARPin K19")) %>%
  mutate(Fill = case_when(
    Partner %in% c("RAF1-RBD", "K27") & Contact ~ "BI1",
    Partner %in% c("DARPin K13", "DARPin K19") & Contact ~ "BI2",
    TRUE ~ NA_character_   # 非接触点不显示
  )) %>%
  filter(!is.na(Fill))  # 只保留有接触的残基

# 设置因子顺序
residue_labels <- plot_df %>% distinct(Position, Residue) %>% arrange(Position) %>%
  mutate(ResLabel = paste0(Residue, Position)) %>% pull(ResLabel)
plot_df$ResLabel <- factor(plot_df$ResLabel, levels = residue_labels)

plot_df$Partner <- factor(plot_df$Partner, 
                          levels = c("RAF1-RBD", "K27", "DARPin K13", "DARPin K19"))
plot_df$Fill <- factor(plot_df$Fill, levels = c("BI1", "BI2"))

# ==============================
# 5. 绘图
# ==============================
ggplot(plot_df, aes(x = ResLabel, y = Partner, fill = Fill)) +
  geom_tile(color = "grey90", width = 1, height = 1) +
  scale_fill_manual(
    values = c(
      "BI1" = alpha("#F4270C", 0.8),
      "BI2" = alpha("#1B38A6", 0.8)
    ),
    name = "Contact Type",
    labels = c("BI1 (RAF1/K27)", "BI2 (K13/K19)")
  ) +
  scale_y_discrete(limits = rev(levels(plot_df$Partner))) + 
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 8) +
  labs(x = "Residue", y = "Binding Partner", fill = "") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )



ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20251108_revise_figures/f3/20251122/Common_and_unique_structural_contacts_with_BI.pdf", 
       width = 10, height = 5.5, dpi = 300)  # 调整高度
