## Reproducible analysis code accompanying the manuscript.
## Run from the repository root (see README.md); paths are resolved from config.R.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')
rm(list=ls())
setwd(file.path(DATA_ROOT, "3.转录组", "1.火山图"))

library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)

# 加载数据
df <- read.table("data.txt", header = 1, check.names = F, sep = "\t")
df$group <- factor(df$group, levels = c("C-L", "C-H", "C-PP", "L-H", "L-PP", "H-PP"))

# 显著性分组
df$group2 <- as.factor(ifelse(df$p_value < 0.05 & abs(df$log2FC) >= 2,
                              ifelse(df$log2FC >= 2, 'up', 'down'), 'NS'))
df$group2 <- factor(df$group2, levels = c("up", "down", "NS"))

# 每组取 p_value < 0.05 且 |log2FC| 最大的前3个基因作为标签
top3_idx <- df %>%
  mutate(idx = row_number()) %>%
  filter(p_value < 0.05) %>%
  group_by(group) %>%
  slice_max(order_by = abs(log2FC), n = 3) %>%
  pull(idx)

df$label <- ''
df$label[top3_idx] <- as.character(df$OTU[top3_idx])

# 提前固定抖动坐标，确保点和标签连线一一对应
set.seed(42)
df$jitter_x <- as.numeric(df$group) + runif(nrow(df), -0.4, 0.4)
df$jitter_y <- df$log2FC

# 每组最大最小值（用于背景柱）
df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(log2FC), min_log2FC = min(log2FC))

# 6个组的填充颜色
fill_colors <- c(
  "C-L" = "#ed7902",
  "C-H" = "#ef5734",
  "C-PP" = "#b5c327",
  "L-H" = "#3a9fbf",
  "L-PP" = "#8e44ad",
  "H-PP" = "#27ae60"
)

# 渐变背景
gradient_colors <- colorRampPalette(c("#e8d5b0", "#c8dcd8"))(100)
gradient_matrix <- matrix(gradient_colors, nrow = 100, ncol = 1)
gradient_raster <- rasterGrob(
  as.raster(gradient_matrix),
  width = unit(1, "npc"),
  height = unit(1, "npc"),
  interpolate = TRUE
)

# 绘图
p_final <- ggplot() +
  # 渐变背景
  annotation_custom(gradient_raster, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  
  # 灰色背景柱
  geom_col(data = df_bg, mapping = aes(group, max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  geom_col(data = df_bg, mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  
  # 数据点（使用预先固定的抖动坐标）
  geom_point(data = df, mapping = aes(x = jitter_x, y = jitter_y, color = group2),
             size = 3, alpha = 0.7) +
  
  # 分组色块
  geom_col(data = df_bg, mapping = aes(x = group, y = 0.4, fill = group), width = 0.8) +
  geom_col(data = df_bg, mapping = aes(x = group, y = -0.4, fill = group), width = 0.8) +
  
  # 分组标签文字
  geom_text(data = df_bg, mapping = aes(x = group, y = 0, label = group),
            size = 4, color = "#dbebfa", fontface = "bold") +
  
  # 颜色映射
  scale_color_manual(values = c("up" = "#e42313", "down" = "#0061d5", "NS" = "#8b8c8d")) +
  scale_fill_manual(values = fill_colors) +
  
  # 基因标签（使用相同的抖动坐标）
  geom_text_repel(data = df[df$label != "", ],
                  mapping = aes(x = jitter_x, y = jitter_y, label = label),
                  max.overlaps = 10000, size = 3,
                  box.padding = unit(0.8, 'lines'),
                  point.padding = unit(0.8, 'lines'),
                  segment.color = 'black',
                  show.legend = FALSE) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +  # 自动选择约10个刻度
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  # 主题
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.8),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.8),
    panel.background = element_rect(fill = NA),
    plot.background = element_rect(fill = "#fafafa")
  ) +
  
  labs(x = "group", y = "Log2FoldChange", fill = NULL, color = NULL) 
  

print(p_final)

ggsave("volcano_plot_6groups.pdf", p_final, width = 14, height = 8)