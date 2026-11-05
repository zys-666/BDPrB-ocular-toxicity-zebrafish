## Reproducible analysis code accompanying the manuscript.
## Run from the repository root (see README.md); paths are resolved from config.R.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')

library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
library(geomtextpath)

# ══════════════════════════════════════════
# 1. 读取数据
# ══════════════════════════════════════════
# 1.1 读取数据
# ══════════════════════════════════════════
input_file <- file.path(DATA_ROOT, "3.转录组", "3.GO", "多组GO前5通路环形图", "Pathway_GO.xlsx")
wb <- loadWorkbook(input_file)

pathway_df <- read.xlsx(wb, "Pathway")
abbre_df   <- read.xlsx(wb, "Abbre")

pathway_df$Pathway <- pathway_df$Pathway %>%
  stringr::str_trim() %>%
  stringr::str_squish()

abbre_df$Pathway <- abbre_df$Pathway %>%
  stringr::str_trim() %>%
  stringr::str_squish()

cat("原始数据行数：\n")
cat("Pathway:", nrow(pathway_df), "\n")
cat("Abbre:", nrow(abbre_df), "\n\n")

# ══════════════════════════════════════════
# 1.2 检查 Abbre 是否有重复 Pathway
# ══════════════════════════════════════════
dup_pathway <- abbre_df %>%
  filter(duplicated(Pathway) | duplicated(Pathway, fromLast = TRUE)) %>%
  arrange(Pathway)

if(nrow(dup_pathway) > 0){
  
  write.table(
    dup_pathway,
    file.path(WORK, "out", "Abbre_duplicate_pathway.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("⚠ 检测到 Abbre 中存在重复 Pathway\n")
  cat("已导出文件: Abbre_duplicate_pathway.txt\n\n")
}

# ══════════════════════════════════════════
# 1.3 去除重复 Pathway（只保留第一个）
# ══════════════════════════════════════════
abbre_unique <- abbre_df %>%
  distinct(Pathway, .keep_all = TRUE)

# ══════════════════════════════════════════
# 1.4 Pathway 匹配 Abbre
# ══════════════════════════════════════════
data <- pathway_df %>%
  left_join(abbre_unique %>% select(Pathway, Abbre), by = "Pathway") %>%
  mutate(individual = Abbre) %>%
  select(-Abbre)

# ══════════════════════════════════════════
# 1.5 检查未匹配 Pathway
# ══════════════════════════════════════════
unmatched <- data %>%
  filter(is.na(individual))

if(nrow(unmatched) > 0){
  
  write.table(
    unmatched,
    file.path(WORK, "out", "unmatched_pathway.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("⚠ 检测到未匹配 Pathway:", nrow(unmatched), "\n")
  cat("已导出文件: unmatched_pathway.txt\n\n")
}else{
  
  cat("✓ 所有 Pathway 均成功匹配 Abbre\n\n")
}

# ══════════════════════════════════════════
# 1.6 结果（检查用）
# ══════════════════════════════════════════
write.table(
  data,
  file.path(WORK, "out", "Pathway_match_check.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("✓ 匹配数据已导出: Pathway_match_check.txt\n\n")

# ══════════════════════════════════════════
# 2. 数据整理
# ══════════════════════════════════════════
data_sorted <- data %>%
  group_by(group, individual, category) %>%
  summarise(value = max(as.numeric(`log10(Pvalue)`)), .groups = "drop") %>%
  mutate(individual = as.character(individual)) %>%
  arrange(group, category, desc(value))          # ★ 先category再desc(value)

datagroup <- unique(data_sorted$group)
n_groups  <- length(datagroup)

gap_size  <- 1

empty_data <- tibble(
  group      = rep(datagroup, each = gap_size),
  individual = paste0("empty_", seq_len(n_groups * gap_size)),
  value      = 0,
  is_empty   = TRUE,
  category   = NA_character_
)

allplotdata <- bind_rows(
  data_sorted %>% mutate(is_empty = FALSE),
  empty_data
) %>%
  arrange(match(group, datagroup), is_empty) %>%  # ★ 去掉desc(value)，保留category排序
  mutate(xid = row_number())

total_bars <- nrow(allplotdata)

allplotdata <- allplotdata %>%
  mutate(
    angle_raw = 90 - 360 * (xid - 0.5) / total_bars,
    hjust = ifelse(angle_raw < -90, 1, 0),
    angle = ifelse(angle_raw < -90, angle_raw + 180, angle_raw)
  )

# ══════════════════════════════════════════
# 3. 核心几何参数
# ══════════════════════════════════════════
inner_r  <- 1
max_val  <- ceiling(max(data_sorted$value, na.rm = TRUE))
y_breaks <- seq(0, max_val, by = 2)
y_scale  <- 1.5
outer_r  <- inner_r + max_val * y_scale

arc_r    <- inner_r - 3      # 组级弧半径
arc_w    <- 2.5

# ★ category弧夹在组弧(arc_r)和柱子底部(inner_r)之间
cat_arc_r <- (arc_r + inner_r) / 2   # 取中间，可微调
cat_arc_w <- 2.0

bar_width    <- 0.95
label_offset <- 0.1

# ══════════════════════════════════════════
# 4. Y轴标签位置（最后一组的 empty 间隔）
# ══════════════════════════════════════════
gap_xids_axis <- allplotdata %>%
  filter(is_empty, group == datagroup[n_groups]) %>%
  pull(xid)

axis_xid <- median(gap_xids_axis)

axis_labels <- tibble(
  xid   = axis_xid,
  y_val = y_breaks * y_scale + inner_r,
  label = as.character(y_breaks)
)

axis_title <- tibble(
  xid   = axis_xid,
  y_val = outer_r + 3,
  label = expression(-log[10](italic(p)~value))
)

# ══════════════════════════════════════════
# 5. 网格线
# ══════════════════════════════════════════
gap_ranges <- allplotdata %>%
  filter(is_empty) %>%
  group_by(group) %>%
  summarise(x_lo = min(xid) - 0.5, x_hi = max(xid) + 0.5, .groups = "drop")

grid_segs <- expand.grid(
  g_idx   = seq_len(nrow(gap_ranges)),
  y_break = y_breaks
) %>%
  mutate(
    x_lo  = gap_ranges$x_lo[g_idx],
    x_hi  = gap_ranges$x_hi[g_idx],
    y_val = y_break * y_scale + inner_r
  )

# ══════════════════════════════════════════
# 6. 彩色圆弧（组级，柱子根部）
# ══════════════════════════════════════════
arc_segs <- allplotdata %>%
  filter(!is_empty) %>%
  group_by(group) %>%
  summarise(x_lo = min(xid) - 0.4, x_hi = max(xid) + 0.4, .groups = "drop") %>%
  rowwise() %>%
  reframe(
    group = group,
    xid   = seq(x_lo, x_hi, length.out = 80),
    y     = arc_r
  )

# ══════════════════════════════════════════
# 6b. ★ 新增：category 彩色圆弧（BP/CC/MF）
# 按 group + category 的连续段分别画弧
# ══════════════════════════════════════════
cat_arc_segs <- allplotdata %>%
  filter(!is_empty) %>%
  # 找每个 group 内各 category 的 xid 范围
  group_by(group, category) %>%
  summarise(x_lo = min(xid) - 0.4, x_hi = max(xid) + 0.4, .groups = "drop") %>%
  rowwise() %>%
  reframe(
    group    = group,
    category = category,
    xid      = seq(x_lo, x_hi, length.out = 60),
    y        = cat_arc_r
  )

# ══════════════════════════════════════════
# 7. 弧形组标签
# ══════════════════════════════════════════
arc_labels <- allplotdata %>%
  filter(!is_empty) %>%
  group_by(group) %>%
  summarise(x_lo = min(xid), x_hi = max(xid), .groups = "drop") %>%
  rowwise() %>%
  reframe(
    group = group,
    label = group,
    xid   = seq(x_lo, x_hi, length.out = 60),
    y     = inner_r - 6
  )

# ══════════════════════════════════════════
# 8. 颜色
# ══════════════════════════════════════════
palette <- c(
  "#D1C7B7", "#BFA6A0", "#7E8AA2", "#7A6FA8",
  "#8FA876", "#B07D7D", "#6C8EBF", "#9C8F6F",
  "#7F9C96", "#A68F8F", "#8B9E70", "#7E8AA2"
)
grp_colors <- setNames(palette[seq_len(n_groups)], datagroup)

# ★ category 颜色（BP/CC/MF 固定配色，可自行调整）
cat_colors <- c(
  "BP" = "#E05C5C",
  "CC" = "#5B9BD5",
  "MF" = "#70AD47"
)

# ══════════════════════════════════════════
# 9. 绘图
# ══════════════════════════════════════════
library(grid)
library(cowplot)

p <- ggplot() +
  
  # 9.1 网格线
  geom_segment(
    data = grid_segs,
    aes(x = x_lo, xend = x_hi, y = y_val, yend = y_val),
    color = "grey65", linewidth = 0.4
  ) +
  
  # 9.2 柱状图
  geom_bar(
    data = allplotdata,
    aes(x = xid, y = value * y_scale + inner_r, fill = group),
    stat = "identity", width = bar_width, color = NA
  ) +
  
  # 9.3 白色底层（遮住柱子根部到inner_r以下）
  geom_bar(
    data = allplotdata,
    aes(x = xid, y = inner_r),
    stat = "identity", width = bar_width,
    fill = "white", color = NA
  ) +
  
  # 9.4 组级彩色圆弧（最内层，arc_r位置）
  geom_path(
    data = arc_segs,
    aes(x = xid, y = y, color = group, group = group),
    linewidth = arc_w
  ) +
  
  # 9.4b ★ category彩色圆弧（arc_r和inner_r之间）
  geom_path(
    data = cat_arc_segs,
    aes(x = xid, y = y,
        color = category,
        group = interaction(group, category)),
    linewidth = cat_arc_w
  ) +
  
  # 9.5 内圆白色填充（遮住组弧内侧）
  geom_rect(
    aes(xmin = 0.5, xmax = total_bars + 0.5,
        ymin = -Inf, ymax = arc_r - arc_w / 2 - 0.3),
    fill = "white", color = NA
  ) +
  
  # 9.6 弧形组标签
  geom_textpath(
    data = arc_labels,
    aes(x = xid, y = y, label = label, color = group),
    size = 2.8, fontface = "bold",
    text_only = TRUE, upright = TRUE
  ) +
  
  # 9.7 individual文字标签
  geom_text(
    data = allplotdata %>% filter(!is_empty),
    aes(
      x     = xid,
      y     = value * y_scale + inner_r + label_offset,
      label = individual,
      angle = angle,
      hjust = hjust
    ),
    color = "grey25", size = 2.3, vjust = 0.5
  ) +
  
  # 9.8 Y轴刻度标签
  geom_text(
    data = axis_labels,
    aes(x = xid, y = y_val, label = label),
    color = "black", size = 2.6,
    angle = 0, fontface = "bold",
    hjust = 0.5, vjust = 0.5
  ) +
  
  # 9.9 Y轴标题
  geom_text(
    data = axis_title,
    aes(x = xid, y = y_val, label = label),
    color = "black", size = 3.2, fontface = "bold",
    angle = 0, hjust = 0.5, parse = TRUE
  ) +
  
  # 9.10 中心标题（只保留GO，图例另外叠加）
  annotate("text",
           x = 1, y = -(abs(arc_r) + arc_w + 11),
           label = "GO",
           size = 7, fontface = "bold", color = "black"
  ) +
  
  scale_fill_manual(values = grp_colors) +
  scale_color_manual(values = c(grp_colors, cat_colors)) +
  scale_x_continuous(limits = c(0.5, total_bars + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-(abs(arc_r) + arc_w + 12), outer_r + 18)) +
  coord_polar(theta = "x", start = 0, clip = "off") +
  theme_void() +
  theme(legend.position = "none")

# ══════════════════════════════════════════
# 9b. 单独构建图例（普通ggplot，不用极坐标）
# ══════════════════════════════════════════
legend_df <- tibble(
  label    = c("BP", "CC", "MF"),
  color    = c("#E05C5C", "#5B9BD5", "#70AD47"),
  y        = c(3, 2.6, 2.2)
)

legend_plot <- ggplot(legend_df) +
  geom_point(aes(x = 0, y = y, color = label),
             shape = 15, size = 3) +
  geom_text(aes(x = 0.3, y = y, label = label, color = label),
            hjust = 0, fontface = "bold", size = 3.5) +
  scale_color_manual(values = c("BP" = "#E05C5C",
                                "CC" = "#5B9BD5",
                                "MF" = "#70AD47")) +
  scale_x_continuous(limits = c(-0.5, 2)) +
  scale_y_continuous(limits = c(0.1, 3.8)) +
  theme_void() +
  theme(legend.position = "none")

# ══════════════════════════════════════════
# 9c. 用cowplot叠加：图例放在圆心位置
# ══════════════════════════════════════════
# x/y/width/height 均为0~1的相对坐标，调整对准圆心
final_plot <- ggdraw(p) +
  draw_plot(legend_plot,
            x = 0.4695,    # 左右位置（0=最左，1=最右）
            y = 0.46,    # 上下位置（0=最下，1=最上）
            width  = 0.08,
            height = 0.12)

final_plot

# ══════════════════════════════════════════
# 10. 保存
# ══════════════════════════════════════════
out_dir <- dirname(input_file)

ggsave(
  file.path(out_dir, "pathway_GO.pdf"), final_plot, width = 10, height = 10)

