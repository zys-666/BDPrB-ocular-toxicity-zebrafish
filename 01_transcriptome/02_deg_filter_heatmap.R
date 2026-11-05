## Reproducible analysis code accompanying the manuscript.
## Run from the repository root (see README.md); paths are resolved from config.R.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')
#!/usr/bin/env Rscript
# ============================================================
# 六组比较 → 四组12样本统一热图
# 每个差异文件只含两组样本，需合并提取全部12样本
# ============================================================

packages <- c("tidyverse", "pheatmap", "RColorBrewer")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# ============================================================
# 1. 参数设置
# ============================================================

FDR_CUTOFF    <- 0.05
LOG2FC_CUTOFF <- 1
TOP_N         <- 50

base_dir <- file.path(DATA_ROOT, "3.转录组", "2.热图")
output_pdf <- file.path(base_dir, "heatmap_four_groups_12samples_50DEGs.pdf")

COLOR_CTL <- "grey80"
COLOR_L   <- "#F4A6A0"
COLOR_H   <- "#6BAED6"
COLOR_T   <- "#74C476"

# 六个比较文件
files <- list(
  C_vs_L = file.path(base_dir, "C_vs_L_all_geneSymbol.xls"),
  C_vs_H = file.path(base_dir, "C_vs_H_all_geneSymbol.xls"),
  C_vs_T = file.path(base_dir, "C_vs_T_all_geneSymbol.xls"),
  L_vs_H = file.path(base_dir, "L_vs_H_all_geneSymbol.xls"),
  L_vs_T = file.path(base_dir, "L_vs_T_all_geneSymbol.xls"),
  H_vs_T = file.path(base_dir, "H_vs_T_all_geneSymbol.xls")
)

# ============================================================
# 2. 读取六组差异结果，筛选显著基因
# ============================================================

read_deg_sig <- function(file_path, comparison_name) {
  if (!file.exists(file_path)) {
    message("⚠️  文件不存在，跳过: ", file_path)
    return(NULL)
  }
  
  df <- read.table(file_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  
  # 相同基因取平均
  df_mean <- df %>%
    group_by(GENE) %>%
    summarise(
      FDR = min(FDR, na.rm = TRUE),
      log2FC = mean(log2FC, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_mean %>%
    mutate(
      significant = (FDR < FDR_CUTOFF & abs(log2FC) > LOG2FC_CUTOFF),
      comparison = comparison_name
    )
}

deg_list <- map(names(files), ~ read_deg_sig(files[[.x]], .x))
deg_all <- bind_rows(deg_list)

# 按综合显著性排序取Top N
gene_scores <- deg_all %>%
  group_by(GENE) %>%
  summarise(
    n_significant = sum(significant, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(score = n_significant * 10 - log10(min_FDR) + max_abs_log2FC) %>%
  arrange(desc(score))

sig_genes <- gene_scores %>%
  slice_head(n = TOP_N) %>%
  pull(GENE)

message("六组比较共筛选到 Top ", length(sig_genes), " 个差异基因")

# ============================================================
# 3. 从六个文件中提取全部12样本的FPKM，合并成完整表达矩阵
# ============================================================

# 定义每个文件包含的样本列
sample_map <- list(
  C_vs_L = c("C1_FPKM", "C2_FPKM", "C3_FPKM", "L1_FPKM", "L2_FPKM", "L3_FPKM"),
  C_vs_H = c("C1_FPKM", "C2_FPKM", "C3_FPKM", "H1_FPKM", "H2_FPKM", "H3_FPKM"),
  C_vs_T = c("C1_FPKM", "C2_FPKM", "C3_FPKM", "T1_FPKM", "T2_FPKM", "T3_FPKM"),
  L_vs_H = c("L1_FPKM", "L2_FPKM", "L3_FPKM", "H1_FPKM", "H2_FPKM", "H3_FPKM"),
  L_vs_T = c("L1_FPKM", "L2_FPKM", "L3_FPKM", "T1_FPKM", "T2_FPKM", "T3_FPKM"),
  H_vs_T = c("H1_FPKM", "H2_FPKM", "H3_FPKM", "T1_FPKM", "T2_FPKM", "T3_FPKM")
)

# 逐个文件读取，提取样本列
expr_parts <- list()

for (comp in names(files)) {
  if (!file.exists(files[[comp]])) next
  
  df <- read.table(files[[comp]], header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  
  # 相同基因取平均
  cols_to_use <- c("GENE", sample_map[[comp]])
  
  df_mean <- df %>%
    select(all_of(cols_to_use)) %>%
    group_by(GENE) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  expr_parts[[comp]] <- df_mean
}

# 合并所有表达量（按GENE连接）
expr_merged <- expr_parts %>%
  reduce(full_join, by = "GENE")

# 处理NA（某些基因可能只在部分比较中出现）
expr_merged <- expr_merged %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# 提取目标基因
mat_all <- expr_merged %>%
  filter(GENE %in% sig_genes) %>%
  column_to_rownames("GENE") %>%
  as.matrix()

# 去重列名（C1_FPKM等可能出现多次），取第一次出现的
mat_all <- mat_all[, !duplicated(colnames(mat_all))]

# 确保列顺序 C-L-H-T
expected_order <- c("C1_FPKM", "C2_FPKM", "C3_FPKM",
                    "L1_FPKM", "L2_FPKM", "L3_FPKM",
                    "H1_FPKM", "H2_FPKM", "H3_FPKM",
                    "T1_FPKM", "T2_FPKM", "T3_FPKM")

mat_all <- mat_all[, expected_order]

message("表达矩阵维度: ", nrow(mat_all), " 基因 × ", ncol(mat_all), " 样本")

# ============================================================
# 4. Z-score 标准化
# ============================================================

z_mat <- t(scale(t(mat_all)))
z_mat[is.nan(z_mat)] <- 0
z_mat[z_mat > 2] <- 2
z_mat[z_mat < -2] <- -2

# ============================================================
# 5. 列注释 & 颜色
# ============================================================

anno_col <- data.frame(
  Group = c(rep("Control", 3), rep("Low", 3), rep("High", 3), rep("Treatment", 3)),
  row.names = colnames(z_mat)
)

anno_colors <- list(
  Group = c(Control = COLOR_CTL, Low = COLOR_L, High = COLOR_H, Treatment = COLOR_T)
)

col_breaks <- seq(-2, 2, length.out = 101)
heatmap_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE",
                                     "white",
                                     "#F4A582", "#D6604D", "#B2182B"))(100)

# ============================================================
# 6. 绘制热图
# ============================================================

pdf(output_pdf, width = 14, height = 12)

pheatmap(
  z_mat,
  color             = heatmap_colors,
  breaks            = col_breaks,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  annotation_col    = anno_col,
  annotation_colors = anno_colors,
  annotation_names_col = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 8,
  fontsize_col      = 10,
  angle_col         = 45,
  main              = paste0("Heatmap of Top ", nrow(z_mat), " DEGs (6 Comparisons)"),
  legend            = TRUE,
  annotation_legend = TRUE,
  border_color      = "white",
  cellwidth         = 22,
  cellheight        = 14
)

dev.off()
message("✅ PDF 已保存至: ", output_pdf)