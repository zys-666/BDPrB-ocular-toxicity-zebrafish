## Reproducible analysis code accompanying the manuscript.
## Run from the repository root (see README.md); paths are resolved from config.R.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')
# ============================================================
# 菌群属水平相对丰度与ELISA指标相关性热图 (Spearman)
# 文件路径: D:\user data\Desktop\相关性
# ============================================================

# 安装所需包（如未安装）
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("dplyr")
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("corrplot")

library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(corrplot)

setwd(file.path(DATA_ROOT, "9.菌群测序", "菌群数量和核心基因表达量的相关性"))

# ============================================================
# 1. 读取菌群数据 (OTU表)
# ============================================================
df <- read.delim("allsample.txt", sep = "\t", header = TRUE, check.names = FALSE)
print(paste("菌群数据形状:", nrow(df), "行 x", ncol(df), "列"))
print(paste("列名:", paste(colnames(df), collapse = ", ")))

# 先定义样本列（排除 OTU ID 和 taxonomy，此时还没有 Genus）
sample_cols <- setdiff(colnames(df), c("OTU ID", "taxonomy"))
cat(paste("\n样本列:", paste(sample_cols, collapse = ", "), "\n"))

# 提取属名函数
extract_genus <- function(taxonomy) {
  parts <- strsplit(taxonomy, "; ")[[1]]
  for (part in parts) {
    if (grepl("^g__", part)) {
      genus <- gsub("g__", "", part)
      if (genus != "" && genus != "unclassified") {
        return(genus)
      }
    }
  }
  return("Unclassified")
}

# 应用提取（此时 sample_cols 已固定，不受影响）
df$Genus <- sapply(df$taxonomy, extract_genus)

# 查看提取结果
cat("\n提取的属名示例:\n")
print(head(df[, c("OTU ID", "Genus")], 10))
cat(paste("\n共有", length(unique(df$Genus)), "个不同的属\n"))

# ❌ 删除下面这行重复定义！
# sample_cols <- setdiff(colnames(df), c("OTU ID", "taxonomy"))

# 按属合并reads数
counts_genus <- df %>%
  group_by(Genus) %>%
  summarise(across(all_of(sample_cols), sum)) %>%
  as.data.frame()
rownames(counts_genus) <- counts_genus$Genus
counts_genus$Genus <- NULL

# 转化为相对丰度（百分比）
genus_rel <- apply(counts_genus, 2, function(x) x / sum(x) * 100)

# 转置：样本为行，菌属为列
genus_rel_t <- as.data.frame(t(genus_rel))

cat(paste("\n相对丰度表形状:", nrow(genus_rel_t), "行 x", ncol(genus_rel_t), "列\n"))
cat("\nTop 20 属（平均相对丰度%）:\n")
top20 <- sort(colMeans(genus_rel_t), decreasing = TRUE)[1:20]
print(round(top20, 3))

# ============================================================
# 2. 读取并处理ELISA数据
# ============================================================
elisa <- read.delim("elisa_data.txt", sep = "\t", header = TRUE, check.names = FALSE)
rownames(elisa) <- elisa[, 1]
elisa <- elisa[, -1]

cat(paste("\nELISA数据形状:", nrow(elisa), "行 x", ncol(elisa), "列\n"))
cat(paste("ELISA指标:", paste(colnames(elisa), collapse = ", "), "\n"))
cat(paste("ELISA样本:", paste(rownames(elisa), collapse = ", "), "\n"))

# 按组取平均
group_means <- data.frame(
  C = colMeans(elisa[grepl("^C", rownames(elisa)), , drop = FALSE]),
  L = colMeans(elisa[grepl("^L", rownames(elisa)), , drop = FALSE]),
  H = colMeans(elisa[grepl("^H", rownames(elisa)), , drop = FALSE]),
  PP = colMeans(elisa[grepl("^PP", rownames(elisa)), , drop = FALSE])
)

cat("\n组平均值:\n")
print(group_means)

# 为每个菌群样本分配对应组的平均值（修复版）
elisa_matched <- as.data.frame(matrix(NA, nrow = length(sample_cols), ncol = ncol(elisa)))
rownames(elisa_matched) <- sample_cols
colnames(elisa_matched) <- colnames(elisa)

for (s in sample_cols) {
  g <- substr(s, 1, 1)
  if (g == "C") elisa_matched[s, ] <- as.numeric(group_means$C)
  else if (g == "L") elisa_matched[s, ] <- as.numeric(group_means$L)
  else if (g == "H") elisa_matched[s, ] <- as.numeric(group_means$H)
  else if (g == "P") elisa_matched[s, ] <- as.numeric(group_means$PP)
}

cat("\n匹配后的ELISA数据:\n")
print(elisa_matched)

# ============================================================
# 3. Z-score标准化
# ============================================================
elisa_zscore <- as.data.frame(scale(elisa_matched))
cat("\nZ-score标准化后的ELISA数据:\n")
print(round(elisa_zscore, 3))
cat("elisa_zscore 列名:", paste(colnames(elisa_zscore), collapse = ", "), "\n")

# ============================================================
# 4. 选择Top 15菌属 + ELISA指标进行相关性分析
# ============================================================
top15_names <- names(sort(colMeans(genus_rel_t), decreasing = TRUE))[1:15]
cat(paste("\nTop 15 菌属:", paste(top15_names, collapse = ", "), "\n"))

# 合并数据
combined <- cbind(genus_rel_t[, top15_names], elisa_zscore)

# 计算Spearman相关性矩阵
cor_matrix <- cor(combined, method = "spearman")

# 提取菌属与ELISA的相关性子矩阵
bacteria_elisa_corr <- cor_matrix[top15_names, colnames(elisa_zscore)]

cat("\n菌属-ELISA Spearman相关性矩阵:\n")
print(round(bacteria_elisa_corr, 3))

# ============================================================
# 5. 绘制相关性热图 (pheatmap)
# ============================================================

# 计算p值矩阵
p_matrix <- matrix(NA, nrow = length(top15_names), ncol = ncol(elisa_zscore))
rownames(p_matrix) <- top15_names
colnames(p_matrix) <- colnames(elisa_zscore)

for (i in 1:length(top15_names)) {
  for (j in 1:ncol(elisa_zscore)) {
    test <- cor.test(combined[, top15_names[i]], combined[, colnames(elisa_zscore)[j]], method = "spearman")
    p_matrix[i, j] <- test$p.value
  }
}

# 标记显著性
sig_matrix <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
sig_matrix[p_matrix < 0.05] <- "*"
sig_matrix[p_matrix < 0.01] <- "**"
sig_matrix[p_matrix < 0.001] <- "***"

# 保存相关性结果到CSV
write.csv(cbind(Genus = rownames(bacteria_elisa_corr), round(bacteria_elisa_corr, 4)), 
          "correlation_results.csv", row.names = FALSE)

# 绘制热图
pdf("correlation_heatmap.pdf", width = 10, height = 12)

# 自定义颜色
breaks <- seq(-1, 1, length.out = 101)
colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                             "#F7F7F7", 
                             "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)

pheatmap(bacteria_elisa_corr,
         color = colors,
         breaks = breaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = round(bacteria_elisa_corr, 2),
         fontsize_number = 10,
         number_color = "black",
         fontsize_row = 11,
         fontsize_col = 11,
         main = "Bacterial Genera vs ELISA Cytokines\n(Spearman Correlation, Z-score normalized)",
         angle_col = 45,
         border_color = "white",
         cellwidth = 60,
         cellheight = 30,
         legend_labels = c("-1", "-0.5", "0", "0.5", "1"),
         legend_breaks = c(-1, -0.5, 0, 0.5, 1))

dev.off()

cat("\nHeatmap written to:", file.path(WORK, "out", "correlation_heatmap.pdf"), "\n")
cat("Correlation table written to:", file.path(WORK, "out", "correlation_results.csv"), "\n")

# ============================================================
# 6. 备选：使用corrplot绘制更美观的热图
# ============================================================
pdf("correlation_heatmap_corrplot.pdf", width = 12, height = 10)

corrplot(bacteria_elisa_corr, 
         method = "color",
         type = "full",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.9,
         number.cex = 0.8,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Bacterial Genera vs ELISA (Spearman r, Z-score normalized)",
         mar = c(0, 0, 2, 0),
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.cex = 1.2,
         pch.col = "black")

dev.off()

cat("corrplot written to:", file.path(WORK, "out", "correlation_heatmap_corrplot.pdf"), "\n")

# ============================================================
# 7. 输出显著相关的结果
# ============================================================
cat("\n========== 显著相关性结果 (p < 0.05) ==========\n")
for (i in 1:nrow(p_matrix)) {
  for (j in 1:ncol(p_matrix)) {
    if (p_matrix[i, j] < 0.05) {
      cat(sprintf("%s vs %s: r = %.3f, p = %.4f %s\n",
                  rownames(p_matrix)[i], colnames(p_matrix)[j],
                  bacteria_elisa_corr[i, j], p_matrix[i, j], sig_matrix[i, j]))
    }
  }
}