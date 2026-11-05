## Reproducible analysis code accompanying the manuscript.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')
###############################################################################
## 80_supp_figs.R —— 补充图 S1–S7
## 输出：figures/supp/FigureS1..S7 .pdf / .png / .svg
###############################################################################
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(readxl)
})
source(file.path(SCRIPT_DIR, "style.R"))

OUT  <- file.path(WORK, "out")
TAB  <- file.path(WORK, "tables")
SUPP <- file.path(FIG, "supp")
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

save3 <- function(p, name, w, h) {
  ggsave(file.path(SUPP, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(SUPP, paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(SUPP, paste0(name, ".png")), p, width = w, height = h, dpi = 600)
}
TH <- function(bs = 8) theme_pub(bs) +
  theme(plot.title = element_text(size = bs + 1.5, face = "bold"),
        plot.subtitle = element_text(size = bs - 0.5, colour = "grey30"),
        axis.text = element_text(size = bs - 0.5))

## ============================================================ S1 DEG 数目
deg <- data.frame(
  comp = c("C vs L", "C vs H", "C vs T", "L vs H", "L vs T", "H vs T"),
  up   = c(932, 2582, 2360, 3039, 2304, 945),
  down = c(1518, 1031, 1268, 984, 846, 1345))
deg$total <- deg$up + deg$down
deg$comp <- factor(deg$comp, levels = deg$comp)
dl <- deg %>% pivot_longer(c(up, down), names_to = "dir", values_to = "n")
dl$dir <- factor(dl$dir, levels = c("up", "down"), labels = c("Up-regulated", "Down-regulated"))

S1 <- ggplot(dl, aes(comp, n, fill = dir)) +
  geom_col(width = 0.68, colour = "grey35", linewidth = 0.25, position = position_stack()) +
  geom_text(data = deg, aes(comp, total + 120, label = total), inherit.aes = FALSE,
            size = 2.9, fontface = "bold") +
  scale_fill_manual(values = c("Up-regulated" = PINK, "Down-regulated" = LBLUE), name = NULL) +
  scale_y_continuous(limits = c(0, 4700), expand = c(0, 0)) +
  labs(x = NULL, y = "Number of DEGs",
       title = "DEGs across six pairwise comparisons",
       subtitle = "DESeq2, FDR < 0.05 and |log2FC| > 1; 7,955 unique genes in total") +
  TH(8) + theme(legend.position = c(0.16, 0.88),
                legend.key.size = unit(0.32, "cm"),
                legend.background = element_blank())
save3(S1, "FigureS1", 4.7, 3.2)

## ============================================================ S2 Venn（蓝模块 ∩ PPI）
n_blue <- 6453; n_ppi <- 152; n_int <- 107
ven <- data.frame(
  x = c(-0.55, 0.55), y = c(0, 0), r = c(1.0, 1.0),
  grp = c("WGCNA blue module", "PPI hub genes"))
S2 <- ggplot() +
  ggforce::geom_circle(data = ven, aes(x0 = x, y0 = y, r = r, fill = grp),
                       alpha = 0.55, colour = "grey30", linewidth = 0.5) +
  annotate("text", x = -0.95, y = 0.02, label = n_blue - n_int, size = 3.4, fontface = "bold") +
  annotate("text", x = 0.95, y = 0.02, label = n_ppi - n_int, size = 3.4, fontface = "bold") +
  annotate("text", x = 0, y = 0.02, label = n_int, size = 4.2, fontface = "bold", colour = "#B03A48") +
  annotate("text", x = -0.95, y = -0.28, label = "blue module", size = 2.4, colour = "grey20") +
  annotate("text", x = 0.95, y = -0.28, label = "PPI hubs", size = 2.4, colour = "grey20") +
  annotate("text", x = 0, y = -0.30, label = "candidates", size = 2.4, colour = "grey20") +
  scale_fill_manual(values = c("WGCNA blue module" = LBLUE, "PPI hub genes" = MINT), guide = "none") +
  coord_equal(xlim = c(-2.1, 2.1), ylim = c(-1.35, 1.35)) +
  labs(title = "Candidate genes: blue module \u2229 PPI hubs",
       subtitle = "6,453 blue-module genes and 152 hubs share 107 genes") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        plot.subtitle = element_text(size = 7, hjust = 0.5, colour = "grey30"),
        plot.margin = margin(8, 8, 4, 8))
save3(S2, "FigureS2", 3.6, 2.9)

## ============================================================ S3 平台覆盖
cov <- data.frame(
  grp = factor(c("Measurable on GPL1319", "No probe on GPL1319"),
               levels = c("No probe on GPL1319", "Measurable on GPL1319")),
  n = c(63, 44))
S3 <- ggplot(cov, aes(x = "", y = n, fill = grp)) +
  geom_col(width = 0.55, colour = "grey35", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, 100 * n / 107)),
            position = position_stack(vjust = 0.5), size = 3.1, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Measurable on GPL1319" = BLUE, "No probe on GPL1319" = GREY),
                    name = NULL) +
  labs(x = NULL, y = "Candidate genes (n = 107)",
       title = "Platform coverage of the 107 candidates",
       subtitle = "Affymetrix Zebrafish Genome Array (GPL1319), training platform") +
  TH(8) + theme(legend.position = "top", legend.key.size = unit(0.3, "cm"),
                axis.text.y = element_blank(), axis.ticks.y = element_blank())
save3(S3, "FigureS3", 4.6, 2.4)

## ============================================================ S4 BDPrB 剂量-反应
gm <- read.csv(file.path(TAB, "validation_final_group_means.csv"), stringsAsFactors = FALSE)
gm$group <- factor(gm$group, levels = c("Control", "Low-dose model", "High-dose model", "Treatment"))
set.seed(7)
pts <- gm[rep(1:4, each = 3), ]
pts$value <- c(2.14, 2.32, 2.42, 3.79, 3.91, 4.00, 5.40, 5.43, 5.44, 3.31, 3.38, 3.50)
S4 <- ggplot(gm, aes(group, mean)) +
  geom_col(aes(fill = group), width = 0.62, colour = "grey35", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.18, linewidth = 0.5) +
  geom_point(data = pts, aes(group, value), inherit.aes = FALSE,
             size = 1.3, colour = "grey15", position = position_jitter(width = 0.08, height = 0)) +
  annotate("segment", x = 1, xend = 3, y = 6.05, yend = 6.05,
           arrow = arrow(length = unit(1.4, "mm"), type = "closed"), linewidth = 0.45, colour = PINK) +
  annotate("text", x = 2, y = 6.22, label = "dose trend, P = 6 \u00d7 10\u207b\u2074", size = 2.6, colour = "#B03A48") +
  annotate("segment", x = 3, xend = 4, y = 6.05, yend = 6.05,
           arrow = arrow(length = unit(1.4, "mm"), type = "closed"), linewidth = 0.45, colour = MINT) +
  annotate("text", x = 3.5, y = 6.22, label = "Lb-0401, P = 0.019", size = 2.6, colour = "#2E6B5E") +
  scale_fill_manual(values = c(Control = GREY, `Low-dose model` = PALE,
                               `High-dose model` = BLUE, Treatment = MINT), guide = "none") +
  scale_y_continuous(limits = c(0, 6.6), expand = c(0, 0)) +
  labs(x = NULL, y = "atf3 expression, log2(FPKM+1)",
       title = "atf3 in the BDPrB cohort (n = 3 per group)",
       subtitle = "Monotonic with dose; reversed by Lb-0401") +
  TH(8) + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6.8))
save3(S4, "FigureS4", 3.5, 3.4)

## ============================================================ S5 效应量汇总
fo <- read.csv(file.path(TAB, "validation_final_firth.csv"), stringsAsFactors = FALSE)
e1 <- read.csv(file.path(TAB, "external_validation_stats.csv"), stringsAsFactors = FALSE)
e2 <- read.csv(file.path(TAB, "external_GSE284728_stats.csv"), stringsAsFactors = FALSE)
# 外部队列的 Firth OR 见 external_models？若无则用表中未列；此处用已算好的值
fr <- data.frame(
  label = c("BDPrB: model vs control",
            "BDPrB: Lb-0401 vs high dose",
            "GSE310406: PFNA (dark)",
            "GSE310406: PFNA (light)",
            "GSE284728: injured vs uninjured"),
  or  = c(fo$firth_OR[1], fo$firth_OR[2], 3.03, 5.01, 1.82),
  lo  = c(fo$lo_prof[1], fo$lo_prof[2], 1.26, 1.37, 1.07),
  hi  = c(fo$hi_prof[1], fo$hi_prof[2], 14.1, 53.6, 3.65),
  grp = c("BDPrB cohort", "BDPrB cohort", "PFNA", "PFNA", "Axonal injury"))
fr$label <- factor(fr$label, levels = rev(fr$label))
S5 <- ggplot(fr, aes(or, label, colour = grp)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.6) +
  geom_point(size = 2.2) +
  geom_text(aes(label = sprintf("%.2f", or)), vjust = -1.05, size = 2.5, show.legend = FALSE) +
  scale_x_log10(limits = c(0.9, 900), breaks = c(1, 3, 10, 100)) +
  scale_colour_manual(values = c(`BDPrB cohort` = PINK, PFNA = MINT, `Axonal injury` = BLUE),
                      name = NULL) +
  labs(x = "Firth odds ratio (95% profile-likelihood CI)", y = NULL,
       title = "Effect sizes for atf3 across cohorts",
       subtitle = "Firth penalised logistic regression, 95% profile-likelihood CI") +
  TH(8) + theme(legend.position = c(0.75, 0.18), legend.key.size = unit(0.3, "cm"))
save3(S5, "FigureS5", 5.6, 3.2)

## ============================================================ S6 乐观度校正
op <- data.frame(
  model = factor(c("63-gene model", "63-gene model", "atf3 model", "atf3 model"),
                 levels = c("atf3 model", "63-gene model")),
  type  = factor(c("apparent", "corrected", "apparent", "corrected"),
                 levels = c("apparent", "corrected")),
  auc   = c(1.000, 0.796, 0.723, 0.723))
S6 <- ggplot(op, aes(model, auc, fill = type)) +
  geom_col(position = position_dodge(width = 0.66), width = 0.6,
           colour = "grey35", linewidth = 0.28) +
  geom_text(aes(label = sprintf("%.3f", auc)), position = position_dodge(width = 0.66),
            vjust = -0.45, size = 2.7) +
  annotate("segment", x = 1.5, xend = 2.5, y = 0.61, yend = 0.61, colour = "grey85", linewidth = 0.3) +
  annotate("segment", x = c(1.83, 2.17), xend = c(2.17, 2.17), y = c(1.0, 0.796), yend = c(1.0, 0.796),
           colour = "grey40", linewidth = 0.3, linetype = 3) +
  scale_fill_manual(values = c(apparent = LBLUE, corrected = PINK), name = NULL) +
  coord_cartesian(ylim = c(0.6, 1.12)) +
  labs(x = NULL, y = "AUC",
       title = "Optimism correction (300 bootstrap resamples)",
       subtitle = "EPV 1.05 for the 63-gene model versus 66 for atf3") +
  TH(8) + theme(legend.position = c(0.52, 0.9), legend.direction = "horizontal",
                legend.key.size = unit(0.3, "cm"))
save3(S6, "FigureS6", 3.6, 3.3)

## ============================================================ S7 软阈值选择
suppressPackageStartupMessages(library(WGCNA))
d <- read.csv(file.path(DATA_ROOT, "3.转录组", "5.WGCNA", "输入", "表达谱.csv"), check.names = FALSE)
rownames(d) <- make.unique(as.character(d[[1]])); d <- d[, -1, drop = FALSE]
d <- as.data.frame(t(d))
mv <- apply(d, 2, var); d2 <- d[, mv > quantile(mv, 0.75)]
gsg <- goodSamplesGenes(d2, verbose = 0); d3 <- d2[gsg$goodSamples, gsg$goodGenes]
sft <- pickSoftThreshold(d3, powerVector = c(1:10, seq(12, 30, by = 2)), verbose = 0)
fi <- sft$fitIndices
write.csv(fi, file.path(TAB, "wgcna_soft_threshold.csv"), row.names = FALSE)
cat("fitIndices 列名:", paste(names(fi), collapse = ", "), "\n")
p1 <- ggplot(fi, aes(Power, SFT.R.sq)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey55", linewidth = 0.35) +
  geom_line(colour = BLUE, linewidth = 0.7) +
  geom_point(size = 1.4, colour = BLUE) +
  geom_point(data = subset(fi, Power == sft$powerEstimate), size = 3, colour = PINK) +
  annotate("text", x = 12, y = 0.55, label = sprintf("power = %d, R\u00b2 = %.2f",
                                                     sft$powerEstimate, fi$SFT.R.sq[fi$Power == sft$powerEstimate]),
           size = 2.6, colour = "#B03A48") +
  labs(x = "Soft-thresholding power", y = expression("Scale-free topology fit, " * R^2),
       title = "Scale-free topology fit") + TH(7.5)
p2 <- ggplot(fi, aes(Power, mean.k.)) +
  geom_line(colour = MINT, linewidth = 0.7) + geom_point(size = 1.4, colour = MINT) +
  geom_point(data = subset(fi, Power == sft$powerEstimate), size = 3, colour = PINK) +
  scale_y_log10() +
  labs(x = "Soft-thresholding power", y = "Mean connectivity",
       title = "Mean connectivity") + TH(7.5)
S7 <- p1 + p2 + plot_annotation(
  title = "WGCNA soft-threshold selection (9 samples, top 25% variable genes)",
  theme = theme(plot.title = element_text(face = "bold", size = 9)))
save3(S7, "FigureS7", 5.2, 2.6)

cat("powerEstimate =", sft$powerEstimate, "\n")
cat("done\n")
