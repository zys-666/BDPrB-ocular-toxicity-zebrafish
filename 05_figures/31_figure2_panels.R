## Reproducible analysis code accompanying the manuscript.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')
###############################################################################
## 31_fig2_final.R — 最终 Fig.2（4×3 = 12 面板）+ 补充材料面板（S1–S5）
##   主图：A 130 组合（带组合名）| B 基因频率 | C 单基因 | D ROC
##         E 校准 | F DCA | G 置换 | H 稳定性
##         I SHAP 蜂巢 | J SHAP 极坐标 | K PFNA | L ONT
##   补充：S1 候选基因流程(107) | S2 平台可测率 | S3 BDPrB 四组 | S4 Firth 森林 | S5 乐观度
###############################################################################
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(pROC); library(patchwork)
})
source(file.path(SCRIPT_DIR, "15_firth_fun.R"))
source(file.path(SCRIPT_DIR, "style.R"))

OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")
FIG <- file.path(WORK, "figures")
PD  <- file.path(FIG, "fig2_panels"); dir.create(PD, showWarnings = FALSE, recursive = TRUE)

auc_safe <- function(y, p) {
  if (length(unique(y)) < 2 || length(unique(p)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
}
## 一次导出 PDF（矢量）/ SVG（可用 AI/Inkscape 直接编辑）/ PNG（预览）
save_all <- function(p, base, w, h) {
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
  ggsave(paste0(base, ".svg"), p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 600)
}
THEME <- function(bs = 7.2) theme_pub(bs) +
  theme(plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 1.5)),
        plot.subtitle = element_text(size = 6.0, margin = margin(b = 1.5)),
        axis.text = element_text(size = 5.8, colour = "grey15"),
        axis.title = element_text(size = 6.2, face = "bold"),
        plot.margin = margin(3, 3, 2, 2))

## ================================================== A 130 组合（带名字）
z <- readRDS(file.path(OUT, "ml130_preds.rds"))
pm <- z$preds %>% group_by(method, fold) %>%
  summarise(auc = auc_safe(y, p), .groups = "drop") %>%
  group_by(method) %>% summarise(m = mean(auc, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()
pm <- pm[order(-pm$m), ]
top <- head(pm, 22)
top$method <- factor(top$method, levels = rev(top$method))
pd <- z$preds %>% filter(method %in% top$method)
pd$method <- factor(pd$method, levels = levels(top$method))
pd <- pd %>% group_by(method) %>% mutate(auc = auc_safe(y, p)) %>% ungroup()

pA <- ggplot(top, aes(m, method)) +
  geom_segment(data = top, aes(x = 0.5, xend = m, yend = method),
               colour = "grey88", linewidth = 0.3) +
  geom_point(data = top, aes(size = "mean"), colour = "#F4A7AE") +
  geom_point(data = pd, aes(x = auc, y = method), size = 0.7, alpha = 0.7,
             colour = "#5878A6", position = position_jitter(height = 0.12)) +
  geom_vline(xintercept = 0.721, linetype = 2, colour = "#F4A7AE", linewidth = 0.4) +
  scale_size_manual(values = c(mean = 1.9), labels = "mean of 3 folds", name = NULL) +
  scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8)) + coord_cartesian(xlim = c(0.45, 0.85)) +
  labs(x = "Cross-validated AUC", y = NULL, title = "130 combinations",
       subtitle = "top 22 of 130\ndashed = atf3 (0.721)") +
  THEME(6.6) + theme(axis.text.y = element_text(size = 4.3),
                     legend.position = "none")
pA <- tag_panel(pA, "A")

## ================================================== B 基因入选频率
gf <- read.csv(file.path(TAB, "ml130_gene_frequency.csv"), stringsAsFactors = FALSE)
g12 <- head(gf[order(-gf$n), ], 12); g12$gene <- factor(g12$gene, levels = rev(g12$gene))
pB <- ggplot(g12, aes(n, gene, fill = gene == "atf3")) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.22) +
  scale_fill_manual(values = c(`TRUE` = "#F4A7AE", `FALSE` = "#A7C7E7"), guide = "none") +
  geom_text(aes(label = n), hjust = -0.15, size = 1.8) +
  scale_x_continuous(limits = c(0, 470), expand = c(0, 0)) +
  labs(x = "Selection count", y = NULL, title = "Gene selection",
       subtitle = "atf3: 388/390 (99.5%)") +
  THEME()
pB <- tag_panel(pB, "B")

## ================================================== C 单基因筛查
uni <- read.csv(file.path(TAB, "univariate_auc.csv"), stringsAsFactors = FALSE)
uni <- uni[order(uni$auc_mean), ]; uni$rank <- seq_len(nrow(uni))
pC <- ggplot(uni, aes(auc_mean, rank)) +
  geom_segment(aes(xend = 0.5, yend = rank), colour = "grey90", linewidth = 0.22) +
  geom_point(aes(colour = gene == "atf3", size = gene == "atf3")) +
  scale_colour_manual(values = c(`TRUE` = "#F4A7AE", `FALSE` = "#A7C7E7"), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.8, `FALSE` = 0.65), guide = "none") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey55", linewidth = 0.35) +
  annotate("text", x = 0.725, y = 60, label = "atf3", size = 2, fontface = "bold",
           colour = "#F4A7AE", hjust = 0) +
  scale_x_continuous(breaks = c(0.5, 0.6, 0.7)) + coord_cartesian(xlim = c(0.48, 0.80)) +
  labs(x = "Univariate CV AUC", y = "63 genes (ranked)",
       title = "Single-gene AUC", subtitle = "atf3 ranks 1st") +
  THEME() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
pC <- tag_panel(pC, "C")

## ================================================== D ROC
am <- readRDS(file.path(OUT, "atf3_model.rds"))
r <- roc(am$y, am$p_loo, quiet = TRUE, direction = "<")
ci <- as.numeric(ci.auc(r))
pD <- ggplot(data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities), aes(fpr, tpr)) +
  geom_abline(linetype = 2, colour = "grey65", linewidth = 0.35) +
  geom_step(colour = "#F4A7AE", linewidth = 1.1) +
  annotate("text", x = 0.22, y = 0.06, hjust = 0, size = 1.8, fontface = "bold",
           colour = "#F4A7AE",
           label = sprintf("AUC %.3f (%.3f-%.3f)\nP = %.3f", as.numeric(auc(r)), ci[1], ci[3], am$p_perm)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, .5, 1)) + scale_y_continuous(breaks = c(0, .5, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "atf3 ROC") +
  THEME()
pD <- tag_panel(pD, "D")

## ================================================== E 校准
cal_bin <- function(y, p, k = 7) {
  q <- unique(quantile(p, probs = seq(0, 1, length.out = k + 1)))
  g <- cut(p, breaks = q, include.lowest = TRUE)
  data.frame(mid = as.numeric(tapply(p, g, mean)), obs = as.numeric(tapply(y, g, mean)))
}
cb1 <- cal_bin(am$y, am$p_loo); cb1$model <- "atf3 (1 gene)"
bestm <- pm$method[1]
pbz <- z$preds %>% filter(method == bestm)
cb2 <- cal_bin(pbz$y, pbz$p); cb2$model <- "best combo"
pE <- ggplot(rbind(cb1, cb2), aes(mid, obs, colour = model)) +
  geom_abline(linetype = 2, colour = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.75) + geom_point(size = 1.2) +
  scale_colour_manual(values = c("#F4A7AE", "#5878A6"), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability", y = "Observed", title = "Calibration") +
  THEME() + theme(legend.position = c(0.66, 0.11), legend.text = element_text(size = 4.6),
                  legend.key.size = unit(0.2, "cm"))
pE <- tag_panel(pE, "E")

## ================================================== F DCA
dca <- read.csv(file.path(TAB, "dca_training.csv"))
pF <- ggplot(dca, aes(threshold)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey55", linewidth = 0.35) +
  geom_line(aes(y = nb_all), linetype = 2, colour = "grey45", linewidth = 0.55) +
  geom_line(aes(y = nb_model), colour = "#F4A7AE", linewidth = 0.85) +
  coord_cartesian(xlim = c(0, 0.85), ylim = c(-0.02, 0.5)) +
  annotate("text", x = 0.55, y = 0.33, label = "atf3 model", colour = "#F4A7AE",
           size = 1.8, fontface = "bold") +
  annotate("text", x = 0.60, y = 0.03, label = "treat all", colour = "grey40", size = 1.8) +
  labs(x = "Threshold probability", y = "Net benefit", title = "Decision curve") +
  THEME()
pF <- tag_panel(pF, "F")

## ================================================== G 置换
pG <- ggplot(data.frame(null = am$null), aes(null)) +
  geom_histogram(bins = 24, fill = "#A7C7E7", colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = as.numeric(auc(r)), colour = "#F4A7AE", linewidth = 0.7) +
  annotate("text", x = as.numeric(auc(r)), y = Inf, vjust = 1.7, hjust = 1.06,
           size = 1.8, colour = "#F4A7AE", fontface = "bold",
           label = sprintf("obs %.3f", as.numeric(auc(r)))) +
  labs(x = "AUC (permuted labels)", y = "Count", title = "Permutation null", subtitle = "B = 1000") +
  THEME()
pG <- tag_panel(pG, "G")

## ================================================== H 稳定性（130 组合 × 20 次重采样）
rep <- readRDS(file.path(OUT, "ml130_preds_rep.rds"))
n_meth <- length(unique(rep$vars$method))
st <- rep$vars %>% group_by(fold, gene) %>% summarise(k = n_distinct(method), .groups = "drop") %>%
  mutate(frac = k / n_meth)
ordg <- st %>% group_by(gene) %>% summarise(m = mean(frac), .groups = "drop") %>%
  arrange(-m) %>% head(12)
st <- st %>% filter(gene %in% ordg$gene)
st$gene <- factor(st$gene, levels = rev(ordg$gene))
st$fold <- factor(st$fold, levels = unique(st$fold))
pH <- ggplot(st, aes(fold, gene, fill = frac)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(colours = GRAD_BWR, limits = c(0, 1), name = "Fraction\nof 130",
                       breaks = c(0, .5, 1)) +
  labs(x = NULL, y = NULL, title = "Selection stability",
       subtitle = "4 repeats x 5 folds\nx 130 combinations") +
  THEME(6.8) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3.4),
                     axis.text.y = element_text(size = 4.4),
                     legend.key.height = unit(0.25, "cm"), legend.key.width = unit(0.16, "cm"))
pH <- tag_panel(pH, "H")

## ================================================== I SHAP 蜂巢
shz <- readRDS(file.path(OUT, "shap_best.rds"))
sh <- shz$shap; Xn <- shz$X; vars <- shz$vars
imp <- data.frame(gene = vars, mean_abs = colMeans(abs(sh)))
imp <- imp[order(-imp$mean_abs), ]
long <- data.frame(id = rep(seq_len(nrow(sh)), times = ncol(sh)),
                   gene = rep(colnames(sh), each = nrow(sh)),
                   shap = as.vector(sh), value = as.vector(as.matrix(Xn)),
                   stringsAsFactors = FALSE)
long$gene <- factor(long$gene, levels = rev(imp$gene))
long$vz <- ave(long$value, long$gene, FUN = function(z) (z - min(z)) / (max(z) - min(z) + 1e-12))
pI <- ggplot(long, aes(shap, gene, colour = vz)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
  geom_jitter(height = 0.2, size = 0.7, alpha = 0.9) +
  scale_colour_gradientn(colours = GRAD_BWR, limits = c(0, 1), name = "Feature\nvalue",
                         breaks = c(0, 1), labels = c("Low", "High")) +
  labs(x = "SHAP value", y = NULL, title = paste0("SHAP (", bestm, ")"),
       subtitle = "atf3 ranks 3rd of 14 genes") +
  THEME(6.8) + theme(axis.text.y = element_text(size = 4.6),
                     legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.18, "cm"))
pI <- tag_panel(pI, "I")

## ================================================== J SHAP 极坐标
n_f <- length(imp$gene)
ang <- data.frame(gene = factor(imp$gene, levels = imp$gene),
                  theta = seq(0, 2*pi - 2*pi/n_f, length.out = n_f))
df <- merge(long, ang, by = "gene")
maxabs <- max(abs(df$shap)); offset <- maxabs * 1.15
df$r <- offset + df$shap
df$theta_j <- df$theta + (as.numeric(ave(df$shap, df$gene, FUN = rank)) - 0.5) *
  (2*pi/n_f*0.55)/ave(df$shap, df$gene, FUN = length)
imp2 <- merge(ang, imp, by = "gene"); imp2$rmean <- imp2$mean_abs / max(imp2$mean_abs) * maxabs * 0.8
wedge <- do.call(rbind, lapply(seq_len(nrow(imp2)), function(i) {
  th <- seq(imp2$theta[i] - 2*pi/n_f*0.40, imp2$theta[i] + 2*pi/n_f*0.40, length.out = 40)
  data.frame(g = i, x = c(0, imp2$rmean[i]*sin(th), 0), y = c(0, imp2$rmean[i]*cos(th), 0))
}))
circ <- function(r, g) data.frame(g = g, x = r*sin(seq(0, 2*pi, length.out = 200)),
                                  y = r*cos(seq(0, 2*pi, length.out = 200)))
lab <- merge(ang, imp, by = "gene"); rr <- offset + maxabs * 1.42
lab$x <- rr*sin(lab$theta); lab$y <- rr*cos(lab$theta)
pJ <- ggplot() +
  geom_polygon(data = wedge, aes(x, y, group = g), fill = "#F4A7AE", alpha = 0.88) +
  geom_path(data = circ(offset, "z"), aes(x, y, group = g), colour = "grey45", linewidth = 0.35) +
  geom_point(data = df, aes(r*sin(theta_j), r*cos(theta_j), colour = vz), size = 0.8, alpha = 0.9) +
  scale_colour_gradientn(colours = GRAD_BWR, limits = c(0, 1), guide = "none") +
  geom_text(data = lab, aes(x, y, label = gene), size = 1.55,
            fontface = ifelse(as.character(lab$gene) == "atf3", "bold", "plain")) +
  coord_equal(clip = "off") +
  labs(title = "SHAP polar", subtitle = "rose = mean |SHAP|; ring = samples") +
  theme_void(base_size = 6.8) +
  theme(plot.title = element_text(size = 7.2, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 5.8, hjust = 0.5),
        plot.margin = margin(3, 3, 2, 3))
pJ <- tag_panel(pJ, "J")

## ================================================== K PFNA
e1 <- read.csv(file.path(TAB, "external_GSE310406_atf3.csv"), stringsAsFactors = FALSE) %>%
  filter(tissue == "Eye")
e1$group <- factor(e1$group, levels = c("Control", "PFNA"))
pK <- ggplot(e1, aes(group, expr_log2, fill = group)) +
  geom_boxplot(width = .55, alpha = .55, outlier.shape = NA, linewidth = 0.32) +
  geom_jitter(width = .09, size = 1, alpha = .9, colour = "grey20") +
  facet_wrap(~light) +
  scale_fill_manual(values = c(Control = "#A7C7E7", PFNA = "#F4A7AE"), guide = "none") +
  labs(x = NULL, y = "atf3 log2(FPKM+1)", title = "GSE310406 (PFNA)",
       subtitle = "AUC 1.000; P = 0.029") +
  THEME() + theme(axis.text.x = element_text(size = 5.2))
pK <- tag_panel(pK, "K")

## ================================================== L ONT
e2 <- read.csv(file.path(TAB, "external_GSE284728_atf3.csv"), stringsAsFactors = FALSE)
e2$group <- factor(e2$group, levels = c("Uninjured", "Injured"))
pL <- ggplot(e2, aes(group, log2cpm, fill = group)) +
  geom_boxplot(width = .55, alpha = .55, outlier.shape = NA, linewidth = 0.32) +
  geom_jitter(width = .09, size = 1, alpha = .9, colour = "grey20") +
  facet_wrap(~tp, nrow = 1) +
  scale_fill_manual(values = c(Uninjured = "#A7C7E7", Injured = "#F4A7AE"), guide = "none") +
  labs(x = NULL, y = "atf3 log2(CPM+1)", title = "GSE284728 (ONT)",
       subtitle = "3.9\u00d7; AUC 0.802; P = 0.032") +
  THEME() + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 4.6))
pL <- tag_panel(pL, "L")

## ================================================== 主图（11 面板：去掉蜂巢图，极坐标占两格）
## A 130组合 | B 基因频率 | C 单基因 | D ROC
## E 校准 | F DCA | G 置换 | H 稳定性
## I PFNA | J ONT | K SHAP 极坐标（跨 2 格）
pK <- tag_panel(pK, "I")
pL <- tag_panel(pL, "J")
pJ <- tag_panel(pJ, "K")
row3 <- pK + pL + pJ + plot_layout(widths = c(1, 1, 2))
main <- (pA | pB | pC | pD) / (pE | pF | pG | pH) / row3 +
  plot_annotation(
    title = "Figure 2 | 130 machine-learning combinations identify atf3 as the core gene",
    theme = theme(plot.title = element_text(face = "bold", size = 7.6)))
save_all(main, file.path(FIG, "Fig2_main"), 7.1, 6.1)

## ================================================== 补充面板 S1–S5
bx <- data.frame(x = c(0.26, 0.26, 0.76), y = c(0.74, 0.26, 0.50), w = 0.44, h = 0.24,
                 lab = c("WGCNA blue\n6,452 genes", "PPI hubs\n(CytoHubba) 152", "107 genes"),
                 fill = c("#A7C7E7", "#A6D0C6", "#FBD9C9"))
S1 <- ggplot() +
  geom_rect(data = bx, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
                           fill = fill), colour = "grey30", linewidth = 0.4) +
  geom_text(data = bx, aes(x, y, label = lab), size = 2.1, fontface = "bold", lineheight = .95) +
  annotate("segment", x = 0.47, xend = 0.53, y = 0.50, yend = 0.50,
           arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.5, colour = "grey30") +
  scale_fill_identity() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(title = "S6 | Candidate genes", subtitle = "WGCNA \u2229 PPI = 107") +
  THEME() + theme(axis.text = element_blank(), axis.title = element_blank(),
                  axis.ticks = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.border = element_blank())

cov <- data.frame(dataset = factor(c("BDPrB (RNA-seq)", "GSE310406 (RNA-seq)",
                                     "GSE284728 (RNA-seq)", "GSE234646 (array)"),
                                   levels = rev(c("BDPrB (RNA-seq)", "GSE310406 (RNA-seq)",
                                                  "GSE284728 (RNA-seq)", "GSE234646 (array)"))),
                  pct = c(100, 96.3, 95.3, 58.9), n = c(107, 103, 102, 63))
S2 <- ggplot(cov, aes(pct, dataset, fill = pct)) +
  geom_col(width = 0.62, colour = "grey30", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%d/107", n)), hjust = -0.15, size = 2) +
  scale_fill_gradientn(colours = GRAD_BWO, limits = c(50, 100), guide = "none") +
  scale_x_continuous(limits = c(0, 125), expand = c(0, 0)) +
  labs(x = "Candidates measurable (%)", y = NULL, title = "S7 | Platform coverage") +
  THEME()

v <- read.csv(file.path(DATA_ROOT, "2.machine_learning", "validation", "expression_norm_log2.csv"),
              row.names = 1, check.names = FALSE)
v <- v[, colnames(v) != "", drop = FALSE]
grp <- c(C1="Control", C2="Control", C3="Control", L1="Low-dose", L2="Low-dose", L3="Low-dose",
         H1="High-dose", H2="High-dose", H3="High-dose",
         PP1="Treatment", PP2="Treatment", PP3="Treatment")
e3 <- data.frame(group = factor(grp[colnames(v)], levels = c("Control", "Low-dose", "High-dose", "Treatment")),
                 val = as.numeric(v["atf3", ]))
S3 <- ggplot(e3, aes(group, val, fill = group)) +
  geom_boxplot(width = .55, alpha = .55, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = .09, size = 1.4, alpha = .95, colour = "grey20") +
  scale_fill_manual(values = c(Control = "#A7C7E7", `Low-dose` = "#CFEAF6",
                               `High-dose` = "#F4A7AE", Treatment = "#A6D0C6"), guide = "none") +
  labs(x = NULL, y = "atf3 log2(FPKM+1)", title = "S8 | BDPrB cohort (n = 3/group)",
       subtitle = "dose P = 6e-4; rescue P = 0.019") +
  THEME() + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5))

forest <- data.frame(
  cohort = factor(c("BDPrB: model vs control", "BDPrB: Lb-0401 vs high dose",
                    "GSE310406: PFNA (dark)", "GSE310406: PFNA (light)",
                    "GSE284728: injured vs uninjured"),
                  levels = rev(c("BDPrB: model vs control", "BDPrB: Lb-0401 vs high dose",
                                 "GSE310406: PFNA (dark)", "GSE310406: PFNA (light)",
                                 "GSE284728: injured vs uninjured"))),
  or = c(9.49, 6.69, 3.03, 5.01, 1.82), lo = c(1.49, 1.29, 1.26, 1.37, 1.07),
  hi = c(541, 117, 14.1, 53.6, 3.65),
  grp = c("BDPrB", "BDPrB", "PFNA", "PFNA", "ONT"))
S4 <- ggplot(forest, aes(or, cohort, colour = grp)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.55) +
  geom_point(size = 2) +
  scale_x_log10(breaks = c(1, 3, 10, 100), limits = c(0.9, 900)) +
  scale_colour_manual(values = c(BDPrB = "#F4A7AE", PFNA = "#A6D0C6", ONT = "#5878A6"), guide = "none") +
  labs(x = "Firth OR (95% CI)", y = NULL, title = "S9 | Effect sizes") +
  THEME()

ss <- data.frame(what = factor(c("63-gene model", "63-gene model", "atf3 model", "atf3 model"),
                               levels = c("atf3 model", "63-gene model")),
                 metric = c("apparent", "corrected", "apparent", "corrected"),
                 auc = c(1.000, 0.796, 0.723, 0.723))
S5 <- ggplot(ss, aes(what, auc, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           colour = "grey30", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", auc)), position = position_dodge(width = 0.7),
            vjust = -0.35, size = 1.9) +
  scale_fill_manual(values = c(apparent = "#A7C7E7", corrected = "#F4A7AE"), name = NULL) +
  coord_cartesian(ylim = c(0.6, 1.12)) +
  labs(x = NULL, y = "AUC", title = "S10 | Optimism correction", subtitle = "EPV 66 vs 1.05") +
  THEME() + theme(legend.position = c(0.5, 0.12), legend.direction = "horizontal",
                  legend.text = element_text(size = 5), legend.key.size = unit(0.22, "cm"))

sup <- (S1 | S2) / (S3 | S4) / (S5 | ggplot() + theme_void()) +
  plot_annotation(title = "Supplementary Figure | Candidate genes, platform coverage and effect sizes",
                  theme = theme(plot.title = element_text(face = "bold", size = 7.6)))
save_all(sup, file.path(FIG, "FigS_S6_S10"), 5.2, 6.0)

## ================================================== 单面板导出（PDF / SVG / PNG）
## 最终图内的面板编号：A..K（I=PFNA, J=ONT, K=SHAP 极坐标）
final <- list(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF, G = pG, H = pH,
              I = pK, J = pL, K = pJ)
for (nm in names(final)) {
  save_all(final[[nm]], file.path(PD, nm), 2.15, 1.9)
}
## 补充材料面板 S1–S5
for (i in 1:5) {
  save_all(get(paste0("S", i)), file.path(PD, paste0("S", i + 5)), 2.2, 2.0)
}
cat("done\n")
