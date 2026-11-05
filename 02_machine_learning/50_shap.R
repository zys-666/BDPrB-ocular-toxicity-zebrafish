source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 50_shap.R — 最优组合模型的 SHAP（蜂巢图 + 极坐标图），突出 atf3
##   最优组合：ml130 三折平均 AUC 最高者（默认 glmBoost+Lasso）
###############################################################################
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(fastshap)
})
source(file.path(SCRIPT_DIR, "mime_functions.R"))
source(file.path(SCRIPT_DIR, "style.R"))
THEME <- function(bs = 8) theme_pub(bs) +
  theme(plot.title = element_text(size = bs, face = "bold"),
        axis.title = element_text(size = bs - 0.5, face = "bold"),
        axis.text = element_text(size = bs - 1.2),
        panel.grid.minor = element_blank())
OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")
FIG <- file.path(WORK, "figures")

d <- readRDS(file.path(OUT, "prep.rds"))
X <- d$X; y <- d$y; genes <- d$genes
classVar <- "Type"

## ---- 最优组合 ----
perf <- read.csv(file.path(TAB, "ml130_combo_performance.csv"), stringsAsFactors = FALSE)
best <- perf$method[which.max(perf$auc_mean)]
cat("最优组合:", best, " AUC =", round(max(perf$auc_mean), 3), "\n")
stages <- strsplit(gsub(" ?\\+ ?", "+", best), "+", fixed = TRUE)[[1]]

## ---- 在全体训练集上重跑该组合 ----
lab <- data.frame(Type = y)
vars <- colnames(X)
for (i in seq_along(stages)) {
  if (i < length(stages)) {
    v <- suppressWarnings(RunML(stages[i], X[, vars, drop = FALSE], lab,
                                mode = "Variable", classVar = classVar))
    vars <- intersect(v, vars)
    cat("阶段", i, stages[i], "保留", length(vars), "个基因\n")
  }
}
fit <- suppressWarnings(RunML(stages[length(stages)], X[, vars, drop = FALSE], lab,
                              mode = "Model", classVar = classVar))
cat("最终基因数:", length(vars), "；含 atf3:", "atf3" %in% vars, "\n")
cat("基因:", paste(sort(vars), collapse = ", "), "\n")
writeLines(sort(vars), file.path(TAB, "best_combo_genes.txt"))

pred_fun <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  as.numeric(CalPredictScore(fit, newdata[, vars, drop = FALSE]))
}

## ---- SHAP ----
set.seed(42)
ns <- min(60, nrow(X))
idx <- sample(nrow(X), ns)
Xbg <- as.data.frame(X)[, vars, drop = FALSE]   # 只用模型实际使用的基因
Xnew <- as.data.frame(X)[idx, vars, drop = FALSE]
sh <- fastshap::explain(object = fit, X = Xbg, newdata = Xnew,
                        pred_wrapper = pred_fun, nsim = 100, adjust = TRUE)
sh <- as.matrix(sh); colnames(sh) <- vars
saveRDS(list(shap = sh, X = Xnew, vars = vars, best = best),
        file.path(OUT, "shap_best.rds"))

imp <- data.frame(gene = vars, mean_abs = colMeans(abs(sh))) %>% arrange(desc(mean_abs))
imp$rank <- seq_len(nrow(imp))
write.csv(imp, file.path(TAB, "shap_importance_best.csv"), row.names = FALSE)
print(head(imp, 10), digits = 3, row.names = FALSE)
cat("\natf3 SHAP 排名:", which(imp$gene == "atf3"), "/", nrow(imp), "\n")

## ---- 蜂巢图 ----
## 直接按位置展开（SHAP 矩阵与 Xnew 的行列一一对应），避免 join 的命名问题
long <- data.frame(
  id    = rep(seq_len(nrow(sh)), times = ncol(sh)),
  gene  = rep(colnames(sh), each = nrow(sh)),
  shap  = as.vector(sh),
  value = as.vector(as.matrix(Xnew)),
  stringsAsFactors = FALSE
)
long$gene <- factor(long$gene, levels = rev(imp$gene))
long$value <- as.numeric(long$value)
long$value_z <- ave(long$value, long$gene, FUN = function(z) (z - min(z)) / (max(z) - min(z) + 1e-12))

p_bee <- ggplot(long, aes(shap, gene, colour = value_z)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_jitter(height = 0.22, size = 0.85, alpha = 0.9) +
  scale_colour_gradientn(colours = GRAD_BWR, limits = c(0, 1), name = "Feature\nvalue",
                         breaks = c(0, 1), labels = c("Low", "High")) +
  labs(x = "SHAP value (impact on predicted probability)", y = NULL,
       title = paste0("SHAP beeswarm — ", best)) +
  THEME(8)
ggsave(file.path(FIG, "fig50_SHAP_beeswarm.pdf"), p_bee, width = 4.6, height = 4.2)
ggsave(file.path(FIG, "fig50_SHAP_beeswarm.svg"), p_bee, width = 4.6, height = 4.2, device = svglite::svglite)
ggsave(file.path(FIG, "fig50_SHAP_beeswarm.png"), p_bee, width = 4.6, height = 4.2, dpi = 500)

## ---- 极坐标蜂巢图（内圈玫瑰 = mean|SHAP|，外环 = 单样本 SHAP）----
df <- long
df$gene <- factor(df$gene, levels = imp$gene)
n_f <- length(levels(df$gene))
ang <- data.frame(gene = factor(levels(df$gene), levels = levels(df$gene)),
                  theta = seq(0, 2*pi - 2*pi/n_f, length.out = n_f))
df <- merge(df, ang, by = "gene")
maxabs <- max(abs(df$shap))
offset <- maxabs * 1.15                       # SHAP=0 所在的圆
df$r <- offset + df$shap
df$theta_j <- df$theta + as.numeric(ave(df$shap, df$gene, FUN = rank)) *
  (2*pi/n_f*0.30)/ave(df$shap, df$gene, FUN = length)
imp2 <- merge(ang, imp, by = "gene")
imp2$rmean <- imp2$mean_abs / max(imp2$mean_abs) * maxabs * 0.85

wedge <- do.call(rbind, lapply(seq_len(nrow(imp2)), function(i) {
  th <- seq(imp2$theta[i] - 2*pi/n_f*0.42, imp2$theta[i] + 2*pi/n_f*0.42, length.out = 40)
  data.frame(g = i, x = c(0, imp2$rmean[i]*sin(th), 0), y = c(0, imp2$rmean[i]*cos(th), 0))
}))
cir <- function(r, g) data.frame(g = g, x = r*sin(seq(0, 2*pi, length.out = 200)),
                                 y = r*cos(seq(0, 2*pi, length.out = 200)))
p_polar <- ggplot() +
  geom_polygon(data = wedge, aes(x, y, group = g), fill = "#F4A7AE", alpha = 0.88) +
  geom_path(data = cir(offset, "z"), aes(x, y, group = g), colour = "grey40", linewidth = 0.4) +
  geom_path(data = cir(offset + maxabs*1.06, "o"), aes(x, y, group = g),
            colour = "grey70", linewidth = 0.3, linetype = 3) +
  geom_point(data = df, aes(r*sin(theta_j), r*cos(theta_j), colour = value_z),
             size = 1.0, alpha = 0.9) +
  scale_colour_gradientn(colours = GRAD_BWR, limits = c(0, 1), name = "Feature\nvalue",
                         breaks = c(0, 1), labels = c("Low", "High")) +
  coord_equal(clip = "off") +
  labs(title = paste0("SHAP polar — ", best),
       subtitle = "centre = mean |SHAP|; ring = per-sample SHAP (grey circle = 0)") +
  theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.1, "cm"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 6.8))
lab_df <- merge(ang, data.frame(gene = imp$gene, rank = imp$rank), by = "gene")
rr <- offset + maxabs * 1.32
lab_df$x <- rr * sin(lab_df$theta); lab_df$y <- rr * cos(lab_df$theta)
lab_df$face <- ifelse(as.character(lab_df$gene) == "atf3", "bold", "plain")
p_polar <- p_polar + geom_text(data = lab_df, aes(x, y, label = gene, fontface = face),
                               size = 2.1)
ggsave(file.path(FIG, "fig50_SHAP_polar.pdf"), p_polar, width = 5.4, height = 5.6)
ggsave(file.path(FIG, "fig50_SHAP_polar.svg"), p_polar, width = 5.4, height = 5.6, device = svglite::svglite)
ggsave(file.path(FIG, "fig50_SHAP_polar.png"), p_polar, width = 5.4, height = 5.6, dpi = 500)
cat("done\n")
