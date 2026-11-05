source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 12_atf3_final.R — atf3 单基因模型的最终评估
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")

d <- readRDS(file.path(OUT, "prep.rds"))
X <- d$X; y <- d$y; meta <- d$meta; genes <- d$genes
auc_safe <- function(y, p) {
  if (length(unique(y)) < 2 || length(unique(p)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
}
g <- "atf3"; x <- X[, g]
fit1 <- function(idx) glm(y ~ x, data = data.frame(y = y[idx], x = x[idx]), family = binomial)

## ------------------------------------------------ 1. 留一 Run
runs <- sort(unique(meta$run))
p_loo <- rep(NA_real_, length(y))
for (r0 in runs) {
  tr <- which(meta$run != r0); te <- which(meta$run == r0)
  fit <- suppressWarnings(fit1(tr))
  p_loo[te] <- predict(fit, newdata = data.frame(x = x[te]), type = "response")
}
r_loo <- roc(y, p_loo, quiet = TRUE, direction = "<")
ci_loo <- as.numeric(ci.auc(r_loo)); auc_loo <- as.numeric(auc(r_loo))

## ------------------------------------------------ 2. 重复 5 折
set.seed(2024); REP <- 20; K <- 5
cv <- numeric(REP); p_cv_last <- NULL
for (r in seq_len(REP)) {
  f <- integer(length(y))
  for (v in c(0, 1)) { i <- which(y == v); f[i] <- sample(rep(1:K, length.out = length(i))) }
  pp <- rep(NA_real_, length(y))
  for (k in 1:K) {
    tr <- which(f != k); te <- which(f == k)
    pp[te] <- predict(suppressWarnings(fit1(tr)), newdata = data.frame(x = x[te]), type = "response")
  }
  cv[r] <- auc_safe(y, pp); if (r == REP) p_cv_last <- pp
}

## ------------------------------------------------ 3. 置换（留一 Run，B=1000）
set.seed(11); B <- 1000; null <- numeric(B)
for (b in seq_len(B)) {
  yb <- sample(y); pp <- rep(NA_real_, length(y))
  for (r0 in runs) {
    tr <- which(meta$run != r0); te <- which(meta$run == r0)
    fb <- suppressWarnings(glm(yb[tr] ~ x[tr], family = binomial))
    pp[te] <- predict(fb, newdata = data.frame(x = x[te]), type = "response")
  }
  null[b] <- auc_safe(yb, pp)
}
p_perm <- (1 + sum(null >= auc_loo, na.rm = TRUE)) / (B + 1)

## ------------------------------------------------ 4. 校准与阈值
lg <- log(p_loo / (1 - p_loo)); cal <- suppressWarnings(glm(y ~ lg, family = binomial))
co <- coords(r_loo, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
brier <- mean((p_loo - y)^2)

## ------------------------------------------------ 5. 分层一致性（两个范式各看一次）
per_ser <- sapply(c("BP", "R"), function(s) {
  i <- which(meta$series == s); xs <- x[i]; ys <- y[i]
  pp <- rep(NA_real_, length(i))
  for (r0 in sort(unique(meta$run[i]))) {
    tr <- which(meta$run[i] != r0); te <- which(meta$run[i] == r0)
    f <- suppressWarnings(glm(ys ~ xs, data = data.frame(ys = ys[tr], xs = xs[tr]),
                              family = binomial))
    pp[te] <- predict(f, newdata = data.frame(xs = xs[te]), type = "response")
  }
  auc_safe(ys, pp)
})
## 单基因逻辑回归的 OR（全训练集）
fit_full <- glm(y ~ x, data = data.frame(y = y, x = x), family = binomial)
or <- exp(coef(fit_full))[2]

## ------------------------------------------------ 6. 验证队列（修正矩阵；Control vs Model）
## 注意：必须读 2.机器学习\验证\ 下的修正矩阵；1.归一化去批次\验证集\ 是旧的错误版本
v <- read.csv(file.path(DATA_ROOT, "2.machine_learning", "validation", "expression_norm_log2.csv"),
              row.names = 1, check.names = FALSE)
v <- v[, colnames(v) != "", drop = FALSE]
grpv <- c(C1="Control", C2="Control", C3="Control",
          L1="Model", L2="Model", L3="Model",
          H1="Model", H2="Model", H3="Model")
selv <- colnames(v) %in% names(grpv)
xv <- as.numeric(v[g, selv]); yv <- as.integer(grpv[colnames(v)[selv]] == "Model")
auc_val <- auc_safe(yv, xv)
cmb <- combn(length(xv), sum(yv == 0))
nullv <- apply(cmb, 2, function(idx) { yy <- integer(length(xv)); yy[idx] <- 1
  auc_safe(yy, xv) })
p_val <- mean(nullv >= auc_val - 1e-12)
cat("\n验证队列（修正矩阵，Control vs Model）：AUC =", round(auc_val, 3),
    " 精确置换 p =", round(p_val, 4), "（", sum(nullv >= auc_val - 1e-12), "/", ncol(cmb), "）\n")

## ------------------------------------------------ 7. 130 组合一致性
z <- readRDS(file.path(OUT, "ml130_preds.rds"))
single <- z$vars %>% dplyr::group_by(fold, method) %>%
  dplyr::summarise(gs = paste(sort(unique(gene)), collapse = ","), n = dplyr::n(), .groups = "drop")
n_single_atf3 <- sum(single$gs == "atf3")
n_with_atf3 <- sum(grepl("(^|,)atf3(,|$)", single$gs))
n_total <- length(unique(z$preds$fold)) * length(unique(z$preds$method))

## ------------------------------------------------ 6b. 63 个基因的单基因模型（留一 Run）
uni_loo <- sapply(genes, function(gg) {
  xg <- X[, gg]; pp <- rep(NA_real_, length(y))
  for (r0 in runs) {
    tr <- which(meta$run != r0); te <- which(meta$run == r0)
    f <- suppressWarnings(glm(y ~ xg, data = data.frame(y = y[tr], xg = xg[tr]),
                              family = binomial))
    pp[te] <- predict(f, newdata = data.frame(xg = xg[te]), type = "response")
  }
  auc_safe(y, pp)
})
uni_tab <- data.frame(gene = genes, auc_leave_run = as.numeric(uni_loo)) %>%
  dplyr::arrange(dplyr::desc(auc_leave_run))
write.csv(uni_tab, file.path(TAB, "univariate_auc_leaveRun.csv"), row.names = FALSE)
cat("\n=== 63 个候选基因的单基因模型（留一 Run）前 8 ===\n"); print(head(uni_tab, 8), digits = 3, row.names = FALSE)
cat("atf3 名次:", which(uni_tab$gene == "atf3"), "/", nrow(uni_tab), "\n")

## ------------------------------------------------ 6c. 时间混杂检查
within_time <- meta %>% dplyr::mutate(v = x) %>% dplyr::group_by(time) %>%
  dplyr::group_modify(function(dd, k) data.frame(auc = auc_safe(dd$y, dd$v),
                                                 .groups = "drop"))
x_z <- ave(x, meta$time, FUN = function(z) (z - mean(z)) / sd(z))
auc_zt <- auc_safe(y, x_z)
cat("\n每个时间点内 atf3 的 AUC 中位数:", round(median(within_time$auc, na.rm = TRUE), 3),
    "; 去时间效应后总体 AUC:", round(auc_zt, 3), "\n")
write.csv(within_time, file.path(TAB, "atf3_within_time_auc.csv"), row.names = FALSE)

res <- data.frame(
  item = c("atf3 单基因 AUC（留一 Run，DeLong 95%CI）",
           "atf3 单基因 AUC（重复 5 折 × 20）",
           "两个范式各自留一 Run 的 AUC（BP / R）",
           "置换 p（留一 Run 全流程，B=1000）",
           "Brier / 校准斜率 / 校准截距",
           "Youden 阈值 / 灵敏度 / 特异度",
           "训练集 atf3 均值（Control / Treat）",
           "单基因 logistic OR（每 +1 log2）",
           "验证队列 atf3 AUC（n=12）",
           "验证队列精确置换 p",
           "130 组合中最终基因集恰为 atf3 的次数",
           "130 组合中含 atf3 的次数（分母为全部运行）",
           "63 基因单基因模型 AUC 排名（atf3）",
           "时间点内 AUC 中位数 / 去时间效应后总体 AUC"),
  value = c(sprintf("%.3f (%.3f-%.3f)", auc_loo, ci_loo[1], ci_loo[3]),
            sprintf("%.3f ± %.3f", mean(cv), sd(cv)),
            sprintf("%.3f / %.3f", per_ser[1], per_ser[2]),
            sprintf("%.4f", p_perm),
            sprintf("%.3f / %.3f / %.3f", brier, coef(cal)[2], coef(cal)[1]),
            sprintf("%.3f / %.2f / %.2f", co$threshold[1], co$sensitivity[1], co$specificity[1]),
            sprintf("%.2f / %.2f", mean(x[y == 0]), mean(x[y == 1])),
            sprintf("%.2f", or),
            sprintf("%.3f", auc_val),
            sprintf("%.4f", p_val),
            sprintf("%d / %d", n_single_atf3, n_total),
            sprintf("%d / %d (%.1f%%)", n_with_atf3, n_total, 100 * n_with_atf3 / n_total),
            sprintf("%d / %d (AUC %.3f)", which(uni_tab$gene == "atf3"), nrow(uni_tab),
                    uni_tab$auc_leave_run[uni_tab$gene == "atf3"]),
            sprintf("%.3f / %.3f", median(within_time$auc, na.rm = TRUE), auc_zt))
)
write.csv(res, file.path(TAB, "atf3_single_gene_model.csv"), row.names = FALSE)
print(res, row.names = FALSE)
saveRDS(list(p_loo = p_loo, y = y, null = null, cv = cv, auc_val = auc_val,
             p_val = p_val, p_perm = p_perm, per_ser = per_ser),
        file.path(OUT, "atf3_model.rds"))
