source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 19_validation_final.R — 用修正后的验证矩阵 + 正确分组（C=对照；L/H=造模；PP=治疗）
##   主对比：Control vs Model(L+H)
##   次要：Control vs H、Control vs L、H vs PP（治疗是否把 atf3 拉回）、Control vs PP
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
source(file.path(SCRIPT_DIR, "15_firth_fun.R"))

OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")
FIG <- file.path(WORK, "figures")

f <- file.path(DATA_ROOT, "2.machine_learning", "validation", "expression_norm_log2.csv")
v <- read.csv(f, row.names = 1, check.names = FALSE)
v <- v[, colnames(v) != "", drop = FALSE]
stopifnot("atf3" %in% rownames(v))
x <- as.numeric(v["atf3", ]); names(x) <- colnames(v)

## 正确分组
grp <- c(C1="Control", C2="Control", C3="Control",
         L1="Low-dose model", L2="Low-dose model", L3="Low-dose model",
         H1="High-dose model", H2="High-dose model", H3="High-dose model",
         PP1="Treatment", PP2="Treatment", PP3="Treatment")
df <- data.frame(id = names(x), group = grp[names(x)], atf3 = as.numeric(x),
                 stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("Control", "Low-dose model",
                                        "High-dose model", "Treatment"))
cat("\n=== 修正后 atf3（log2 FPKM）===\n")
print(df, row.names = FALSE, digits = 3)

sm <- df %>% group_by(group) %>%
  summarise(n = n(), mean = mean(atf3), sd = sd(atf3), min = min(atf3), max = max(atf3),
            .groups = "drop")
cat("\n=== 分组均值 ===\n"); print(sm, digits = 3, row.names = FALSE)
write.csv(sm, file.path(TAB, "validation_final_group_means.csv"), row.names = FALSE)

auc1 <- function(a, b, label) {
  y <- c(rep(1, length(a)), rep(0, length(b)))
  data.frame(contrast = label, n_case = length(a), n_ctrl = length(b),
             AUC = as.numeric(pROC::auc(pROC::roc(y, c(a, b), quiet = TRUE, direction = "<"))),
             stringsAsFactors = FALSE)
}
g <- function(k) df$atf3[df$group == k]
res <- bind_rows(
  auc1(c(g("Low-dose model"), g("High-dose model")), g("Control"), "Control vs Model (L+H)"),
  auc1(g("High-dose model"), g("Control"), "Control vs High-dose model"),
  auc1(g("Low-dose model"), g("Control"), "Control vs Low-dose model"),
  auc1(g("Treatment"), g("Control"), "Control vs Treatment (PP)"),
  auc1(g("High-dose model"), g("Treatment"), "Treatment vs High-dose model")
)
res$方向 <- c("model 更高", "high 更高", "low 更高", "PP 更高", "PP 更低（治疗回落）")
cat("\n=== 各对比 AUC ===\n"); print(res, digits = 3, row.names = FALSE)
write.csv(res, file.path(TAB, "validation_final_auc.csv"), row.names = FALSE)

## 精确置换检验（Control 3 vs Model 6 -> C(9,3)=84 种分配）
a <- c(g("Low-dose model"), g("High-dose model")); b <- g("Control")
obs <- as.numeric(pROC::auc(pROC::roc(c(rep(1, length(a)), rep(0, length(b))),
                                     c(a, b), quiet = TRUE, direction = "<")))
cmb <- combn(length(c(a, b)), length(b))
null <- apply(cmb, 2, function(i) {
  yy <- integer(length(c(a, b))); yy[i] <- 1
  as.numeric(pROC::auc(pROC::roc(yy, c(a, b), quiet = TRUE, direction = "<")))
})
p_exact <- mean(null >= obs - 1e-12)
cat(sprintf("\nControl vs Model: AUC = %.3f，精确置换 p = %d/%d = %.4f\n",
            obs, sum(null >= obs - 1e-12), ncol(cmb), p_exact))

## 剂量趋势（C < L < H）：精确置换 —— 三组共 9 个样本，C(9,3,3,3)=1680 种分配
dose <- rep(c(1, 2, 3), each = 3)
vals <- c(g("Control"), g("Low-dose model"), g("High-dose model"))
ord_obs <- all(max(vals[dose == 1]) < min(vals[dose == 2])) &&
           all(max(vals[dose == 2]) < min(vals[dose == 3]))
## 计算所有 1680 种分组方式中"完全单调"的次数
cnt <- 0
idx_all <- 1:9
for (i in combn(idx_all, 3, simplify = FALSE)) {
  rest <- setdiff(idx_all, i)
  for (j in combn(rest, 3, simplify = FALSE)) {
    k <- setdiff(rest, j)
    if (max(vals[i]) < min(vals[j]) && max(vals[j]) < min(vals[k])) cnt <- cnt + 1
  }
}
cat(sprintf("\n剂量趋势（C<L<H）精确置换 p = %d/1680 = %.5f\n", cnt, cnt / 1680))

## Firth
f1 <- firth_report(a, b, "Control vs Model (L+H)")
f2 <- firth_report(g("High-dose model"), g("Treatment"), "High-dose vs Treatment")
cat("\n=== Firth ===\n"); print(bind_rows(f1, f2), digits = 3, row.names = FALSE)
write.csv(bind_rows(f1, f2), file.path(TAB, "validation_final_firth.csv"), row.names = FALSE)

## 与发现集方向
d <- readRDS(file.path(OUT, "prep.rds"))
tm <- tapply(d$X[, "atf3"], d$y, mean)
cat(sprintf("\n发现集 GSE234646: %.2f -> %.2f (%.2f 倍, 上升)\n", tm[1], tm[2], 2^(tm[2] - tm[1])))
cat(sprintf("验证集 Control vs Model: %.2f -> %.2f (%.2f 倍, %s)\n",
            mean(b), mean(a), 2^(mean(a) - mean(b)),
            ifelse(mean(a) > mean(b), "上升 → 与发现集同向", "下降 → 方向相反")))

## 图
library(ggplot2)
p <- ggplot(df, aes(group, atf3, fill = group)) +
  geom_boxplot(width = .5, alpha = .45, outlier.shape = NA) +
  geom_jitter(width = .1, size = 2, alpha = .85) +
  scale_fill_manual(values = c("#2166AC", "#F4A582", "#B2182B", "#1B7837"), guide = "none") +
  labs(x = NULL, y = "atf3 (log2 FPKM)",
       title = "atf3 in the BDPrB cohort (corrected matrix)",
       subtitle = "C = control; L / H = BDPrB model (low / high dose); PP = high dose + Lb-0401") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "fig19_validation_atf3_groups.pdf"), p, width = 6.6, height = 4.6)
ggsave(file.path(FIG, "fig19_validation_atf3_groups.png"), p, width = 6.6, height = 4.6, dpi = 300)
cat("\n完成\n")
