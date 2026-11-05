source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 23_GSE284728.R — 视神经切断（ONT）损伤模型：atf3 的独立外部验证
##   I_* = injured（1/5/7 dpi），U_* = uninjured
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(ggplot2) })
source(file.path(SCRIPT_DIR, "15_firth_fun.R"))
TAB <- file.path(WORK, "tables")
FIG <- file.path(WORK, "figures")
EXT <- file.path(DATA_ROOT, "ext")

auc_safe <- function(y, p) {
  if (length(unique(y)) < 2 || length(unique(p)) < 2) return(NA_real_)
  a <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))); max(a, 1 - a)
}
exact_perm <- function(case, ctrl, max_comb = 2e5) {
  v <- c(case, ctrl); n <- length(v); k <- length(ctrl)
  obs <- auc_safe(c(rep(1, length(case)), rep(0, k)), v)
  tot <- choose(n, k)
  if (tot <= max_comb) {
    cmb <- combn(n, k)
    nul <- apply(cmb, 2, function(i) { yy <- rep(1L, n); yy[i] <- 0L; auc_safe(yy, v) })
    p <- mean(nul >= obs - 1e-12); np <- ncol(cmb)
  } else {
    set.seed(1); B <- 20000
    nul <- replicate(B, { yy <- sample(c(rep(1L, length(case)), rep(0L, k)))
                          auc_safe(yy, v) })
    p <- (1 + sum(nul >= obs - 1e-12)) / (B + 1); np <- B
  }
  c(auc = obs, p = p, n_perm = np)
}

d <- read.delim(gzfile(file.path(EXT, "GSE284728_BulkRNAseq_raw_counts.tsv.gz")),
                check.names = FALSE)
cat("dim:", dim(d), "\n")
cnt_cols <- grep("^[IU]_[0-9]|^[IU]_[0-9]+DPI", colnames(d), value = TRUE, ignore.case = TRUE)
cnt_cols <- grep("^[IU]_", colnames(d), value = TRUE)
cnt_cols <- setdiff(cnt_cols, grep("Mean|minus", cnt_cols, value = TRUE))
cat("计数列数:", length(cnt_cols), "\n")
cat(paste(cnt_cols, collapse = " | "), "\n")

i <- which(toupper(d$Gene_name) == "ATF3")
cat("atf3 行:", length(i), "\n")
stopifnot(length(i) == 1)
cnt <- as.numeric(d[i, cnt_cols])
names(cnt) <- cnt_cols
tot <- colSums(as.matrix(d[, cnt_cols]))
cpm <- cnt / tot * 1e6
dd <- data.frame(sample = names(cpm), cpm = as.numeric(cpm),
                 log2cpm = log2(as.numeric(cpm) + 1), stringsAsFactors = FALSE)
dd$group <- ifelse(grepl("^I_", dd$sample), "Injured", "Uninjured")
dd$tp <- sub("^[IU]_([0-9]+)DP[iI].*$", "\\1", dd$sample)
dd$tp <- paste0(dd$tp, " dpi")
print(dd %>% arrange(tp, group, sample), digits = 3, row.names = FALSE)

res <- list()
a <- dd$log2cpm[dd$group == "Injured"]; b <- dd$log2cpm[dd$group == "Uninjured"]
e <- exact_perm(a, b)
f <- tryCatch(firth_report(a, b, "GSE284728 ONT injured vs uninjured"), error = function(x) NULL)
res[["all"]] <- data.frame(dataset = "GSE284728", subset = "all time points",
                           n_case = length(a), n_ctrl = length(b),
                           mean_case = mean(a), mean_ctrl = mean(b),
                           log2FC = mean(a) - mean(b), fold = 2^(mean(a) - mean(b)),
                           AUC = e["auc"], p_exact = e["p"], n_perm = e["n_perm"],
                           direction = ifelse(mean(a) > mean(b), "injury higher", "injury lower"))
for (t0 in sort(unique(dd$tp))) {
  s <- dd[dd$tp == t0, ]
  aa <- s$log2cpm[s$group == "Injured"]; bb <- s$log2cpm[s$group == "Uninjured"]
  if (length(aa) < 2 || length(bb) < 2) next
  ee <- exact_perm(aa, bb)
  res[[t0]] <- data.frame(dataset = "GSE284728", subset = paste0("ONT_", t0),
                          n_case = length(aa), n_ctrl = length(bb),
                          mean_case = mean(aa), mean_ctrl = mean(bb),
                          log2FC = mean(aa) - mean(bb), fold = 2^(mean(aa) - mean(bb)),
                          AUC = ee["auc"], p_exact = ee["p"], n_perm = ee["n_perm"],
                          direction = ifelse(mean(aa) > mean(bb), "injury higher", "injury lower"))
}
tab <- bind_rows(res)
cat("\n=== 结果 ===\n"); print(tab, digits = 3, row.names = FALSE)
write.csv(tab, file.path(TAB, "external_GSE284728_stats.csv"), row.names = FALSE)
write.csv(dd, file.path(TAB, "external_GSE284728_atf3.csv"), row.names = FALSE)
if (!is.null(f)) { cat("\n=== Firth（9 vs 9）===\n"); print(f, digits = 3, row.names = FALSE) }

p <- ggplot(dd, aes(group, log2cpm, fill = group)) +
  geom_boxplot(width = .5, alpha = .45, outlier.shape = NA) +
  geom_jitter(width = .12, size = 1.8, alpha = .85) +
  facet_wrap(~tp, nrow = 1) +
  scale_fill_manual(values = c(Uninjured = "#2166AC", Injured = "#B2182B"), guide = "none") +
  labs(x = NULL, y = "atf3 log2(CPM+1)",
       title = "GSE284728 — optic nerve transection (RGCs)",
       subtitle = "I = injured, U = uninjured; 1 / 5 / 7 days post injury") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "fig23_GSE284728_atf3.pdf"), p, width = 9, height = 4.4)
ggsave(file.path(FIG, "fig23_GSE284728_atf3.png"), p, width = 9, height = 4.4, dpi = 300)
cat("\n图: figures/fig23_GSE284728_atf3.png\n")
