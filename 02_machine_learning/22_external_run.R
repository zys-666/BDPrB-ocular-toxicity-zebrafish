source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 22_external_run.R — 实际跑 GSE310406（PFNA 眼）与 GSE313277（光损伤视网膜）
##   注意：两个文件的样本名格式与 21 里的假设不同，这里按实际列名分组。
###############################################################################
suppressPackageStartupMessages({
  library(dplyr); library(pROC); library(readxl); library(AnnotationDbi); library(org.Dr.eg.db)
})
source(file.path(SCRIPT_DIR, "15_firth_fun.R"))
TAB <- file.path(WORK, "tables")
FIG <- file.path(WORK, "figures")
EXT <- file.path(DATA_ROOT, "ext")

auc_safe <- function(y, p) {
  if (length(unique(y)) < 2 || length(unique(p)) < 2) return(NA_real_)
  a <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
  max(a, 1 - a)
}
exact_perm <- function(case, ctrl) {
  v <- c(case, ctrl); n <- length(v); k <- length(ctrl)
  obs <- auc_safe(c(rep(1, length(case)), rep(0, k)), v)
  cmb <- combn(n, k)
  nul <- apply(cmb, 2, function(i) { yy <- rep(1L, n); yy[i] <- 0L; auc_safe(yy, v) })
  c(auc = obs, p = mean(nul >= obs - 1e-12), n_perm = ncol(cmb))
}

## ================================================== 1. GSE310406 (FPKM)
cat("\n########## GSE310406 (PFNA eye vs control) ##########\n")
d1 <- as.data.frame(read_excel(file.path(EXT, "GSE310406_PFNA-fpkm.xlsx"), sheet = 1))
cat("dim:", dim(d1), "\n")
cat("gene_id 前 3:", paste(head(d1$gene_id, 3), collapse = ", "), "\n")
is_symbol <- any(toupper(d1$gene_id) == "ATF3")
if (!is_symbol) {
  cat("gene_id 不是符号，尝试 Ensembl 映射…\n")
  map <- AnnotationDbi::select(org.Dr.eg.db, keys = d1$gene_id, keytype = "ENSEMBL",
                               columns = "SYMBOL")
  d1$gene_id <- map$SYMBOL[match(d1$gene_id, map$ENSEMBL)]
}
i <- which(toupper(d1$gene_id) == "ATF3")
cat("atf3 行数:", length(i), "\n")
X <- d1[i[1], -1, drop = FALSE]
v <- as.numeric(X)
names(v) <- colnames(X)
## 只保留标准命名的纯眼样本：Cntrl_Eye_DD_1_fpkm / PFNA_Eye_LD_3_fpkm 等
## （数据里混入了两个 "Cntrl_Brain_LDvsCntrl_Eye_LD_*" 的混合样本，必须排除）
keep <- grepl("^(Cntrl|PFNA)_Eye_(DD|LD)_[0-9]+_fpkm$", names(v))
cat("排除的非标准列:", paste(names(v)[!keep], collapse = ", "), "\n")
v <- v[keep]
dd <- data.frame(sample = names(v), expr = v, stringsAsFactors = FALSE)
dd$tissue <- ifelse(grepl("Eye", dd$sample), "Eye",
             ifelse(grepl("Brain", dd$sample), "Brain", NA))
dd$light  <- ifelse(grepl("_DD", dd$sample), "DD", ifelse(grepl("_LD", dd$sample), "LD", NA))
dd$group  <- ifelse(grepl("^PFNA", dd$sample), "PFNA", ifelse(grepl("^Cntrl", dd$sample), "Control", NA))
dd$expr_log2 <- log2(dd$expr + 1)
cat("\n样本表（仅眼）:\n"); print(dd %>% filter(tissue == "Eye") %>% arrange(light, group, sample), digits = 3, row.names = FALSE)

res1 <- list()
for (lt in c("DD", "LD")) {
  s <- dd %>% filter(tissue == "Eye", light == lt)
  a <- s$expr_log2[s$group == "PFNA"]; b <- s$expr_log2[s$group == "Control"]
  if (length(a) < 2 || length(b) < 2) next
  e <- exact_perm(a, b)
  f <- tryCatch(firth_report(a, b, paste0("GSE310406 eye PFNA vs control (", lt, ")")),
                error = function(x) NULL)
  res1[[lt]] <- data.frame(dataset = "GSE310406", subset = paste0("Eye_", lt),
                           n_case = length(a), n_ctrl = length(b),
                           mean_case = mean(a), mean_ctrl = mean(b),
                           log2FC = mean(a) - mean(b), fold = 2^(mean(a) - mean(b)),
                           AUC = e["auc"], p_exact = e["p"], n_perm = e["n_perm"],
                           direction = ifelse(mean(a) > mean(b), "PFNA higher", "PFNA lower"))
  res1[[paste0(lt, "_firth")]] <- f
}
tab1 <- bind_rows(res1[!grepl("_firth", names(res1))])
cat("\n=== GSE310406 结果 ===\n"); print(tab1, digits = 3, row.names = FALSE)

## ================================================== 2. GSE313277 (counts)
cat("\n########## GSE313277 (light-induced retinal degeneration) ##########\n")
d2 <- read.csv(gzfile(file.path(EXT, "GSE313277_count_matrix.csv.gz")), row.names = 1,
               check.names = FALSE)
cat("dim:", dim(d2), "  前 3 行名:", paste(head(rownames(d2), 3), collapse = ", "), "\n")
map2 <- AnnotationDbi::select(org.Dr.eg.db, keys = rownames(d2), keytype = "ENSEMBL",
                              columns = "SYMBOL")
sym <- map2$SYMBOL[match(rownames(d2), map2$ENSEMBL)]
j <- which(toupper(sym) == "ATF3")
cat("映射到 atf3 的行:", length(j), "->",
    if (length(j)) rownames(d2)[j] else NA, "\n")
if (length(j)) {
  cnt <- as.numeric(d2[j[1], ])
  names(cnt) <- colnames(d2)
  cpm <- cnt / colSums(d2) * 1e6
  dd2 <- data.frame(sample = names(cpm), cpm = as.numeric(cpm),
                    log2cpm = log2(as.numeric(cpm) + 1), stringsAsFactors = FALSE)
  dd2$group <- ifelse(grepl("LIRD", dd2$sample, ignore.case = TRUE), "LIRD",
               ifelse(grepl("Ctrl|CTRL", dd2$sample), "Control", NA))
  dd2 <- dd2[!is.na(dd2$group), ]
  print(dd2 %>% arrange(group, sample), digits = 3, row.names = FALSE)
  a <- dd2$log2cpm[dd2$group == "LIRD"]; b <- dd2$log2cpm[dd2$group == "Control"]
  cat("\nLIRD n =", length(a), " Control n =", length(b), "\n")
  e2 <- exact_perm(a, b)
  f2 <- tryCatch(firth_report(a, b, "GSE313277 LIRD vs control"), error = function(x) NULL)
  tab2 <- data.frame(dataset = "GSE313277", subset = "Retina_LIRD",
                     n_case = length(a), n_ctrl = length(b),
                     mean_case = mean(a), mean_ctrl = mean(b),
                     log2FC = mean(a) - mean(b), fold = 2^(mean(a) - mean(b)),
                     AUC = e2["auc"], p_exact = e2["p"], n_perm = e2["n_perm"],
                     direction = ifelse(mean(a) > mean(b), "LIRD higher", "LIRD lower"))
  print(tab2, digits = 3, row.names = FALSE)
  write.csv(dd2, file.path(TAB, "external_GSE313277_atf3.csv"), row.names = FALSE)
}

out <- bind_rows(tab1, if (exists("tab2")) tab2 else NULL)
write.csv(out, file.path(TAB, "external_validation_stats.csv"), row.names = FALSE)
write.csv(dd, file.path(TAB, "external_GSE310406_atf3.csv"), row.names = FALSE)
cat("\n===== 汇总表 =====\n"); print(out, digits = 3, row.names = FALSE)

## ================================================== 3. 图
suppressPackageStartupMessages(library(ggplot2))
p1 <- dd %>% filter(tissue == "Eye") %>%
  ggplot(aes(group, expr_log2, fill = group)) +
  geom_boxplot(width = .55, alpha = .45, outlier.shape = NA) +
  geom_jitter(width = .12, size = 1.8, alpha = .85) +
  facet_wrap(~light) +
  scale_fill_manual(values = c(Control = "#2166AC", PFNA = "#B2182B"), guide = "none") +
  labs(x = NULL, y = "atf3 log2(FPKM+1)", title = "GSE310406 — PFNA vs control (eye)",
       subtitle = "DD / LD = constant dark vs 14L:10D") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "fig22_GSE310406_atf3.pdf"), p1, width = 7, height = 4.2)
ggsave(file.path(FIG, "fig22_GSE310406_atf3.png"), p1, width = 7, height = 4.2, dpi = 300)

if (exists("dd2")) {
  p2 <- ggplot(dd2, aes(group, log2cpm, fill = group)) +
    geom_boxplot(width = .45, alpha = .45, outlier.shape = NA) +
    geom_jitter(width = .12, size = 1.8, alpha = .85) +
    scale_fill_manual(values = c(Control = "#2166AC", LIRD = "#B2182B"), guide = "none") +
    labs(x = NULL, y = "atf3 log2(CPM+1)", title = "GSE313277 — light-induced retinal degeneration") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG, "fig22_GSE313277_atf3.pdf"), p2, width = 5.4, height = 4.2)
  ggsave(file.path(FIG, "fig22_GSE313277_atf3.png"), p2, width = 5.4, height = 4.2, dpi = 300)
}
cat("\n完成\n")
