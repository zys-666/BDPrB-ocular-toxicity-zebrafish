source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 40_build_107.R — 按"blue 模块 ∩ PPI = 107"重建候选基因面板，并做单基因筛查
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(readxl) })
ROOT <- DATA_ROOT
OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")

## ---- 107 基因 ----
blue <- read.csv(file.path(ROOT, "2.\u673a\u5668\u5b66\u4e60/\u8bad\u7ec3/blue.csv"),
                 check.names = FALSE)
ppi <- as.data.frame(read_excel(file.path(ROOT, "2.\u673a\u5668\u5b66\u4e60/\u8bad\u7ec3/HC\u548cLC\u7684ppi\u7684\u4ea4\u96c6\u57fa\u56e0.xls")))
bg <- unique(as.character(blue[[1]])); pg <- unique(as.character(ppi[[1]]))
cand <- sort(unique(intersect(bg, pg)))
cat("blue 模块:", length(bg), " PPI:", length(pg), " 交集:", length(cand), "\n")

## ---- 表达矩阵 ----
norm <- read.delim(file.path(ROOT, "2.\u673a\u5668\u5b66\u4e60/\u8bad\u7ec3/\u8f93\u5165/normalize.txt"),
                   row.names = 1, check.names = FALSE)
info <- read.csv(file.path(ROOT, "1.\u5f52\u4e00\u5316\u53bb\u6279\u6b21/\u8bad\u7ec3\u96c6/sample_info.csv"),
                 stringsAsFactors = FALSE)
present <- cand[cand %in% rownames(norm)]
cat("107 中在表达矩阵里的:", length(present), "；缺失:", paste(setdiff(cand, present), collapse = ", "), "\n")
cat("atf3 是否在面板中:", "atf3" %in% present, "\n")

mi <- match(sub("_.*$", "", colnames(norm)), info$GSM_ID)
meta <- data.frame(col = colnames(norm), series = ifelse(grepl("^BP_", sub("^GSM[0-9]+_", "",
                  info$FileName[mi])), "BP", "R"),
                  run = paste0("Run", sub("^R(un)?", "", info$Replicate[mi])),
                  time = info$Time[mi], group = info$Group[mi], stringsAsFactors = FALSE)
meta$y <- ifelse(meta$group == "Control", 0L, ifelse(meta$group == "Treatment", 1L, NA_integer_))
keep <- !is.na(meta$y)
X <- t(as.matrix(norm[present, keep])); storage.mode(X) <- "double"
meta <- meta[keep, ]
cat("样本:", nrow(X), " 特征:", ncol(X), " 类别:", paste(names(table(meta$y)), table(meta$y), collapse = "/"), "\n")

saveRDS(list(X = X, y = meta$y, meta = meta, genes = present), file.path(OUT, "prep107.rds"))
writeLines(present, file.path(TAB, "candidates_107.txt"))

## ---- 单基因筛查（留一 Run）----
auc_safe <- function(y, p) {
  if (length(unique(y)) < 2 || length(unique(p)) < 2) return(NA_real_)
  a <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))); max(a, 1 - a)
}
runs <- sort(unique(meta$run))
uni <- data.frame(gene = present, auc = NA_real_)
for (g in present) {
  pp <- rep(NA_real_, length(meta$y))
  xg <- X[, g]
  for (r0 in runs) {
    tr <- which(meta$run != r0); te <- which(meta$run == r0)
    f <- suppressWarnings(glm(y ~ xg, data = data.frame(y = meta$y[tr], xg = xg[tr]), family = binomial))
    pp[te] <- predict(f, newdata = data.frame(xg = xg[te]), type = "response")
  }
  uni$auc[uni$gene == g] <- auc_safe(meta$y, pp)
}
uni <- uni[order(-uni$auc), ]
rownames(uni) <- NULL
write.csv(uni, file.path(TAB, "univariate_auc_107.csv"), row.names = FALSE)
cat("\n=== 107 基因单基因 AUC 前 12 ===\n"); print(head(uni, 12), digits = 3, row.names = FALSE)
cat("atf3 名次:", which(uni$gene == "atf3"), "/", nrow(uni),
    "  AUC =", round(uni$auc[uni$gene == "atf3"], 3), "\n")
