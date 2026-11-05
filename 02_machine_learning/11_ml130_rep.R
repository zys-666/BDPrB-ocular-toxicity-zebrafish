source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 11_ml130_rep.R — 130 组合 × 重复分层 5 折（带选择步骤缓存）
##   与 10_ml130.R 相同的规则：特征选择在折内、允许单基因、三阶段全执行。
##   缓存：同一折内相同的"阶段 + 输入列集合"只计算一次（原模板对每个组合重复计算，
##   是 130 组合跑一小时的主因）。
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
source(file.path(SCRIPT_DIR, "mime_functions.R"))

OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")
REF <- file.path(DATA_ROOT, "130机器学习组合", "refer.txt")
classVar <- "Type"; MAXVAR <- 20; MINVAR <- 1
K <- as.integer(Sys.getenv("ML130_K", "5")); NREP <- as.integer(Sys.getenv("ML130_REP", "1"))

d <- readRDS(file.path(OUT, "prep.rds"))
X <- d$X; y <- d$y
methods <- read.delim(REF, header = TRUE, stringsAsFactors = FALSE)$Model
methods <- trimws(gsub("\\s+", " ", methods))
methods <- gsub(" ?\\+ ?", "+", methods)
methods <- gsub("\\s*\\[\\s*", "[", methods); methods <- gsub("\\s*\\]\\s*", "]", methods)
methods <- unique(methods[nzchar(methods)])
cat("组合数:", length(methods), "\n")

## ---- 带缓存的选择 ----
cache <- new.env(parent = emptyenv())
get_sel <- function(fid, stage, Xs, lab_tr, tag) {
  key <- paste(fid, stage, tag, sep = "||")
  if (!is.null(cache[[key]])) return(cache[[key]])
  v <- tryCatch(suppressWarnings(RunML(stage, Xs, lab_tr, mode = "Variable", classVar = classVar)),
                error = function(e) character(0))
  v <- intersect(v, colnames(Xs))
  if (length(v) > MAXVAR) v <- v[seq_len(MAXVAR)]
  cache[[key]] <- v
  v
}

run_combo <- function(fid, stages, Xtr, ytr, Xte) {
  lab_tr <- data.frame(Type = ytr)
  vars <- colnames(Xtr); tag <- "root"
  for (i in seq_along(stages)) {
    st <- stages[i]; Xs <- Xtr[, vars, drop = FALSE]
    if (i < length(stages)) {
      v <- get_sel(fid, st, Xs, lab_tr, tag)
      if (length(v) < MINVAR) return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      vars <- v; tag <- paste0(tag, ">", st, ":", length(v))
    } else {
      if (grepl("^plsRglm", st) && length(vars) < 3)
        return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      fit <- tryCatch(suppressWarnings(RunML(st, Xs, lab_tr, mode = "Model", classVar = classVar)),
                      error = function(e) NULL)
      if (is.null(fit)) return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      p <- tryCatch(CalPredictScore(fit, Xte[, vars, drop = FALSE]),
                    error = function(e) rep(NA_real_, nrow(Xte)))
      if (length(p) != nrow(Xte)) return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      return(list(p = as.numeric(p), vars = vars))
    }
  }
  list(p = rep(NA_real_, nrow(Xte)), vars = character(0))
}

set.seed(99)
folds <- list()
for (r in seq_len(NREP)) {
  f <- integer(length(y))
  for (v in c(0, 1)) { i <- which(y == v); f[i] <- sample(rep(seq_len(K), length.out = length(i))) }
  for (k in seq_len(K)) folds[[paste0("rep", r, "_f", k)]] <- list(train = which(f != k), test = which(f == k))
}

pred_rec <- list(); var_rec <- list(); t0 <- Sys.time()
for (fn in names(folds)) {
  tr <- folds[[fn]]$train; te <- folds[[fn]]$test
  mu <- colMeans(X[tr, , drop = FALSE]); sdv <- apply(X[tr, , drop = FALSE], 2, sd); sdv[sdv == 0] <- 1
  Xtr <- scale(X[tr, , drop = FALSE], center = mu, scale = sdv)
  Xte <- scale(X[te, , drop = FALSE], center = mu, scale = sdv)
  ytr <- y[tr]
  for (m in methods) {
    out <- run_combo(fn, strsplit(m, "+", fixed = TRUE)[[1]], Xtr, ytr, Xte)
    pred_rec[[length(pred_rec) + 1]] <- data.frame(fold = fn, method = m, idx = te,
                                                   y = y[te], p = out$p, stringsAsFactors = FALSE)
    if (length(out$vars)) var_rec[[length(var_rec) + 1]] <- data.frame(
      fold = fn, method = m, gene = out$vars, stringsAsFactors = FALSE)
  }
  cat(sprintf("%-8s done, %.1f min\n", fn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

preds <- dplyr::bind_rows(pred_rec); vars <- dplyr::bind_rows(var_rec)
saveRDS(list(preds = preds, vars = vars, methods = methods, folds = folds),
        file.path(OUT, "ml130_preds_rep.rds"))

auc_safe <- function(y, p) {
  ok <- is.finite(p)
  if (sum(ok) < 10 || length(unique(y[ok])) < 2 || length(unique(p[ok])) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y[ok], p[ok], quiet = TRUE, direction = "<")))
}
perf <- preds %>% dplyr::group_by(method) %>% dplyr::group_modify(function(dd, k) {
  fa <- dd %>% dplyr::group_by(fold) %>% dplyr::summarise(a = auc_safe(y, p), .groups = "drop")
  data.frame(auc_mean = mean(fa$a, na.rm = TRUE), auc_sd = sd(fa$a, na.rm = TRUE),
             auc_pooled = auc_safe(dd$y, dd$p), n_eval = sum(is.finite(dd$p)),
             n_fold_ok = sum(is.finite(fa$a)))
}) %>% dplyr::arrange(dplyr::desc(auc_mean))
write.csv(perf, file.path(TAB, "ml130rep_combo_performance.csv"), row.names = FALSE)

gfreq <- vars %>% dplyr::count(gene, name = "n") %>%
  dplyr::mutate(freq_overall = n / (length(methods) * length(folds))) %>%
  dplyr::arrange(dplyr::desc(n))
write.csv(gfreq, file.path(TAB, "ml130rep_gene_frequency.csv"), row.names = FALSE)

cat("\n=== 重复 5 折：组合 AUC 前 12 ===\n"); print(head(perf, 12), digits = 3, row.names = FALSE)
cat("\n=== 基因入选频率 前 12 ===\n"); print(head(gfreq, 12), digits = 3, row.names = FALSE)
cat("\natf3 名次:", which(gfreq$gene == "atf3"), " 入选", gfreq$n[gfreq$gene == "atf3"],
    "/", length(methods) * length(folds), "\n")
cat("总耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
