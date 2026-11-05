source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 10_ml130.R — 130 种机器学习组合（Mime 框架）在修正后的数据上重跑
##
## 与原模板的关键差别（都是为了消除泄漏 / 让单基因模型能留下来）：
##   1) 特征选择放在每一折内部（原模板在全体训练集上一次做完，再在外部测试集评估）
##   2) 三阶段组合按顺序执行全部分阶段，最后一阶段才是学习器
##      （原模板对 3 段组合只用了第 1 段选基因 + 第 2 段当模型，最后的学习器被丢掉）
##   3) min.selected.var = 0 → 允许"只留下 1 个基因"的模型（原模板 length<=1 直接 skip，
##      这正是 atf3 单基因模型永远出不来的原因）
##   4) 标准化只用训练折的均值/方差，测试折不参与
##   5) 交叉验证按 Run 分块（留一 Run），不用随机划分
##
## 输出：out/ml130_*.rds、tables/ml130_*.csv
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
source(file.path(SCRIPT_DIR, "mime_functions.R"))

OUT <- file.path(WORK, "out")
TAB <- file.path(WORK, "tables")
REF <- file.path(DATA_ROOT, "130机器学习组合", "refer.txt")

classVar   <- "Type"
MAXVAR     <- 20      # 与原模板一致
MINVAR     <- 1       # 允许单基因模型（原模板为 1 且用 <= 跳过）
NFOLD_TEST <- as.integer(Sys.getenv("ML130_TESTFOLD", "0"))

d <- readRDS(file.path(OUT, "prep.rds"))
X <- d$X; y <- d$y; meta <- d$meta; genes <- d$genes

## ---- 方法列表 ----
methods <- read.delim(REF, header = TRUE, stringsAsFactors = FALSE)$Model
methods <- trimws(gsub("\\s+", " ", methods))
methods <- gsub(" ?\\+ ?", "+", methods)          # 统一成 A+B
methods <- gsub("\\s*\\[\\s*", "[", methods)      # Stepglm [both] -> Stepglm[both]
methods <- gsub("\\s*\\]\\s*", "]", methods)
methods <- unique(methods[nzchar(methods)])
cat("组合数:", length(methods), "\n")

split_stages <- function(m) strsplit(m, "\\+", fixed = FALSE)[[1]]

## ---- 折 ----
runs <- sort(unique(meta$run))
folds <- lapply(runs, function(r0) list(train = which(meta$run != r0), test = which(meta$run == r0)))
names(folds) <- runs
if (NFOLD_TEST > 0) folds <- folds[seq_len(min(NFOLD_TEST, length(folds)))]

## ---- 单组合执行 ----
run_combo <- function(stages, Xtr, ytr, Xte) {
  vars <- colnames(Xtr)
  lab_tr <- data.frame(Type = ytr)
  for (i in seq_along(stages)) {
    st <- stages[i]
    Xs <- Xtr[, vars, drop = FALSE]
    if (i < length(stages)) {
      v <- tryCatch(suppressWarnings(
             RunML(st, Xs, lab_tr, mode = "Variable", classVar = classVar)),
             error = function(e) character(0))
      v <- intersect(v, vars)
      if (length(v) < MINVAR) return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      if (length(v) > MAXVAR) v <- v[seq_len(MAXVAR)]
      vars <- v
    } else {
      ## plsRglm 在只剩 1 个变量时会挂起（cv.plsRglm 死循环），按原框架语义跳过该组合
      if (grepl("^plsRglm", st) && length(vars) < 3)
        return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      fit <- tryCatch(RunML(st, Xs, lab_tr, mode = "Model", classVar = classVar),
                      error = function(e) NULL)
      if (is.null(fit)) return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      p <- tryCatch(CalPredictScore(fit, Xte[, vars, drop = FALSE]),
                    error = function(e) rep(NA_real_, nrow(Xte)))
      p <- as.numeric(p)
      if (length(p) != nrow(Xte)) return(list(p = rep(NA_real_, nrow(Xte)), vars = character(0)))
      return(list(p = p, vars = vars))
    }
  }
  list(p = rep(NA_real_, nrow(Xte)), vars = character(0))
}

## ---- 交叉验证主循环 ----
pred_rec <- list(); var_rec <- list()
t0 <- Sys.time()
for (fn in names(folds)) {
  tr <- folds[[fn]]$train; te <- folds[[fn]]$test
  ## 只用训练折标准化
  mu <- colMeans(X[tr, , drop = FALSE]); sdv <- apply(X[tr, , drop = FALSE], 2, sd)
  sdv[sdv == 0] <- 1
  Xtr <- scale(X[tr, , drop = FALSE], center = mu, scale = sdv)
  Xte <- scale(X[te, , drop = FALSE], center = mu, scale = sdv)
  ytr <- y[tr]
  t1 <- Sys.time()
  for (m in methods) {
    st <- split_stages(m)
    out <- run_combo(st, Xtr, ytr, Xte)
    pred_rec[[length(pred_rec) + 1]] <- data.frame(
      fold = fn, method = m, idx = te, y = y[te], p = out$p, stringsAsFactors = FALSE)
    if (length(out$vars)) {
      var_rec[[length(var_rec) + 1]] <- data.frame(
        fold = fn, method = m, gene = out$vars, stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("fold %-3s (%d 训练 / %d 测试) 完成，用时 %.1f s，累计 %.1f min\n",
              fn, length(tr), length(te),
              as.numeric(difftime(Sys.time(), t1, units = "secs")),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

preds <- bind_rows(pred_rec); vars <- bind_rows(var_rec)
saveRDS(list(preds = preds, vars = vars, methods = methods), file.path(OUT, "ml130_preds.rds"))

auc_safe <- function(y, p) {
  ok <- is.finite(p)
  if (sum(ok) < 10 || length(unique(y[ok])) < 2 || length(unique(p[ok])) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y[ok], p[ok], quiet = TRUE, direction = "<")))
}

## ---- 每个组合的性能 ----
perf <- preds %>% dplyr::group_by(method) %>% dplyr::group_modify(function(dd, k) {
  foldau <- dd %>% dplyr::group_by(fold) %>%
    dplyr::summarise(a = auc_safe(y, p), n_ok = sum(is.finite(p)), .groups = "drop")
  data.frame(auc_mean = mean(foldau$a, na.rm = TRUE),
             auc_sd = sd(foldau$a, na.rm = TRUE),
             auc_pooled = auc_safe(dd$y, dd$p),
             n_eval = sum(is.finite(dd$p)),
             n_genes_mean = mean(table(vars$method)[k$method], na.rm = TRUE))
}) %>% dplyr::arrange(dplyr::desc(auc_mean))
write.csv(perf, file.path(TAB, "ml130_combo_performance.csv"), row.names = FALSE)

## ---- 基因入选频率 ----
gfreq <- vars %>% dplyr::count(gene, name = "n") %>%
  dplyr::mutate(freq_method = n / length(methods),
                freq_all = n / (length(methods) * length(folds))) %>%
  dplyr::arrange(dplyr::desc(n))
write.csv(gfreq, file.path(TAB, "ml130_gene_frequency.csv"), row.names = FALSE)

cat("\n=== 每个组合的测试集 AUC（前 20）===\n")
print(head(perf[, c("method","auc_mean","auc_sd","auc_pooled","n_eval")], 20), digits = 3, row.names = FALSE)
cat("\n=== 基因入选次数（前 20，共", length(folds), "折 ×", length(methods), "组合）===\n")
print(head(gfreq, 20), digits = 3, row.names = FALSE)
cat("\natf3 名次:", which(gfreq$gene == "atf3"), "/", nrow(gfreq), "\n")
cat("\n总耗时:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
