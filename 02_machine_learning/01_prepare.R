source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## 01_prepare.R — 修正后的数据准备
##   修正点：
##   1) 6 个 0 hr 基线样本被原脚本标为 Treat  → 本脚本剔除（另存一份含基线的敏感性数据）
##   2) BP_RunX 与 RX 是同一 GEO 系列里的两套命名体系（两个消融范式）
##      → 记录 Series / Run / Time，供分组交叉验证使用
##   3) 特征面板固定为已发表的 63 基因（不做任何再次筛选），保证与原文可比
## 输出：out/prep.rds
###############################################################################
suppressPackageStartupMessages({ library(dplyr) })

ROOT <- DATA_ROOT
OUT <- file.path(WORK, "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

norm <- read.delim(file.path(ROOT, "2.机器学习/训练/输入/normalize.txt"),
                   row.names = 1, check.names = FALSE)
info <- read.csv(file.path(ROOT, "1.归一化去批次/训练集/sample_info.csv"), stringsAsFactors = FALSE)

ig <- read.delim(file.path(ROOT, "2.机器学习/训练/输入/interGenes.txt"), header = FALSE)[[1]]
ig <- unique(ig)
present <- ig[ig %in% rownames(norm)]

## ---- 样本元信息 ----
mi <- match(sub("_.*$", "", colnames(norm)), info$GSM_ID)
meta <- data.frame(
  col    = colnames(norm),
  GSM    = sub("_.*$", "", colnames(norm)),
  series = ifelse(grepl("^BP_", sub("^GSM[0-9]+_", "", info$FileName[mi])), "BP", "R"),
  run    = paste0("Run", sub("^R(un)?", "", info$Replicate[mi])),
  time   = info$Time[mi],
  trt    = info$Treatment[mi],
  group  = info$Group[mi],
  stringsAsFactors = FALSE
)
meta$y <- ifelse(meta$group == "Control", 0L, ifelse(meta$group == "Treatment", 1L, NA_integer_))
meta$stratum <- paste(meta$series, meta$run, sep = "_")   # 6 个交叉验证块
meta$cond    <- paste(meta$time, meta$trt, meta$series, sep = "|")

cat("总样本:", nrow(meta), " 基线:", sum(is.na(meta$y)), "\n")
cat("特征面板: interGenes", length(ig), "，在表达矩阵中", length(present), "，缺失:",
    paste(setdiff(ig, present), collapse = ","), "\n")

X <- t(as.matrix(norm[present, , drop = FALSE]))
storage.mode(X) <- "double"
X <- X[, present, drop = FALSE]
cat("X dim:", dim(X), "\n")

keep <- !is.na(meta$y)
dat <- list(
  X = X[keep, , drop = FALSE],
  y = meta$y[keep],
  meta = meta[keep, ],
  genes = present,
  n_all = nrow(meta)
)

## 类别/块平衡度
cat("\n[类别平衡] Control =", sum(dat$y == 0), " Treat =", sum(dat$y == 1), "\n")
cat("\n[交叉验证块 x 类别]\n"); print(table(dat$meta$stratum, dat$y))
cat("\n[时间点 x 类别]\n"); print(table(dat$meta$time, dat$y))
cat("\n[系列 x 类别]\n"); print(table(dat$meta$series, dat$y))
cat("\n[时间点 x 系列]\n"); print(table(dat$meta$time, dat$meta$series))

## 零方差检查
v <- apply(dat$X, 2, var)
cat("\n零方差特征:", sum(v == 0), " 最小方差:", signif(min(v), 3), "\n")

saveRDS(dat, file.path(OUT, "prep.rds"))
cat("\n已保存:", file.path(OUT, "prep.rds"), "\n")
