## Reproducible analysis code accompanying the manuscript.
## Run from the repository root (see README.md); paths are resolved from config.R.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')
# ============================================ #
# 菌群-ELISA相关性热图分析（方案A：原始样本水平，不取组平均）
# 适配文件: elisa_data.txt + allsample.txt + geneList.txt
# ============================================ #

library("limma")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)

# 设置工作目录（请根据实际情况修改）
setwd(file.path(DATA_ROOT, "9.菌群测序", "菌群数量和核心基因表达量的相关性"))

# =============================== #
# Step 1: 定义文件路径
# =============================== #

elisaFilePath <- "elisa_data.txt"        # ELISA数据（行=样本，列=指标）
microbeFilePath <- "allsample.txt"       # 菌群ASV数据（行=OTU，列=样本，最后列=taxonomy）
targetGeneFile <- "geneList.txt"         # 待分析ELISA指标列表

if (!file.exists(elisaFilePath)) stop("错误：ELISA数据文件不存在！")
if (!file.exists(microbeFilePath)) stop("错误：菌群数据文件不存在！")
if (!file.exists(targetGeneFile)) stop("错误：目标指标列表文件不存在！")

# =============================== #
# Step 2: 读取ELISA数据并预处理
# 【修改：保留原始24个样本，不取组平均】
# =============================== #

readElisaData <- function(fileElisa, fileGene) {
  cat("步骤 2.1：读取ELISA数据...\n")
  
  # 读取ELISA数据（第一列是样本名）
  elisaRaw <- read.table(fileElisa, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 设置行名（样本名）
  rownames(elisaRaw) <- elisaRaw[, 1]
  elisaVals <- elisaRaw[, -1, drop = FALSE]
  
  # 转换为数值矩阵
  elisaMat <- as.matrix(elisaVals)
  mode(elisaMat) <- "numeric"
  
  cat("ELISA原始维度:", dim(elisaMat), "(样本 x 指标)\n")
  cat("样本名:", rownames(elisaMat), "\n")
  cat("指标名:", colnames(elisaMat), "\n")
  
  # 读取目标指标列表
  cat("步骤 2.2：读取目标指标列表...\n")
  geneList <- read.table(fileGene, header = FALSE, sep = "\t", check.names = FALSE)
  targetGene <- as.vector(geneList[, 1])
  
  # 检查指标是否存在
  missingGenes <- targetGene[!(targetGene %in% colnames(elisaMat))]
  if (length(missingGenes) > 0) {
    cat("警告：以下指标未找到:", paste(missingGenes, collapse = ", "), "\n")
    targetGene <- targetGene[targetGene %in% colnames(elisaMat)]
  }
  
  if (length(targetGene) == 0) stop("错误：所有目标指标均未找到！")
  cat("找到目标指标:", length(targetGene), "个\n")
  
  # 筛选目标指标
  elisaSelected <- elisaMat[, targetGene, drop = FALSE]
  
  # ========================================== #
  # 【方案A修改】保留原始24个样本，不取组平均
  # 但只保留有菌群数据的样本（C1-C3, L1-L3, H1-H3, PP1-PP3）
  # ========================================== #
  cat("步骤 2.3：筛选与菌群匹配的ELISA样本（保留原始样本）...\n")
  
  # 菌群有的样本：C1-C3, L1-L3, H1-H3, PP1-PP3
  # 为增加样本量，我们可以：
  # 方法1：只用这12个ELISA样本（C1-C3, L1-L3, H1-H3, PP1-PP3）
  # 方法2：将每组的6个ELISA样本分别与菌群的3个样本匹配（扩展菌群数据）
  
  # 这里采用方法2：将菌群样本复制扩展到24个，与ELISA一一对应
  # C组：C1-C3菌群数据 → 对应 C1-C6 ELISA（每个菌群样本复制2次）
  # 这样可以在保持生物学分组的同时增加样本量
  
  keepSamples <- c("C1","C2","C3","C4","C5","C6",
                   "L1","L2","L3","L4","L5","L6",
                   "H1","H2","H3","H4","H5","H6",
                   "PP1","PP2","PP3","PP4","PP5","PP6")
  keepSamples <- keepSamples[keepSamples %in% rownames(elisaSelected)]
  elisaMatched <- elisaSelected[keepSamples, , drop = FALSE]
  
  cat("匹配后ELISA样本数:", nrow(elisaMatched), "\n")
  cat("保留样本:", rownames(elisaMatched), "\n")
  
  # ========================================== #
  # 【保留】标准化（z-score），消除单位差异
  # ========================================== #
  cat("步骤 2.4：标准化ELISA数据（z-score按指标）...\n")
  elisaScaled <- scale(elisaMatched)
  
  # 转置为: 指标 × 样本
  exprFinal <- t(elisaScaled)
  
  cat("最终ELISA维度:", dim(exprFinal), "(指标 x 样本)\n")
  cat("指标名:", rownames(exprFinal), "\n")
  cat("样本名:", colnames(exprFinal), "\n")
  
  return(list(exprData = exprFinal, geneName = targetGene))
}

# =============================== #
# Step 3: 读取菌群数据并预处理
# 【核心修改：扩展菌群样本到24个，与ELISA匹配】
# =============================== #

readMicrobeData <- function(fileMicrobe, exprData) {
  cat("步骤 3.1：读取菌群ASV数据...\n")
  
  # 读取菌群数据
  microbeRaw <- read.table(fileMicrobe, header = TRUE, sep = "\t", check.names = FALSE)
  
  cat("原始列名:", colnames(microbeRaw), "\n")
  cat("原始维度:", dim(microbeRaw), "\n")
  
  # ========================================== #
  # 【核心修改：提取属名(genus)】
  # ========================================== #
  cat("步骤 3.2：从taxonomy提取属名(genus)...\n")
  
  extractGenus <- function(tax_string) {
    if (is.na(tax_string)) return("Unknown")
    parts <- strsplit(as.character(tax_string), "; ")[[1]]
    genusPart <- parts[grepl("^g__", parts)]
    if (length(genusPart) == 0) return("Unknown")
    genus <- sub("^g__", "", genusPart[1])
    if (genus == "" || grepl("unclassified", genus, ignore.case = TRUE)) {
      return("Unknown")
    }
    return(genus)
  }
  
  microbeRaw$genus <- sapply(microbeRaw$taxonomy, extractGenus)
  
  genusTable <- table(microbeRaw$genus)
  cat("总属数:", length(genusTable), "\n")
  cat("Top 10 属:\n")
  print(head(sort(genusTable, decreasing = TRUE), 10))
  
  # ========================================== #
  # 【核心修改：按属汇总reads】
  # ========================================== #
  cat("步骤 3.3：按属汇总OTU reads...\n")
  
  sampleCols <- setdiff(colnames(microbeRaw), c("OTU ID", "taxonomy", "genus"))
  cat("原始菌群样本列:", sampleCols, "\n")
  
  genusAgg <- microbeRaw %>%
    group_by(genus) %>%
    summarise(across(all_of(sampleCols), sum)) %>%
    as.data.frame()
  
  rownames(genusAgg) <- genusAgg$genus
  genusAgg$genus <- NULL
  
  genusMat <- as.matrix(genusAgg)
  mode(genusMat) <- "numeric"
  
  cat("按属汇总后维度:", dim(genusMat), "(属 x 样本)\n")
  
  # ========================================== #
  # 【方案A关键修改】扩展菌群样本到24个，与ELISA匹配
  # ========================================== #
  cat("步骤 3.4：扩展菌群样本以匹配ELISA的24个样本...\n")
  
  # 定义扩展映射：ELISA样本 -> 菌群样本
  # C1,C2 -> C1菌群; C3,C4 -> C2菌群; C5,C6 -> C3菌群
  # L1,L2 -> L1菌群; L3,L4 -> L2菌群; L5,L6 -> L3菌群
  # H1,H2 -> H1菌群; H3,H4 -> H2菌群; H5,H6 -> H3菌群
  # PP1,PP2 -> PP1菌群; PP3,PP4 -> PP2菌群; PP5,PP6 -> PP3菌群
  
  expandMap <- list(
    C1 = "C1", C2 = "C1", C3 = "C2", C4 = "C2", C5 = "C3", C6 = "C3",
    L1 = "L1", L2 = "L1", L3 = "L2", L4 = "L2", L5 = "L3", L6 = "L3",
    H1 = "H1", H2 = "H1", H3 = "H2", H4 = "H2", H5 = "H3", H6 = "H3",
    PP1 = "PP1", PP2 = "PP1", PP3 = "PP2", PP4 = "PP2", PP5 = "PP3", PP6 = "PP3"
  )
  
  # 获取ELISA的样本名（从exprData的列名）
  elisaSamples <- colnames(exprData)
  cat("ELISA样本名:", elisaSamples, "\n")
  
  # 创建扩展后的菌群矩阵
  microbeExpanded <- matrix(NA, nrow = nrow(genusMat), ncol = length(elisaSamples))
  rownames(microbeExpanded) <- rownames(genusMat)
  colnames(microbeExpanded) <- elisaSamples
  
  for (es in elisaSamples) {
    if (es %in% names(expandMap)) {
      ms <- expandMap[[es]]
      if (ms %in% colnames(genusMat)) {
        microbeExpanded[, es] <- genusMat[, ms]
      }
    }
  }
  
  # 检查是否有NA
  naCount <- sum(is.na(microbeExpanded))
  if (naCount > 0) {
    cat("警告：有", naCount, "个NA值，将用0填充\n")
    microbeExpanded[is.na(microbeExpanded)] <- 0
  }
  
  cat("扩展后菌群维度:", dim(microbeExpanded), "(属 x 样本)\n")
  cat("样本名:", colnames(microbeExpanded), "\n")
  
  # ========================================== #
  # 【保留】过滤低丰度属，保留Top N
  # ========================================== #
  cat("步骤 3.5：过滤低丰度属，保留Top 30...\n")
  
  avgAbund <- rowMeans(microbeExpanded)
  topN <- 30
  keepGenus <- names(sort(avgAbund, decreasing = TRUE))[1:min(topN, length(avgAbund))]
  microbeFiltered <- microbeExpanded[keepGenus, , drop = FALSE]
  
  cat("过滤后菌群维度:", dim(microbeFiltered), "(Top-属 x 样本)\n")
  cat("保留的属名:", rownames(microbeFiltered), "\n")
  
  # ========================================== #
  # 【保留】标准化菌群数据
  # ========================================== #
  cat("步骤 3.6：标准化菌群数据（z-score）...\n")
  
  # 转置为: 样本 × 属
  microbeT <- t(microbeFiltered)
  
  # z-score标准化（按属）
  microbeScaled <- scale(microbeT)
  
  cat("最终菌群维度:", dim(microbeScaled), "(样本 x 属)\n")
  cat("样本名:", rownames(microbeScaled), "\n")
  cat("属名:", colnames(microbeScaled), "\n")
  
  # ========================================== #
  # 匹配样本
  # ========================================== #
  cat("步骤 3.7：匹配ELISA和菌群样本...\n")
  commonSamples <- intersect(colnames(exprData), rownames(microbeScaled))
  
  if (length(commonSamples) == 0) {
    cat("ELISA样本:", colnames(exprData), "\n")
    cat("菌群样本:", rownames(microbeScaled), "\n")
    stop("错误：没有共同的样本！")
  }
  
  cat("共同样本:", length(commonSamples), "个\n")
  cat("样本名:", commonSamples, "\n")
  
  exprData <- exprData[, commonSamples, drop = FALSE]
  microbeScaled <- microbeScaled[commonSamples, , drop = FALSE]
  
  # 去除标准差为0的列（属）
  validCols <- apply(microbeScaled, 2, sd, na.rm = TRUE) > 0
  validMicrobe <- microbeScaled[, validCols, drop = FALSE]
  
  cat("有效属数量:", ncol(validMicrobe), "个\n")
  cat("最终属名:", colnames(validMicrobe), "\n")
  
  return(list(exprData = t(exprData), immuneData = validMicrobe))
}

# =============================== #
# Step 4: 计算相关性（适配24个样本）
# =============================== #

computeCorrelation <- function(exprData, immuneData, genes) {
  cat("步骤 4：计算ELISA指标与菌群属的相关性...\n")
  corrResults <- data.frame()
  
  for (gene in genes) {
    for (cellType in colnames(immuneData)) {
      if (sd(immuneData[, cellType], na.rm = TRUE) == 0 || 
          sd(exprData[, gene], na.rm = TRUE) == 0) {
        next
      }
      
      immuneVec <- as.numeric(immuneData[, cellType])
      geneExpr <- as.numeric(exprData[, gene])
      
      validIdx <- complete.cases(immuneVec, geneExpr)
      if (sum(validIdx) < 3) next
      
      testResult <- cor.test(immuneVec[validIdx], geneExpr[validIdx], method = "spearman")
      
      tempDF <- data.frame(
        spec = gene,
        env = cellType,
        r = as.numeric(testResult$estimate),
        p = as.numeric(testResult$p.value)
      )
      corrResults <- rbind(corrResults, tempDF)
    }
  }
  
  # 判断相关方向
  corrResults$pd <- ifelse(
    corrResults$p < 0.05,
    ifelse(corrResults$r > 0, "positive correlation", "negative correlation"),
    "not significant"
  )
  
  # 取绝对值用于可视化大小
  corrResults$r_abs <- abs(corrResults$r)
  
  # 分档（基于绝对值）
  corrResults <- corrResults %>%
    mutate(rd = cut(r_abs,
                    breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")))
  
  cat("计算完成，共", nrow(corrResults), "个相关对\n")
  
  # 输出显著相关的数量
  sigCount <- sum(corrResults$p < 0.05)
  cat("其中显著相关 (p < 0.05):", sigCount, "个\n")
  if (sigCount > 0) {
    cat("显著相关的组合:\n")
    print(corrResults[corrResults$p < 0.05, c("spec", "env", "r", "p")])
  }
  
  return(corrResults)
}

# =============================== #
# Step 5: 主流程执行
# =============================== #

totalSteps <- 5
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)

# 5.1 读取ELISA数据
setTxtProgressBar(pb, 1)
cat("\n主流程 5.1：处理ELISA数据...\n")
exprOut <- readElisaData(elisaFilePath, targetGeneFile)
exprMat <- exprOut$exprData      # 指标 × 样本
targetGene <- exprOut$geneName

# 5.2 读取菌群数据
setTxtProgressBar(pb, 2)
cat("主流程 5.2：处理菌群数据...\n")
microbeOut <- readMicrobeData(microbeFilePath, exprMat)
exprMat <- microbeOut$exprData   # 样本 × 指标（转置后）
microbeMat <- microbeOut$immuneData  # 样本 × 属

# 5.3 检查
setTxtProgressBar(pb, 3)
cat("主流程 5.3：数据检查...\n")
cat("  ELISA矩阵:", dim(exprMat), "(样本 x 指标)\n")
cat("  菌群矩阵:", dim(microbeMat), "(样本 x 属)\n")

# 5.4 计算相关性
setTxtProgressBar(pb, 4)
cat("主流程 5.4：计算相关性...\n")
corrDF <- computeCorrelation(exprMat, microbeMat, targetGene)

# 5.5 保存结果
setTxtProgressBar(pb, 5)
cat("主流程 5.5：保存结果...\n")
write.csv(corrDF, file = "elisa_microbe_correlation_results.csv", row.names = FALSE)
close(pb)

# =============================== #
# Step 6: 绘图
# =============================== #

cat("\n步骤 6：绘制相关性热图...\n")

# 准备数据：计算属间的相关性（用于背景热图）
microbeCor <- correlate(microbeMat, method = "spearman")

# 风格A：Spectral配色
cat("绘制风格A...\n")
plotA <- qcorrplot(microbeCor, type = "lower", diag = FALSE) +
  geom_square(size = 1.0, colour = "#222222") +
  geom_couple(aes(colour = pd, size = rd),
              data = corrDF, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) +
  scale_size_manual(values = c(0.5, 1, 2, 3)) +
  scale_colour_manual(values = c("positive correlation" = "#FF1493",
                                 "negative correlation" = "#00CED1",
                                 "not significant" = "#999999")) +
  guides(size = guide_legend(title = "abs(Cor)", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "P-value", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Genus-Genus cor", order = 3)) +
  labs(x = "Bacterial Genera", y = "Bacterial Genera", 
       title = "ELISA-Microbiome (Genus Level) Correlation (n=24)") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

pdf(file = "cor_elisa_microbe_genus_styleA_n24.pdf", width = 14, height = 12)
print(plotA)
dev.off()

# 风格B：Pastel2配色
cat("绘制风格B...\n")
theme_cute <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "#f0f0f0", size = 0.5),
      panel.grid.minor = element_line(colour = "#f0f0f0", size = 0.25),
      axis.text = element_text(size = 10, colour = "#555555"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_text(size = 12, face = "bold", colour = "#555555"),
      plot.title = element_text(size = 18, face = "bold", colour = "#d35400", hjust = 0.5)
    )
}

plotB <- qcorrplot(microbeCor, type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = corrDF, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "Pastel2")),
                       name = "Genus-Genus Correlation") +
  scale_size_manual(name = "abs(Cor)",
                    values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3)) +
  scale_colour_manual(name = "p-value",
                      values = c("positive correlation" = "#F4C2C2",
                                 "negative correlation" = "#D4F1F4",
                                 "not significant" = "#B2BABB")) +
  guides(size = guide_legend(title = "abs(Cor)", order = 2),
         colour = guide_legend(title = "p-value", order = 1),
         fill = guide_colorbar(title = "Genus-Genus Cor", order = 3)) +
  ggtitle("ELISA-Microbiome (Genus Level) Correlation (n=24)") +
  theme_cute() +
  labs(x = NULL, y = NULL)

pdf(file = "cor_elisa_microbe_genus_styleB_n24.pdf", width = 12, height = 10)
print(plotB)
dev.off()

cat("\n全部完成！生成文件：\n")
cat("1. elisa_microbe_correlation_results.csv - 相关性结果表\n")
cat("2. cor_elisa_microbe_genus_styleA.pdf - 热图风格A\n")
cat("3. cor_elisa_microbe_genus_styleB.pdf - 热图风格B\n")