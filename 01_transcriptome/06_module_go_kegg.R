## Reproducible analysis code accompanying the manuscript.
## Run from the repository root (see README.md); paths are resolved from config.R.
if (!file.exists('config.R')) stop('Please run this script from the repository root.')
source('config.R')

#加载
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
library(stringr)
R.utils::setOption("clusterProfiler.download.method","auto")
#其他物种
#https://www.genome.jp/kegg/catalog/org_list.html
#https://www.jianshu.com/p/84e70566a6c6

#工作目录
setwd(file.path(DATA_ROOT, "3.转录组", "5.WGCNA"))

#阈值
pvalueFilter=0.05       #p值过滤条件
qvalueFilter=1       #矫正后的p值过滤条件

#颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")

#读入
rt=read.table("id.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#转换因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
#去除NA
gene=entrezIDs[entrezIDs!="NA"]
#可能存在多个entrezIDs
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

#保存结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#GO的数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="GObarplot.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GObubble.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

#################GO圈图整理。可以从由开始的id.txt生成GO.txt开始，也可以使用整理好的GO圈图.txt从这里开始###################
# ========== 加载包 ==========
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(stringr)

# ========== 设置工作路径 ==========
path <- file.path(DATA_ROOT, "3.转录组", "5.WGCNA")

# ========== 读入数据 ==========
GO <- read.delim(file.path(path, "GO圈图.txt"), 
                 sep = "\t", 
                 header = TRUE, 
                 stringsAsFactors = FALSE, 
                 check.names = FALSE)

# ========== 参数 ==========
ontology.col <- c("#E05C5C", "#5B9BD5", "#70AD47")

# ========== 数据整理 ==========
data <- GO[order(GO$pvalue),]
datasig <- data[data$pvalue < 0.05, , drop = FALSE]

BP <- datasig[datasig$ONTOLOGY == "BP", , drop = FALSE]
CC <- datasig[datasig$ONTOLOGY == "CC", , drop = FALSE]
MF <- datasig[datasig$ONTOLOGY == "MF", , drop = FALSE]

BP <- head(BP, 9)
CC <- head(CC, 9)
MF <- head(MF, 9)

data <- rbind(BP, CC, MF)
main.col <- ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio <- Gene / BgGene
logpvalue <- -log(data$pvalue, 10)

logpvalue.col <- brewer.pal(n = 8, name = "Reds")
f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10, 15, 20), colors = logpvalue.col)
BgGene.col <- f(logpvalue)

# 创建Description查找表（关键！用ID作为键）
desc_lookup <- setNames(data$Description, data$ID)

# ========== 构建circos数据框（用ID作为sector标识）==========
df <- data.frame(GO = data$ID, start = 1, end = max(BgGene))
rownames(df) <- df$GO

bed2 <- data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 <- data.frame(GO = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 <- data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5

# ========== 画图 ==========
pdf(file.path(path, "GO.circlize.pdf"), width = 14, height = 14)
par(omi = c(0.1, 0.1, 0.1, 1.5))

circos.genomicInitialize(df, plotType = "none")

# ---------- 第1轨道：GO ID标签 + 背景色 ----------
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, 
              cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.5, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 50), labels.facing = "clockwise")
}

# ---------- 第2轨道：背景基因数 ----------
circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, 
                                         col = value[, 2], border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], 
                                         adj = 0, cex = 0.7, ...)
                    })

# ---------- 第3轨道：选中基因数 ----------
circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, 
                                         col = '#7A6FA8', border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], 
                                         cex = 0.8, adj = 0, ...)
                    })

# ---------- 第4轨道（最内圈）：rich factor + Description ----------
circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, 
                    bg.border = "white", bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      
                      cell.xlim <- get.cell.meta.data("cell.xlim")
                      cell.ylim <- get.cell.meta.data("cell.ylim")
                      
                      # 画网格线
                      for (j in 1:9) {
                        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      
                      # 画rich factor条形
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[, 1], 
                                         col = value[, 2], border = NA, ...)
                      
                      # ========== 新增：显示Description（参考代码逻辑）==========
                      sector.name <- get.cell.meta.data("sector.index")  # 获取GO ID
                      desc <- desc_lookup[sector.name]                    # 查找Description
                      desc <- str_wrap(desc, width = 25)                  # 自动换行
                      
                      xcenter <- (region[[1]] + region[[2]]) / 2          # sector中心
                      
                      circos.text(x = xcenter,
                                  y = 0,                                 # 放在轨道底部
                                  labels = desc,
                                  facing = "clockwise",                  # 沿圆周方向
                                  niceFacing = TRUE,
                                  adj = c(0, 0.5),
                                  cex = 1.0,      # ← 改这里：调大字体（如0.8, 1.0）
                                  font = 2,       # ← 加这一行：2=粗体
                                  lineheight = 0.85)
                    })

circos.clear()

# ---------- 图例 ----------
middle.legend <- Legend(
  labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
  type = "points", pch = c(15, 15, 17),
  legend_gp = gpar(col = c('pink', '#7A6FA8', ontology.col[1])),
  title = "", nrow = 3, size = unit(3, "mm")
)
circle_size <- unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend <- Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),
  type = "points", pch = 15,
  legend_gp = gpar(col = ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3, size = unit(3, "mm"),
  grid_height = unit(5, "mm"), grid_width = unit(5, "mm")
)
logp.legend <- Legend(
  labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
  type = "points", pch = 16,
  legend_gp = gpar(col = logpvalue.col), title = "-log10(Pvalue)",
  title_position = "topcenter", grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm"), size = unit(3, "mm")
)
lgd <- packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")

dev.off()

#KEGG
#读入
rt=read.table("id.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#转换为id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
#去除NA
gene=entrezIDs[entrezIDs!="NA"]
#可能存在多个entrezIDs
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

#保存结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#显示通路
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="KEGGbarplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
pdf(file="KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

#################KEGG圈图整理。可以从由开始的id.txt生成KEGG.txt开始，也可以使用整理好的KEGG圈图.txt从这里开始###################
#KEGG圈图
# ========== 加载包 ==========
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(stringr)

# ========== 设置工作路径 ==========
path <- file.path(DATA_ROOT, "3.转录组", "5.WGCNA")

# ========== 读入已整理好的KEGG圈图数据 ==========
KEGG <- read.delim(file.path(path, "KEGG圈图.txt"), 
                   sep = "\t", 
                   header = TRUE, 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)

# ========== 参数 ==========
Pathway.col <- c("#7E8AA2", "#BFA6A0", "#D1C7B7")

# ========== 数据整理 ==========
showNum <- 30
data <- KEGG[order(KEGG$p.adjust), ]

if (nrow(KEGG) > showNum) {
  data <- data[1:showNum, ]
}

data$Pathway <- "KEGG"
main.col <- Pathway.col[as.numeric(as.factor(data$Pathway))]

BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio <- Gene / BgGene
logpvalue <- -log(data$pvalue, 10)

logpvalue.col <- brewer.pal(n = 8, name = "Reds")
f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10, 15, 20), colors = logpvalue.col)
BgGene.col <- f(logpvalue)

# 创建Description查找表（关键！用ID作为键）
desc_lookup <- setNames(data$Description, data$ID)

# ========== 构建circos数据框（用ID作为sector标识）==========
df <- data.frame(KEGG = data$ID, start = 1, end = max(BgGene))
rownames(df) <- df$KEGG

bed2 <- data.frame(KEGG = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 <- data.frame(KEGG = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 <- data.frame(KEGG = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5

# ========== 画图 ==========
pdf(file.path(path, "KEGG.circlize.pdf"), width = 14, height = 14)
par(omi = c(0.1, 0.1, 0.1, 1.5))

circos.genomicInitialize(df, plotType = "none")

# ---------- 第1轨道：KEGG ID标签 + 背景色 ----------
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, 
              cex = 1.0, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.5, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 50), labels.facing = "clockwise")
}

# ---------- 第2轨道：背景基因数 ----------
circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, 
                                         col = value[, 2], border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], 
                                         adj = 0, cex = 0.7, ...)
                    })

# ---------- 第3轨道：选中基因数 ----------
circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, 
                                         col = '#7A6FA8', border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], 
                                         cex = 0.8, adj = 0, ...)
                    })

# ---------- 第4轨道（最内圈）：rich factor + Description ----------
circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, 
                    bg.border = "white", bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      
                      cell.xlim <- get.cell.meta.data("cell.xlim")
                      cell.ylim <- get.cell.meta.data("cell.ylim")
                      
                      # 画网格线
                      for (j in 1:9) {
                        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      
                      # 画rich factor条形
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[, 1], 
                                         col = "#8FA876", border = NA, ...)
                      
                      # ========== 显示Description（同GO圈图逻辑）==========
                      sector.name <- get.cell.meta.data("sector.index")  # 获取KEGG ID
                      desc <- desc_lookup[sector.name]                    # 查找Description
                      desc <- str_wrap(desc, width = 25)                  # 自动换行
                      
                      xcenter <- (region[[1]] + region[[2]]) / 2          # sector中心
                      
                      circos.text(x = xcenter,
                                  y = 0,                                 # 轨道底部
                                  labels = desc,
                                  facing = "clockwise",                  # 沿圆周方向
                                  niceFacing = TRUE,
                                  adj = c(0, 0.5),
                                  cex = 1.0,                             # 字体大小
                                  font = 2,                              # 粗体
                                  lineheight = 0.85)
                    })

circos.clear()

# ---------- 图例 ----------
middle.legend <- Legend(
  labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
  type = "points", pch = c(15, 15, 17),
  legend_gp = gpar(col = c('pink', '#7A6FA8', "#8FA876")),
  title = "", nrow = 3, size = unit(3, "mm")
)
circle_size <- unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend <- Legend(
  labels = c("KEGG"),
  type = "points", pch = 15,
  legend_gp = gpar(col = Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3, size = unit(3, "mm"),
  grid_height = unit(5, "mm"), grid_width = unit(5, "mm")
)
logp.legend <- Legend(
  labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
  type = "points", pch = 16,
  legend_gp = gpar(col = logpvalue.col), title = "-log10(Pvalue)",
  title_position = "topcenter", grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm"), size = unit(3, "mm")
)
lgd <- packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")

dev.off()