source("config.R")  ## 统一路径配置，见 config.R`r`n###############################################################################
## style.R — 统一作图风格（使用用户指定的 8 色板）
##   #F7FAFC  #A7C7E7  #CFEAF6  #F4A7AE  #FBD9C9  #5878A6  #A6D0C6  #DADDE3
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

BG    <- "#F7FAFC"   # 背景
LBLUE <- "#A7C7E7"   # 浅蓝
PALE  <- "#CFEAF6"   # 淡青
PINK  <- "#F4A7AE"   # 粉（重点 / 核心基因）
PEACH <- "#FBD9C9"   # 蜜桃
BLUE  <- "#5878A6"   # 主蓝（处理组 / 强调）
MINT  <- "#A6D0C6"   # 薄荷（治疗组）
GREY  <- "#DADDE3"   # 浅灰（对照）

PAL8 <- c(BG, LBLUE, PALE, PINK, PEACH, BLUE, MINT, GREY)

## 连续色阶：蓝 -> 浅蓝 -> 淡青 -> 白 -> 蜜桃 -> 粉
GRAD_BWR <- c(BLUE, LBLUE, PALE, "#FFFFFF", PEACH, PINK)
## 热图（覆盖率等）：浅灰 -> 淡青 -> 浅蓝 -> 蓝
GRAD_BWO <- c(GREY, PALE, LBLUE, BLUE)

## 分组配色（全篇统一）
COL_GROUP <- c(Control = GREY, `Low-dose` = PALE, `Low-dose model` = PALE,
               `High-dose` = BLUE, `High-dose model` = BLUE,
               Treatment = MINT, Model = BLUE,
               Uninjured = GREY, Injured = BLUE,
               PFNA = BLUE, DMSO = GREY, case = BLUE, control = GREY)

theme_pub <- function(base_size = 9, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
      panel.border = element_rect(colour = "grey55", linewidth = 0.5),
      axis.text = element_text(colour = "grey15"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "grey30", hjust = 0),
      strip.background = element_rect(fill = GREY, colour = "grey55", linewidth = 0.5),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.title = element_text(face = "bold")
    )
}

tag_panel <- function(p, lab, size = 13) {
  p + labs(tag = lab) +
    theme(plot.tag = element_text(face = "bold", size = size),
          plot.tag.position = c(0, 1))
}
save_panel <- function(p, name, w = 1.9, h = 1.6) {
  ggsave(file.path(PANEL_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(PANEL_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 400)
}
