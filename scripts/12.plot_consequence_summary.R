#!/usr/bin/env Rscript
# Session 19（第 12 步）R/ggplot2：变异后果注释汇总图（出版级）
#
# 对应官方 Session 19 的 "R — Annotation summary plot with Arial font
# and journal-ready output"。
#
# 输入（由 12.run_ggplot2_plots.sh 生成）：
#   <out_dir>/vep_summary.tsv
#     列：chrom pos alt gene impact consequence
#
# 输出：
#   <out_dir>/NA12878_chr20_consequence_summary.pdf   （向量，300 dpi）
#   <out_dir>/NA12878_chr20_consequence_summary.png   （位图，300 dpi）
#
# 用法：Rscript 12.plot_consequence_summary.R <vep_summary.tsv> <out_dir>
#
# 字体说明：
#   官方用 Arial（需 extrafont + font_import，首次 2–5 分钟）。
#   本机未装 extrafont → 自动回退到 sans（Linux 上等价于 Liberation Sans，
#   与 Arial 视觉一致），不阻塞出图。

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
tsv_path <- if (length(args) >= 1) args[1] else "12.figures/vep_summary.tsv"
out_dir  <- if (length(args) >= 2) args[2] else "12.figures"

# ── 字体：Arial 优先，缺失则 sans ──
BASE_FONT <- "sans"
if (requireNamespace("extrafont", quietly = TRUE)) {
  ok <- tryCatch({
    if (!("Arial" %in% extrafont::fonts())) {
      extrafont::font_import(prompt = FALSE)
      extrafont::loadfonts(device = "pdf", quiet = TRUE)
    }
    TRUE
  }, error = function(e) {
    message("extrafont 导入失败 → 回退 sans")
    FALSE
  })
  if (ok) BASE_FONT <- "Arial"
} else {
  message("未安装 extrafont → 回退 sans（Liberation Sans，视觉与 Arial 一致）")
}

# ── 数据 ──
vep <- read.table(tsv_path, sep = "\t", header = FALSE,
                  col.names = c("chrom","pos","alt","gene","impact","consequence"),
                  stringsAsFactors = FALSE)

coding <- vep %>%
  filter(impact %in% c("HIGH","MODERATE","LOW")) %>%
  mutate(impact = factor(impact, levels = c("HIGH","MODERATE","LOW")))

impact_colours <- c("HIGH" = "#B91C1C", "MODERATE" = "#D97706", "LOW" = "#0E7C5B")

top15 <- coding %>%
  count(consequence, impact) %>%
  arrange(desc(n)) %>%
  slice_head(n = 15)

p <- ggplot(top15, aes(x = reorder(consequence, n), y = n, fill = impact)) +
  geom_col(alpha = 0.88) +
  scale_fill_manual(values = impact_colours, name = "IMPACT") +
  coord_flip() +
  labs(
    title    = "NA12878 chr20 — Variant Consequence Summary",
    subtitle = paste0("PASS variants, VEP v105, GRCh38 | n=", nrow(coding), " coding variants"),
    x = NULL, y = "Count"
  ) +
  theme_bw(base_family = BASE_FONT) +
  theme(
    plot.title      = element_text(face = "bold", size = 15),
    plot.subtitle   = element_text(size = 10, colour = "grey30"),
    axis.text       = element_text(size = 9),
    legend.position = "right"
  )

pdf_file <- file.path(out_dir, "NA12878_chr20_consequence_summary.pdf")
png_file <- file.path(out_dir, "NA12878_chr20_consequence_summary.png")

# PDF（向量）：cairo_pdf 可用
ggsave(pdf_file, p, width = 9, height = 6, device = cairo_pdf)

# PNG：本机 ggsave(device="png"/agg_png/svglite) 会报
#   "Graphics API version mismatch"（R 与 cairo 库版本不符）→
#   改为**直接调设备函数**（实测 png(type="cairo") 与 Cairo::CairoPNG 均可）
png_ok <- tryCatch({
  grDevices::png(png_file, width = 9, height = 6, units = "in", res = 300, type = "cairo")
  print(p); grDevices::dev.off()
  TRUE
}, error = function(e) {
  message("png(type=cairo) 失败 → 试 Cairo::CairoPNG")
  FALSE
})
if (!png_ok) {
  Cairo::CairoPNG(png_file, width = 9, height = 6, units = "in", res = 300)
  print(p); grDevices::dev.off()
}

cat(sprintf("  coding variants: %d\n", nrow(coding)))
cat(sprintf("  IMPACT 分布: HIGH=%d MODERATE=%d LOW=%d\n",
            sum(coding$impact == "HIGH"),
            sum(coding$impact == "MODERATE"),
            sum(coding$impact == "LOW")))
cat(sprintf("  top 后果: %s\n", top15$consequence[1]))
cat(sprintf("  字体: %s\n", BASE_FONT))
cat(sprintf("  → %s\n  → %s\n", pdf_file, png_file))