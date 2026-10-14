#!/usr/bin/env bash
# Session 05：FastQC + MultiQC 质控
# 教程 Session 05 主流程：对 chr20 paired FASTQ 跑 FastQC，再用 MultiQC 聚合
# 流程第 1 步（数据下载除外）：原始 FASTQ 质控
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

FASTQ_DIR=$WS/rawData/fastq
FASTQC_OUT=$ROOT/01.qc/fastqc
MULTIQC_OUT=$ROOT/01.qc/multiqc
LOG_DIR=$WS/logs
mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT" "$LOG_DIR"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

cd "$FASTQ_DIR"

echo "[$(date)] === 01. FastQC + MultiQC ===" | tee "$LOG_DIR/01.fastqc.log"
echo "输入: $FASTQ_DIR/${SAMPLE}_R*.fastq.gz" | tee -a "$LOG_DIR/01.fastqc.log"
echo "FastQC 输出: $FASTQC_OUT" | tee -a "$LOG_DIR/01.fastqc.log"
echo "MultiQC 输出: $MULTIQC_OUT" | tee -a "$LOG_DIR/01.fastqc.log"
echo | tee -a "$LOG_DIR/01.fastqc.log"

# 1) FastQC：教程默认单线程，本机 64 核用 -t 48 加速
echo "[$(date)] FastQC 开始..." | tee -a "$LOG_DIR/01.fastqc.log"
fastqc "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz" \
  -o "$FASTQC_OUT" \
  -t 48 \
  2>>"$LOG_DIR/01.fastqc.log"

echo "[$(date)] FastQC 完成" | tee -a "$LOG_DIR/01.fastqc.log"
ls -lh "$FASTQC_OUT" | tee -a "$LOG_DIR/01.fastqc.log"
echo | tee -a "$LOG_DIR/01.fastqc.log"

# 2) MultiQC 聚合
echo "[$(date)] MultiQC 开始..." | tee -a "$LOG_DIR/01.fastqc.log"
multiqc "$FASTQC_OUT" \
  -o "$MULTIQC_OUT" \
  --filename multiqc_report.html \
  --force \
  2>>"$LOG_DIR/01.fastqc.log"

echo "[$(date)] MultiQC 完成" | tee -a "$LOG_DIR/01.fastqc.log"
ls -lh "$MULTIQC_OUT" | tee -a "$LOG_DIR/01.fastqc.log"

echo | tee -a "$LOG_DIR/01.fastqc.log"
echo "报告路径:" | tee -a "$LOG_DIR/01.fastqc.log"
echo "  FastQC HTML: $FASTQC_OUT/${SAMPLE}_R1_fastqc.html" | tee -a "$LOG_DIR/01.fastqc.log"
echo "  MultiQC HTML: $MULTIQC_OUT/multiqc_report.html" | tee -a "$LOG_DIR/01.fastqc.log"