#!/usr/bin/env bash
# 自建 BWA-MEM2 索引
# 教程 Session 03: Building BWA-MEM2 index (~45 min, ~28 GB RAM peak)
# 本机: 2 TB RAM, 64 CPU → 估计 < 30 分钟
# 输出 5 个文件与 Broad 预构建版等价:
#   .0123, .amb, .ann, .bwt.2bit.64, .pac
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
RAW=$WS/rawData
LOGDIR=$WS/logs
LOG=$LOGDIR/bwa_mem2_index.log
mkdir -p "$RAW/ref" "$LOGDIR"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

cd "$RAW/ref"

FASTA="$RAW/ref/Homo_sapiens_assembly38.fasta"
[ -f "$FASTA" ] || { echo "ERROR: FASTA not found: $FASTA"; exit 1; }

# 教程原命令:
# nohup bwa-mem2 index Homo_sapiens_assembly38.fasta   2> ~/wgs_tutorial/logs/bwa_mem2_index.log &
# 教程原注释:
#   - 默认参数（-a 是 BWA-MEM2 默认的 BWT 算法 + 2-bit 编码，与 Broad 预构建版一致）
#   - 输出位置：FASTA 同目录下
#   - 增量构建：已存在文件会被跳过（如 .bwt.2bit.64 存在则不重建）

echo "[$(date)] Start bwa-mem2 index" | tee "$LOG"
echo "FASTA: $FASTA" | tee -a "$LOG"
echo "Log: $LOG" | tee -a "$LOG"
echo "Output files (5):" | tee -a "$LOG"
echo "  Homo_sapiens_assembly38.fasta.0123" | tee -a "$LOG"
echo "  Homo_sapiens_assembly38.fasta.amb" | tee -a "$LOG"
echo "  Homo_sapiens_assembly38.fasta.ann" | tee -a "$LOG"
echo "  Homo_sapiens_assembly38.fasta.bwt.2bit.64" | tee -a "$LOG"
echo "  Homo_sapiens_assembly38.fasta.pac" | tee -a "$LOG"
echo | tee -a "$LOG"

# 启动（前台运行，写到日志，主进程退出即完成）
bwa-mem2 index "$FASTA" 2>&1 | tee -a "$LOG"

echo "[$(date)] bwa-mem2 index complete" | tee -a "$LOG"
ls -lh "$RAW/ref"/Homo_sapiens_assembly38.fasta.* | tee -a "$LOG"