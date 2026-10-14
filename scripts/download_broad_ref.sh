#!/usr/bin/env bash
# 只下载教程指定的 Broad GATK GRCh38 参考基因组
# 源：AWS S3 broad-references 桶
# 不下 BWA 索引、不下资源包
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
RAW=$WS/rawData
LOGDIR=$WS/logs
LOG=$LOGDIR/broad_ref_dl.log
mkdir -p "$RAW/ref" "$LOGDIR"

export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999

BASE="https://s3.amazonaws.com/broad-references/hg38/v0"

cd "$RAW/ref"

# 参考基因组 FASTA + 配套索引（samtools .fai + GATK .dict）
for f in \
  "Homo_sapiens_assembly38.fasta" \
  "Homo_sapiens_assembly38.fasta.fai" \
  "Homo_sapiens_assembly38.dict"; do
  if [ -f "$f" ] && [ -s "$f" ]; then
    echo "[$(date)] 已存在，跳过: $f" >> "$LOG"
    continue
  fi
  echo "[$(date)] 开始下载: $f" >> "$LOG"
  wget -c --tries=10 --retry-connrefused --waitretry=10 -a "$LOG" \
    "${BASE}/${f}" -O "$f"
  echo "[$(date)] 完成: $f ($(ls -lh $f | awk '{print $5}'))" >> "$LOG"
done

echo "[$(date)] 参考基因组下载完成" >> "$LOG"
ls -lh "$RAW/ref" | tee -a "$LOG"