#!/usr/bin/env bash
# 下教程 Session 10/13 需要的 GATK 资源包（除 dbSNP 外）
# dbSNP138 在 S3 只有未压缩版 10 GB，按"不浪费流量"原则跳过
# Mills/1000G_phase1/HapMap/1000G_omni2.5 都是压缩版，体量合理
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
RAW=$WS/rawData
LOGDIR=$WS/logs
LOG=$LOGDIR/resources_dl_v2.log
mkdir -p "$RAW/resources" "$LOGDIR"

export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999

BASE="https://s3.amazonaws.com/broad-references/hg38/v0"

cd "$RAW/resources"

for f in \
  "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
  "1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
  "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" \
  "hapmap_3.3.hg38.vcf.gz" \
  "hapmap_3.3.hg38.vcf.gz.tbi" \
  "1000G_omni2.5.hg38.vcf.gz" \
  "1000G_omni2.5.hg38.vcf.gz.tbi"; do
  if [ -f "$f" ] && [ -s "$f" ]; then
    echo "[$(date)] 已存在，跳过: $f" >> "$LOG"
    continue
  fi
  echo "[$(date)] 开始下载: $f" >> "$LOG"
  wget -c --tries=10 --retry-connrefused --waitretry=10 -a "$LOG" \
    "${BASE}/${f}" -O "$f"
  echo "[$(date)] 完成: $f ($(ls -lh $f | awk '{print $5}'))" >> "$LOG"
done

echo "[$(date)] 资源包（除 dbSNP）下载完成" >> "$LOG"
ls -lh "$RAW/resources" | tee -a "$LOG"