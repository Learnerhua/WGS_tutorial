#!/usr/bin/env bash
# 纯下载 NA12878 CRAM（教程 Session 02 下载源 ENA ERR3239334）
# 只下，不抽 FASTQ，不做 samtools 调用，不联动后续流程
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
RAW=$WS/rawData
LOGDIR=$WS/logs
LOG=$LOGDIR/download_cram.log
mkdir -p "$RAW/fastq" "$LOGDIR"

export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999

EBI="https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334"

cd "$RAW/fastq"
wget -c --tries=10 --retry-connrefused --waitretry=10 -o "$LOG" \
  "${EBI}/NA12878.final.cram" -O NA12878.final.cram
wget -c -o "$LOG" \
  "${EBI}/NA12878.final.cram.crai" -O NA12878.final.cram.crai

ls -lh "$RAW/fastq"