#!/usr/bin/env bash
# 教程 Session 02 Option A：从 CRAM 抽取 chr20 → paired FASTQ
# 教程原文命令（已按 rawData/ 路径适配）：
#   samtools view -b -T Homo_sapiens_assembly38.fasta -@ 8 NA12878.final.cram chr20 | \
#     samtools sort -n -@ 8 | \
#     samtools fastq -@ 8 \
#       -1 NA12878_chr20_R1.fastq.gz -2 NA12878_chr20_R2.fastq.gz \
#       -0 /dev/null -s /dev/null -n -
#
# 参数说明：
#   view -b   : 输出 BAM 格式（管道流）
#   view -T   : 参考 FASTA（CRAM 解码必需）
#   view -@ 8 : 8 线程解码
#   view chr20: 只抽取 chr20（教程用 chr20 做示例，chr20 占基因组 ~2%）
#   sort -n   : 按 read 名排序（paired FASTQ 输出要求 R1/R2 相邻）
#   sort -@ 8 : 8 线程排序
#   fastq -1/-2 : R1/R2 输出文件名
#   fastq -0 : 单端 read → /dev/null（丢弃）
#   fastq -s : singleton → /dev/null（丢弃）
#   fastq -n : 保留 read 原始名称（不加 /1 /2 后缀）
#   fastq -@ 8: 8 线程压缩
#
# 估计耗时：~25–40 分钟（CRAM 解码 + gzip 写盘是瓶颈）
# 磁盘：流式不写中间 BAM，峰值 ~19 GB（CRAM 15 GB + 输出 4 GB）
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
RAW=$WS/rawData
LOG=$WS/logs/extract_chr20_fastq.log
mkdir -p "$WS/logs"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH
export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999

CRAM="$RAW/fastq/NA12878.final.cram"
CRAI="$RAW/fastq/NA12878.final.cram.crai"
REF="$RAW/ref/Homo_sapiens_assembly38.fasta"

# 前置检查
[ -f "$CRAM" ] || { echo "ERROR: CRAM not found: $CRAM"; exit 1; }
[ -f "$CRAI" ] || { echo "ERROR: CRAI not found: $CRAI"; exit 1; }
[ -f "$REF" ]  || { echo "ERROR: REF not found: $REF"; exit 1; }

cd "$RAW/fastq"

echo "[$(date)] === Start extract chr20 from CRAM ===" | tee "$LOG"
echo "CRAM: $CRAM ($(du -h $CRAM | cut -f1))" | tee -a "$LOG"
echo "REF:  $REF ($(du -h $REF | cut -f1))" | tee -a "$LOG"
echo "CRAI: $CRAI" | tee -a "$LOG"
echo "Threads: 8 (教程推荐值，本机 64 核用8 核足够)" | tee -a "$LOG"
echo | tee -a "$LOG"

# 1) samtools quickcheck 验证 CRAM 完整性
echo "[$(date)] samtools quickcheck CRAM..." | tee -a "$LOG"
if samtools quickcheck "$CRAM"; then
  echo "  ✓ CRAM OK" | tee -a "$LOG"
else
  echo "  ✗ CRAM corrupted" | tee -a "$LOG"
  exit 1
fi

# 2) 流式管道：CRAM → chr20 BAM → name-sort → paired FASTQ
echo "[$(date)] Start pipeline: view → sort -n → fastq" | tee -a "$LOG"
samtools view -b -T "$REF" -@ 8 "$CRAM" chr20 2>>"$LOG" | \
  samtools sort -n -@ 8 2>>"$LOG" | \
  samtools fastq -@ 8 \
    -1 NA12878_chr20_R1.fastq.gz \
    -2 NA12878_chr20_R2.fastq.gz \
    -0 /dev/null \
    -s /dev/null \
    -n - 2>>"$LOG"

echo "[$(date)] Pipeline done" | tee -a "$LOG"
ls -lh NA12878_chr20_R*.fastq.gz | tee -a "$LOG"

# 3) 校验 R1/R2 read 数必须完全相等
echo "[$(date)] Verify R1/R2 read counts..." | tee -a "$LOG"
R1=$(zcat NA12878_chr20_R1.fastq.gz | wc -l | awk '{print int($1/4)}')
R2=$(zcat NA12878_chr20_R2.fastq.gz | wc -l | awk '{print int($1/4)}')
echo "  R1 reads: $R1" | tee -a "$LOG"
echo "  R2 reads: $R2" | tee -a "$LOG"
if [ "$R1" -eq "$R2" ] && [ "$R1" -gt 0 ]; then
  echo "  ✓ Paired: $R1 reads each" | tee -a "$LOG"
else
  echo "  ✗ MISMATCH: R1=$R1 R2=$R2 — files corrupted" | tee -a "$LOG"
  exit 1
fi

echo | tee -a "$LOG"
echo "[$(date)] === Done ===" | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
ls -lh NA12878_chr20_R1.fastq.gz NA12878_chr20_R2.fastq.gz | tee -a "$LOG"
echo | tee -a "$LOG"
echo "CRAM 仍保留（你选择 A：不删 CRAM）。如需删除释放 15 GB：" | tee -a "$LOG"
echo "  rm $CRAM $CRAI" | tee -a "$LOG"