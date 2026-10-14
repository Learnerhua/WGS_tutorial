#!/usr/bin/env bash
# Session 09：比对后 QC — 覆盖度、插入片段、WGS 指标
# 流程第 4 步：比对后验证质量，再进入 Session 10 的 MarkDuplicates + BQSR
# 跑 mosdepth（覆盖度）+ Picard 3 个工具（InsertSize/AlignmentSummary/WGSMetrics）+ MultiQC 聚合
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

BAM=$ROOT/03.aligned/${SAMPLE}_sorted.bam
REF=$WS/rawData/ref/Homo_sapiens_assembly38.fasta
DICT=$WS/rawData/ref/Homo_sapiens_assembly38.dict
OUT_ROOT=$ROOT/04.post_alignment_qc
MOSD=$OUT_ROOT/mosdepth
PICARD=$OUT_ROOT/picard
MULTIQC=$OUT_ROOT/multiqc
LOG_DIR=$WS/logs
LOG=$LOG_DIR/04.post_alignment_qc.log
mkdir -p "$OUT_ROOT" "$MOSD" "$PICARD" "$MULTIQC" "$LOG_DIR"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

echo "[$(date)] === 04. Post-Alignment QC ===" | tee "$LOG"
echo "BAM: $BAM" | tee -a "$LOG"
echo "Output: $OUT_ROOT/" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 1) mosdepth：覆盖度与 breadth（流式算法，快） ───────────────────────
# --threads 48：BAM 读取是 IO 密集型；64 核机器用 48 线程
# --quantize 0:5:10:20:30:50:100:：深度分箱，用于算 % bases at ≥Nx（breadth）
# --by 500：500 bp 窗口覆盖度，用于均匀性可视化
# --mapq 20：排除 MAPQ<20 的 read（多映射/重复区）→ 与 GATK 变异调用的"有效覆盖"对齐
# --no-per-base：跳过 per-base BED（文件太大；窗口覆盖已足够 QC）
# --thresholds 1,5,10,20,30,50：在 summary 里直接报告各阈值的 % bases 占比
echo "[$(date)] mosdepth..." | tee -a "$LOG"
mosdepth \
  --threads 48 \
  --quantize 0:5:10:20:30:50:100: \
  --by 500 \
  --mapq 20 \
  --no-per-base \
  --thresholds 1,5,10,20,30,50 \
  "$MOSD/${SAMPLE}" \
  "$BAM" \
  2>>"$LOG"

echo "[$(date)] mosdepth done" | tee -a "$LOG"
ls -lh "$MOSD" | tee -a "$LOG"
echo | tee -a "$LOG"

# mosdepth summary 行格式：chrom  length  bases  mean  min  max
echo "=== mosdepth summary (chr20) ===" | tee -a "$LOG"
awk -F'\t' '$1=="chr20"{printf "Mean: %.2fx  Min: %sx  Max: %sx  Length: %s bp\n", $4, $5, $6, $2}' \
  "$MOSD/${SAMPLE}.mosdepth.summary.txt" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 2) Picard CollectInsertSizeMetrics ─────────────────────────────────────
# 反映文库片段化效果；PCR-free 预期为 350-420 bp 的钟形分布
echo "[$(date)] Picard CollectInsertSizeMetrics..." | tee -a "$LOG"
gatk CollectInsertSizeMetrics \
  -I "$BAM" \
  -O "$PICARD/${SAMPLE}_insert_size_metrics.txt" \
  -H "$PICARD/${SAMPLE}_insert_size_histogram.pdf" \
  --INCLUDE_DUPLICATES false \
  --MINIMUM_PCT 0.0 \
  2>>"$LOG"
echo | tee -a "$LOG"

# ── 3) Picard CollectAlignmentSummaryMetrics ──────────────────────────────
# 关键字段：PCT_PF_READS_ALIGNED（≥0.985）、STRAND_BALANCE（~0.50）、PF_MISMATCH_RATE（<0.005）
echo "[$(date)] Picard CollectAlignmentSummaryMetrics..." | tee -a "$LOG"
gatk CollectAlignmentSummaryMetrics \
  -R "$REF" \
  -I "$BAM" \
  -O "$PICARD/${SAMPLE}_alignment_summary.txt" \
  2>>"$LOG"
echo | tee -a "$LOG"

# ── 4) Picard CollectWgsMetrics（最全面的 WGS QC 工具） ───────────────────
# Fold-80 base penalty（均匀性指标）+ 覆盖度分布
# --COVERAGE_CAP 250：忽略极端深度（着丝粒/端粒），防止 mean 偏倚
# --MINIMUM_MAPPING_QUALITY 20：与 mosdepth --mapq 20 一致
# --MINIMUM_BASE_QUALITY 20：覆盖度计算的最低碱基质量
# --INTERVALS：限制到 chr20（否则 GENOME_TERRITORY=全基因组，mean 被稀释成 ~0.8x）
#   注意：GATK 4.6 不接受教程里的 `-L chr20` 简写，且要求 interval_list
#   必须带完整序列字典（与 BAM 的 @SQ 一致），故从 .dict + .fai 生成
ILIST=$OUT_ROOT/chr20.interval_list
{
  grep '^@HD' "$DICT" 2>/dev/null || printf '@HD\tVN:1.6\tSO:coordinate\n'
  grep '^@SQ' "$DICT" 2>/dev/null
  awk -F'\t' '$1=="chr20"{printf "%s\t1\t%s\t+\t%s\n", $1, $2, $1}' "$REF.fai"
} > "$ILIST"

echo "[$(date)] Picard CollectWgsMetrics..." | tee -a "$LOG"
gatk CollectWgsMetrics \
  -R "$REF" \
  -I "$BAM" \
  -O "$PICARD/${SAMPLE}_wgs_metrics.txt" \
  --COVERAGE_CAP 250 \
  --COUNT_UNPAIRED false \
  --MINIMUM_MAPPING_QUALITY 20 \
  --MINIMUM_BASE_QUALITY 20 \
  --INTERVALS "$ILIST" \
  2>>"$LOG"
echo | tee -a "$LOG"

# ── 5) 解析 Picard 指标 + QC gate（按列名取值，不依赖固定列号） ──────────
echo "=== Picard metrics ===" | tee -a "$LOG"
python3 - "$PICARD" "$SAMPLE" <<'PYEOF' | tee -a "$LOG"
import os
import sys

picard_dir, sample = sys.argv[1], sys.argv[2]


def read_table(path, header_key):
    """读取 Picard 指标文件；返回 (表头列表, 数据行列表)。"""
    header, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(header_key):
                header = line.rstrip('\n').split('\t')
            elif header and line.strip() and not line.startswith('#'):
                rows.append(line.rstrip('\n').split('\t'))
    return header, rows


def get(header, row, name, default='NA'):
    try:
        return row[header.index(name)]
    except (ValueError, TypeError, IndexError):
        return default


# ── 插入片段 ──
h, rows = read_table(os.path.join(picard_dir, f'{sample}_insert_size_metrics.txt'),
                     'MEDIAN_INSERT_SIZE')
r = rows[0]
print(f"Insert size:  median={get(h,r,'MEDIAN_INSERT_SIZE')} bp  "
      f"mean={float(get(h,r,'MEAN_INSERT_SIZE',0)):.1f} bp  "
      f"SD={float(get(h,r,'STANDARD_DEVIATION',0)):.1f} bp  "
      f"read_pairs={get(h,r,'READ_PAIRS')}")

# ── 比对汇总（PAIR 行） ──
h, rows = read_table(os.path.join(picard_dir, f'{sample}_alignment_summary.txt'),
                     'CATEGORY')
r = next((row for row in rows if row[0] == 'PAIR'), None)
if r:
    print(f"Alignment:    PCT_PF_READS_ALIGNED={float(get(h,r,'PCT_PF_READS_ALIGNED',0)):.4f}  "
          f"PF_MISMATCH_RATE={float(get(h,r,'PF_MISMATCH_RATE',0)):.5f}  "
          f"STRAND_BALANCE={float(get(h,r,'STRAND_BALANCE',0)):.4f}")

# ── WGS 指标 ──
h, rows = read_table(os.path.join(picard_dir, f'{sample}_wgs_metrics.txt'),
                     'GENOME_TERRITORY')
r = rows[0]
mean_depth = float(get(h, r, 'MEAN_COVERAGE', 0))
fold80 = float(get(h, r, 'FOLD_80_BASE_PENALTY', 99))
pct_20x = float(get(h, r, 'PCT_20X', 0)) * 100
pct_1x = float(get(h, r, 'PCT_1X', 0)) * 100
pct_exc_dupe = float(get(h, r, 'PCT_EXC_DUPE', 0)) * 100
print(f"WGS metrics:  MEAN_COVERAGE={mean_depth:.2f}x  "
      f"FOLD_80_BASE_PENALTY={fold80:.2f}  "
      f"PCT_20X={pct_20x:.1f}%  PCT_1X={pct_1x:.1f}%  "
      f"PCT_EXC_DUPE={pct_exc_dupe:.1f}%")

# ── QC gate（官方阈值） ──
print()
print("=== QC GATE (Post-Alignment) ===")
checks = [
    (f"Mean depth >=20x", mean_depth >= 20, f"{mean_depth:.1f}x"),
    (f"Fold-80 <3.0", fold80 < 3.0, f"{fold80:.2f}"),
    (f"% bases >=20x >80%", pct_20x > 80, f"{pct_20x:.1f}%"),
]
all_pass = True
for label, ok, val in checks:
    print(f"  {'PASS' if ok else 'FAIL'}: {label}  (实测 {val})")
    if not ok:
        all_pass = False
print()
print("QC GATE:", "PASS — proceed to Session 10 (MarkDuplicates + BQSR)" if all_pass
      else "FAIL — investigate before proceeding")
PYEOF

# ── 6) samtools flagstat（04 步自留版本，避免被混进 logs 触发样本名错位） ─
# 用法：把 flagstat 输出到 04 步目录，让 samtools 模块的样本名 = Picard 一致
echo | tee -a "$LOG"
echo "[$(date)] samtools flagstat..." | tee -a "$LOG"
samtools flagstat "$BAM" > "$OUT_ROOT/${SAMPLE}.flagstat.txt"
echo "[$(date)] flagstat done" | tee -a "$LOG"

# ── 7) MultiQC：聚合所有 QC 报告 ───────────────────────────────────────
# 只扫描 OUT_ROOT（不扫 LOG_DIR）——LOG_DIR 里的 03 步日志会被 samtools/flagstat
# 模块按内容匹配捕获，导致样本名错位（03.bwa 与 NA12878_chr20_sorted 不同行）
echo | tee -a "$LOG"
echo "[$(date)] MultiQC..." | tee -a "$LOG"
multiqc \
  "$OUT_ROOT" \
  -o "$MULTIQC" \
  --filename "${SAMPLE_SHORT}_PostAlignment_QC" \
  --title "${SAMPLE_LABEL} — Post-Alignment QC Report" \
  --comment "mosdepth + Picard + samtools flagstat | WGS Session 09 | Multiome Academy" \
  --force \
  2>>"$LOG"
echo "[$(date)] MultiQC done" | tee -a "$LOG"

echo | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
echo "  mosdepth:  $MOSD/${SAMPLE}*" | tee -a "$LOG"
echo "  Picard:    $PICARD/${SAMPLE}_{insert_size_metrics,alignment_summary,wgs_metrics}.txt" | tee -a "$LOG"
echo "  MultiQC:   $MULTIQC/NA12878_PostAlignment_QC.html" | tee -a "$LOG"