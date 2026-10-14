#!/usr/bin/env bash
# Session 10：MarkDuplicates + BaseRecalibrator + ApplyBQSR
# 流程第 5 步：去除 PCR/光学重复 → 重校准碱基质量分数 → 输出分析就绪的 BAM
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

BAM_IN=$ROOT/03.aligned/${SAMPLE}_sorted.bam
REF=$WS/rawData/ref/Homo_sapiens_assembly38.fasta
DICT=$WS/rawData/ref/Homo_sapiens_assembly38.dict
RES=$WS/rawData/resources
OUT_ROOT=$ROOT/05.bqsr
LOG=$WS/logs/05.bqsr.log
BAM_MD=$OUT_ROOT/${SAMPLE}_markdup_sorted.bam
METRICS=$OUT_ROOT/${SAMPLE}_duplication_metrics.txt
RECAL=$OUT_ROOT/${SAMPLE}_recal_data.table
BAM_BQSR=$OUT_ROOT/${SAMPLE}_BQSR.bam

mkdir -p "$OUT_ROOT" "$WS/logs"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

echo "[$(date)] === 05. MarkDuplicates + BQSR ===" | tee "$LOG"
echo "Input: $BAM_IN" | tee -a "$LOG"
echo "Output: $OUT_ROOT/" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 0) 准备 interval_list（GATK 4.6 的 -L 简写不成立，必须用 interval_list） ─
# 详见坑 10：必须含完整 3366 条 @SQ（与 BAM header 一致）
ILIST=$OUT_ROOT/chr20.interval_list
{
  grep '^@HD' "$DICT" 2>/dev/null || printf '@HD\tVN:1.6\tSO:coordinate\n'
  grep '^@SQ' "$DICT" 2>/dev/null
  awk -F'\t' '$1=="chr20"{printf "%s\t1\t%s\t+\t%s\n", $1, $2, $1}' "$REF.fai"
} > "$ILIST"
echo "interval_list: $ILIST ($(wc -l < "$ILIST") 行)" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 1) BAM 完整性 + header 检查（防御 BWA @PG 旧问题） ──────────────────
echo "[$(date)] BAM quickcheck..." | tee -a "$LOG"
samtools quickcheck -v "$BAM_IN" || {
  echo "ERROR: BAM 完整性检查失败" | tee -a "$LOG"
  exit 1
}
echo | tee -a "$LOG"

# ── 2) MarkDuplicates：标记 PCR/光学重复（FLAG bit 1024） ─────────────────
# 断点重入：如果 markdup BAM 已生成且 quickcheck 通过，跳过这步
if [ ! -f "$BAM_MD" ] || ! samtools quickcheck "$BAM_MD" 2>/dev/null; then
  echo "[$(date)] MarkDuplicates..." | tee -a "$LOG"
  gatk MarkDuplicates \
    -I "$BAM_IN" \
    -O "$BAM_MD" \
    -M "$METRICS" \
    --CREATE_INDEX true \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --TMP_DIR /tmp \
    --VALIDATION_STRINGENCY SILENT \
    2>>"$LOG"
  echo "[$(date)] MarkDuplicates done" | tee -a "$LOG"
else
  echo "[$(date)] MarkDuplicates skipped (BAM exists and OK)" | tee -a "$LOG"
fi

# 提取关键统计
echo | tee -a "$LOG"
echo "=== MarkDuplicates 关键指标 ===" | tee -a "$LOG"
samtools flagstat "$BAM_MD" | grep -E "in total|duplicate" | tee -a "$LOG"
python3 - "$METRICS" <<'PYEOF' | tee -a "$LOG"
import sys
hdr, row = None, None
for line in open(sys.argv[1]):
    if line.startswith('LIBRARY\tUNPAIRED_READS_EXAMINED'):
        hdr = line.rstrip('\n').split('\t')
    elif hdr and line.strip() and not line.startswith('#'):
        row = line.rstrip('\n').split('\t')
        break

def g(name):
    try: return row[hdr.index(name)]
    except ValueError: return 'NA'

print(f"  PERCENT_DUPLICATION:   {g('PERCENT_DUPLICATION')}")
print(f"  ESTIMATED_LIBRARY_SIZE: {g('ESTIMATED_LIBRARY_SIZE')}")
print(f"  READ_PAIR_DUPLICATES:   {g('READ_PAIR_DUPLICATES')}")
print(f"  READ_PAIR_OPTICAL_DUPLICATES: {g('READ_PAIR_OPTICAL_DUPLICATES')}")
PYEOF
echo | tee -a "$LOG"

# ── 3) BaseRecalibrator：建 BQSR 重校准模型（不修改 BAM） ──────────────────
# --known-sites × 3：Mills + 1000G_phase1 + dbSNP138（教程用 dbsnp138；我们没下，
#                   用 Mills + 1000G_phase1 替代——见 2.6 节说明）
# -O 输出 recal_data.table（重校准模型文件，不修改 BAM）
echo "[$(date)] BaseRecalibrator (build model)..." | tee -a "$LOG"
gatk BaseRecalibrator \
  -I "$BAM_MD" \
  -R "$REF" \
  --known-sites "$RES/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  --known-sites "$RES/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
  --known-sites "$RES/hapmap_3.3.hg38.vcf.gz" \
  -O "$RECAL" \
  -L chr20 \
  2>>"$LOG"
echo "[$(date)] BaseRecalibrator done" | tee -a "$LOG"
ls -lh "$RECAL" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 4) ApplyBQSR：把重校准模型应用到 BAM（生成分析就绪 BAM） ─────────────
# --bqsr-recal-file       上一步生成的模型
# --emit-original-quals true
#                         把原始 Q 写到 OQ tag（调试 + 对比前/后质量分布）
# --static-quantized-quals N N N
#                         把 Q 压成 4 个值（<10, ~10, ~20, ~30）
#                         BAM 体积减小 ~15%，对 GATK 工具精度影响可忽略
echo "[$(date)] ApplyBQSR (apply to BAM)..." | tee -a "$LOG"
gatk ApplyBQSR \
  -R "$REF" \
  -I "$BAM_MD" \
  --bqsr-recal-file "$RECAL" \
  -O "$BAM_BQSR" \
  --emit-original-quals true \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 \
  -L chr20 \
  2>>"$LOG"
echo "[$(date)] ApplyBQSR done" | tee -a "$LOG"

echo | tee -a "$LOG"
samtools index "$BAM_BQSR"
ls -lh "$BAM_BQSR"* | tee -a "$LOG"

# ── 5) 最终 BAM 完整性 + header 检查 ──────────────────────────────────────
echo | tee -a "$LOG"
echo "[$(date)] Final BAM verification..." | tee -a "$LOG"
samtools quickcheck -v "$BAM_BQSR" && echo "quickcheck OK" | tee -a "$LOG"
echo "@RG（必需）:" | tee -a "$LOG"
samtools view -H --no-PG "$BAM_BQSR" | grep '^@RG' | tee -a "$LOG"
echo "@HD sort order:" | tee -a "$LOG"
samtools view -H --no-PG "$BAM_BQSR" | grep '^@HD' | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 6) QC GATE ──────────────────────────────────────────────────────────────
echo "=== QC GATE (Post-MarkDuplicates + BQSR) ===" | tee -a "$LOG"
python3 - "$METRICS" <<'PYEOF' | tee -a "$LOG"
import sys
hdr, row = None, None
for line in open(sys.argv[1]):
    if line.startswith('LIBRARY\tUNPAIRED_READS_EXAMINED'):
        hdr = line.rstrip('\n').split('\t')
    elif hdr and line.strip() and not line.startswith('#'):
        row = line.rstrip('\n').split('\t')
        break

def g(name):
    try: return row[hdr.index(name)]
    except ValueError: return None

pct_dup = float(g('PERCENT_DUPLICATION')) * 100
est_lib = float(g('ESTIMATED_LIBRARY_SIZE'))
optical = float(g('READ_PAIR_OPTICAL_DUPLICATES'))
read_pairs = float(g('READ_PAIRS_EXAMINED'))

print(f"Pct duplicated:      {pct_dup:.2f}%        (PCR-free 文库预期 <10%)")
print(f"Estimated library:   {est_lib:,.0f} molecules")
print(f"Optical duplicates:  {optical:,.0f} pairs  ({optical/read_pairs*100:.2f}% of examined)")

# PCR-free 文库预期 <10%（教程预期）；FFPE 可达 40–60%
ok = pct_dup < 10
print()
print("QC GATE:", "PASS — proceed to Session 11 (HaplotypeCaller)" if ok
      else "FAIL — investigate library prep")
PYEOF

echo | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
echo "  MarkDup BAM:  $BAM_MD" | tee -a "$LOG"
echo "  BQSR BAM:     $BAM_BQSR" | tee -a "$LOG"
echo "  Metrics:      $METRICS" | tee -a "$LOG"
echo "  Recal table:  $RECAL" | tee -a "$LOG"
