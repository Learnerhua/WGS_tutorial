#!/usr/bin/env bash
# Session 15（Step 1/4）：hap.py F1 基准测试
#
# 目的：拿 4.8 步过滤后的 NA12878 chr20 PASS 变异 vs GIAB 真值集，
#       算 precision / recall / F1，客观检验前 8 步流程的准确性。
#
# 教程预期成绩（NA12878 GATK Best Practices，chr20）：
#   SNV  Precision ≥ 99.5%  Recall ≥ 99.0%  F1 ≥ 99.2%
#   INDEL Precision ≥ 97.0%  Recall ≥ 95.0%  F1 ≥ 96.0%
#
# 用法：bash scripts/09.run_happy_benchmark.sh
# 预计耗时：3–8 分钟

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

LOG="$WS/logs/09.happy_benchmark.log"
mkdir -p "$ROOT/09.benchmark" "$WS/logs"

# 工具路径（WGS env 里的 bcftools；hap.py 在 py27 env）
export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

# ── 输入 ──
VCF_IN="$ROOT/07.var.filtered/${SAMPLE}_filtered.vcf.gz"   # 4.8 步过滤后
REF="$WS/rawData/ref/Homo_sapiens_assembly38.fasta"
TRUTH="$WS/rawData/databases/giab/HG001_truth.vcf.gz"
CONF_BED="$WS/rawData/databases/giab/HG001_highconf.bed"
PY27=/oldhome/ouyjh/miniforge3/envs/py27/bin
mkdir -p "$ROOT/09.benchmark"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 15 step 1/4: hap.py F1 benchmark"
echo "=============================================="
echo

# ── 0) Pre-flight ──
log "Pre-flight checks..."
for f in "$VCF_IN" "$REF" "$TRUTH" "$CONF_BED"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing $f" | tee -a "$LOG"
    exit 1
  fi
done
[ -x "$PY27/hap.py" ] || { echo "ERROR: hap.py not found at $PY27/hap.py" | tee -a "$LOG"; exit 1; }

# 单样本提取（hap.py 期望单样本 VCF）
# 注：教程是 trio VCF（NA12891/92/78），用 -s NA12878 抽子代；
# 本项目是 NA12878_chr20 单样本（4.3 步起就是 chr20 抽取），无需 -s 抽取
SAMPLE_NAME=$(bcftools view -h "$VCF_IN" 2>/dev/null | awk -F'\t' '/^#CHROM/{print $10; exit}')
log "检测到样本名：${SAMPLE_NAME:-<none>}"

SINGLE="$ROOT/09.benchmark/${SAMPLE_NAME}_pass_singletons.vcf.gz"
if [ ! -f "$SINGLE" ]; then
  if [ "$SAMPLE_NAME" = "$SAMPLE" ]; then      # 配置样本 == VCF 样本 → 已是单样本
    # 本项目：直接是单样本 VCF，不需要 -s
    log "已是单样本 VCF，直接复制 + 索引..."
    cp "$VCF_IN" "$SINGLE" && bcftools index -t "$SINGLE" 2>>"$LOG"
  else
    # 通用路径：trio 抽 NA12878
    log "trio VCF 抽 ${SAMPLE_SHORT} 子样本..."
    bcftools view -s "$SAMPLE_SHORT" "$VCF_IN" \
      | bcftools view -f PASS -O z -o "$SINGLE" 2>>"$LOG"
    bcftools index -t "$SINGLE" 2>>"$LOG"
  fi
else
  log "SKIP  单样本 VCF 已存在"
fi

# high-conf BED 覆盖范围（教程 bio-box 要求）
n_chr20=$(grep -c "^chr20" "$CONF_BED" 2>/dev/null || echo 0)
log "high-confidence intervals on chr20: $n_chr20"

# ── 1) 跑 hap.py（vcfeval 引擎，更精确）──
OUT="$ROOT/09.benchmark/hapy_${SAMPLE}"
if [ ! -f "$OUT.summary.csv" ]; then
  log "Running hap.py (vcfeval, 4 threads)..."
  log "  truth: $TRUTH"
  log "  query: $SINGLE"
  log "  bed:   $CONF_BED"
  log "  ref:   $REF"
  log "  region:chr20"
  "$PY27/hap.py" \
    "$TRUTH" "$SINGLE" \
    -f "$CONF_BED" \
    -r "$REF" \
    -o "$OUT" \
    --engine=xcmp \
    -l chr20 \
    --threads 4 >>"$LOG" 2>&1
  log "hap.py finished, exit=$?"
else
  log "SKIP  hap.py output 已存在"
fi

# ── 2) 解析 summary.csv ──
log "=== hap.py summary ==="
if [ ! -f "$OUT.summary.csv" ]; then
log "ERROR: $OUT.summary.csv 不存在"
exit 1
fi
cat "$OUT.summary.csv" | tee -a "$LOG"
echo | tee -a "$LOG"

# 直接读 hap.py 官方指标列（不自算）：
#   0=Type 1=Filter 10=METRIC.Recall 11=METRIC.Precision 13=METRIC.F1_Score
# 注：hap.py 的 Precision 对 FP 分类有内部处理，与 TP/(TP+FP) 略有差异（INDEL 上差 0.02pp），
#     故以 hap.py 官方列为准。
python3 - "$OUT.summary.csv" "$ROOT/09.benchmark/hapy_metrics.tsv" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, csv
src, dst = sys.argv[1], sys.argv[2]
with open(src) as f:
    rows = list(csv.reader(f))
with open(dst, 'w') as f:
    f.write('Type\tFilter\tRecall\tPrecision\tF1\n')
    for r in rows[1:]:
        if len(r) > 13 and r[0] in ('SNP', 'INDEL', 'MNP'):
            try:
                recall, prec, f1 = float(r[10]), float(r[11]), float(r[13])
                f.write(f'{r[0]}\t{r[1]}\t{recall*100:.2f}%\t{prec*100:.2f}%\t{f1*100:.2f}%\n')
                print(f'  {r[0]:5s} {r[1]:5s}  Recall={recall*100:6.2f}%  Precision={prec*100:6.2f}%  F1={f1*100:6.2f}%')
            except (ValueError, IndexError) as e:
                print(f'  {r[0]} 解析失败: {e}')
PYEOF
echo | tee -a "$LOG"

# ── 3) 产出一览 ──
echo
log "=== 产出文件 ==="
ls -lh "$ROOT/09.benchmark/" | tail -n +2 | tee -a "$LOG"
echo
log "QC GATE: F1 应 ≥ 99% SNV / ≥ 96% INDEL"

# ── 4) 教程预期对照 ──
cat <<EOF | tee -a "$LOG"

[expected from tutorial] ${SAMPLE_LABEL} GATK Best Practices:
  SNV  Precision ≥ 99.5%   Recall ≥ 99.0%   F1 ≥ 99.2%
  INDEL Precision ≥ 97.0%  Recall ≥ 95.0%   F1 ≥ 96.0%

[interpreting results]
  F1 ≥ 预期值：流程准确，无需调整
  precision 低：过多 false positive → 4.8 步阈值偏松
  recall 低：   漏检真变异 → 4.6 步 HC 调参 / 4.5 步 BQSR 重做
  SNV/INDEL 偏差大：先看 INDEL（通常较难）
EOF

log "Done. Next: bash scripts/09.run_database_crossref.sh"
