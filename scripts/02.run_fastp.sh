#!/usr/bin/env bash
# Session 06：fastp 接头修剪与预处理 + 修剪前后 MultiQC 对比
# 流程第 2 步：对 FastQC 质控通过的数据做接头去除 + 质量修剪
# 官方教程 fastp profile（30× PCR-free WGS, insert 350-450bp）
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

FASTQ_DIR=$WS/rawData/fastq
OUT_ROOT=$ROOT/02.trimmed
TRIM_DIR=$OUT_ROOT
FASP_QC=$OUT_ROOT/fastp
FASTQC_PRE=$ROOT/01.qc/fastqc          # 01 步的修剪前 FastQC（只读引用）
FASTQC_POST=$OUT_ROOT/fastqc_post    # 修剪后 FastQC
MULTIQC_OUT=$OUT_ROOT/multiqc        # pre/post 对比报告
LOG_DIR=$WS/logs
mkdir -p "$OUT_ROOT" "$OUT_ROOT/config" "$FASP_QC" "$FASTQC_POST" "$MULTIQC_OUT" "$LOG_DIR"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

echo "[$(date)] === 02. fastp trimming + pre/post MultiQC ===" | tee "$LOG_DIR/02.fastp.log"
echo "输入: $FASTQ_DIR/${SAMPLE}_R{1,2}.fastq.gz" | tee -a "$LOG_DIR/02.fastp.log"
echo "输出: $OUT_ROOT/" | tee -a "$LOG_DIR/02.fastp.log"
echo | tee -a "$LOG_DIR/02.fastp.log"

# ── 1) fastp 核心参数（官方 WGS profile，逐 flag 见教程） ─────────────
# --detect_adapter_for_pe: PE 数据 overlap 检测接头（最准）
# --cut_tail: 3' 端滑窗质量修剪（窗口 4 bp，均值 <Q15 才剪——保守，保留 Q15-30 给 BQSR）
# --length_required 36: BWA-MEM2 最小有效比对长度
# --qualified_quality_phred 15 + --unqualified_percent_limit 40: 过滤 >40% 碱基 <Q15 的 read
# --n_base_limit 5: 过滤 >5 个 N 的 read
# --correction: overlap 区 R1/R2 不一致时取高质量碱基
echo "[$(date)] fastp trimming 开始..." | tee -a "$LOG_DIR/02.fastp.log"
fastp \
  --in1 "$FASTQ_DIR/${SAMPLE}_R1.fastq.gz" \
  --in2 "$FASTQ_DIR/${SAMPLE}_R2.fastq.gz" \
  --out1 "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
  --out2 "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
  --json "$FASP_QC/${SAMPLE}_fastp.json" \
  --html "$FASP_QC/${SAMPLE}_fastp.html" \
  --thread 48 \
  --detect_adapter_for_pe \
  --cut_tail \
  --cut_tail_window_size 4 \
  --cut_tail_mean_quality 15 \
  --length_required 36 \
  --qualified_quality_phred 15 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --correction \
  2>>"$LOG_DIR/02.fastp.log"

echo "[$(date)] fastp 完成" | tee -a "$LOG_DIR/02.fastp.log"
echo | tee -a "$LOG_DIR/02.fastp.log"

# ── 2) 验证 R1/R2 read 数一致 ─────────────────────────────────────────
R1=$(seqkit stats -T "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" | awk 'NR==2{print $4}')
R2=$(seqkit stats -T "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" | awk 'NR==2{print $4}')
if [ "$R1" = "$R2" ]; then
  echo "✓ Paired reads match: $R1 pairs" | tee -a "$LOG_DIR/02.fastp.log"
else
  echo "✗ MISMATCH: $R1 vs $R2 — re-run fastp" | tee -a "$LOG_DIR/02.fastp.log"
fi
echo | tee -a "$LOG_DIR/02.fastp.log"

# ── 3) 修剪后 FastQC（对比用） ────────────────────────────────────────
echo "[$(date)] FastQC (post-trim) 开始..." | tee -a "$LOG_DIR/02.fastp.log"
fastqc \
  "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
  "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
  -t 48 \
  --extract \
  -o "$FASTQC_POST" \
  2>>"$LOG_DIR/02.fastp.log"
echo "[$(date)] FastQC (post-trim) 完成" | tee -a "$LOG_DIR/02.fastp.log"
echo | tee -a "$LOG_DIR/02.fastp.log"

# ── 4) MultiQC：合并 修剪前 FastQC + 修剪后 FastQC + fastp 报告 ──────
# --config: 保留 _trimmed 后缀（否则 pre/post 样本名冲突，post 被跳过）
echo "[$(date)] MultiQC pre/post 对比 开始..." | tee -a "$LOG_DIR/02.fastp.log"

# 生成 MultiQC 配置（保留 _trimmed 后缀；沙箱模式下 OUT_ROOT/config 为空，须自建）
MQ_CONFIG="$OUT_ROOT/config/multiqc_config.yaml"
if [ ! -f "$MQ_CONFIG" ]; then
  cat > "$MQ_CONFIG" <<'MQEOF'
# MultiQC 配置：保留样本名中的 _trimmed 等后缀
# 原因：MultiQC 默认清洗样本名会剥离 _trimmed，导致修剪前(NA12878_chr20_R1)
# 与修剪后(NA12878_chr20_R1_trimmed)样本名冲突，后者被当重复跳过，
# pre/post 对比报告中丢失修剪后 FastQC 数据。
fn_clean_trim: []
fn_clean_exts: []
MQEOF
  echo "[$(date)] 已生成 $MQ_CONFIG" | tee -a "$LOG_DIR/02.fastp.log"
fi

multiqc \
  "$FASTQC_PRE" \
  "$FASTQC_POST" \
  "$FASP_QC" \
  -o "$MULTIQC_OUT" \
  --config "$OUT_ROOT/config/multiqc_config.yaml" \
  --filename "${SAMPLE}_pre_post_trimming" \
  --title "Pre vs Post Trimming — ${SAMPLE_LABEL}" \
  --comment "fastp 0.24.0 | WGS Session 06 | Multiome Academy" \
  --force \
  2>>"$LOG_DIR/02.fastp.log"
echo "[$(date)] MultiQC 完成" | tee -a "$LOG_DIR/02.fastp.log"
echo "对比报告: $MULTIQC_OUT/${SAMPLE}_pre_post_trimming.html" | tee -a "$LOG_DIR/02.fastp.log"
echo | tee -a "$LOG_DIR/02.fastp.log"

# ── 5) fastp JSON 关键 QC 指标自动检查（官方阈值） ───────────────────
echo "=== QC gate ===" | tee -a "$LOG_DIR/02.fastp.log"
if [ -f "$FASP_QC/${SAMPLE}_fastp.json" ]; then
  python3 - "$FASP_QC/${SAMPLE}_fastp.json" <<'PYEOF' | tee -a /Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial/logs/02.fastp.log
import json, sys
with open(sys.argv[1]) as f:
    d = json.load(f)
s = d["summary"]
before, after = s["before_filtering"], s["after_filtering"]
pass_rate = after["total_reads"] / before["total_reads"] * 100
q30_after = after["q30_rate"] * 100
adapter = d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
ins = d.get("insert_size", {}).get("peak", "N/A")

print(f"Input reads:      {before['total_reads']:,}")
print(f"Output reads:     {after['total_reads']:,}")
print(f"Reads passing:    {pass_rate:.2f}%  (需 >97%)")
print(f"Q30 rate after:   {q30_after:.1f}%  (需 >=80%)")
print(f"Adapter trimmed:  {adapter:,}")
print(f"Insert size peak: {ins} bp")
ok = pass_rate > 97 and q30_after >= 80
print("QC GATE:", "PASS — proceed to alignment" if ok else "FAIL — investigate")
PYEOF
else
  echo "ERROR: fastp JSON 未生成" | tee -a "$LOG_DIR/02.fastp.log"
fi

echo | tee -a "$LOG_DIR/02.fastp.log"
echo "输出文件:" | tee -a "$LOG_DIR/02.fastp.log"
echo "  Trimmed FASTQ: $OUT_ROOT/${SAMPLE}_R{1,2}_trimmed.fastq.gz" | tee -a "$LOG_DIR/02.fastp.log"
echo "  fastp 报告: $FASP_QC/${SAMPLE}_fastp.html" | tee -a "$LOG_DIR/02.fastp.log"
echo "  对比报告: $MULTIQC_OUT/NA12878_pre_post_trimming.html" | tee -a "$LOG_DIR/02.fastp.log"