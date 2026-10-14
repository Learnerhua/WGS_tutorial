#!/usr/bin/env bash
# Session 08：BWA-MEM2 比对 + samtools sort + 索引 + 验证
# 流程第 3 步：将修剪后 FASTQ 比对到参考基因组，输出 sorted BAM + BAI
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

REF=$WS/rawData/ref/Homo_sapiens_assembly38.fasta
TRIM_DIR=$ROOT/02.trimmed
OUT_DIR=$ROOT/03.aligned
TMP_DIR=/tmp/wgs_sort  # 本地 SSD 临时目录（samtools sort 临时文件）
LOG_DIR=$WS/logs
R1=$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz
R2=$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz
BAM=$OUT_DIR/${SAMPLE}_sorted.bam
BAI=$OUT_DIR/${SAMPLE}_sorted.bam.bai

mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

echo "[$(date)] === 03. BWA-MEM2 alignment + sort ===" | tee "$LOG_DIR/03.bwa.log"
echo "R1: $R1" | tee -a "$LOG_DIR/03.bwa.log"
echo "R2: $R2" | tee -a "$LOG_DIR/03.bwa.log"
echo "REF: $REF" | tee -a "$LOG_DIR/03.bwa.log"
echo "OUT: $BAM" | tee -a "$LOG_DIR/03.bwa.log"
echo | tee -a "$LOG_DIR/03.bwa.log"

# Read Group（官方要求 4 个必备字段 ID/SM/PL/LB + PU，GATK 强依赖）
# -t 必须从 FASTQ header 读（这里硬编码教程示例值；实际样本需替换）
# \t 是字面 tab，printf 生成
RG_FLOWCELL="${RG_FLOWCELL:-HMVFKDRXX.2}"        # 测序 flowcell.lane（RG ID 后缀）
RG_PU="${RG_PU:-${RG_FLOWCELL}.ACGTACGT}"          # platform unit
RG=$(printf '@RG\tID:%s.%s\tSM:%s\tPL:ILLUMINA\tLB:LIB-%s-PCRfree\tPU:%s' \
  "$SAMPLE" "$RG_FLOWCELL" "$SAMPLE" "$SAMPLE_SHORT" "$RG_PU")

echo "[$(date)] BWA-MEM2 + samtools sort pipeline starting..." | tee -a "$LOG_DIR/03.bwa.log"
echo "(32 threads BWA + 16 threads sort; bwa-mem2 outputs SAM to stdout, piped into samtools sort)" | tee -a "$LOG_DIR/03.bwa.log"
echo | tee -a "$LOG_DIR/03.bwa.log"

# 流式管道：BWA-MEM2 输出 SAM → stdin 进 samtools sort
# -K 10000000 (-k): 写 tmp 文件前需累积的 read 数（官方默认，本机大 RAM 用不上；此处不设）
bwa-mem2 mem \
  -t 32 \
  -R "$RG" \
  "$REF" \
  "$R1" "$R2" \
  2>>"$LOG_DIR/03.bwa.log" | \
samtools sort \
  -@ 16 \
  -m 4G \
  -T "$TMP_DIR/${SAMPLE}_sort" \
  -o "$BAM" \
  -

echo "[$(date)] Alignment + sort complete" | tee -a "$LOG_DIR/03.bwa.log"

# 建索引（必须 .bai 随机访问 BAM）
echo "[$(date)] samtools index..." | tee -a "$LOG_DIR/03.bwa.log"
samtools index "$BAM"
echo "[$(date)] Index complete" | tee -a "$LOG_DIR/03.bwa.log"

# ── 3) 清理 BAM header：修复 BWA-MEM2 @PG 行未转义 tab ──────────────────
# BWA-MEM2 把 -R 的 RG 字符串原样写进 @PG 的 CL 字段（tab 未转义）→
# htsjdk（GATK/Picard）严格解析会报错：
#   SAMFormatException: Problem parsing @PG key:value pair ID:xxx clashes with ID:bwa-mem2
# samtools 只警告，但 GATK 步骤会全部失败。用脚本把 CL 内的 tab 换成空格后 reheader 写回。
echo "[$(date)] Sanitizing @PG header..." | tee -a "$LOG_DIR/03.bwa.log"
python3 "$WS/scripts/fix_bam_header.py" "$BAM" "$BAM.sanitized"
mv "$BAM.sanitized" "$BAM"
samtools index "$BAM"          # header 变了，重建索引
echo "[$(date)] Header sanitized" | tee -a "$LOG_DIR/03.bwa.log"

ls -lh "$BAM" "$BAI" | tee -a "$LOG_DIR/03.bwa.log"
echo | tee -a "$LOG_DIR/03.bwa.log"

# ── 4) BAM 验证 + QC（官方 8 项检查）───────────────────────────────────
echo "=== Alignment Statistics ===" | tee -a "$LOG_DIR/03.bwa.log"
samtools flagstat "$BAM" | tee -a "$LOG_DIR/03.bwa.log"

echo | tee -a "$LOG_DIR/03.bwa.log"
echo "=== Read Groups (@RG, required by GATK) ===" | tee -a "$LOG_DIR/03.bwa.log"
samtools view -H "$BAM" | grep '^@RG' | tee -a "$LOG_DIR/03.bwa.log"

echo | tee -a "$LOG_DIR/03.bwa.log"
echo "=== Sort order (must be coordinate) ===" | tee -a "$LOG_DIR/03.bwa.log"
samtools view -H "$BAM" | grep '^@HD' | grep -o 'SO:[a-z]*' | tee -a "$LOG_DIR/03.bwa.log"

echo | tee -a "$LOG_DIR/03.bwa.log"
echo "=== chr20 mean depth (quick estimate) ===" | tee -a "$LOG_DIR/03.bwa.log"
samtools depth -a "$BAM" -r chr20 | \
  awk 'BEGIN{s=0;n=0}{s+=$3;n++}END{printf "chr20 mean depth: %.1fx (over %d positions)\n", s/n, n}' | \
  tee -a "$LOG_DIR/03.bwa.log"

echo | tee -a "$LOG_DIR/03.bwa.log"
echo "=== MAPQ distribution (top 5) ===" | tee -a "$LOG_DIR/03.bwa.log"
samtools view "$BAM" | awk '{print $5}' | sort | uniq -c | sort -rn | head -5 | tee -a "$LOG_DIR/03.bwa.log"

echo | tee -a "$LOG_DIR/03.bwa.log"
echo "=== QC GATE (expected thresholds) ===" | tee -a "$LOG_DIR/03.bwa.log"
# 从 flagstat 括号内提取百分比（如 "(99.97% : N/A)" → 99.97）
MAPPED_PCT=$(samtools flagstat "$BAM" | awk '/^[0-9]+ \+ [0-9]+ mapped \(/{gsub(/[()%]/,"",$5); print $5; exit}')
PAIRED_PCT=$(samtools flagstat "$BAM" | awk '/properly paired \(/{gsub(/[()%]/,"",$6); print $6; exit}')
SINGLE_PCT=$(samtools flagstat "$BAM" | awk '/singletons \(/{gsub(/[()%]/,"",$5); print $5; exit}')
echo "Mapped: ${MAPPED_PCT}% (need >=98.5%)" | tee -a "$LOG_DIR/03.bwa.log"
echo "Properly paired: ${PAIRED_PCT}% (need >=97%)" | tee -a "$LOG_DIR/03.bwa.log"
echo "Singletons: ${SINGLE_PCT}% (need <1%)" | tee -a "$LOG_DIR/03.bwa.log"

# 自动判定（浮点比较）
GATE=$(awk -v m="$MAPPED_PCT" -v p="$PAIRED_PCT" -v s="$SINGLE_PCT" \
  'BEGIN{print (m>=98.5 && p>=97 && s<1) ? "PASS" : "FAIL"}')

echo | tee -a "$LOG_DIR/03.bwa.log"
if [ "$GATE" = "PASS" ]; then
  echo "QC GATE: PASS — proceed to MarkDuplicates (Session 10)" | tee -a "$LOG_DIR/03.bwa.log"
else
  echo "QC GATE: FAIL — investigate" | tee -a "$LOG_DIR/03.bwa.log"
fi

echo | tee -a "$LOG_DIR/03.bwa.log"
echo "Output files:" | tee -a "$LOG_DIR/03.bwa.log"
echo "  Sorted BAM: $BAM" | tee -a "$LOG_DIR/03.bwa.log"
echo "  BAI index:  $BAI" | tee -a "$LOG_DIR/03.bwa.log"