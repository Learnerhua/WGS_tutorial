#!/usr/bin/env bash
# Session 11：HaplotypeCaller 变异调用（GVCF 模式）
# 流程第 6 步：对 BQSR BAM 做变异检测，输出 GVCF（每位置都记录，含 REF 块）
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

BAM=$ROOT/05.bqsr/${SAMPLE}_BQSR.bam
REF=$WS/rawData/ref/Homo_sapiens_assembly38.fasta
OUT_ROOT=$ROOT/06.gatk
LOG=$WS/logs/06.gatk.log
GVCF=$OUT_ROOT/${SAMPLE}.g.vcf.gz

mkdir -p "$OUT_ROOT" "$WS/logs"

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

echo "[$(date)] === 06. HaplotypeCaller (GVCF) ===" | tee "$LOG"
echo "Input:  $BAM" | tee -a "$LOG"
echo "Output: $GVCF" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 0) Pre-flight 检查（防御 Session 10 BAM header 不合法 / GATK 未装） ─────
echo "[$(date)] Pre-flight checks..." | tee -a "$LOG"
samtools quickcheck -v "$BAM" || { echo "ERROR: BAM quickcheck 失败"; exit 1; }
[ -f "${BAM}.bai" ] || { echo "ERROR: BAM 索引缺失"; exit 1; }
samtools view -H "$BAM" | grep -q "^@RG" || { echo "ERROR: 无 @RG 标签"; exit 1; }
gatk --version 2>&1 | grep -q "4\." || { echo "ERROR: GATK 4.x 未找到"; exit 1; }
echo "Pre-flight OK" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 1) HaplotypeCaller：核心变异检测（断点重入：GVCF 已存在则跳过） ────────
if [ ! -f "$GVCF" ]; then
  # -ERC GVCF                输出 GVCF 模式（每位点都记录，含 REF 块）
  #                          → Session 12 Joint Genotyping 用
  # --native-pair-hmm-threads 8
  #                          PairHMM 计算最耗时（CPU 重活）；8 线程 64 核机器
  #                          教程用 4，生产用 8（4→8 仅 1.5×加速但仍显著）
  # -G StandardAnnotation    标准注释（DP, AD, FS, SOR, MQ 等）
  # -G AS_StandardAnnotation 等位基因特异注释（Session 13 VQSR 必需；硬过滤不用）
  # -L chr20                 限制 chr20 子集
  # 注：教程用 --dbsnp 加 dbSNP138 注释；我们 dbSNP138 未下 → 注释略少但变异调用不受影响
  echo "[$(date)] HaplotypeCaller..." | tee -a "$LOG"
  gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "$GVCF" \
    -ERC GVCF \
    --native-pair-hmm-threads 8 \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -L chr20 \
    2>>"$LOG"
  echo "[$(date)] HaplotypeCaller done" | tee -a "$LOG"
  ls -lh "$GVCF"* | tee -a "$LOG"
  echo | tee -a "$LOG"
else
  echo "[$(date)] HaplotypeCaller skipped (GVCF exists, will re-validate only)" | tee -a "$LOG"
  echo | tee -a "$LOG"
fi

# ── 2) 验证 GVCF：结构 + 记录数（用 bcftools stats，权威且快） ─────
echo "[$(date)] Validation..." | tee -a "$LOG"
echo "Header（关键 INFO 字段）:" | tee -a "$LOG"
bcftools view -h "$GVCF" 2>/dev/null | grep "^##source\|^##GATKCommand\|^##reference\|^##fileformat\|^##contig=<ID=chr20" | head -5 | tee -a "$LOG"
echo | tee -a "$LOG"

echo "=== 记录统计 ===" | tee -a "$LOG"
bcftools stats "$GVCF" 2>/dev/null | awk -F'\t' '
/^SN/ {
  # 字段 3 是 "number of ...:"；字段 4 是数字
  if ($3 ~ /number of SNPs/) snv = $4
  if ($3 ~ /number of indels/) indel = $4
  if ($3 ~ /number of MNPs/) mnp = $4
  if ($3 ~ /number of records/) total = $4
  if ($3 ~ /number of no-ALTs/) no_alt = $4
  if ($3 ~ /number of multiallelic sites/) multi = $4
}
END {
  gsub("[ \t]", "", snv); gsub("[ \t]", "", indel); gsub("[ \t]", "", mnp)
  gsub("[ \t]", "", total); gsub("[ \t]", "", no_alt); gsub("[ \t]", "", multi)
  printf "  Total records:         %s\n", total
  printf "  REF confidence blocks: %s  (no-ALTs, chr20 完整 REF 覆盖)\n", no_alt
  printf "  SNV (SNPs):            %s\n", snv
  printf "  Indel:                 %s\n", indel
  printf "  MNP:                   %s\n", mnp
  printf "  Multiallelic sites:    %s\n", multi
  printf "  SNV+INDEL total:       %d\n", snv+indel
  print ""
  # NA12878 chr20 30× 预期：SNV+INDEL 120k-160k
  exp_min = 120000; exp_max = 160000
  ok = (snv+indel >= exp_min && snv+indel <= exp_max)
  printf "  SNV+INDEL 预期 120k-160k, 实测 %d  → %s\n", snv+indel, (ok ? "OK" : "CHECK")
}' | tee -a "$LOG"

echo | tee -a "$LOG"
echo "Output: $GVCF" | tee -a "$LOG"
echo "Tutorial next step: Session 12 (Joint Genotyping, **skipped - missing trio data**)" | tee -a "$LOG"