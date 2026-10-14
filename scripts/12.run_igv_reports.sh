#!/usr/bin/env bash
# Session 19（第 12 步）igv-reports：协作式 HTML 变异报告
#
# 对应官方 Session 19 的 "BASH — Install and run igv-reports"。
# 产出**自包含 HTML**（--standalone：所有数据内嵌，离线可分享）。
# 注意：igv-reports 1.17 用 `--output`（教程写的 `-o` 在 1.17 上不被识别）。
#
# 两份报告：
#   1) top50 SNV 报告     —— 按 QUAL 取前 50 个 PASS SNV（官方做法）
#   2) SV 断点报告        —— 用 4.13 的 SURVIVOR 共识 VCF，--flanking 2000
#                            （官方提示：看断点两侧的 discordant read pairs）
#
# 用法：
#   bash scripts/12.run_igv_reports.sh            # 两份都做
#   bash scripts/12.run_igv_reports.sh snv|sv

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

OUT="$ROOT/12.figures"
LOG="$WS/logs/12.run_igv_reports.log"
VCF="$ROOT/07.var.filtered/${SAMPLE}_filtered.vcf.gz"
SV_VCF="$ROOT/11.sv/merged/${SAMPLE}_SVs_merged.vcf.gz"
BAM="$ROOT/05.bqsr/${SAMPLE}_BQSR.bam"
REF="$WS/rawData/ref/Homo_sapiens_assembly38.fasta"

ONLY="${1:-all}"
mkdir -p "$OUT" "$WS/logs"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
want() { [ "$ONLY" = "all" ] || [ "$ONLY" = "$1" ]; }

echo "=============================================="
echo " Session 19: igv-reports HTML reports"
echo "=============================================="
echo " SNV VCF : $VCF"
echo " SV VCF  : $SV_VCF"
echo " BAM     : $BAM"
echo " Out     : $OUT"
echo " Run     : $ONLY"
echo

# ── 0) 工具检查/安装 ──
log "检查 igv-reports..."
if ! python3 -c "import igv_reports" 2>/dev/null; then
  log "  未安装 → pip 安装"
  python3 -m pip install --quiet igv-reports 2>&1 | tail -3 | tee -a "$LOG"
fi
python3 -c "import igv_reports; print('  igv-reports', getattr(igv_reports,'__version__','?'))" 2>&1 | tee -a "$LOG"
CREATEPORT=$(command -v create_report || echo "")
[ -n "$CREATEPORT" ] || CREATEPORT="python3 -m igv_reports"
log "  create_report: $CREATEPORT"
echo

# ── pre-flight ──
for f in "$BAM" "$REF"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" | tee -a "$LOG"; exit 1; }
done
[ -f "$BAM.bai" ] || { echo "ERROR: missing BAM index" | tee -a "$LOG"; exit 1; }
echo

# ══════════════════════════════════════════════════════
# 报告 1：top 50 PASS SNV
# ══════════════════════════════════════════════════════
if want snv; then
HTML="$OUT/${SAMPLE_SHORT}_variant_report.html"
TOP="$OUT/top50_snvs.vcf.gz"
if [ -f "$HTML" ]; then
  log "SKIP  SNV 报告已存在"
else
log "[1/2] 取 top 50 PASS SNV（按 QUAL）..."
HDR="$OUT/.top50.header"
bcftools view -h "$VCF" > "$HDR" 2>/dev/null
bcftools view -H -f PASS -v snps "$VCF" 2>/dev/null \
  | sort -t$'\t' -k6,6gr | head -50 > "$OUT/.top50.body"
cat "$HDR" "$OUT/.top50.body" | bgzip -c > "$TOP"
bcftools index -f "$TOP" 2>/dev/null
rm -f "$HDR" "$OUT/.top50.body"
log "  top50 SNV：$(bcftools view -H "$TOP" 2>/dev/null | wc -l) 个"

log "  生成 HTML（--standalone 自包含）..."
$CREATEPORT "$TOP" \
  --fasta "$REF" \
  --tracks "$BAM" \
  --info-columns "DP" \
  --sample-columns "GT,AD,DP" \
  --standalone \
  --output "$HTML" 2>&1 | tail -5 | tee -a "$LOG"
fi
echo
fi

# ══════════════════════════════════════════════════════
# 报告 2：SV 断点（--flanking 2000）
# ══════════════════════════════════════════════════════
if want sv; then
HTML_SV="$OUT/${SAMPLE_SHORT}_sv_breakpoint_report.html"
if [ -f "$HTML_SV" ]; then
  log "SKIP  SV 报告已存在"
elif [ ! -f "$SV_VCF" ]; then
  log "SKIP  SV VCF 不存在（需先跑 11.run_sv.sh merge）"
else
log "[2/2] SV 断点报告（--flanking 2000）..."
# 上限 50 个（igv-reports 超 200 个会变慢）
SV_TOP="$OUT/.sv_top.vcf.gz"
bcftools view -h "$SV_VCF" > "$OUT/.sv.hdr" 2>/dev/null
bcftools view -H "$SV_VCF" 2>/dev/null | head -50 > "$OUT/.sv.body"
cat "$OUT/.sv.hdr" "$OUT/.sv.body" | bgzip -c > "$SV_TOP"
bcftools index -f "$SV_TOP" 2>/dev/null
rm -f "$OUT/.sv.hdr" "$OUT/.sv.body"
log "  SV 取 $(bcftools view -H "$SV_TOP" 2>/dev/null | wc -l) 个断点"

$CREATEPORT "$SV_TOP" \
  --fasta "$REF" \
  --tracks "$BAM" "$SV_VCF" \
  --flanking 2000 \
  --standalone \
  --output "$HTML_SV" 2>&1 | tail -5 | tee -a "$LOG"
rm -f "$SV_TOP" "$SV_TOP.tbi"
fi
echo
fi

# ── 小结 ──
log "=== 产出 ==="
for f in "$OUT"/${SAMPLE_SHORT}_variant_report.html "$OUT"/${SAMPLE_SHORT}_sv_breakpoint_report.html "$OUT"/top50_snvs.vcf.gz; do
  [ -f "$f" ] && printf "  %-52s %9s B\n" "$(basename "$f")" "$(stat -c%s "$f")" | tee -a "$LOG"
done 2>/dev/null
echo
log "Done."