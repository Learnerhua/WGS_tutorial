#!/usr/bin/env bash
# Session 19（第 12 步）R/ggplot2 出版图：包装脚本
#
# 两步：
#   1) 从 VEP 注释 VCF 提取「变异 → 后果 / IMPACT / 基因」表（vep_summary.tsv）
#   2) 调 Rscript 画出版级汇总图（PDF + 300 dpi PNG）
#
# 官方用 awk 从 CSQ 列取第一段；本项目改用 Python 取**最严重**的那段
# （VEP 一条变异可有多个逗号分隔的 CSQ，取最严重更准确）。
#
# 用法：bash scripts/12.run_ggplot2_plots.sh

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
LOG="$WS/logs/12.run_ggplot2_plots.log"
VEP="$ROOT/08.annotation/${SAMPLE}_VEP.vcf.gz"
TSV="$OUT/vep_summary.tsv"
RSCRIPT="$WS/scripts/12.plot_consequence_summary.R"

mkdir -p "$OUT" "$WS/logs"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 19: R/ggplot2 publication figure"
echo "=============================================="
echo " VEP VCF : $VEP"
echo " Out     : $OUT"
echo

# ── Pre-flight ──
log "Pre-flight checks..."
[ -f "$VEP" ] || { echo "ERROR: missing $VEP" | tee -a "$LOG"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found" | tee -a "$LOG"; exit 1; }
Rscript -e 'suppressPackageStartupMessages({library(ggplot2);library(dplyr)})' 2>/dev/null \
  || { echo "ERROR: R 包 ggplot2 / dplyr 缺失" | tee -a "$LOG"; exit 1; }
log "R $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
log "ggplot2 $(Rscript -e 'cat(as.character(packageVersion("ggplot2")))' 2>/dev/null)"
echo

# ── 1) 提取 vep_summary.tsv ──
if [ -f "$TSV" ]; then
  log "[1/2] SKIP  vep_summary.tsv 已存在（$(wc -l < "$TSV") 行）"
else
log "[1/2] 从 VEP VCF 提取后果表（取最严重 CSQ）..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$VEP" 2>/dev/null \
  | python3 -c '
import sys
# CSQ 格式：Allele|Consequence|IMPACT|SYMBOL|...
# 一条变异可有多个逗号分隔的 CSQ 段 → 取 IMPACT 最严重的一段
RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
out = []
for line in sys.stdin:
    p = line.rstrip("\n").split("\t")
    if len(p) < 5 or not p[4] or p[4] == ".":
        continue
    chrom, pos, ref, alt, csq = p[0], p[1], p[2], p[3], p[4]
    best = None
    for seg in csq.split(","):
        f = seg.split("|")
        if len(f) < 4:
            continue
        cons, imp, sym = f[1], f[2], f[3]
        if not imp:
            continue
        if best is None or RANK.get(imp, 9) < RANK.get(best[1], 9):
            best = (cons, imp, sym)
    if best:
        # 列顺序与 R 脚本一致：chrom pos alt gene impact consequence
        out.append((chrom, pos, alt, best[2] or "-", best[1], best[0]))
w = sys.stdout.write
for r in out:
    w("\t".join(r) + "\n")
' > "$TSV"
log "  提取 $(wc -l < "$TSV") 条变异"
fi
echo

# ── 2) R 出图 ──
log "[2/2] Rscript 出图（PDF + 300 dpi PNG）..."
Rscript "$RSCRIPT" "$TSV" "$OUT" 2>&1 | tee -a "$LOG"
echo

# ── 小结 ──
log "=== 产出 ==="
for f in "$OUT/${SAMPLE}_consequence_summary.pdf" "$OUT/${SAMPLE}_consequence_summary.png" "$TSV"; do
  [ -f "$f" ] && printf "  %-58s %8s B\n" "$(basename "$f")" "$(stat -c%s "$f")" | tee -a "$LOG"
done
echo
log "Done."