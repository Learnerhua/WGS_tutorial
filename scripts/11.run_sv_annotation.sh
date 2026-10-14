#!/usr/bin/env bash
# Session 18（第 11 步）后处理：SV 基因注释 + 临床 review 表
#
# 输入：
#   11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz        SURVIVOR 2-caller 共识
#   11.sv/merged/NA12878_chr20_SVs_svdb_annotated.tsv    gnomAD-SV v4.1 已知/novel 标签
#   rawData/resources/Homo_sapiens_GRCh38_genes.bed      基因 BED（chr20）
#
# 输出：
#   NA12878_sv_calls.bed           共识 SV → BED
#   NA12878_sv_gene_annotated.tsv  SV × 基因
#   NA12878_sv_clinical_review.tsv 最终表（含 gnomAD 标签 + 优先级）
#
# 用法：bash scripts/11.run_sv_annotation.sh

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:/oldhome/ouyjh/miniforge3/envs/py27/bin:$PATH

OUT="$ROOT/11.sv"
LOG="$WS/logs/11.run_sv_annotation.log"
MERGED_VCF="$OUT/merged/${SAMPLE}_SVs_merged.vcf.gz"
ANNO_TSV="$OUT/merged/${SAMPLE}_SVs_svdb_annotated.tsv"
GENE_BED="$WS/rawData/resources/Homo_sapiens_GRCh38_genes.bed"

CALLS_BED="$OUT/merged/${SAMPLE_SHORT}_sv_calls.bed"
GENE_TSV="$OUT/merged/${SAMPLE_SHORT}_sv_gene_annotated.tsv"
FINAL_TSV="$OUT/merged/${SAMPLE_SHORT}_sv_clinical_review.tsv"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 18 post: SV gene annotation + clinical review"
echo "=============================================="
echo

# ── Pre-flight ──
log "Pre-flight checks..."
for f in "$MERGED_VCF" "$GENE_BED"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" | tee -a "$LOG"; exit 1; }
done
[ -f "$ANNO_TSV" ] || log "  WARN: $ANNO_TSV 不存在（先跑 11.run_sv.sh annotate）——gnomAD 标签将留空"
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found" | tee -a "$LOG"; exit 1; }
log "bedtools: $(bedtools --version 2>&1 | head -1)"
log "基因 BED 染色体命名：$(awk 'NR==1{print $1}' "$GENE_BED")  （$(wc -l < "$GENE_BED") 个基因）"
echo

# ── 1) 共识 SV → BED ──
log "[1/3] 共识 SV → BED..."
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$MERGED_VCF" 2>/dev/null \
  | awk 'BEGIN{OFS="\t"} $3!="" && $3+0>0 {print $1, $2-1, $3, $4, $5}' > "$CALLS_BED"
log "  共识 SV：$(wc -l < "$CALLS_BED") 个"
echo

# ── 2) bedtools 基因注释 ──
log "[2/3] bedtools 基因注释..."
bedtools intersect -a "$CALLS_BED" -b "$GENE_BED" -wa -wb 2>>"$LOG" \
  | awk 'BEGIN{OFS="\t"} {print $1, $2+1, $3, $4, $5, $9}' > "$GENE_TSV"
n_rows=$(wc -l < "$GENE_TSV")
n_sv=$(awk -F'\t' '{print $1":"$2":"$3}' "$GENE_TSV" | sort -u | wc -l)
log "  SV × 基因：$n_rows 行；有基因重叠的 SV：$n_sv 个"
echo

# ── 3) 合并 gnomAD 标签 → 临床 review 表 ──
log "[3/3] 合并 gnomAD-SV 标签 → 临床 review 表..."
ANNO_TSV="$ANNO_TSV" GENE_TSV="$GENE_TSV" CALLS_BED="$CALLS_BED" FINAL_TSV="$FINAL_TSV" \
python3 - <<'PYEOF' 2>&1 | tee -a "$LOG"
import os

anno_p = os.environ['ANNO_TSV']
gene_p = os.environ['GENE_TSV']
calls_p = os.environ['CALLS_BED']
out = os.environ['FINAL_TSV']

# gnomAD 标签： key = chrom:start:end  （start 为 0-based BED 坐标）
gnomad = {}
if os.path.exists(anno_p):
    with open(anno_p) as f:
        next(f, None)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 6:
                gnomad[(p[0], p[1], p[2])] = (p[5], p[6] if len(p) > 6 else '')

# 基因： key = chrom:start:end（0-based）→ 基因名列表
genes = {}
with open(gene_p) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 6:
            key = (p[0], str(int(p[1]) - 1), p[2])   # 还原为 BED 0-based
            genes.setdefault(key, []).append(p[5])

rows = []
n_known = n_novel = 0
with open(calls_p) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 4:
            continue
        chrom, start, end, svtype = p[0], p[1], p[2], p[3]
        svlen = p[4] if len(p) > 4 else '.'
        ing, gid = gnomad.get((chrom, start, end), ('NO', ''))
        if ing == 'YES':
            n_known += 1
            pri = 'LOW_known_population_SV'
        else:
            n_novel += 1
            pri = 'HIGH_novel_needs_clinical_review'
        gl = ';'.join(sorted(set(genes.get((chrom, start, end), [])))) or '-'
        rows.append((chrom, start, end, svtype, svlen, ing, gl, pri))

with open(out, 'w') as f:
    f.write('chrom\tstart\tend\tsvtype\tsvlen\tin_gnomad_sv\tgenes\tpriority\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')

print(f'  合计 {len(rows)} 个共识 SV')
print(f'    gnomAD-SV 已知（LOW）：{n_known}')
print(f'    novel（HIGH）        ：{n_novel}')
print()
print('  === 有基因重叠的 SV（按优先级）===')
shown = [r for r in rows if r[6] != '-']
for r in sorted(shown, key=lambda x: x[7]):
    print(f'    {r[0]}:{r[1]}-{r[2]}  {r[3]:4s}  {r[4]:>8s} bp  [{r[7][:4]}]  {r[6][:70]}')
if not shown:
    print('    （无）')
PYEOF
echo

# ── 小结 ──
log "=== 产出 ==="
for f in "$CALLS_BED" "$GENE_TSV" "$FINAL_TSV"; do
  [ -f "$f" ] && printf "  %-58s %8s B\n" "${f#$WS/}" "$(stat -c%s "$f")" | tee -a "$LOG"
done
echo
log "Done."