#!/usr/bin/env bash
# Session 17（第 10 步）后处理：CNV 基因注释 + 临床相关性过滤
#
# 三步：
#   1) _called.cns → BED   由 cn 列推导 DEL/DUP（教程门槛：size ≥ 10 kb）
#   2) bedtools            与基因 BED 求交（CNV 覆盖了哪些基因）
#   3) 过滤                教程四项门槛，输出临床 review 表
#
# ⚠️ CNVkit 版本差异：教程（0.9.10）的 call 输出末列是 `DEL`/`DUP` 文字；
#   本项目 0.9.14 输出的是 **`cn` 数值列**（0/1=DEL，2=neutral，3+=DUP），
#   故本脚本由 cn 推导调用类型。列序：
#     chromosome start end gene log2 cn depth probes weight
#
# 教程的过滤门槛：
#   (1) size ≥ 10 kb            —— 30× 深度下更小的多为噪声
#   (2) log2 ≥ +0.2 (DUP) 或 ≤ -0.3 (DEL) —— 幅度过小在噪声范围内
#   (3) probes ≥ 3              —— 单 bin 调用不可靠
#   (4) 不在已知常见 CNP 区      —— 需 DGV（本项目未下载，跳过并在日志标注）
#
# 用法：bash scripts/10.run_cnv_annotation.sh

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

CNV_BIN=/oldhome/ouyjh/miniforge3/envs/cnv_env/bin
export PATH="$CNV_BIN:/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH"

OUT="$ROOT/10.cnv/results"
LOG="$WS/logs/10.run_cnv_annotation.log"
CALLED="$OUT/${SAMPLE_SHORT}_called.cns"
GENE_BED="$WS/rawData/resources/Homo_sapiens_GRCh38_genes.bed"

CALLS_BED="$OUT/${SAMPLE_SHORT}_cnv_calls.bed"
ANNO_TXT="$OUT/${SAMPLE_SHORT}_cnv_annotated.txt"
FILT_TSV="$OUT/${SAMPLE_SHORT}_cnv_filtered.tsv"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 17 post: CNV annotation + filtering"
echo "=============================================="
echo

# ── Pre-flight ──
log "Pre-flight checks..."
[ -f "$CALLED" ] || { echo "ERROR: missing $CALLED（先跑 10.run_cnvkit.sh）" | tee -a "$LOG"; exit 1; }
command -v bedtools >/dev/null || {
  echo "ERROR: bedtools not found. Run: bash scripts/10.setup_cnv_env.sh" | tee -a "$LOG"; exit 1; }
[ -f "$GENE_BED" ] || { echo "ERROR: missing $GENE_BED" | tee -a "$LOG"; exit 1; }
log "bedtools: $(bedtools --version 2>&1 | head -1)"
log "gene BED: $GENE_BED（$(wc -l < "$GENE_BED") 个基因）"

# 染色体命名一致性检查
CNV_CHR=$(awk 'NR==2{print $1}' "$CALLED")
GENE_CHR=$(awk 'NR==1{print $1}' "$GENE_BED")
log "chr 命名：CNV='$CNV_CHR'  基因 BED='$GENE_CHR'"
if [ "$CNV_CHR" != "$GENE_CHR" ]; then
  log "WARN 命名不一致 → 给基因 BED 加 chr 前缀"
  awk 'BEGIN{OFS="\t"} {if($1!~/^chr/)$1="chr"$1; print}' "$GENE_BED" \
    > "$WS/rawData/resources/Homo_sapiens_GRCh38_genes_chr.bed"
  GENE_BED="$WS/rawData/resources/Homo_sapiens_GRCh38_genes_chr.bed"
fi
echo

# ── Step 1: _called.cns → BED（由 cn 推导 DEL/DUP）──
log "[1/4] 由 cn 列推导 DEL/DUP → BED..."
if [ ! -f "$CALLS_BED" ]; then
python3 - "$CALLED" "$CALLS_BED" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, csv
src, dst = sys.argv[1], sys.argv[2]
n_all = n_kept = 0
with open(src) as fi, open(dst, 'w') as fo:
    rd = csv.DictReader(fi, delimiter='\t')
    for r in rd:
        n_all += 1
        try:
            cn = int(float(r['cn']))
        except (KeyError, TypeError, ValueError):
            continue
        call = 'DEL' if cn < 2 else ('DUP' if cn > 2 else 'neutral')
        if call == 'neutral':
            continue
        start, end = int(float(r['start'])), int(float(r['end']))
        size = end - start
        if size < 10000:          # 门槛 (1) size ≥ 10 kb
            continue
        n_kept += 1
        # BED: chrom start end call log2 size
        fo.write(f"{r['chromosome']}\t{start}\t{end}\t{call}\t{r.get('log2','')}\t{size}\n")
print(f"  总分割段 {n_all} 个 → DEL/DUP 且 ≥10kb 的 {n_kept} 个")
PYEOF
else
  log "  SKIP BED 已存在"
fi
echo

# ── Step 2: bedtools 基因注释 ──
log "[2/4] bedtools 基因注释..."
if [ ! -f "$ANNO_TXT" ]; then
  # -wa -wb：输出 A 记录 + B 记录；B 的第 4 列是基因名（故输出 $10）
  bedtools intersect -a "$CALLS_BED" -b "$GENE_BED" -wa -wb 2>>"$LOG" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$10}' \
  > "$ANNO_TXT"
  log "  注释完成：$(wc -l < "$ANNO_TXT") 行（CNV × 基因 对）"
else
  log "  SKIP 注释已存在"
fi
echo

# ── Step 3: 过滤 + 汇总 ──
log "[3/4] 应用教程过滤门槛..."
python3 - "$CALLS_BED" "$ANNO_TXT" "$CALLED" "$FILT_TSV" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, csv, collections
calls_bed, anno_txt, called_cns, out = sys.argv[1:5]

# 基因注释：按 (chrom,start,end) 聚合
gene_of = collections.defaultdict(list)
with open(anno_txt) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 7: continue
        gene_of[(p[0], p[1], p[2])].append(p[6])

# probes 数（门槛 3）：从 _called.cns 取
probes_of = {}
with open(called_cns) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        probes_of[(r['chromosome'], str(int(float(r['start']))), str(int(float(r['end']))))] = r.get('probes', '')

rows = []
with open(calls_bed) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 6: continue
        chrom, start, end, call, log2s, siz = p[:6]
        try: l2 = float(log2s)
        except ValueError: continue
        # 门槛 (2)：DUP 需 ≥ +0.2；DEL 需 ≤ -0.3
        if call == 'DUP' and l2 < 0.2:  continue
        if call == 'DEL' and l2 > -0.3: continue
        # 门槛 (3)：probes ≥ 3
        pr = probes_of.get((chrom, start, end), '')
        try:
            if int(float(pr)) < 3: continue
        except (TypeError, ValueError):
            pass
        genes = ';'.join(sorted(set(gene_of.get((chrom, start, end), [])))) or '-'
        rows.append((chrom, start, end, int(siz), call, f'{l2:.3f}', pr, genes))

with open(out, 'w') as f:
    f.write('chrom\tstart\tend\tsize_bp\tcall\tlog2\tprobes\tgenes\n')
    for r in rows:
        f.write('\t'.join(str(x) for x in r) + '\n')

print(f'  通过门槛 (1)(2)(3) 的调用：{len(rows)} 个')
print()
if rows:
    print('  === 临床 review 表 ===')
    for r in rows[:20]:
        print(f'    {r[0]}:{r[1]}-{r[2]}  {r[3]:>9} bp  {r[4]:3}  log2={r[5]:>7}  probes={r[6]:>4}  {r[7][:60]}')
else:
    print('  （无通过门槛的调用）')
PYEOF
echo

# ── Step 4: 门槛 4 —— DGV 过滤（剔除已知良性拷贝数多态）──
log "[4/4] 门槛 4：DGV 过滤（已知良性 CNP）..."
DGV_BED="$WS/rawData/databases/dgv/dgvGold.chr20.bed"
FILT_DGV="$OUT/${SAMPLE_SHORT}_cnv_final.tsv"

if [ ! -f "$DGV_BED" ]; then
  log "  SKIP —— DGV BED 不存在（先跑 bash scripts/10.download_dgv.sh）"
  log "  最终表 = 门槛 1–3 结果（$FILT_TSV）"
elif [ ! -s "$FILT_TSV" ] || [ "$(wc -l < "$FILT_TSV")" -le 1 ]; then
  log "  SKIP —— 门槛 1–3 后无剩余调用"
  cp "$FILT_TSV" "$FILT_DGV"
else
  log "  DGV chr20 条目：$(wc -l < "$DGV_BED")"
  # 取 CNV 调用（BED 前 3 列）
  awk 'NR>1{OFS="\t"; print $1,$2,$3}' "$FILT_TSV" > "$OUT/.tmp_calls.bed"
  # 与 DGV 求交：-f 0.5 = CNV 调用至少 50% 落在 DGV 条目内 → 判为已知良性
  bedtools intersect -a "$OUT/.tmp_calls.bed" -b "$DGV_BED" -f 0.5 -u 2>>"$LOG" \
    > "$OUT/.tmp_dgv_hit.bed"
  n_hit=$(wc -l < "$OUT/.tmp_dgv_hit.bed")
  log "  与 DGV 重叠 ≥50% 的调用：$n_hit 个（将被剔除）"

python3 - "$FILT_TSV" "$OUT/.tmp_dgv_hit.bed" "$FILT_DGV" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys
filt, hit, out = sys.argv[1], sys.argv[2], sys.argv[3]
bad = set()
with open(hit) as f:
    for line in f:
        p = line.split('\t')
        if len(p) >= 3:
            bad.add((p[0], p[1], p[2]))
kept, removed = [], []
with open(filt) as f:
    hdr = f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 3: continue
        (removed if (p[0], p[1], p[2]) in bad else kept).append(p)
with open(out, 'w') as f:
    f.write(hdr)
    for r in kept:
        f.write('\t'.join(r) + '\n')
print(f'  剔除 {len(removed)} 个已知良性多态 → 剩余 {len(kept)} 个候选')
if kept:
    for r in kept:
        print(f'    {r[0]}:{r[1]}-{r[2]}  {r[3]} bp  {r[4]}  log2={r[5]}  {r[7][:60]}')
else:
    print('    （无剩余候选 —— 全部为已知良性拷贝数多态）')
PYEOF
  rm -f "$OUT/.tmp_calls.bed" "$OUT/.tmp_dgv_hit.bed"
fi
echo

# ── 小结 ──
log "=== 产出 ==="
for f in "$CALLS_BED" "$ANNO_TXT" "$FILT_TSV" "$FILT_DGV"; do
  [ -f "$f" ] && printf "  %-56s %8s B\n" "$(basename "$f")" "$(stat -c%s "$f")" | tee -a "$LOG"
done
echo
log "Done."
