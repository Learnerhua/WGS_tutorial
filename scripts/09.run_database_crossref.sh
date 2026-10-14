#!/usr/bin/env bash
# Session 15（Step 2–4/4）：变异数据库交叉查询
#
# 4 个子任务：
#   1. ClinVar isec       — 找与你的 VCF 交集的 ClinVar 收录变异
#   2. ClinVar annotate   — 把 CLNSIG/CLNREVSTAT/CLNDN 加进 VCF INFO
#   3. gnomAD annotate    — 把 AF/各人种频率加进 VCF INFO
#   4. OMIM + HPO 查询    — 给 43 个候选基因映射疾病与表型
#
# 用法：bash scripts/09.run_database_crossref.sh
# 预计耗时：5–10 分钟（gnomAD 11.8 GB 注释占大头）
#
# 注意：本脚本所有 heredoc 与 `if` 平级（不缩进在 if 里），
#       避免 bash 5.2 对嵌套 heredoc 的解析误报。

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

LOG="$WS/logs/09.database_crossref.log"
mkdir -p "$ROOT/09.databases" "$WS/logs"
export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

VCF_IN="$ROOT/07.var.filtered/${SAMPLE}_filtered.vcf.gz"
PRIO="$ROOT/08.annotation/${SAMPLE}_VEP_prioritised.tsv"

CLINVAR="$WS/rawData/databases/clinvar/clinvar.vcf.gz"
GNOMAD="$WS/rawData/databases/gnomad/gnomad.genomes.v4.1.sites.chr20.vcf.bgz"
CONSTRAINT="$WS/rawData/databases/gnomad/gnomad.v4.1.constraint_metrics.tsv"
OMIM="$WS/rawData/databases/omim/mim2gene.txt"
HPO_PHENO2GENE="$WS/rawData/databases/hpo/hpo_phenotype_to_genes.txt"
HPO_GENE2PHENO="$WS/rawData/databases/hpo/hpo_genes_to_phenotype.txt"

OUT="$ROOT/09.databases"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 15 step 2-4/4: database crossref"
echo "=============================================="
echo

# ── Pre-flight ──
log "Pre-flight checks..."
for f in "$VCF_IN" "$PRIO" "$CLINVAR" "$GNOMAD" "$CONSTRAINT" "$OMIM" \
         "$HPO_PHENO2GENE" "$HPO_GENE2PHENO"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing $f" | tee -a "$LOG"
    exit 1
  fi
done

# ── Step 0: 染色体名修正（ClinVar 1→chr1）──
CLINVAR_CHR="$OUT/clinvar.chr_prefix.vcf.gz"
if [ ! -f "$CLINVAR_CHR" ]; then
log "Step 0: 生成 ClinVar chr-prefix 版本..."
printf '1\tchr1\n2\tchr2\n3\tchr3\n4\tchr4\n5\tchr5\n6\tchr6\n7\tchr7\n8\tchr8\n9\tchr9\n10\tchr10\n11\tchr11\n12\tchr12\n13\tchr13\n14\tchr14\n15\tchr15\n16\tchr16\n17\tchr17\n18\tchr18\n19\tchr19\n20\tchr20\n21\tchr21\n22\tchr22\nX\tchrX\nY\tchrY\nMT\tchrMT\n' > "$OUT/clinvar_chr_map.txt"
bcftools annotate --rename-chrs "$OUT/clinvar_chr_map.txt" "$CLINVAR" -Oz -o "$CLINVAR_CHR" 2>>"$LOG"
bcftools index -t "$CLINVAR_CHR" 2>>"$LOG"
else
log "Step 0: SKIP ClinVar chr-prefix 版本已存在"
fi
CLINVAR="$CLINVAR_CHR"

# ── Step 1: ClinVar isec ──
log "[1/4] ClinVar isec..."
if [ ! -f "$OUT/clinvar_isec/0000.vcf" ]; then
bcftools isec \
  "$VCF_IN" "$CLINVAR" \
  -p "$OUT/clinvar_isec" \
  -n =2 \
  --regions chr20 2>>"$LOG"
fi
# -n=2 时：0000.vcf = 交集（我们 VCF 的记录），0001.vcf = 交集（ClinVar 的记录）
n_both=$(grep -vc '^#' "$OUT/clinvar_isec/0000.vcf" 2>/dev/null || echo 0)
log "  ClinVar ∩ 你的 VCF（chr20）：$n_both 个"
if [ ! -f "$OUT/clinvar_isec/0000.vcf.gz" ] && [ -s "$OUT/clinvar_isec/0000.vcf" ]; then
  bgzip -c "$OUT/clinvar_isec/0000.vcf" > "$OUT/clinvar_isec/0000.vcf.gz"
  tabix -p vcf "$OUT/clinvar_isec/0000.vcf.gz" 2>>"$LOG"
fi

# ── Step 2: ClinVar annotate ──
log "[2/4] ClinVar annotate..."
CLINVAR_ANN="$OUT/${SAMPLE}_clinvar_annotated.vcf.gz"
if [ ! -f "$CLINVAR_ANN" ]; then
bcftools annotate \
  -a "$CLINVAR" \
  -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN \
  "$VCF_IN" -Oz -o "$CLINVAR_ANN" 2>>"$LOG"
bcftools index -t "$CLINVAR_ANN" 2>>"$LOG"
else
log "[2/4] SKIP ClinVar annotate 已存在"
fi

# 分类统计 + 精确提取致病。
# ⚠️ 不用 bcftools 正则：`CLNSIG ~ "pathogenic"` 大小写敏感（漏 Pathogenic），
#    而 `~ "[Pp]athogenic"` 会误匹配 Conflicting_classifications_of_pathogenicity。
#    故改用 Python 按 , / | 拆词元做精确匹配。
CLS="$OUT/clinvar_classification.tsv"
if [ ! -f "$CLS" ]; then
python3 - "$CLINVAR_ANN" "$CLS" <<'PYEOF' 2>&1 | tee -a "$LOG"
import subprocess, collections, sys
vcf, out = sys.argv[1], sys.argv[2]
r = subprocess.run(['bcftools','query',
    '-f','%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n', vcf],
    capture_output=True, text=True)
n = 0; cats = collections.Counter(); patho = []
for line in r.stdout.splitlines():
    p = line.split('\t')
    if len(p) < 5: continue
    sig = p[4]
    if sig in ('', '.'): continue
    n += 1
    toks = [t.strip() for t in sig.replace('|', ',').split(',') if t.strip()]
    for t in toks: cats[t] += 1
    if any(t in ('Pathogenic', 'Likely_pathogenic') for t in toks):
        patho.append(p)
with open(out, 'w') as f:
    f.write('Category\tCount\n')
    for k, v in cats.most_common():
        f.write(f'{k}\t{v}\n')
print(f'  有 CLNSIG 的变异：{n}')
print('  CLNSIG 分类分布：')
for k, v in cats.most_common():
    print(f'    {k:45} {v}')
print(f'  真·致病（Pathogenic / Likely_pathogenic）：{len(patho)}')
if patho:
    print('  致病变异明细：')
    for p in patho:
        print(f'    {p[0]}:{p[1]} {p[2]}>{p[3]}  {p[4]}  {p[6] if len(p)>6 else ""}')
PYEOF
fi

# ── Step 3: gnomAD annotate + 频率 tier 统计 ──
log "[3/4] gnomAD annotate..."
GNOMAD_ANN="$OUT/${SAMPLE}_gnomad_annotated.vcf.gz"
if [ ! -f "$GNOMAD_ANN" ]; then
# ⚠️ 关键：我们的 VCF 已带 GATK 的 INFO/AF，与 gnomAD 的 AF 撞名。
# bcftools annotate 遇到已存在的 INFO 字段会中止（曾只写出 37,030/114,156 条）。
# 解法：先 -x INFO/AF 删掉 GATK 的 AF，再注释 gnomAD 的 AF。
log "  注意：先删除 GATK 自带 INFO/AF 以避免与 gnomAD AF 撞名"
bcftools annotate -x INFO/AF -Oz -o "$OUT/.noAF.vcf.gz" "$VCF_IN" 2>>"$LOG"
bcftools index -t "$OUT/.noAF.vcf.gz" 2>>"$LOG"
bcftools annotate \
  -a "$GNOMAD" \
  -c INFO/AF,INFO/AF_afr,INFO/AF_eas,INFO/AF_nfe,INFO/AF_sas \
  -r chr20 \
  "$OUT/.noAF.vcf.gz" 2>>"$LOG" \
| bcftools view -f PASS -O z -o "$GNOMAD_ANN" 2>>"$LOG"
bcftools index -t "$GNOMAD_ANN" 2>>"$LOG"
rm -f "$OUT/.noAF.vcf.gz" "$OUT/.noAF.vcf.gz.tbi"
n_ann=$(bcftools view -H "$GNOMAD_ANN" 2>/dev/null | wc -l)
log "  gnomAD 注释后 PASS 记录数：$n_ann（应 ≈114,156）"
else
log "[3/4] SKIP gnomAD annotate 已存在"
fi

# 频率 tier 统计（独立于上面的断点，缺失即重算）
TIER="$OUT/gnomad_freq_tier.tsv"
if [ ! -f "$TIER" ] && [ -f "$GNOMAD_ANN" ]; then
log "  频率 tier 统计..."
python3 - "$GNOMAD_ANN" "$TIER" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, subprocess
vcf, out = sys.argv[1], sys.argv[2]
res = subprocess.run(['bcftools','query','-f','%INFO/AF\n', vcf],
                     capture_output=True, text=True)
novel = ultra = rare = common = 0
for line in res.stdout.splitlines():
    v = line.strip()
    if v == '' or v == '.':
        novel += 1
    else:
        try:
            f = float(v)
            if   f == 0:   novel += 1
            elif f < 0.001: ultra += 1
            elif f < 0.01:  rare  += 1
            else:           common += 1
        except ValueError:
            novel += 1
total = novel + ultra + rare + common
with open(out, 'w') as f:
    f.write('Tier\tCount\tPercent\n')
    f.write(f'novel (AF=0/缺失)\t{novel}\t{novel/total*100:.1f}%\n')
    f.write(f'ultra-rare (AF<0.001)\t{ultra}\t{ultra/total*100:.1f}%\n')
    f.write(f'rare (0.001<=AF<0.01)\t{rare}\t{rare/total*100:.1f}%\n')
    f.write(f'common (AF>=0.01)\t{common}\t{common/total*100:.1f}%\n')
print(f'  novel         : {novel:>5} ({novel/total*100:.1f}%)')
print(f'  ultra-rare    : {ultra:>5} ({ultra/total*100:.1f}%)')
print(f'  rare          : {rare:>5} ({rare/total*100:.1f}%)')
print(f'  common        : {common:>5} ({common/total*100:.1f}%)')
print(f'  total         : {total}')
PYEOF
fi

# ── Step 4: OMIM + HPO + gnomAD constraint（43 个候选基因）──
log "[4/4] OMIM + HPO + gnomAD constraint 查询..."
ANNO="$OUT/${SAMPLE}_43candidates_clinically_annotated.tsv"
if [ ! -f "$ANNO" ]; then
python3 - "$PRIO" "$ANNO" "$OMIM" "$HPO_PHENO2GENE" "$HPO_GENE2PHENO" "$CONSTRAINT" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, csv
prio, out, omim, hpo_p2g, hpo_g2p, constraint = sys.argv[1:7]

# 读优先级候选（按基因聚合，保留变异位点）
gene_vars = {}
with open(prio) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row.get('SYMBOL') or row.get('Gene', '')
        if not g: continue
        key = f"{row.get('CHROM','')}:{row.get('POS','')}"
        gene_vars.setdefault(g, []).append(key)
genes = list(gene_vars.keys())
n_var = sum(len(v) for v in gene_vars.values())
print(f'  候选基因：{len(genes)} 个（对应 {n_var} 个变异）')

# OMIM mim2gene.txt
omim_gene2mim = {}
with open(omim) as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 4: continue
        mim, type_, entrez, symbol = p[0], p[1], p[2], p[3]
        if type_ in ('gene','phenotype','both','predominantly phenotypes'):
            omim_gene2mim.setdefault(symbol, []).append((mim, type_))
print(f'  OMIM 已索引 {len(omim_gene2mim)} 个基因')

# HPO phenotype_to_genes.txt
hpo_p2g_map = {}
with open(hpo_p2g) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row.get('gene_symbol', row.get('#gene_id', ''))
        t = row.get('hpo_id', row.get('#hpo_id', ''))
        if g and t: hpo_p2g_map.setdefault(g, []).append(t)
print(f'  HPO gene→phenotype 已索引 {len(hpo_p2g_map)} 个基因')

# gnomAD constraint（列名是 lof.pLI / lof.oe_ci.upper(LOEUF) / mis.z_score；
# 每基因多转录本 → 只取 canonical==true 那一行）
constraint_map = {}
with open(constraint) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row.get('gene', '')
        if not g: continue
        if row.get('canonical') != 'true' and g in constraint_map: continue
        constraint_map[g] = {
            'pLI':          row.get('lof.pLI', ''),
            'oe_lof_upper': row.get('lof.oe_ci.upper', ''),   # LOEUF
            'mis_z':        row.get('mis.z_score', ''),
        }
print(f'  gnomAD constraint 已索引 {len(constraint_map)} 个基因')

# 写输出
with open(out, 'w') as f:
    f.write('Gene\tN_var\tVariants\tOMIM_MIM\tHPO_terms\tgnomAD_pLI\tgnomAD_LOEUF\tgnomAD_mis_z\n')
    for g in genes:
        omim_list = omim_gene2mim.get(g, [])
        omim_str  = ';'.join(m for m,_ in omim_list) if omim_list else ''
        hpo_terms = sorted(set(hpo_p2g_map.get(g, [])))
        hpo_str   = ';'.join(hpo_terms) if hpo_terms else ''
        c = constraint_map.get(g, {})
        vs = gene_vars.get(g, [])
        f.write(f'{g}\t{len(vs)}\t{";".join(vs)}\t{omim_str}\t{hpo_str}\t'
                f'{c.get("pLI","")}\t{c.get("oe_lof_upper","")}\t{c.get("mis_z","")}\n')
print(f'  wrote: {out}')
print()
print('  候选基因中有 OMIM 关联的：',
      sum(1 for g in genes if g in omim_gene2mim), '/', len(genes))
print('  候选基因中有 HPO 关联的：',
      sum(1 for g in genes if g in hpo_p2g_map), '/', len(genes))
print('  候选基因中有 gnomAD constraint 数据的：',
      sum(1 for g in genes if g in constraint_map), '/', len(genes))
PYEOF
else
log "[4/4] SKIP 临床注释已存在"
fi

echo
log "=== 产出文件 ==="
ls -lh "$OUT/" | tail -n +2 | tee -a "$LOG"
echo
log "All Session 15 crossref tasks done."
log "Next: 把所有结果（hap.py F1 / ClinVar 命中 / gnomAD 频率分层 / 43 个候选临床注释）直接写进教程 4.10.4 节"
