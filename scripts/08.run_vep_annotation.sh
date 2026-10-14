#!/usr/bin/env bash
# Session 14：Variant Annotation（**仅 VEP**）
# 流程第 8 步：给 Session 13 过滤后的变异加功能注释
#
# 本脚本只跑 VEP。第 08 步的另一个平行工具 SnpEff 用独立脚本：
#   scripts/08.run_snpeff_annotation.sh   （库小 627 MB，可先跑）
# 两者是平行备选（各读同一输入 VCF、各出各的输出），不是先后关系。
# 本脚本依赖已下载的 VEP cache（14 GB），由 scripts/download_vep_cache.sh 负责。
#
# 与官方教程的差异（必须注明）：
#   1. 教程用 VEP v113，我们 conda 装的是 VEP 105 → cache 必须用 105
#   2. 教程写 --af_gnomade / --af_gnomadg（gnomAD v4），VEP 105 无这两个选项，
#      对应的是 --af_gnomad（已确认 VEP 105 模块里只有 af_gnomad / af_1kg / af_exac / af_pops）
#   3. 教程用 `vep` 直接调用；本环境需显式指定 perl 与模块路径
#
# 输入：07.var.filtered/NA12878_chr20_filtered.vcf.gz（Session 13 输出，124,135 变异）
# 输出：08.annotation/
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

VCF_IN=$ROOT/07.var.filtered/${SAMPLE}_filtered.vcf.gz
OUT_ROOT=$ROOT/08.annotation
LOG=$WS/logs/08.annotation.log
# VEP 环境（conda vep_env，VEP 105）
VEP_PERL=/oldhome/ouyjh/miniforge3/envs/vep_env/bin/perl
VEP_LIB1=/oldhome/ouyjh/miniforge3/envs/vep_env/lib/site_perl/5.26.2
VEP_LIB2=/oldhome/ouyjh/miniforge3/envs/vep_env/lib/5.26.2
VEP_BIN=/oldhome/ouyjh/miniforge3/envs/vep_env/share/ensembl-vep-105.0-0/vep
CACHE_DIR=$WS/rawData/vep_cache
# 参考 FASTA：离线模式算 HGVS 坐标（--hgvs）**必须**提供（VEP 原报错：
# "Cannot generate HGVS coordinates (--hgvs) in offline mode without a FASTA file"）
# 需同目录有 .fai 索引（本项目已有）
REF=$WS/rawData/ref/Homo_sapiens_assembly38.fasta

VCF_OUT=$OUT_ROOT/${SAMPLE}_VEP.vcf.gz
TSV_ALL=$OUT_ROOT/${SAMPLE}_VEP_all.tsv
TSV_PRIO=$OUT_ROOT/${SAMPLE}_VEP_prioritised.tsv
STATS=$OUT_ROOT/${SAMPLE}_VEP_stats.html

mkdir -p "$OUT_ROOT" "$WS/logs"

# 工具路径（WGS env 里的 bcftools / tabix；VEP 用 vep_env 的 perl）
export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

echo "[$(date)] === 08. Variant Annotation (VEP 105) ===" | tee "$LOG"
echo "Input:  $VCF_IN" | tee -a "$LOG"
echo "Output: $OUT_ROOT/" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 0) Pre-flight ─────────────────────────────────────────────────
echo "[$(date)] Pre-flight checks..." | tee -a "$LOG"

if [ ! -f "$VCF_IN" ]; then
  echo "ERROR: input VCF not found: $VCF_IN" | tee -a "$LOG"
  exit 1
fi

if [ ! -f "$CACHE_DIR/homo_sapiens/105_GRCh38/info.txt" ]; then
  echo "ERROR: VEP cache missing. Run scripts/download_vep_cache.sh first." | tee -a "$LOG"
  echo "       expected: $CACHE_DIR/homo_sapiens/105_GRCh38/info.txt" | tee -a "$LOG"
  exit 1
fi

# 参考 FASTA + .fai（--hgvs 离线模式必须；缺失则 VEP 会立即报
# "Cannot generate HGVS coordinates (--hgvs) in offline mode without a FASTA file"）
if [ ! -f "$REF" ]; then
  echo "ERROR: reference FASTA not found: $REF" | tee -a "$LOG"
  exit 1
fi
if [ ! -f "${REF}.fai" ]; then
  echo "ERROR: reference FASTA index missing: ${REF}.fai" | tee -a "$LOG"
  echo "       fix: samtools faidx $REF" | tee -a "$LOG"
  exit 1
fi

N_IN=$(bcftools view -H "$VCF_IN" 2>/dev/null | wc -l)
echo "Input variants: $N_IN" | tee -a "$LOG"
echo "VEP cache:      $CACHE_DIR/homo_sapiens/105_GRCh38" | tee -a "$LOG"
echo "Reference:      $REF" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 1) VEP 注释（断点重入：输出已存在则跳过） ──────────────────────
if [ ! -f "$VCF_OUT" ]; then
  echo "[$(date)] Running VEP annotation (this may take 10-30 min)..." | tee -a "$LOG"

  # 各参数含义：
  #   --cache --offline : 用本地 cache，不发网络请求
  #   --dir_cache       : cache 目录（内含 homo_sapiens/105_GRCh38）
  #   --fasta           : 参考 FASTA（离线 + --hgvs 时**必需**，否则报错退出）
  #   --assembly GRCh38 : 参考基因组版本
  #   --sift b          : SIFT 预测 + 分数（b = both）
  #   --polyphen b      : PolyPhen-2 预测 + 分数
  #   --ccds            : 标注 CCDS 金标准编码转录本
  #   --hgvs            : HGVS 命名（如 p.Arg123Cys、c.367C>T）
  #   --symbol          : 基因符号（HGNC）
  #   --numbers         : 外显子/内含子编号
  #   --canonical       : 标记 canonical 转录本
  #   --biotype         : 转录本类型（protein_coding / lncRNA ...）
  #   --af_gnomad       : gnomAD 等位基因频率（VEP 105 的写法；105 无 af_gnomade/af_gnomadg）
  #   --pubmed          : 关联文献 PubMed ID（经 ClinVar）
  #   --uniprot         : 受累蛋白的 UniProt 编号
  #   --vcf             : 输出注释后的 VCF（注释在 INFO 的 CSQ 字段）
  #   --compress_output bgzip : BGZF 压缩（可建 tabix 索引）
  #   --fork 8          : 8 进程并行
  #   --stats_file      : HTML 统计报告
  $VEP_PERL \
    -I $VEP_LIB1 -I $VEP_LIB2 \
    $VEP_BIN \
    --input_file "$VCF_IN" \
    --output_file "$VCF_OUT" \
    --cache --offline \
    --dir_cache "$CACHE_DIR" \
    --fasta "$REF" \
    --assembly GRCh38 \
    --sift b \
    --polyphen b \
    --ccds \
    --hgvs \
    --symbol \
    --numbers \
    --canonical \
    --biotype \
    --af_gnomad \
    --pubmed \
    --uniprot \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --fork 8 \
    --stats_file "$STATS" \
    2>>"$LOG"

  echo "[$(date)] VEP done" | tee -a "$LOG"
else
  echo "[$(date)] VEP output exists, skip annotation" | tee -a "$LOG"
fi

# ── 2) 索引（tabix） ──────────────────────────────────────────────
if [ ! -f "${VCF_OUT}.tbi" ]; then
  echo "[$(date)] Indexing VEP output..." | tee -a "$LOG"
  tabix -p vcf "$VCF_OUT" 2>>"$LOG"
fi
ls -lh "$VCF_OUT"* | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 3) 解析 CSQ：导出全表 + 优先级筛选表 ──────────────────────────
echo "[$(date)] Parsing CSQ into tables..." | tee -a "$LOG"
python3 - "$VCF_OUT" "$TSV_ALL" "$TSV_PRIO" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, gzip, re

vcf_path, out_all, out_prio = sys.argv[1], sys.argv[2], sys.argv[3]

# 从 VCF 头部读出 CSQ 的字段顺序（不同版本字段不同，必须按表头解析）
csq_fields = None
opener = gzip.open if vcf_path.endswith('.gz') else open
with opener(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('##INFO=<ID=CSQ'):
            m = re.search(r'Format:\s*([^"]+)"', line)
            if m:
                csq_fields = m.group(1).split('|')
            break

if csq_fields is None:
    print("ERROR: CSQ header not found — VEP may have failed")
    sys.exit(1)

print(f"CSQ fields ({len(csq_fields)}): {','.join(csq_fields[:12])} ...")

# 关注的字段
want = ['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature','BIOTYPE','HGVSc','HGVSp',
        'SIFT','PolyPhen','CLIN_SIG','gnomAD_AF','Existing_variation','PUBMED']

# HIGH/MODERATE = 蛋白层面有影响的后果
def is_impactful(impact):
    return impact in ('HIGH', 'MODERATE')

def af_ok(af_str):
    """gnomAD AF < 0.01 视为罕见；缺失视为未知（保留）"""
    if not af_str or af_str == '':
        return True
    try:
        return float(af_str) < 0.01
    except ValueError:
        return True

n_total = 0
n_pass = 0
n_anno = 0
n_prio = 0
impact_count = {}

with opener(vcf_path, 'rt') as fin, \
     open(out_all, 'w') as fa, open(out_prio, 'w') as fp:
    hdr = ['CHROM','POS','REF','ALT','FILTER','Consequence','IMPACT','SYMBOL','Gene',
           'HGVSc','HGVSp','SIFT','PolyPhen','CLIN_SIG','gnomAD_AF','Existing_variation']
    fa.write('\t'.join(hdr) + '\n')
    fp.write('\t'.join(hdr) + '\n')

    for line in fin:
        if line.startswith('#'):
            continue
        n_total += 1
        c = line.rstrip('\n').split('\t')
        if len(c) < 8:
            continue
        chrom, pos, _id, ref, alt, _qual, filt, info = c[:8]
        if filt != 'PASS':
            continue
        n_pass += 1

        m = re.search(r'(?:^|;)CSQ=([^;]+)', info)
        if not m:
            continue

        # 一个位点可能有多个 CSQ（每个转录本一条）；取最严重的后果
        entries = m.group(1).split(',')
        best = None
        best_rank = 99
        rank = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
        for e in entries:
            vals = e.split('|')
            d = dict(zip(csq_fields, vals))
            r = rank.get(d.get('IMPACT', 'MODIFIER'), 3)
            if r < best_rank:
                best_rank = r
                best = d
        if best is None:
            continue
        n_anno += 1
        impact_count[best.get('IMPACT','?')] = impact_count.get(best.get('IMPACT','?'), 0) + 1

        row = [chrom, pos, ref, alt, filt] + [best.get(k, '') for k in
               ['Consequence','IMPACT','SYMBOL','Gene','HGVSc','HGVSp','SIFT',
                'PolyPhen','CLIN_SIG','gnomAD_AF','Existing_variation']]
        fa.write('\t'.join(row) + '\n')

        # 优先级筛选：蛋白层面有影响 且 gnomAD 频率罕见
        if is_impactful(best.get('IMPACT','')) and af_ok(best.get('gnomAD_AF','')):
            n_prio += 1
            fp.write('\t'.join(row) + '\n')

print()
print("=== Annotation summary ===")
print(f"Input records:        {n_total:,}")
print(f"PASS variants:        {n_pass:,}")
print(f"Annotated (has CSQ):  {n_anno:,}")
print(f"Prioritised (HIGH/MODERATE & gnomAD AF<0.01): {n_prio:,}")
print()
print("IMPACT distribution (PASS variants):")
for k in ('HIGH','MODERATE','LOW','MODIFIER'):
    if k in impact_count:
        print(f"  {k:9s}: {impact_count[k]:,}")
print()
print("Tables written:")
print(f"  {out_all}")
print(f"  {out_prio}")
PYEOF

echo | tee -a "$LOG"

# ── 4) QC GATE ────────────────────────────────────────────────────
echo "=== QC GATE (Post-Annotation) ===" | tee -a "$LOG"
echo "Check: VEP output non-empty, CSQ present, prioritised table non-empty" | tee -a "$LOG"
if [ -s "$TSV_ALL" ] && [ -s "$TSV_PRIO" ]; then
  echo "QC GATE: PASS — proceed to Session 15 (Variant databases)" | tee -a "$LOG"
else
  echo "QC GATE: FAIL — check VEP log" | tee -a "$LOG"
fi
echo | tee -a "$LOG"

echo "Output files:" | tee -a "$LOG"
ls -lh "$OUT_ROOT" | tee -a "$LOG"
echo "Next: Session 15 (Variant databases — ClinVar / gnomAD / OMIM-HPO)" | tee -a "$LOG"