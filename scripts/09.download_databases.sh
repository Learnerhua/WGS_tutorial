#!/usr/bin/env bash
# Session 15（第 09 步）数据准备：下载变异注释数据库
#
# 下载 5 组数据库到 rawData/databases/（~12.4 GB）：
#   1. ClinVar GRCh38      ~200 MB  临床致病性（五级分类 + 星级）
#   2. gnomAD v4.1 chr20   ~11.8 GB 人群频率（最大）
#   3. gnomAD constraint   ~96 MB   基因约束指标 pLI/LOEUF
#   4. GIAB HG001 v4.2.1   ~143 MB  真值集（hap.py 基准测试用）
#   5. HPO                 ~88 MB   基因-表型关联
#
# 用法：
#   bash scripts/09.download_databases.sh            # 全部
#   bash scripts/09.download_databases.sh giab       # 只下 GIAB（139 MB，跑 hap.py 用）
#   bash scripts/09.download_databases.sh giab clinvar hpo   # 可选组合
#   可选组名：clinvar gnomad gnomad_constraint giab hpo
#
# 说明：
#   - OMIM（mim2gene.txt）不在本脚本内：需在 omim.org 注册账号后手动下载
#   - NCBI 的 ftp-trace 主机经代理**不稳定**（同一 URL 会随机 404），故每个文件
#     都带「重试循环 + 字节数校验」，失败会自动重试而非中断整个下载
#   - 断点续传：wget -c；已下完的文件（字节数匹配）直接跳过
#   - 教程两处与实测不符（脚本已按实测修正）：
#       * gnomAD chr20 实际 11.8 GB（教程写 8 GB）
#       * GIAB high-conf BED 实际文件名 `..._benchmark.bed`（教程写 `_noinconsistent.bed`，不存在）

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
DB="$WS/rawData/databases"
LOG="$WS/logs/09.download_databases.log"
PROXY=http://localhost:9999
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"
RETRIES=10          # 每个文件的下载尝试次数（应对 ftp-trace 随机 404）

mkdir -p "$DB"/{clinvar,gnomad,giab,hpo} "$WS/logs"

# ── 数据集 URL 与预期字节数 ──
CLINVAR_BASE="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1"
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38"
HPO_BASE="https://purl.obolibrary.org/obo/hp/hpoa"

# 目标选择（默认全部）
WANT=("$@")
[ ${#WANT[@]} -eq 0 ] && WANT=(clinvar gnomad gnomad_constraint giab hpo)
want() { local x; for x in "${WANT[@]}"; do [ "$x" = "$1" ] && return 0; done; return 1; }

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# ── 核心：带重试 + 字节校验的下载 ──
# 参数：url  dest  expected_bytes  label
fetch() {
  local url="$1" dest="$2" exp="$3" label="$4"
  local cur=0

  if [ -f "$dest" ]; then
    cur=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$cur" = "$exp" ]; then
      log "SKIP  $label (already complete: $exp bytes)"
      return 0
    fi
    log "RESUME $label (have $cur / $exp bytes)"
  fi

  local n
  for n in $(seq 1 "$RETRIES"); do
    log "attempt $n/$RETRIES: $label"
    env http_proxy="$PROXY" https_proxy="$PROXY" all_proxy="$PROXY" \
      wget -c --tries=3 --retry-connrefused --waitretry=5 \
           --timeout=60 --read-timeout=120 --dns-timeout=30 \
           -U "$UA" -q -O "$dest" "$url" 2>>"$LOG" || true

    cur=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$cur" = "$exp" ]; then
      log "OK    $label ($exp bytes)"
      return 0
    fi
    log "  partial: $cur / $exp bytes"
    sleep 5
  done

  log "FAIL  $label — after $RETRIES attempts got $cur / $exp bytes"
  return 1
}

echo "=============================================="
echo " Session 15 (step 09) — database downloads"
echo " Target dir : $DB"
echo " Datasets   : ${WANT[*]}"
echo "=============================================="
echo

FAILED=0

# ── 1) ClinVar ──
if want clinvar; then
  echo "── ClinVar GRCh38 (~200 MB) ──"
  fetch "$CLINVAR_BASE/clinvar.vcf.gz"     "$DB/clinvar/clinvar.vcf.gz"     200145964 "ClinVar VCF" || FAILED=1
  fetch "$CLINVAR_BASE/clinvar.vcf.gz.tbi" "$DB/clinvar/clinvar.vcf.gz.tbi"    611775 "ClinVar tbi" || FAILED=1
  echo
fi

# ── 2) gnomAD chr20 genome VCF（最大，11.8 GB）──
if want gnomad; then
  echo "── gnomAD v4.1 chr20 genomes (~11.8 GB) ──"
  fetch "$GNOMAD_BASE/vcf/genomes/gnomad.genomes.v4.1.sites.chr20.vcf.bgz" \
        "$DB/gnomad/gnomad.genomes.v4.1.sites.chr20.vcf.bgz" 11838569571 "gnomAD chr20 VCF" || FAILED=1
  fetch "$GNOMAD_BASE/vcf/genomes/gnomad.genomes.v4.1.sites.chr20.vcf.bgz.tbi" \
        "$DB/gnomad/gnomad.genomes.v4.1.sites.chr20.vcf.bgz.tbi" 58550 "gnomAD chr20 tbi" || FAILED=1
  echo
fi

# ── 3) gnomAD constraint metrics ──
if want gnomad_constraint; then
  echo "── gnomAD v4.1 constraint metrics (~96 MB) ──"
  fetch "$GNOMAD_BASE/constraint/gnomad.v4.1.constraint_metrics.tsv" \
        "$DB/gnomad/gnomad.v4.1.constraint_metrics.tsv" 95546041 "gnomAD constraint TSV" || FAILED=1
  echo
fi

# ── 4) GIAB HG001 v4.2.1（hap.py 基准测试）──
if want giab; then
  echo "── GIAB HG001 GRCh38 v4.2.1 (~143 MB) ──"
  fetch "$GIAB_BASE/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
        "$DB/giab/HG001_truth.vcf.gz" 125932193 "GIAB truth VCF" || FAILED=1
  fetch "$GIAB_BASE/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
        "$DB/giab/HG001_truth.vcf.gz.tbi" 1589208 "GIAB truth tbi" || FAILED=1
  # 注意：教程写的 `..._noinconsistent.bed` 不存在，实际是 `..._benchmark.bed`
  fetch "$GIAB_BASE/HG001_GRCh38_1_22_v4.2.1_benchmark.bed" \
        "$DB/giab/HG001_highconf.bed" 15479939 "GIAB high-conf BED" || FAILED=1
  echo
fi

# ── 5) HPO ──
if want hpo; then
  echo "── HPO gene-phenotype files (~88 MB) ──"
  fetch "$HPO_BASE/phenotype_to_genes.txt" "$DB/hpo/hpo_phenotype_to_genes.txt" 67246520 "HPO phenotype_to_genes" || FAILED=1
  fetch "$HPO_BASE/genes_to_phenotype.txt" "$DB/hpo/hpo_genes_to_phenotype.txt" 20821704 "HPO genes_to_phenotype" || FAILED=1
  echo
fi

# ── 汇总 ──
echo "=============================================="
echo " Download summary"
echo "=============================================="
du -sh "$DB"/* 2>/dev/null | while read -r sz name; do printf "  %-10s %s\n" "$sz" "$name"; done
echo
if [ "$FAILED" -eq 0 ]; then
  echo "All requested datasets OK."
  echo "Next: bash scripts/10.run_happy_benchmark.sh   (hap.py F1 benchmark)"
else
  echo "WARNING: some downloads failed — re-run this script to resume."
fi

echo
echo "Not included (manual, needs account): OMIM mim2gene.txt -> https://omim.org/downloads"
[ "$FAILED" -eq 0 ] || exit 1
