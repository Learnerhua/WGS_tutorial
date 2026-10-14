#!/usr/bin/env bash
# Session 17 / 18 后续数据准备（约 215 MB）
#
# 4 项（教程给的 5 项里 1 项不存在，自动替换）：
#   1. Ensembl GRCh38 基因 BED  ← 教程 URL 不存在；改用 GTF 提取
#   2. DELLY hg38 排除区域 TSV
#   3. 1000G Phase 3 SV (hg19) VCF + tbi
#   4. gnomAD SV v2.1 (hg38) VCF + tbi
#
# 用法：bash scripts/09.download_databases.sh session17 session18 [任选]
#       不传参数默认全部下载

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
LOG="$WS/logs/09.download_extras.log"
PROXY=http://localhost:9999
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"
RETRIES=10

mkdir -p "$WS/rawData/databases/sv" "$WS/rawData/resources"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# ── 期望字节数（实测）──
declare -A EXP=(
  [Homo_sapiens.GRCh38.111.chr.gtf.gz]=54346520
  [human.hg38.excl.tsv]=84394
  [ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz]=18298662
  [ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi]=666732
  [gnomad_v2.1_sv.sites.vcf.gz]=140885330
  [gnomad_v2.1_sv.sites.vcf.gz.tbi]=452300
)
declare -A URL=(
  [Homo_sapiens.GRCh38.111.chr.gtf.gz]="https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.chr.gtf.gz"
  [human.hg38.excl.tsv]="https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv"
  [ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz]="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"
  [ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi]="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi"
  [gnomad_v2.1_sv.sites.vcf.gz]="https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz"
  [gnomad_v2.1_sv.sites.vcf.gz.tbi]="https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz.tbi"
)
declare -A DEST=(
  [Homo_sapiens.GRCh38.111.chr.gtf.gz]="$WS/rawData/resources/Homo_sapiens.GRCh38.111.chr.gtf.gz"
  [human.hg38.excl.tsv]="$WS/rawData/resources/human.hg38.excl.tsv"
  [ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz]="$WS/rawData/databases/sv/1000G_SVs_hg19.vcf.gz"
  [ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi]="$WS/rawData/databases/sv/1000G_SVs_hg19.vcf.gz.tbi"
  [gnomad_v2.1_sv.sites.vcf.gz]="$WS/rawData/databases/sv/gnomad_sv_v2.1.vcf.gz"
  [gnomad_v2.1_sv.sites.vcf.gz.tbi]="$WS/rawData/databases/sv/gnomad_sv_v2.1.vcf.gz.tbi"
)

# 组选择
WANT=("$@")
[ ${#WANT[@]} -eq 0 ] && WANT=(session17 session18)
want() { for x in "${WANT[@]}"; do [ "$x" = "$1" ] && return 0; done; return 1; }

# ── 通用下载：带重试 + 字节校验 + 断点续传 ──
fetch() {
  local key="$1"
  local url="${URL[$key]}" dest="${DEST[$key]}" exp="${EXP[$key]}"
  local cur=0

  if [ -f "$dest" ]; then
    cur=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$cur" = "$exp" ]; then
      log "SKIP  $key (already complete: $exp bytes)"
      return 0
    fi
    log "RESUME $key (have $cur / $exp bytes)"
  fi

  local n
  for n in $(seq 1 "$RETRIES"); do
    log "attempt $n/$RETRIES: $key"
    env http_proxy="$PROXY" https_proxy="$PROXY" all_proxy="$PROXY" \
      wget -c --tries=3 --retry-connrefused --waitretry=5 \
           --timeout=60 --read-timeout=120 --dns-timeout=30 \
           -U "$UA" -q -O "$dest" "$url" 2>>"$LOG" || true

    cur=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$cur" = "$exp" ]; then
      log "OK    $key ($exp bytes)"
      return 0
    fi
    log "  partial: $cur / $exp bytes"
    sleep 5
  done

  log "FAIL  $key — got $cur / $exp bytes after $RETRIES attempts"
  return 1
}

echo "=============================================="
echo " Session 17/18 extra databases (~215 MB)"
echo " Datasets: ${WANT[*]}"
echo "=============================================="
echo

FAILED=0

# Session 17：Ensembl 基因 BED（用 GTF 替代，下载完后提取）
if want session17; then
  echo "── Session 17: Ensembl GRCh38 gene annotation ──"
  fetch "Homo_sapiens.GRCh38.111.chr.gtf.gz" || FAILED=1
  echo
  if [ -f "${DEST[Homo_sapiens.GRCh38.111.chr.gtf.gz]}" ] \
     && [ ! -f "$WS/rawData/resources/Homo_sapiens_GRCh38_genes.bed" ]; then
    log "提取 GRCh38 基因 BED（仅 chr20）..."
    zcat "${DEST[Homo_sapiens.GRCh38.111.chr.gtf.gz]}" \
      | awk -F'\t' '$3=="gene" && ($1=="20" || $1=="chr20") {
          # split 第 9 列 attributes；找 "gene_name" 后面的引号内字段
          n = split($9, a, /[; \t]+/);
          for (i=1; i<=n; i++) {
            if (a[i] == "gene_name") {
              gsub(/"/, "", a[i+1]);
              chrom = ($1 ~ /^chr/) ? $1 : "chr"$1;
              print chrom, $4-1, $5, a[i+1];
              break;
            }
          }
        }' OFS='\t' \
      | sort -k1,1 -k2,2n -u \
      > "$WS/rawData/resources/Homo_sapiens_GRCh38_genes.bed"
    n=$(wc -l < "$WS/rawData/resources/Homo_sapiens_GRCh38_genes.bed" 2>/dev/null || echo 0)
    log "BED 提取完成：$n 行 chr20 基因"
  fi
  echo
fi

# Session 18：DELLY 排除 + 1000G SV + gnomAD SV
if want session18; then
  echo "── Session 18: SV population databases + DELLY excl ──"
  fetch "human.hg38.excl.tsv" || FAILED=1
  fetch "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"     || FAILED=1
  fetch "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi" || FAILED=1
  fetch "gnomad_v2.1_sv.sites.vcf.gz"     || FAILED=1
  fetch "gnomad_v2.1_sv.sites.vcf.gz.tbi" || FAILED=1
  echo
fi

echo "=============================================="
echo " Summary"
echo "=============================================="
du -sh "$WS/rawData/databases/sv" 2>/dev/null \
  | while read -r sz name; do printf "  %-10s %s\n" "$sz" "$name"; done
du -sh "$WS/rawData/resources" 2>/dev/null \
  | while read -r sz name; do printf "  %-10s %s\n" "$sz" "$name"; done
echo
if [ "$FAILED" -eq 0 ]; then
  echo "All requested datasets OK."
else
  echo "WARNING: some downloads failed — re-run to resume."
fi

[ "$FAILED" -eq 0 ] || exit 1
