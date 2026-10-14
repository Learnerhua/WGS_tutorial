#!/usr/bin/env bash
# Session 18（第 11 步）主流程：Manta + DELLY SV calling + SURVIVOR 合并
#
# 工具组合（互补策略）：
#   - Manta   (1.6.0)  ：配对读段 + 分裂读段联合打分 → diploidSV.vcf.gz
#   - DELLY   (1.2.6)  ：discordant pair 聚类 + 分裂读段精修
#   - SURVIVOR (1.0.7) ：多 caller 合并 → 取 2-caller 共识
#
# 教程数据差异：教程用 CEPH trio（NA12878 + NA12891 + NA12892）
#   本项目只有 NA12878（同 4.12 CNV 步）→ DELLY cohort genotyping 无法做
#   → DELLY 单样本 call，跳过 cohort joint genotyping，改用 gnomAD-SV 作人群频率参考
#
# 关键坑（已实测）：
#   1) Manta 预编译包用相对 import → 必须 cd 到其 bin/ 目录调用
#   2) runWorkflow.py 是 configManta.py **运行后生成在 runDir 内**的，不在 PATH
#   3) DELLY 1.2.6 只支持 4 类（DEL/DUP/INV/INS），无 TRA
#   4) bcftools index 对 BCF 必须用 -f（CSI），-t 会报 "TBI not supported for BCF"
#   5) delly merged 输出是 **sites-only**（无样本列）→ delly filter -f germline 不适用
#   6) SURVIVOR 1.0.7 **不能读 .vcf.gz**（当纯文本解析 → 输出乱码/0 条）→ 必须先解压为 .vcf
#
# 用法：
#   bash scripts/11.run_sv.sh                # 全流程
#   bash scripts/11.run_sv.sh manta|delly|merge|annotate
#
# 产物（11.sv/）见脚本末尾

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
CHROM="${CHROM:-chr20}"                                     # 染色体
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

# Manta：预编译包（Illumina GitHub release），放 ~/Download；用相对 import，须 cwd=bin/
MANTA=/oldhome/ouyjh/Download/manta-1.6.0.centos6_x86_64
export PATH="$MANTA/bin:/oldhome/ouyjh/miniforge3/envs/WGS/bin:/oldhome/ouyjh/miniforge3/envs/py27/bin:$PATH"

OUT="$ROOT/11.sv"
LOG="$WS/logs/11.run_sv.log"
mkdir -p "$OUT"/{manta,delly,merged} "$WS/logs"

REF="$WS/rawData/ref/Homo_sapiens_assembly38.fasta"
BAM="$ROOT/05.bqsr/${SAMPLE_SHORT}_${CHROM}_BQSR.bam"
EXCL="$WS/rawData/resources/human.hg38.excl.tsv"
# gnomAD-SV v4.1（GRCh38, chr 前缀）—— 教程给的 v2.1 实为 GRCh37，见 11.download_sv_db.sh
GNOMAD_SV="$WS/rawData/databases/sv/gnomad.v4.1.sv.sites.vcf.gz"

MANTA_PASS="$OUT/manta/${SAMPLE_SHORT}_manta_pass.vcf.gz"
DELLY_MERGED="$OUT/delly/${SAMPLE_SHORT}_merged.bcf"
DELLY_PASS="$OUT/delly/${SAMPLE_SHORT}_delly_pass.vcf.gz"
MERGED_VCF="$OUT/merged/NA12878_chr20_SVs_merged.vcf.gz"
ANNO_TSV="$OUT/merged/NA12878_chr20_SVs_svdb_annotated.tsv"

ONLY="${1:-all}"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
want() { [ "$ONLY" = "all" ] || [ "$ONLY" = "$1" ]; }

echo "=============================================="
echo " Session 18: Manta + DELLY + SURVIVOR SV calling"
echo "=============================================="
echo " Sample : $SAMPLE"
echo " BAM    : $BAM"
echo " Out    : $OUT"
echo " Run    : $ONLY"
echo

# ── Pre-flight ──
log "Pre-flight checks..."
for f in "$BAM" "$BAM.bai" "$REF" "$REF.fai" "$EXCL" "$GNOMAD_SV" "$GNOMAD_SV.tbi"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" | tee -a "$LOG"; exit 1; }
done
for t in delly SURVIVOR bcftools bgzip tabix bedtools; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" | tee -a "$LOG"; exit 1; }
done
[ -x "$MANTA/bin/configManta.py" ] || { echo "ERROR: Manta not at $MANTA" | tee -a "$LOG"; exit 1; }
log "manta:    $(cd "$MANTA/bin" && ./configManta.py --version 2>&1 | head -1)"
log "delly:    $(delly --version 2>&1 | head -1)"
log "SURVIVOR: $(SURVIVOR 2>&1 | grep -i version | head -1)"

# ── 1) Manta ──
if want manta; then
log "[1/4] Manta germline single-sample..."
CHR_REGIONS="$OUT/manta/chr20_regions.bed.gz"
if [ ! -f "$OUT/manta/run/results/variants/diploidSV.vcf.gz" ]; then
  log "  生成 chr20 regions BED（从 .fai 取长度，避免硬编码）..."
  awk '$1=="chr20"{print $1"\t0\t"$2}' "$REF.fai" | bgzip > "$CHR_REGIONS"
  tabix -p bed "$CHR_REGIONS"

  log "  configManta..."
  rm -rf "$OUT/manta/run"
  ( cd "$MANTA/bin" && ./configManta.py \
      --bam "$BAM" \
      --referenceFasta "$REF" \
      --callRegions "$CHR_REGIONS" \
      --runDir "$OUT/manta/run" ) >>"$LOG" 2>&1

  log "  runWorkflow（本地 8 核；chr20 约 1.5 分钟）..."
  "$OUT/manta/run/runWorkflow.py" --mode local --jobs 8 >>"$LOG" 2>&1
else
  log "  SKIP Manta 已有结果"
fi

if [ ! -f "$MANTA_PASS" ]; then
  log "  过滤：PASS + 排除 IMPRECISE-only..."
  RAW="$OUT/manta/run/results/variants/diploidSV.vcf.gz"
  bcftools view -f PASS "$RAW" 2>>"$LOG" \
    | bcftools view -e 'INFO/IMPRECISE=1' 2>>"$LOG" \
    | bcftools view -i 'FILTER="PASS"' -Oz -o "$MANTA_PASS" 2>>"$LOG"
  bcftools index -f "$MANTA_PASS"
fi
log "  Manta PASS：$(bcftools view -H "$MANTA_PASS" | wc -l) 个"
bcftools query -f '%SVTYPE\n' "$MANTA_PASS" 2>/dev/null | sort | uniq -c | sort -rn | sed 's/^/    /' | tee -a "$LOG"
echo
fi

# ── 2) DELLY ──
if want delly; then
log "[2/4] DELLY germline（4 类分跑）..."
if [ ! -f "$DELLY_MERGED" ]; then
  # delly 1.2.6 只支持 DEL/DUP/INV/INS（无 TRA）
  for SV in DEL DUP INV INS; do
    log "  delly call -t $SV..."
    delly call -t "$SV" -g "$REF" -x "$EXCL" \
      -o "$OUT/delly/${SAMPLE_SHORT}_${SV}.bcf" "$BAM" >>"$LOG" 2>&1
  done
  log "  delly merge 4 类..."
  delly merge \
    "$OUT/delly/${SAMPLE_SHORT}_DEL.bcf" \
    "$OUT/delly/${SAMPLE_SHORT}_DUP.bcf" \
    "$OUT/delly/${SAMPLE_SHORT}_INV.bcf" \
    "$OUT/delly/${SAMPLE_SHORT}_INS.bcf" \
    -o "$DELLY_MERGED" >>"$LOG" 2>&1
  bcftools index -f "$DELLY_MERGED"
else
  log "  SKIP delly merge 已有"
fi

log "  注：cohort genotyping 跳过（无 NA12891/92 BAM）→ 用 gnomAD-SV 替代 cohort 信息"

if [ ! -f "$DELLY_PASS" ]; then
  # merged 是 sites-only（无样本列）→ delly filter -f germline 不适用；直接筛 PASS
  log "  单样本模式：bcftools view -f PASS..."
  bcftools view -f PASS "$DELLY_MERGED" 2>>"$LOG" \
    | bcftools view -f PASS -Oz -o "$DELLY_PASS" 2>>"$LOG"
  bcftools index -f "$DELLY_PASS"
fi
log "  DELLY PASS：$(bcftools view -H "$DELLY_PASS" | wc -l) 个"
bcftools query -f '%SVTYPE\n' "$DELLY_PASS" 2>/dev/null | sort | uniq -c | sort -rn | sed 's/^/    /' | tee -a "$LOG"
echo
fi

# ── 3) SURVIVOR 合并 ──
if want merge; then
log "[3/4] SURVIVOR merge（Manta ∩ DELLY，2-caller 共识）..."
if [ ! -f "$MERGED_VCF" ]; then
  # 坑 6：SURVIVOR 不能读 .vcf.gz → 先解压为纯 VCF
  log "  解压为纯 VCF（SURVIVOR 不支持 bgzip 输入）..."
  PLAIN_MANTA="$OUT/merged/.manta_plain.vcf"
  PLAIN_DELLY="$OUT/merged/.delly_plain.vcf"
  bcftools view "$MANTA_PASS" -Ov -o "$PLAIN_MANTA" 2>>"$LOG"
  bcftools view "$DELLY_PASS" -Ov -o "$PLAIN_DELLY" 2>>"$LOG"

  printf "%s\n%s\n" "$PLAIN_MANTA" "$PLAIN_DELLY" > "$OUT/merged/sv_vcf_list.txt"
  log "  SURVIVOR merge 500 2 1 0 0 50..."
  SURVIVOR merge "$OUT/merged/sv_vcf_list.txt" 500 2 1 0 0 50 \
    "$OUT/merged/.merged_plain.vcf" >>"$LOG" 2>&1

  # SURVIVOR 输出**未按坐标排序** → 必须先 sort：
  #   否则 (a) 无法建索引，(b) IGV 拒绝加载，(c) bcftools -r 区域查询静默返回空
  log "  bcftools sort（SURVIVOR 输出未排序）..."
  bcftools sort "$OUT/merged/.merged_plain.vcf" -Ov \
    -o "$OUT/merged/.merged_sorted.vcf" 2>>"$LOG"

  bgzip -c "$OUT/merged/.merged_sorted.vcf" > "$MERGED_VCF"
  # 用 -t 建 .tbi（IGV 偏好 tbi；-f 会建 csi）
  bcftools index -t -f "$MERGED_VCF" 2>>"$LOG"
  rm -f "$PLAIN_MANTA" "$PLAIN_DELLY" "$OUT/merged/.merged_plain.vcf" "$OUT/merged/.merged_sorted.vcf"
fi
n_merged=$(bcftools view -H "$MERGED_VCF" | wc -l)
log "  SURVIVOR 共识（≥2 caller）：$n_merged 个"
bcftools query -f '%SVTYPE\n' "$MERGED_VCF" 2>/dev/null | sort | uniq -c | sort -rn | sed 's/^/    /' | tee -a "$LOG"
echo
fi

# ── 4) SV 数据库（gnomAD-SV v4.1, GRCh38）人群频率注释 ──
if want annotate; then
log "[4/4] gnomAD-SV v4.1 人群频率注释（按位置重叠 ≥50%）..."

if [ ! -f "$ANNO_TSV" ]; then
  # 我们的 SV → BED
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$MERGED_VCF" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} $3!="" && $3+0>0 {print $1, $2-1, $3, $4, $5}' > "$OUT/merged/.ours.bed"
  # gnomAD-SV v4.1 chr20 → BED（v4.1 已是 chr 前缀 + hg38，无需改名 / liftOver）
  bcftools query -r chr20 -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' "$GNOMAD_SV" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} $3!="" && $3+0>0 {print $1, $2-1, $3, $4}' > "$OUT/merged/.gnomad.bed"

  log "  我们的 SV: $(wc -l < "$OUT/merged/.ours.bed") 个 / gnomAD-SV v4.1 chr20: $(wc -l < "$OUT/merged/.gnomad.bed") 个"

  # 重叠 ≥50% **且 SVTYPE 一致** 判为「已知群体 SV」
  #   -wa -wb 输出 A(1-5) + B(6-9)；B 的 SVTYPE 在第 9 列
  bedtools intersect -a "$OUT/merged/.ours.bed" -b "$OUT/merged/.gnomad.bed" -f 0.5 -wa -wb 2>>"$LOG" \
    | awk 'BEGIN{OFS="\t"} $4==$9 {print $1, $2, $3}' | sort -u > "$OUT/merged/.hit.bed"

  {
    printf "chrom\tstart\tend\tsvtype\tsvlen\tin_gnomad_sv\tpriority\n"
    awk 'BEGIN{OFS="\t"}
      NR==FNR {hit[$1":"$2":"$3]=1; next}
      { key=$1":"$2":"$3
        ing = (key in hit) ? "YES" : "NO"
        pri = (ing=="YES") ? "LOW(known population SV)" : "HIGH(novel, needs clinical review)"
        print $1, $2, $3, $4, $5, ing, pri }' \
      "$OUT/merged/.hit.bed" "$OUT/merged/.ours.bed"
  } > "$ANNO_TSV"

  rm -f "$OUT/merged/.ours.bed" "$OUT/merged/.gnomad.bed" "$OUT/merged/.hit.bed"
fi

n_all=$(($(wc -l < "$ANNO_TSV") - 1))
n_known=$(awk -F'\t' 'NR>1 && $6=="YES"' "$ANNO_TSV" | wc -l)
n_novel=$(awk -F'\t' 'NR>1 && $6=="NO"' "$ANNO_TSV" | wc -l)
log "  合计 $n_all 个：gnomAD-SV 已知 $n_known / novel $n_novel"
echo
log "  === 各类型分布 ==="
awk -F'\t' 'NR>1{c[$4"\t"$6]++} END{for(k in c){split(k,a,"\t"); printf "    %-10s %-4s %d\n", a[1], (a[2]=="YES"?"known":"novel"), c[k]}}' "$ANNO_TSV" | sort | tee -a "$LOG"
echo
fi

# ── 小结 ──
log "=== 产出 ==="
for f in "$MANTA_PASS" "$DELLY_PASS" "$MERGED_VCF" "$ANNO_TSV"; do
  [ -f "$f" ] && printf "  %-60s %9s B\n" "${f#$WS/}" "$(stat -c%s "$f")" | tee -a "$LOG"
done
echo
log "Done. Next: bash scripts/11.run_sv_annotation.sh  （基因注释 + 临床 review）"