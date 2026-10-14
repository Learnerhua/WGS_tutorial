#!/usr/bin/env bash
# Session 14（第 08 步）备选工具：SnpEff 快速本地注释
#
# 本脚本只跑 SnpEff。第 08 步的另一个平行工具 VEP 用独立脚本：
#   scripts/download_vep_cache.sh          （先下 cache 14 GB）
#   scripts/08.run_vep_annotation.sh       （再跑 VEP 注释）
# 两者是平行备选（各读同一输入 VCF、各出各的输出），不是先后关系。
#
# 背景 / 为什么不用 `snpEff download`：
#   snpEff 5.2-3 自动下载只会请求仓库路径 v5_2/snpEff_v5_2_GRCh38.105.zip，
#   但该文件在 v5_2 下**不存在**（v5_2 直接跳到 GRCh38.110）→ 404。
#   实际存在的路径是 v5_1/snpEff_v5_1_GRCh38.105.zip（627 MB）。
#   因此本脚本改为手动下载 + 解压到项目自管目录，再用 -dataDir 指过去。
#
# 与 VEP 的关系：两者都属于第 08 步（Session 14），是**平行的备选工具**，
#   不是先后关系。VEP 是教程主线（cache 约 14 GB，下载较慢）；
#   SnpEff 数据库仅 627 MB，可先跑起来做快速注释。
#
# 输入：07.var.filtered/NA12878_chr20_filtered.vcf.gz
# 输出：08.annotation/snpeff/
set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

VCF_IN=$ROOT/07.var.filtered/${SAMPLE}_filtered.vcf.gz
OUT_ROOT=$ROOT/08.annotation/snpeff
LOG=$WS/logs/08.snpeff.log
# SnpEff 数据库
SNPEFF_BASE=$WS/rawData/snpeff_db          # 解压根目录（zip 内是 data/<genome>/...）
DB_ZIP=/tmp/snpEff_v5_1_GRCh38.105.zip
DB_URL="https://snpeff-public.s3.amazonaws.com/databases/v5_1/snpEff_v5_1_GRCh38.105.zip"
DB_SIZE=627276234                          # 官方字节数（用于完整性校验）
GENOME=GRCh38.105                          # Ensembl 105 注释，与 VEP 105 版本对齐

VCF_OUT=$OUT_ROOT/${SAMPLE}_snpeff.vcf.gz
STATS_HTML=$OUT_ROOT/${SAMPLE}_snpeff_summary.html
STATS_CSV=$OUT_ROOT/${SAMPLE}_snpeff_summary.csv

# 代理（VPN）
export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

MAX_RETRY=15

mkdir -p "$OUT_ROOT" "$SNPEFF_BASE" "$WS/logs"

echo "[$(date)] === 08. SnpEff annotation (第 08 步备选工具) ===" | tee "$LOG"
echo "Input:  $VCF_IN" | tee -a "$LOG"
echo "Output: $OUT_ROOT/" | tee -a "$LOG"
echo "Genome: $GENOME (Ensembl 105)" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 0) Pre-flight ─────────────────────────────────────────────────
if [ ! -f "$VCF_IN" ]; then
  echo "ERROR: input VCF not found: $VCF_IN" | tee -a "$LOG"; exit 1
fi
if ! command -v snpEff >/dev/null 2>&1; then
  echo "ERROR: snpEff not found in PATH" | tee -a "$LOG"; exit 1
fi
echo "SnpEff: $(snpEff -version 2>&1 | head -1)" | tee -a "$LOG"
echo "Input variants: $(bcftools view -H "$VCF_IN" 2>/dev/null | wc -l)" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 1) 下载数据库（断点续传 + 重试 + 字节校验）──────────────────────
DB_OK=0
if [ -f "$SNPEFF_BASE/data/$GENOME/snpEffectPredictor.bin" ] || [ -f "$SNPEFF_BASE/$GENOME/snpEffectPredictor.bin" ]; then
  echo "[$(date)] SnpEff database already installed, skip download" | tee -a "$LOG"
  DB_OK=1
else
  CUR=0
  [ -f "$DB_ZIP" ] && CUR=$(stat -c%s "$DB_ZIP")
  echo "[$(date)] Downloading SnpEff DB ($GENOME, $DB_SIZE bytes). Resume point: $CUR" | tee -a "$LOG"

  attempt=0
  while [ "$attempt" -lt "$MAX_RETRY" ]; do
    attempt=$((attempt + 1))
    echo "[$(date)] attempt $attempt/$MAX_RETRY ..." | tee -a "$LOG"
    wget -c --tries=3 --waitretry=10 --timeout=60 --read-timeout=120 \
      "$DB_URL" -O "$DB_ZIP" || true
    CUR=$(stat -c%s "$DB_ZIP" 2>/dev/null || echo 0)
    echo "[$(date)] got $CUR / $DB_SIZE bytes" | tee -a "$LOG"
    [ "$CUR" -ge "$DB_SIZE" ] && break
    sleep 5
  done

  if [ "$CUR" -ne "$DB_SIZE" ]; then
    echo "ERROR: download incomplete ($CUR / $DB_SIZE). Re-run to resume." | tee -a "$LOG"
    exit 1
  fi
  echo "Size check OK: $CUR bytes" | tee -a "$LOG"
fi

# ── 2) 解压（zip 内含 data/<genome>/...）──────────────────────────
if [ "$DB_OK" -eq 0 ]; then
  echo "[$(date)] Extracting to $SNPEFF_BASE ..." | tee -a "$LOG"
  unzip -o -q "$DB_ZIP" -d "$SNPEFF_BASE" 2>>"$LOG"
fi

# ── 3) 自动定位 dataDir（-dataDir 覆盖 config 的 data.dir）─────────
if [ -f "$SNPEFF_BASE/data/$GENOME/snpEffectPredictor.bin" ]; then
  DATA_DIR="$SNPEFF_BASE/data"
elif [ -f "$SNPEFF_BASE/$GENOME/snpEffectPredictor.bin" ]; then
  DATA_DIR="$SNPEFF_BASE"
else
  echo "ERROR: snpEffectPredictor.bin not found under $SNPEFF_BASE" | tee -a "$LOG"
  echo "  expected: $SNPEFF_BASE/data/$GENOME/ or $SNPEFF_BASE/$GENOME/" | tee -a "$LOG"
  exit 1
fi
echo "Data dir resolved: $DATA_DIR" | tee -a "$LOG"
echo "  database size: $(du -sh "$DATA_DIR/$GENOME" 2>/dev/null | cut -f1)" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 4) 运行 SnpEff 注释（断点重入）─────────────────────────────────
if [ ! -f "$VCF_OUT" ]; then
  echo "[$(date)] Running SnpEff annotation ..." | tee -a "$LOG"

  # -nodownload  : 禁止自动下载（仓库 v5_2 路径 404，必须显式用本地库）
  # -dataDir     : 指向含 <genome>/ 子目录的数据目录
  # -stats/-csvStats : 生成 HTML / CSV 统计报告
  # -v           : 详细日志
  # -Xmx16g      : Java heap 16 GB。GRCh38.105 数据库加载阶段（转录本/蛋白模型
  #               全装进内存）默认 8 GB 会 OOM，实测需 ≥ 16 GB 稳定运行
  # 说明：默认注释全部转录本；如需只保留 canonical 转录本可加 -canon
  snpEff -nodownload -dataDir "$DATA_DIR" \
    -Xmx16g \
    -stats "$STATS_HTML" -csvStats "$STATS_CSV" \
    -v "$GENOME" "$VCF_IN" > "$VCF_OUT.raw" 2>>"$LOG"

  echo "[$(date)] SnpEff done, compressing ..." | tee -a "$LOG"
  bgzip -c "$VCF_OUT.raw" > "$VCF_OUT"
  rm -f "$VCF_OUT.raw"
  tabix -p vcf "$VCF_OUT" 2>>"$LOG"
else
  echo "[$(date)] SnpEff output exists, skip annotation" | tee -a "$LOG"
fi
ls -lh "$VCF_OUT"* | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 5) 解析 ANN 字段：统计 consequence / IMPACT ────────────────────
echo "[$(date)] Parsing ANN field ..." | tee -a "$LOG"
python3 - "$VCF_OUT" "$OUT_ROOT/${SAMPLE}_snpeff_impact.tsv" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, gzip, re
from collections import Counter

vcf, out_tsv = sys.argv[1], sys.argv[2]
opener = gzip.open if vcf.endswith('.gz') else open

n_total = n_pass = n_ann = 0
impact = Counter()
consq = Counter()
rows = []

with opener(vcf, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        n_total += 1
        c = line.rstrip('\n').split('\t')
        if len(c) < 8:
            continue
        chrom, pos, _id, ref, alt, _q, filt, info = c[:8]
        if filt != 'PASS':
            continue
        n_pass += 1

        # SnpEff: ANN=Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|...
        m = re.search(r'(?:^|;)ANN=([^;]+)', info)
        if not m:
            continue
        n_ann += 1

        best = None
        rank = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
        br = 99
        for e in m.group(1).split(','):
            v = e.split('|')
            if len(v) < 4:
                continue
            imp = v[2]
            r = rank.get(imp, 3)
            if r < br:
                br, best = r, v
        if best is None:
            continue
        impact[best[2]] += 1
        # 一条注释可能含多个 consequence，用 & 连接
        for cc in best[1].split('&'):
            consq[cc] += 1
        rows.append((chrom, pos, ref, alt, best[1], best[2], best[3] if len(best) > 3 else '',
                     best[9] if len(best) > 9 else '', best[10] if len(best) > 10 else ''))

print()
print("=== SnpEff annotation summary ===")
print(f"Input records:       {n_total:,}")
print(f"PASS variants:       {n_pass:,}")
print(f"Annotated (has ANN): {n_ann:,}")
print()
print("IMPACT distribution:")
for k in ('HIGH', 'MODERATE', 'LOW', 'MODIFIER'):
    if k in impact:
        print(f"  {k:9s}: {impact[k]:,}")
print()
print("Top consequences:")
for k, v in consq.most_common(12):
    print(f"  {k:34s}: {v:,}")

with open(out_tsv, 'w') as o:
    o.write("CHROM\tPOS\tREF\tALT\tConsequence\tIMPACT\tGene\tHGVSc\tHGVSp\n")
    for r in rows:
        o.write('\t'.join(str(x) for x in r) + '\n')
print()
print(f"Table written: {out_tsv}")
PYEOF

echo | tee -a "$LOG"

# ── 6) QC GATE ────────────────────────────────────────────────────
echo "=== QC GATE (SnpEff) ===" | tee -a "$LOG"
if [ -s "$VCF_OUT" ] && [ -s "$OUT_ROOT/${SAMPLE}_snpeff_impact.tsv" ]; then
  echo "QC GATE: PASS — SnpEff annotation available" | tee -a "$LOG"
else
  echo "QC GATE: FAIL — check $LOG" | tee -a "$LOG"
fi
echo | tee -a "$LOG"
echo "Outputs:" | tee -a "$LOG"
ls -lh "$OUT_ROOT" | tee -a "$LOG"