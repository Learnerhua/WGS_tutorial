#!/usr/bin/env bash
# Session 17（第 10 步）主流程：CNVkit 拷贝数变异检测
#
# 教程 5 步（+2 步可视化）：
#   0) bins      — 生成 50 kb bin 坐标（chr20 全长）
#   1) coverage  — 计算每个 bin 的覆盖深度
#   2) reference — 构建参考（多样本，或单样本 flat reference）
#   3) fix       — 样本相对参考归一化 + GC 校正
#   4) segment   — CBS 分割
#   5) call      — 赋予拷贝数（cn 数值）
#   6) scatter   — 单染色体散点图
#   7) diagram   — 全基因组示意图
#
# ⚠️ CNVkit 版本差异（重要）
#   教程基于 **cnvkit 0.9.10**；本项目装的是 **0.9.14**（0.9.8 与 pandas 2.x 不兼容，
#   见 2.3.4）。0.9.14 的 CLI 有三处变化，本脚本已按 0.9.14 实测通过的写法：
#
#   | 步骤 | 教程 0.9.10 | 0.9.14 实际（本脚本）|
#   |---|---|---|
#   | coverage | `coverage BAM --method wgs --targets /dev/null` | `coverage BAM <bins.bed> -f REF`（无 --method；bin 由 BED 给出）|
#   | reference | `reference n1.cnn n2.cnn -f REF` | 同左；**flat 模式**为 `reference -f REF -t <bins.bed> --no-edge` |
#   | fix | `fix bins.cnn /dev/null ref.cnn` | `fix bins.cnn -r ref.cnn --no-edge`（antitarget 直接省略）|
#   | call 输出 | 末列 `DEL`/`DUP` 文字 | **`cn` 数值列**（0/1=DEL，2=neutral，3+=DUP）|
#
#   `--method wgs` 在 0.9.14 移到了 `batch` 子命令；单步流程改用显式 bins BED。
#
# ⚠️ 本项目数据差异（与教程）
#   教程用 CEPH trio：NA12878 待测 + NA12891/92 作参考正常样本。
#   本项目**只有 NA12878**，故 reference 自动降级为 **flat reference**
#   （CNVkit 官方支持的单样本模式：所有 bin 期望 log2=0，仅做 GC 校正）。
#   代价：无多样本平均去噪；实测 log2 仍以 0 为中心（p05=-0.078 / p95=+0.078），
#   能检出具大 CNV，但难以区分系统偏差与真实小 CNV。
#   若日后补齐 NA12891/92 BAM，脚本自动切回教程的多样本参考模式。
#
# 用法：
#   bash scripts/10.run_cnvkit.sh              # 全流程
#   bash scripts/10.run_cnvkit.sh coverage     # 只跑某一段
#   （段名：bins|coverage|reference|fix|segment|call|viz）
#
# 环境变量：
#   BIN_SIZE=50000                  bin 大小，默认 50 kb
#   SEG_METHOD=cbs|hmm|haar|none    分割算法，默认 cbs（需 R+DNAcopy）

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

OUT="$ROOT/10.cnv"
LOG="$WS/logs/10.run_cnvkit.log"
mkdir -p "$OUT/bins" "$OUT/ref" "$OUT/results" "$WS/logs"

# ── 输入 ──
REF="$WS/rawData/ref/Homo_sapiens_assembly38.fasta"
CHROM=chr20
BAM="$ROOT/05.bqsr/${SAMPLE_SHORT}_${CHROM}_BQSR.bam"
# 参考正常样本（教程用 NA12891/92；本项目若补齐则自动启用）
NORMAL_BAMS=("$ROOT/05.bqsr/NA12891_chr20_BQSR.bam" "$ROOT/05.bqsr/NA12892_chr20_BQSR.bam")

BIN_SIZE="${BIN_SIZE:-50000}"
SEG_METHOD="${SEG_METHOD:-cbs}"
ONLY="${1:-all}"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
want() { [ "$ONLY" = "all" ] || [ "$ONLY" = "$1" ]; }

echo "=============================================="
echo " Session 17: CNVkit copy-number calling"
echo "=============================================="
echo " Sample     : $SAMPLE ($CHROM)"
echo " BAM        : $BAM"
echo " Bin size   : $BIN_SIZE"
echo " Seg method : $SEG_METHOD"
echo " Run        : $ONLY"
echo

# ── Pre-flight ──
log "Pre-flight checks..."
command -v cnvkit.py >/dev/null || {
  echo "ERROR: cnvkit.py not found. Run: bash scripts/10.setup_cnv_env.sh" | tee -a "$LOG"; exit 1; }
for f in "$BAM" "$BAM.bai" "$REF" "$REF.fai"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" | tee -a "$LOG"; exit 1; }
done
log "cnvkit: $(cnvkit.py version 2>&1 | head -1)"

# 检测可用的正常样本（决定 reference 模式）
USABLE_NORMALS=()
for nb in "${NORMAL_BAMS[@]}"; do
  [ -f "$nb" ] && USABLE_NORMALS+=("$nb")
done
if [ ${#USABLE_NORMALS[@]} -gt 0 ]; then
  log "参考模式：multiple normals (${#USABLE_NORMALS[@]} 个 —— $(basename -a "${USABLE_NORMALS[@]}" | tr '\n' ' '))"
else
  log "参考模式：FLAT reference（本项目无 NA12891/92 —— 单样本模式，见脚本头注释）"
fi
echo

# ── Step 0: bins ──
BINS_BED="$OUT/bins/${CHROM}_${BIN_SIZE}bp.bed"
if want bins || want coverage; then
log "[0/7] bins — 生成 ${CHROM} 的 ${BIN_SIZE} bp bin 坐标..."
if [ ! -f "$BINS_BED" ]; then
python3 - "$REF.fai" "$CHROM" "$BIN_SIZE" "$BINS_BED" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys
fai_p, chrom, binsize, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
length = None
with open(fai_p) as f:
    for line in f:
        p = line.split('\t')
        if p[0] == chrom:
            length = int(p[1]); break
if length is None:
    print(f'  ERROR: {chrom} not found in {fai_p}'); sys.exit(1)
n = 0
with open(out, 'w') as f:
    for s in range(0, length, binsize):
        e = min(s + binsize, length)
        if e - s < binsize // 2:   # 末尾残片丢弃
            continue
        f.write(f'{chrom}\t{s}\t{e}\t-\n')
        n += 1
print(f'  {chrom} 长度 {length:,} bp → {n} 个 {binsize:,} bp bins')
PYEOF
else
  log "  SKIP bins BED 已存在"
fi
echo
fi

# ── Step 1: coverage ──
if want coverage; then
log "[1/7] coverage — 计算 bin 覆盖深度（约 2–3 分钟）..."
for b in "$BAM" "${USABLE_NORMALS[@]:-}"; do
  [ -z "$b" ] && continue
  s=$(basename "$b" | sed "s/_${CHROM}_BQSR.bam//")
  cnn="$OUT/bins/${s}_${CHROM}.cnn"
  if [ -f "$cnn" ]; then log "  SKIP ${s}（已存在）"; continue; fi
  log "  ${s} ..."
  cnvkit.py coverage "$b" "$BINS_BED" -f "$REF" -o "$cnn" >>"$LOG" 2>&1 \
    || log "  FAIL ${s} coverage"
done
ls -lh "$OUT/bins/" | tail -n +2 | tee -a "$LOG"
echo
fi

# ── Step 2: reference ──
if want reference; then
log "[2/7] reference — 构建参考..."
if [ ${#USABLE_NORMALS[@]} -gt 0 ]; then
  REF_CNN="$OUT/ref/CEPH_normal_reference.cnn"
  if [ ! -f "$REF_CNN" ]; then
    NORM_CNNS=()
    for b in "${USABLE_NORMALS[@]}"; do
      s=$(basename "$b" | sed "s/_${CHROM}_BQSR.bam//")
      NORM_CNNS+=("$OUT/bins/${s}_${CHROM}.cnn")
    done
    cnvkit.py reference "${NORM_CNNS[@]}" -f "$REF" --no-edge -o "$REF_CNN" >>"$LOG" 2>&1 \
      && log "  多样本参考完成（${#NORM_CNNS[@]} 个 normal）" || log "  FAIL reference"
  else
    log "  SKIP 多样本参考已存在"
  fi
else
  # 平坦参考（单样本模式）：必须用 -t 指定 bins（0.9.14 行为）
  REF_CNN="$OUT/ref/flat_reference.cnn"
  if [ ! -f "$REF_CNN" ]; then
    cnvkit.py reference -f "$REF" -t "$BINS_BED" --no-edge -o "$REF_CNN" >>"$LOG" 2>&1 \
      && log "  平坦参考完成（flat reference，无 normal 样本）" || log "  FAIL reference"
  else
    log "  SKIP 平坦参考已存在"
  fi
fi
wc -l "$REF_CNN" 2>/dev/null | tee -a "$LOG"
echo
fi

# 统一 REF_CNN（供后续段使用）
if [ ${#USABLE_NORMALS[@]} -gt 0 ]; then
  REF_CNN="$OUT/ref/CEPH_normal_reference.cnn"
else
  REF_CNN="$OUT/ref/flat_reference.cnn"
fi

# ── Step 3: fix ──
if want fix; then
log "[3/7] fix — 归一化 + GC 校正..."
BINS_CNN="$OUT/bins/${SAMPLE_SHORT}_${CHROM}.cnn"
CNR="$OUT/results/${SAMPLE_SHORT}.cnr"
if [ ! -f "$CNR" ]; then
  # 0.9.14：antitarget 直接省略（WGS 无 antitarget bins），参考用 -r 传
  cnvkit.py fix "$BINS_CNN" -r "$REF_CNN" -i "$SAMPLE" --no-edge -o "$CNR" >>"$LOG" 2>&1 \
    && log "  fix 完成" || log "  FAIL fix"
else
  log "  SKIP .cnr 已存在"
fi
wc -l "$CNR" 2>/dev/null | tee -a "$LOG"
echo
fi

# ── Step 4: segment ──
if want segment; then
log "[4/7] segment — $SEG_METHOD 分割..."
CNR="$OUT/results/${SAMPLE_SHORT}.cnr"
CNS="$OUT/results/${SAMPLE_SHORT}.cns"
if [ ! -f "$CNS" ]; then
  if [ "$SEG_METHOD" = "cbs" ]; then
    cnvkit.py segment "$CNR" -m cbs -t 0.0001 -o "$CNS" >>"$LOG" 2>&1; rc=$?
    if [ $rc -ne 0 ]; then
      log "  CBS 失败（R/DNAcopy 缺失？）→ 自动降级为 hmm"
      SEG_METHOD=hmm
      cnvkit.py segment "$CNR" -m hmm -o "$CNS" >>"$LOG" 2>&1; rc=$?
    fi
  else
    cnvkit.py segment "$CNR" -m "$SEG_METHOD" -o "$CNS" >>"$LOG" 2>&1; rc=$?
  fi
  [ $rc -eq 0 ] && log "  分割完成（method=$SEG_METHOD）" || log "  FAIL segment"
else
  log "  SKIP .cns 已存在"
fi
log "  分割段数：$(( $(wc -l < "$CNS" 2>/dev/null || echo 1) - 1 ))"
echo
fi

# ── Step 5: call ──
if want call; then
log "[5/7] call — 赋予拷贝数状态..."
CNS="$OUT/results/${SAMPLE_SHORT}.cns"
CALLED="$OUT/results/${SAMPLE_SHORT}_called.cns"
if [ ! -f "$CALLED" ]; then
  cnvkit.py call "$CNS" --ploidy 2 --purity 1.0 -o "$CALLED" >>"$LOG" 2>&1 \
    && log "  call 完成" || log "  FAIL call"
else
  log "  SKIP _called.cns 已存在"
fi
# 注意：0.9.14 输出 cn 数值列（非 DEL/DUP 文字）——DEL=cn<2，DUP=cn>2
awk 'NR>1{c=$6; s=(c<2?"DEL":(c>2?"DUP":"neutral")); print "  " $1 ":" $2 "-" $3 "\tcn=" c "\t" s "\tlog2=" $5 "\tprobes=" $8}' \
  "$CALLED" 2>/dev/null | tee -a "$LOG"
echo
fi

# ── Step 6-7: visualization ──
if want viz; then
CNR="$OUT/results/${SAMPLE_SHORT}.cnr"
CALLED="$OUT/results/${SAMPLE_SHORT}_called.cns"

log "[6/7] scatter — 单染色体散点图..."
[ -f "$OUT/results/${SAMPLE_SHORT}_${CHROM}_scatter.png" ] || \
  cnvkit.py scatter "$CNR" -s "$CALLED" -c "$CHROM" \
    -o "$OUT/results/${SAMPLE_SHORT}_${CHROM}_scatter.png" >>"$LOG" 2>&1

log "[7/7] diagram — 全基因组示意图..."
[ -f "$OUT/results/${SAMPLE_SHORT}_diagram.pdf" ] || \
  cnvkit.py diagram "$CNR" -s "$CALLED" \
    -o "$OUT/results/${SAMPLE_SHORT}_diagram.pdf" >>"$LOG" 2>&1

# 自定义 log2 剖面图（Python + matplotlib，标注 ±1 / +0.58 参考线）
PNG="$OUT/results/${SAMPLE_SHORT}_${CHROM}_log2_profile.png"
if [ ! -f "$PNG" ]; then
log "  自定义 log2 剖面图..."
SML="$SAMPLE_LABEL" python3 - "$CNR" "$CALLED" "$PNG" <<'PYEOF' 2>&1 | tee -a "$LOG"
import sys, pandas as pd, matplotlib
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt

cnr_p, cns_p, out_png = sys.argv[1], sys.argv[2], sys.argv[3]
cnr = pd.read_csv(cnr_p, sep='\t', comment='#')
cns = pd.read_csv(cns_p, sep='\t', comment='#')

fig, ax = plt.subplots(figsize=(14, 5))
ax.scatter(cnr['start']/1e6, cnr['log2'], s=1.5, color='#888888', alpha=0.5, label='bins')
for _, s in cns.iterrows():
    ax.hlines(s['log2'], s['start']/1e6, s['end']/1e6, colors='#C0392B', lw=1.8)
ax.axhline(0,    color='#0E7C5B', lw=0.8, ls='--', label='log2=0 (neutral)')
ax.axhline(-1.0, color='#1F6FB2', lw=0.8, ls=':',  label='log2=-1 (het DEL)')
ax.axhline(0.58, color='#B2793A', lw=0.8, ls=':',  label='log2=+0.58 (het DUP)')
ax.set_xlabel('chr20 position (Mb)', fontsize=11)
ax.set_ylabel('log2 ratio (sample/reference)', fontsize=11)
ax.set_title(f'{os.environ.get("SML", "NA12878 chr20")} CNV Profile - CNVkit', fontsize=12, fontweight='bold')
ax.set_ylim(-2, 2)
ax.legend(fontsize=8, loc='upper right', ncol=2)
plt.tight_layout()
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'  写出 {out_png}  ({len(cnr)} bins, {len(cns)} segments)')
PYEOF
fi

ls -lh "$OUT/results/" | tail -n +2 | tee -a "$LOG"
echo
fi

# ── 小结 ──
log "=== 产出 ==="
log "  bins/    : $(ls "$OUT/bins" 2>/dev/null | wc -l) 个文件"
log "  ref/     : $(ls "$OUT/ref" 2>/dev/null | wc -l) 个文件"
log "  results/ : $(ls "$OUT/results" 2>/dev/null | wc -l) 个文件"
echo
log "Next: bash scripts/10.run_cnv_annotation.sh   （基因注释 + 过滤）"
log "注：GATK GermlineCNVCaller 需 >=30 样本（本项目 1 个）——教程定位为 cohort 演示，本项目不跑。"
