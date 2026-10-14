#!/usr/bin/env bash
# Session 17（第 10 步）数据下载：DGV 金标准 CNV 目录（用于剔除良性拷贝数多态）
#
# 背景：教程 4.12 门槛 4「不在已知常见 CNP 区」需要 DGV，但教程未给下载地址。
#       DGV 原站 dgv.tcag.ca 已下线；UCSC 提供 hg38 版轨道：
#         https://hgdownload.soe.ucsc.edu/gbdb/hg38/dgv/
#   取 dgvGold.bb（金标准：健康对照高置信 CNV，31.5 MB）
#   格式为 bigBed（.bb 二进制）→ 需 bigBedToBed 转 BED
#
# 用法：bash scripts/10.download_dgv.sh

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
CNV_BIN=/oldhome/ouyjh/miniforge3/envs/cnv_env/bin
export PATH="$CNV_BIN:/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH"

DEST="$WS/rawData/databases/dgv"
LOG="$WS/logs/10.download_dgv.log"
mkdir -p "$DEST" "$WS/logs"

BASE="https://hgdownload.soe.ucsc.edu/gbdb/hg38/dgv"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 17 data: DGV gold-standard CNVs"
echo "=============================================="
echo

# ── 0) 工具检查/安装 ──
log "检查 bigBedToBed..."
if ! command -v bigBedToBed >/dev/null; then
  log "  未安装 → conda 装 ucsc-bigbedtobed 到 cnv_env"
  /oldhome/ouyjh/miniforge3/bin/conda install -n cnv_env -y \
    -c conda-forge -c bioconda ucsc-bigbedtobed 2>&1 | tail -6 | tee -a "$LOG"
fi
command -v bigBedToBed >/dev/null && log "  bigBedToBed OK" || { echo "ERROR: bigBedToBed 不可用" | tee -a "$LOG"; exit 1; }

# ── 1) 下载 dgvGold.bb ──
BB="$DEST/dgvGold.bb"
if [ ! -f "$BB" ]; then
  log "[1/3] 下载 dgvGold.bb（31.5 MB）..."
  # 第一步：显式**绕过**代理直连（UCSC 直连可用；不绕过会被失效代理卡住）
  if ! env -u http_proxy -u https_proxy -u all_proxy \
       wget -c -q --tries=3 --timeout=60 -O "$BB" "$BASE/dgvGold.bb" 2>>"$LOG"; then
    log "  直连失败 → 走代理重试"
    http_proxy=http://localhost:9999 https_proxy=http://localhost:9999 \
      wget -c -q --tries=3 --timeout=60 -O "$BB" "$BASE/dgvGold.bb" 2>>"$LOG"
  fi
else
  log "[1/3] SKIP dgvGold.bb 已存在"
fi
log "  大小: $(stat -c%s "$BB" 2>/dev/null) B"
echo

# ── 2) bigBed → BED ──
BED_ALL="$DEST/dgvGold.hg38.bed"
if [ ! -f "$BED_ALL" ]; then
  log "[2/3] bigBedToBed 转换..."
  bigBedToBed "$BB" "$BED_ALL" 2>>"$LOG" && log "  完成：$(wc -l < "$BED_ALL") 行"
else
  log "[2/3] SKIP BED 已存在"
fi
echo

# ── 3) 提取 chr20（与 CNV 结果同染色体）──
BED20="$DEST/dgvGold.chr20.bed"
if [ ! -f "$BED20" ]; then
  log "[3/3] 提取 chr20..."
  awk -F'\t' '$1=="chr20"{print $1"\t"$2"\t"$3"\t"$4}' "$BED_ALL" > "$BED20"
  log "  chr20 条目：$(wc -l < "$BED20") 行"
else
  log "[3/3] SKIP chr20 BED 已存在"
fi
echo

log "=== 产出 ==="
for f in "$BB" "$BED_ALL" "$BED20"; do
  [ -f "$f" ] && printf "  %-40s %10s B\n" "$(basename "$f")" "$(stat -c%s "$f")" | tee -a "$LOG"
done
echo
log "Done. Next: bash scripts/10.run_cnv_annotation.sh  （含 DGV 过滤）"
