#!/usr/bin/env bash
# Session 18（第 11 步）数据下载：gnomAD-SV v4.1（GRCh38）群体结构变异数据库
#
# 为什么不用教程的 v2.1：
#   教程给的是 gnomad_v2.1_sv.sites.vcf.gz —— 实测 contig 长度 249,250,621 = **GRCh37**，
#   与本项目 GRCh38 数据**不同基因组版本**，按位置比对无效。
#   → 改用 **v4.1**（hg38，chr 前缀），与项目数据同版本，无需 liftOver / 染色体改名。
#
# 为什么用「分块下载」而不是 wget -c：
#   实测 Google Storage 大文件直连**速度剧烈波动**（11.7 MB/s ↔ 300 KB/s），
#   wget/curl 会中途断开；而 `wget -c` 续传在本环境会**错位**（文件超出真实大小 → 静默损坏）。
#   → 采用**显式 Range 分块**（每块独立校验）+ 拼接后 gzip 校验，可确定性重试。
#
# 用法：bash scripts/11.download_sv_db.sh

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
DEST="$WS/rawData/databases/sv"
LOG="$WS/logs/11.download_sv_db.log"
mkdir -p "$DEST" "$WS/logs"

BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv"
FN="gnomad.v4.1.sv.sites.vcf.gz"
CHUNK=50000000          # 50 MB / 块

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
fetch_head() { env -u http_proxy -u https_proxy curl -sS -m 30 -I "$1" 2>/dev/null; }

echo "=============================================="
echo " Session 18 data: gnomAD-SV v4.1 (GRCh38)"
echo "=============================================="
echo

cd "$DEST"

# ── 1) 主 VCF（分块下载）──
TOTAL=$(fetch_head "$BASE/$FN" | awk 'BEGIN{IGNORECASE=1} /^content-length/{print $2}' | tr -d '\r')
[ -n "$TOTAL" ] || { log "✗ 无法取到远端大小（网络问题？）"; exit 1; }
log "目标 $FN：$TOTAL B"

if [ -f "$FN" ] && [ "$(stat -c%s "$FN")" = "$TOTAL" ]; then
  log "SKIP 已完整（$TOTAL B）"
else
  n=$(( (TOTAL + CHUNK - 1) / CHUNK ))
  log "分块下载：$n 块 × $CHUNK B"
  mkdir -p .chunks
  fail=0
  for i in $(seq 0 $((n-1))); do
    s=$(( i * CHUNK )); e=$(( s + CHUNK - 1 ))
    [ $e -ge $TOTAL ] && e=$(( TOTAL - 1 ))
    want=$(( e - s + 1 ))
    f=".chunks/part_$(printf '%03d' $i)"
    [ "$(stat -c%s "$f" 2>/dev/null || echo 0)" = "$want" ] && continue

    ok=0
    for try in 1 2 3 4 5 6; do
      env -u http_proxy -u https_proxy \
        curl -sS -m 180 --retry 2 --retry-delay 2 -r "${s}-${e}" -o "$f" "$BASE/$FN" 2>>"$LOG"
      if [ "$(stat -c%s "$f" 2>/dev/null || echo 0)" = "$want" ]; then ok=1; break; fi
      log "  块 $i 重试 $try（$(stat -c%s "$f" 2>/dev/null || echo 0)/$want）"
    done
    [ $ok = 0 ] && { log "  ✗ 块 $i 最终失败"; fail=1; }
    [ $(( i % 5 )) = 0 ] && log "  进度 $i/$n"
  done
  [ $fail = 0 ] || { log "✗ 有块失败，退出（.chunks 保留以便重跑）"; exit 1; }

  log "拼接..."
  cat .chunks/part_* > "$FN"
  sz=$(stat -c%s "$FN")
  if [ "$sz" != "$TOTAL" ]; then
    log "✗ 拼接后大小不符（$sz / $TOTAL）——删除以防误用"
    rm -f "$FN"; exit 1
  fi
  log "✓ 大小精确匹配（$sz B）"
  rm -rf .chunks
fi

# ── 2) 索引（小文件，直接下）──
if [ ! -f "$FN.tbi" ]; then
  log "下载 $FN.tbi..."
  env -u http_proxy -u https_proxy \
    curl -sS -m 120 --retry 3 -o "$FN.tbi" "$BASE/$FN.tbi" 2>>"$LOG"
fi
WANT_TBI=$(fetch_head "$BASE/$FN.tbi" | awk 'BEGIN{IGNORECASE=1} /^content-length/{print $2}' | tr -d '\r')
have_tbi=$(stat -c%s "$FN.tbi" 2>/dev/null || echo 0)
[ -n "$WANT_TBI" ] && [ "$have_tbi" != "$WANT_TBI" ] && log "⚠ tbi 大小不符（$have_tbi / $WANT_TBI）" || log "✓ tbi OK（$have_tbi B）"

# ── 3) 校验 ──
log "gzip 完整性检查..."
if gzip -t "$FN" 2>>"$LOG"; then log "✓ gzip OK"; else log "✗ gzip 损坏"; exit 1; fi

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH
asm=$(bcftools view -h "$FN" 2>/dev/null | grep -m1 "ID=chr1," | grep -oE "assembly=[a-z0-9]+")
log "校验：chr1 版本=$asm  chr20 记录=$(bcftools view -H -r chr20 "$FN" 2>/dev/null | wc -l)"
echo
log "=== 产出 ==="
ls -lh "$DEST" | grep "v4.1" | awk '{printf "  %-42s %s\n", $9, $5}' | tee -a "$LOG"
echo
log "Done."