#!/usr/bin/env bash
# Session 14 前置步骤：下载 VEP cache（GRCh38, Ensembl 105）
#
# 背景：
#   - 官方教程用 VEP v113 cache；我们的 conda vep_env 装的是 VEP 105
#     所以必须用 105 的 cache（VEP 与 cache 版本必须严格 1:1 匹配）
#   - vep_install 工具在本环境报 "DBI module not found"（其内部 @INC 有问题），
#     因此本脚本绕过 vep_install，直接 wget 官方 tar.gz 再解压
#   - EBI (ftp.ensembl.org) 单连接限速 ~400 KB/s，且 SSL 约 30 分钟会断一次，
#     所以用 wget -c + 循环重试（已验证 -c -O 会发 Range 请求正确续传）
#
# 用法：
#   bash scripts/download_vep_cache.sh                  # 直接下载（可反复运行续传）
#   bash scripts/download_vep_cache.sh /path/to.tar.gz  # 使用已手动下载好的 tar.gz
#
# 说明：浏览器/下载器下载更快时，把下好的 tar.gz 路径作为参数传入即可。

set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
CACHE_DIR=$WS/rawData/vep_cache
TARBALL=/tmp/homo_sapiens_vep_105_GRCh38.tar.gz
URL="https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_vep_105_GRCh38.tar.gz"

# 官方文件精确字节数（用于完整性校验）
EXPECTED_SIZE=15428379257

# 代理（VPN）
export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999

# 最大重试次数（应对 EBI 定期断连）
MAX_RETRY=20

mkdir -p "$CACHE_DIR"

echo "=== VEP cache download (GRCh38, Ensembl 105) ==="
echo "Target dir:   $CACHE_DIR"
echo "Cache version: 105 (must match VEP 105)"
echo "Expected size: $EXPECTED_SIZE bytes (14.37 GiB)"
echo

# ── 已有 cache 则直接跳过（断点重入） ──────────────────────────────
if [ -f "$CACHE_DIR/homo_sapiens/105_GRCh38/info.txt" ]; then
  echo "[$(date)] Cache already installed, nothing to do."
  du -sh "$CACHE_DIR/homo_sapiens" 2>/dev/null
  exit 0
fi

# ── 情况 A：调用方传入了一个已下载好的 tar.gz ──────────────────────
if [ "$#" -ge 1 ] && [ -n "${1:-}" ]; then
  SRC="$1"
  if [ ! -f "$SRC" ]; then
    echo "ERROR: provided file not found: $SRC"
    exit 1
  fi
  echo "Using provided tarball: $SRC"
  cp -f "$SRC" "$TARBALL"
fi

# ── 情况 B：从 EBI 下载 / 续传 ─────────────────────────────────────
CUR_SIZE=0
[ -f "$TARBALL" ] && CUR_SIZE=$(stat -c%s "$TARBALL")
if [ "$CUR_SIZE" -ge "$EXPECTED_SIZE" ]; then
  echo "Existing tarball already complete ($CUR_SIZE bytes), skip download."
else
  echo "[$(date)] Download resume point: $CUR_SIZE bytes"
  echo

  attempt=0
  while [ "$attempt" -lt "$MAX_RETRY" ]; do
    attempt=$((attempt + 1))
    echo "[$(date)] Download attempt $attempt/$MAX_RETRY ..."

    # wget -c 断点续传；失败不退出（|| true），由外层循环判断
    wget -c --tries=3 --waitretry=10 --timeout=60 --read-timeout=120 \
      "$URL" -O "$TARBALL" || true

    CUR_SIZE=$(stat -c%s "$TARBALL" 2>/dev/null || echo 0)
    echo "[$(date)] Got $CUR_SIZE / $EXPECTED_SIZE bytes"

    if [ "$CUR_SIZE" -ge "$EXPECTED_SIZE" ]; then
      echo "[$(date)] Download complete."
      break
    fi

    sleep 10
  done

  if [ "$CUR_SIZE" -lt "$EXPECTED_SIZE" ]; then
    echo "ERROR: download incomplete after $MAX_RETRY attempts ($CUR_SIZE / $EXPECTED_SIZE)."
    echo "Re-run this script to resume, or pass a manually downloaded tarball:"
    echo "  bash scripts/download_vep_cache.sh /path/to/homo_sapiens_vep_105_GRCh38.tar.gz"
    exit 1
  fi
fi

# ── 完整性校验（字节数） ──────────────────────────────────────────
FINAL_SIZE=$(stat -c%s "$TARBALL")
if [ "$FINAL_SIZE" -ne "$EXPECTED_SIZE" ]; then
  echo "WARNING: size mismatch — got $FINAL_SIZE, expected $EXPECTED_SIZE."
  echo "File may be corrupt. Remove it and re-download."
  exit 1
fi
echo "Size check OK: $FINAL_SIZE bytes"

# ── 解压 ──────────────────────────────────────────────────────────
echo "[$(date)] Extracting to $CACHE_DIR/homo_sapiens ..."
mkdir -p "$CACHE_DIR/homo_sapiens"
tar -xzf "$TARBALL" -C "$CACHE_DIR/homo_sapiens/"

# ── 解压后校验 ────────────────────────────────────────────────────
if [ ! -f "$CACHE_DIR/homo_sapiens/105_GRCh38/info.txt" ]; then
  echo "ERROR: extraction finished but info.txt not found — tarball may be corrupt."
  exit 1
fi

echo "[$(date)] Cache installed OK:"
ls -lh "$CACHE_DIR/homo_sapiens/105_GRCh38/info.txt"
du -sh "$CACHE_DIR/homo_sapiens"

# 解压成功后才删 tarball（失败时保留以便续传/重试）
rm -f "$TARBALL"
echo
echo "Done. Cache ready at: $CACHE_DIR/homo_sapiens/105_GRCh38/"
echo "Next step: run scripts/08.run_vep_annotation.sh"