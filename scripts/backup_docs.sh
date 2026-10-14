#!/usr/bin/env bash
# 自动备份教程和开发文档
# 防止大型修改（如删错、批量替换）造成不可逆数据丢失
# 使用：bash scripts/backup_docs.sh
#
# 备份策略：
#   - 每次运行创建时间戳快照
#   - 保留最近 5 个备份
#   - 备份目录：项目根目录的 backups/

set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
BACKUP_DIR=$WS/backups
TIMESTAMP=$(date +%Y%m%d-%H%M%S)

mkdir -p "$BACKUP_DIR"

# 备份目标文件
FILES_TO_BACKUP=(
  "WGS_tutorial.md"
  "DEVELOPMENT.md"
  "scripts/01.run_fastqc.sh"
  "scripts/02.run_fastp.sh"
  "scripts/03.run_bwa_mem2.sh"
  "scripts/04.run_post_alignment_qc.sh"
  "scripts/05.run_markdup_bqsr.sh"
  "scripts/06.run_haplotypecaller.sh"
)

echo "=== 自动备份 ==="
echo "目标目录: $BACKUP_DIR"
echo "时间戳: $TIMESTAMP"
echo

# 创建快照备份
SNAPSHOT_DIR="$BACKUP_DIR/snapshot-$TIMESTAMP"
mkdir -p "$SNAPSHOT_DIR"

for file in "${FILES_TO_BACKUP[@]}"; do
  if [ -f "$WS/$file" ]; then
    mkdir -p "$SNAPSHOT_DIR/$(dirname $file)"
    cp "$WS/$file" "$SNAPSHOT_DIR/$file"
    echo "  ✓ $file ($(wc -c < "$WS/$file") bytes)"
  else
    echo "  ⚠ $file (不存在)"
  fi
done

echo
echo "=== 快照位置 ==="
echo "$SNAPSHOT_DIR"
echo

# 清理旧快照：保留最近 5 个
echo "=== 清理旧快照（保留最近 5 个） ==="
SNAPSHOTS=($(ls -d "$BACKUP_DIR"/snapshot-* 2>/dev/null | sort -r))
if [ ${#SNAPSHOTS[@]} -gt 5 ]; then
  for old in "${SNAPSHOTS[@]:5}"; do
    rm -rf "$old"
    echo "  删除: $(basename $old)"
  done
else
  echo "  无需清理（当前 ${#SNAPSHOTS[@]} 个快照）"
fi

echo
echo "=== 当前快照列表 ==="
ls -lt "$BACKUP_DIR"/snapshot-* 2>/dev/null | awk '{print "  "$NF}' | head -10

echo
echo "=== 恢复方法 ==="
echo "如需恢复某个文件："
echo "  cp $BACKUP_DIR/snapshot-<时间戳>/<文件路径> <原位置>"
echo
echo "示例（恢复教程到最新快照）："
LATEST=$(ls -dt "$BACKUP_DIR"/snapshot-* 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then
  echo "  cp $LATEST/WGS_tutorial.md $WS/WGS_tutorial.md"
fi