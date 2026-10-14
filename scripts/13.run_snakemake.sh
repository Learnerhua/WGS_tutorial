#!/usr/bin/env bash
# 13.run_snakemake.sh — Session 20 启动脚本
#
# 职责：
#   1. 校验 config.yaml 里所有路径都存在
#   2. 默认 --dry-run（**实际不跑**，只生成执行计划）
#   3. 显示 DAG 图、详细命令清单
#   4. 真跑用 --run 参数
#
# 用法：
#   bash scripts/13.run_snakemake.sh                  # 默认 dry-run
#   bash scripts/13.run_snakemake.sh --run            # 真跑
#   bash scripts/13.run_snakemake.sh --dag            # 生成 DAG 图
#   bash scripts/13.run_snakemake.sh --rerun          # --rerun-incomplete
#   bash scripts/13.run_snakemake.sh --until 4.7      # 跑到 4.7 步就停
#   bash scripts/13.run_snakemake.sh --clean          # 清理产物

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
cd "$WS" || { echo "❌ 工作目录不存在: $WS"; exit 1; }

# 本机 conda 26.5.2 的 activate 不可靠 → 直接把 WGS 环境加进 PATH（绝对路径原则）
WGS_BIN=/oldhome/ouyjh/miniforge3/envs/WGS/bin
export PATH="$WGS_BIN:$PATH"
export SNAKEMAKE_OUTPUT_CACHE=/dev/null  # 不启用缓存

SNAKEFILE=scripts/13.snakefile
CONFIG=scripts/13.config.yaml
LOG_DIR=logs
mkdir -p "$LOG_DIR"

# ─────────────────────────────────────────────────
# 颜色
# ─────────────────────────────────────────────────
RED='\033[0;31m'; GRN='\033[0;32m'; YEL='\033[1;33m'; BLU='\033[0;34m'; NC='\033[0m'

# ─────────────────────────────────────────────────
# 0. snakemake 在 PATH 里？
# ─────────────────────────────────────────────────
if ! command -v snakemake >/dev/null 2>&1; then
  echo -e "${RED}❌ snakemake 不在 PATH${NC}"
  echo "  激活 WGS 环境: source /oldhome/ouyjh/miniforge3/etc/profile.d/conda.sh"
  echo "  conda activate WGS"
  exit 1
fi
SNK_VER=$(snakemake --version 2>&1)
echo -e "${GRN}✅ Snakemake ${SNK_VER}${NC}"

# ─────────────────────────────────────────────────
# 1. 路径校验
# ─────────────────────────────────────────────────
echo
echo "=========================================="
echo " 路径校验（config.yaml → 实际文件）"
echo "=========================================="

PYTHON=/oldhome/ouyjh/miniforge3/envs/WGS/bin/python3
if [ ! -f "$PYTHON" ]; then
  echo -e "${RED}❌ 找不到 Python 解释器${NC}"; exit 1
fi

if ! "$PYTHON" scripts/13.validate_paths.py; then
  echo
  echo -e "${RED}❌ 路径校验失败，请先在 config.yaml 中修正${NC}"
  exit 1
fi
echo -e "${GRN}✅ 全部路径就绪${NC}"

# ─────────────────────────────────────────────────
# 2. Snakemake 运行
# ─────────────────────────────────────────────────
CORES=${CORES:-32}

case "${1:-}" in
  --dag)
    echo
    echo "=========================================="
    echo " 生成 DAG 图"
    echo "=========================================="
    # dot 优先用 WGS 环境的（conda 装）
    DOT=/oldhome/ouyjh/miniforge3/envs/WGS/bin/dot
    [ -x "$DOT" ] || DOT=$(command -v dot || true)
    if [ -z "$DOT" ]; then
      echo -e "${RED}❌ 找不到 dot（graphviz）。安装：conda install -n WGS -y graphviz${NC}"
      exit 1
    fi
    mkdir -p 13.snakemake
    # 规则级 DAG（清晰，推荐）
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" --rulegraph 2>/dev/null \
      | "$DOT" -Tpng -Gdpi=150 > 13.snakemake/rulegraph.png
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" --rulegraph 2>/dev/null \
      | "$DOT" -Tpdf > 13.snakemake/rulegraph.pdf
    # 文件级 DAG（完整依赖）
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" --dag 2>/dev/null \
      | "$DOT" -Tpng -Gdpi=120 > 13.snakemake/dag.png
    echo "  → 13.snakemake/rulegraph.png  (规则级，推荐)"
    echo "  → 13.snakemake/rulegraph.pdf  (矢量)"
    echo "  → 13.snakemake/dag.png        (文件级)"
    ;;

  --clean)
    echo
    echo -e "${YEL}⚠ 清理所有产物（慎用！）${NC}"
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" \
              --cores 1 clean
    ;;

  --rerun)
    shift
    echo
    echo "=========================================="
    echo " 增量重跑（--rerun-incomplete）"
    echo "=========================================="
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" \
              --cores "$CORES" --rerun-incomplete "$@"
    ;;

  --until)
    shift
    target="$1"; shift
    echo
    echo "=========================================="
    echo " 跑到 target = $target 后停"
    echo "=========================================="
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" \
              --cores "$CORES" --until "$target" "$@"
    ;;

  --run)
    shift
    echo
    echo "=========================================="
    echo " 全流程执行（实际跑）"
    echo "=========================================="
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" \
              --cores "$CORES" "$@"
    ;;

  --help|-h)
    cat <<EOF
用法：$0 [选项]
  (无参数)       路径校验 + dry-run，显示执行计划
  --run          实际跑全流程
  --dag          生成 DAG 图（需 graphviz 的 dot）
  --rerun        重跑（--rerun-incomplete）
  --until X      跑到 X 步骤后停
  --clean        清理所有产物
  --help         帮助

环境变量：
  CORES=32  设定 Snakemake 使用的核心数
EOF
    ;;

  *)
    echo
    echo "=========================================="
    echo " Dry-run（不实际执行）"
    echo "=========================================="
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG" \
              --cores "$CORES" -n -p
    echo
    echo "=========================================="
    echo " 小结"
    echo "=========================================="
    echo "  • 跑全流程:  bash scripts/13.run_snakemake.sh --run"
    echo "  • 跑到某步:  bash scripts/13.run_snakemake.sh --until 4.7"
    echo "  • 看依赖图:  bash scripts/13.run_snakemake.sh --dag"
    echo "  • 清理:      bash scripts/13.run_snakemake.sh --clean"
    ;;
esac
