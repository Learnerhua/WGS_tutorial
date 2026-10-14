#!/usr/bin/env bash
# Session 17（第 10 步）环境准备：CNVkit + bedtools
#
# 为什么建独立环境：
#   现有 WGS env 里的 cnvkit 0.9.8 与 pandas 2.3.3 **不兼容**
#   （pandas 2.0 移除了 Int64Index，而 cnvkit 0.9.8 的 skgenome 仍在 import 它）：
#     ImportError: cannot import name 'Int64Index' from 'pandas'
#   直接升级 WGS env 的 cnvkit 可能牵动 pandas 版本，影响其他已跑通的步骤
#   → 故建独立 env（与本项目 vep_env / py27 的隔离做法一致）
#
# 安装策略：**分批**（一次性求解易卡死，见项目经验）
#   批 1：cnvkit + bedtools            （核心）
#   批 2：bioconductor-dnacopy         （CBS 分割所需，会带入 R，较大）
#
# CBS 说明：教程用 `cnvkit.py segment --method cbs`，该算法调用 R 的 DNAcopy 包。
#   若批 2 安装失败，可改用纯 Python 的 hmm 方法（见 10.run_cnvkit.sh 的 SEG_METHOD）。
#
# 用法：bash scripts/10.setup_cnv_env.sh

set -uo pipefail

CONDA=/oldhome/ouyjh/miniforge3/bin/conda
ENV_NAME=cnv_env
LOG_DIR=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial/logs
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/10.setup_cnv_env.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

echo "=============================================="
echo " Session 17 setup: cnvkit + bedtools env"
echo "=============================================="
echo

# ── 0) 已存在则检查可用性 ──
if "$CONDA" env list 2>/dev/null | grep -qE "^${ENV_NAME}\s"; then
  log "env '$ENV_NAME' 已存在，验证工具可用性..."
  if "$CONDA" run -n "$ENV_NAME" cnvkit.py version 2>>"$LOG"; then
    log "cnvkit 可用 ✅ —— 跳过安装"
    "$CONDA" run -n "$ENV_NAME" bedtools --version 2>>"$LOG" || log "WARN bedtools 缺失"
    exit 0
  else
    log "env 存在但 cnvkit 不可用 → 重建"
  fi
fi

# ── 1) 创建 env（只装 python）──
log "[1/3] 创建 env '$ENV_NAME'（python 3.11）..."
"$CONDA" create -n "$ENV_NAME" -y python=3.11 2>&1 | tail -5 | tee -a "$LOG"

# ── 2) 批 1：cnvkit + bedtools ──
log "[2/3] 安装 cnvkit + bedtools（批 1）..."
"$CONDA" install -n "$ENV_NAME" -y \
  -c conda-forge -c bioconda \
  "cnvkit>=0.9.13" bedtools 2>&1 | tail -15 | tee -a "$LOG"

# ── 3) 批 2：CBS 依赖（R + DNAcopy）──
log "[3/3] 安装 bioconductor-dnacopy（批 2，CBS 分割用，会带入 R）..."
if "$CONDA" install -n "$ENV_NAME" -y -c conda-forge -c bioconda \
     bioconductor-dnacopy 2>&1 | tail -15 | tee -a "$LOG"; then
  log "DNAcopy 安装完成"
else
  log "WARN DNAcopy 安装失败 —— CBS 不可用，请改用 SEG_METHOD=hmm"
fi

# ── 4) 验证 ──
echo
log "=== 验证 ==="
"$CONDA" run -n "$ENV_NAME" cnvkit.py version 2>&1 | tee -a "$LOG"
"$CONDA" run -n "$ENV_NAME" bedtools --version 2>&1 | tee -a "$LOG"
"$CONDA" run -n "$ENV_NAME" python3 -c "import cnvlib; print('cnvlib', cnvlib.__version__)" 2>&1 | tee -a "$LOG"
"$CONDA" run -n "$ENV_NAME" python3 -c "import pandas; print('pandas', pandas.__version__)" 2>&1 | tee -a "$LOG"
"$CONDA" run -n "$ENV_NAME" R --no-save -q -e 'if(require("DNAcopy",quietly=TRUE)) cat("DNAcopy OK\n") else cat("DNAcopy MISSING\n")' 2>&1 | tail -2 | tee -a "$LOG"

echo
log "Done. Next: bash scripts/10.run_cnvkit.sh"
