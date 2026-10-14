#!/usr/bin/env bash
# 构建 WGS_tutorial_git（位于项目根目录内）
#
# 真实拷贝（4 类）：
#   ① 根文档  WGS_tutorial.md / DEVELOPMENT.md / SKILL.md
#   ② scripts/ 全部脚本
#   ③ 12.figures 可视化文件（png/pdf/html + screenshots）
#   ④ 额外白名单 CORE_FILES（教程内嵌图 + 核心结论表 + 核心数据）
# 占位（touch）：其余全部
#
# 排除项（不镜像）：
#   • 镜像自身（避免递归嵌套）
#   • .snakemake/（运行期状态，含机器绝对路径）
set -euo pipefail

SRC=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
MIRROR_NAME=WGS_tutorial_git
DST="$SRC/$MIRROR_NAME"

# ── 额外白名单：教程内嵌图 + 核心结论表 + 核心数据 ──────────────
CORE_FILES="
04.post_alignment_qc/picard/NA12878_chr20_insert_size_histogram_fixed.png
13.snakemake/rulegraph.png
13.snakemake/rulegraph.pdf
13.snakemake/dag.png
04.post_alignment_qc/NA12878_chr20.flagstat.txt
04.post_alignment_qc/mosdepth/NA12878_chr20.mosdepth.summary.txt
04.post_alignment_qc/picard/NA12878_chr20_insert_size_metrics.txt
05.bqsr/NA12878_chr20_duplication_metrics.txt
05.bqsr/NA12878_chr20_recal_data.table
02.trimmed/fastp/NA12878_chr20_fastp.json
09.benchmark/hapy_NA12878_chr20.summary.csv
10.cnv/results/NA12878_cnv_final.tsv
10.cnv/results/NA12878_cnv_calls.bed
11.sv/merged/NA12878_sv_clinical_review.tsv
11.sv/merged/NA12878_sv_gene_annotated.tsv
11.sv/merged/NA12878_sv_calls.bed
09.databases/NA12878_chr20_43candidates_clinically_annotated.tsv
09.databases/gnomad_freq_tier.tsv
09.databases/clinvar_classification.tsv
07.var.filtered/NA12878_chr20_filtered.vcf.gz
07.var.filtered/NA12878_chr20_filtered.vcf.gz.tbi
08.annotation/NA12878_chr20_VEP.vcf.gz
08.annotation/NA12878_chr20_VEP.vcf.gz.tbi
11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz
11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz.tbi
11.sv/manta/NA12878_manta_pass.vcf.gz
11.sv/delly/NA12878_delly_pass.vcf.gz
12.figures/vep_summary.tsv
12.figures/top50_snvs.vcf.gz
12.figures/.sv_top.vcf.gz.csi
"

is_core() {
  local f=$1
  printf '%s\n' "$CORE_FILES" | grep -qxF "$f"
}

# ★ 脱敏：镜像里的文档要去掉本机隐私（绝对路径 / 用户名 / 代理端口）★
# 源项目保留真实路径，镜像版替换为模板占位符。
desensitize() {
  python3 - "$1" <<'PYEOF'
import sys, re
p = sys.argv[1]
t = open(p).read()

# ① §2.2 代理段整体替换（去掉本机代理端口与函数名）
old = """### 2.2 网络代理

本机访问境外数据源需走 VPN 代理（`~/.zshrc` 中 `vpn()` 函数）：

```bash
export http_proxy=http://localhost:9999
export https_proxy=http://localhost:9999
export all_proxy=http://localhost:9999
```

所有 wget/curl 命令前需先设置上述环境变量（脚本内已包含）。"""
new = """### 2.2 网络代理

如果访问境外数据源受限，需配置 HTTP 代理：

- 设置 `http_proxy` / `https_proxy` / `all_proxy` 环境变量指向本地代理端口
- 在 shell rc 文件中封装一个切换函数以快速启停
- 所有 wget / curl 命令前自动 source 该脚本

**国内镜像**：本项目实测使用清华源（`https://mirrors.tuna.tsinghua.edu.cn/anaconda/...`）。"""
if old in t:
    t = t.replace(old, new, 1)

# ② 绝对路径 → 模板占位符
t = t.replace('/oldhome/ouyjh/miniforge3/envs/', '<conda_root>/envs/')
t = t.replace('/oldhome/ouyjh/miniforge3/etc/',  '<conda_root>/etc/')
t = t.replace('/oldhome/ouyjh/miniforge3',      '<conda_root>')
t = t.replace('/oldhome/ouyjh/Download',        '<下载根目录>')
t = t.replace('/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial', '<项目根>/WGS_tutorial')
t = t.replace('/Nas3Data/users/ouyjh/',         '<项目根>/')
t = t.replace('/Users/oyjh/Workspace/Development/WGS_tutorial', '<项目根>/WGS_tutorial')

# ③ 用户名
t = t.replace('ouyjh', '<你的用户名>')
t = t.replace('oyjh',  '<你的用户名>')

# ④ 代理端口（只匹配 4 位端口，不误伤 99.9999% 这类测序质量数字）
t = re.sub(r'localhost:\d{4}', 'localhost:<port>', t)

open(p, 'w').write(t)
PYEOF
}

# ★ 镜像专属文件（源分析项目里没有，重建时必须保留）★
# README 只属于 GitHub 仓库，不属于分析项目本身。
# 这里在清空镜像前备份、重建后恢复，避免重建把它们丢掉。
KEEP="README.md README_CN.md"
KEEPBAK=$(mktemp -d)
for k in $KEEP; do
  [ -s "$DST/$k" ] && cp "$DST/$k" "$KEEPBAK/$k" && echo "  [保留] $k"
done

rm -rf "$DST"
mkdir -p "$DST"

real_files=0
ph_files=0
dirs=0

cd "$SRC"

# ① 目录结构
while IFS= read -r d; do
  rel="${d#./}"
  [ "$rel" = "." ] && continue
  [ "$rel" = "$MIRROR_NAME" ] && continue
  case "$rel" in .snakemake|.snakemake/*) continue ;; esac
  case "$rel" in logs|logs/*|backups|backups/*) continue ;; esac
  mkdir -p "$DST/$rel"
  dirs=$((dirs+1))
done < <(find . -path "./$MIRROR_NAME" -prune -o -type d -print | sort)

# ② 文件
while IFS= read -r f; do
  rel="${f#./}"
  case "$rel" in "$MIRROR_NAME"/*) continue ;; esac
  case "$rel" in .snakemake/*|.snakemake) continue ;; esac
  case "$rel" in logs/*|backups/*) continue ;; esac
  case "$rel" in DEVELOPMENT.md) continue ;; esac   # 开发文档不对外发布

  real=0
  case "$rel" in
    WGS_tutorial.md|SKILL.md)                                   real=1 ;;
    scripts/*)                                                  real=1 ;;
    12.figures/*.png|12.figures/*.pdf|12.figures/*.html)        real=1 ;;
    12.figures/screenshots/*)                                   real=1 ;;
  esac
  if [ "$real" = "0" ] && is_core "$rel"; then real=1; fi

  if [ "$real" = "1" ]; then
    cp -p "$SRC/$rel" "$DST/$rel"
    # 文档需脱敏（去掉本机绝对路径 / 用户名 / 代理端口）
    case "$rel" in
      WGS_tutorial.md|SKILL.md) desensitize "$DST/$rel" ;;
    esac
    real_files=$((real_files+1))
  else
    : > "$DST/$rel"
    ph_files=$((ph_files+1))
  fi
done < <(find . -path "./$MIRROR_NAME" -prune -o -type f -print | sort)

# 恢复镜像专属文件
for k in $KEEP; do
  if [ -s "$KEEPBAK/$k" ]; then
    cp "$KEEPBAK/$k" "$DST/$k"
    echo "  [恢复] $k"
    real_files=$((real_files+1))   # 计入真实文件统计
  fi
done
rm -rf "$KEEPBAK"

echo "═══════════════════════════════════"
echo " 目标: $DST"
echo " 目录:   $dirs 个"
echo " 真实:   $real_files 个文件"
echo " 占位:   $ph_files 个文件"
echo " 总计:   $((real_files+ph_files)) 个文件"
echo "═══════════════════════════════════"
