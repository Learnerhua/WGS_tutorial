#!/usr/bin/env bash
# Session 14 可选步骤：下载 Funcotator 数据源（Broad/GATK 原生注释）
#
# ⚠️ 重要前提（本机当前不满足，故脚本会先做前置检查并明确报错）：
#   broad-public-datasets 是 **Requester Pays（请求方付费）** 桶。
#   Google 要求请求方提供「计费的 GCP 项目 ID」并完成认证，否则返回：
#     400 UserProjectMissing            （没带 userProject）
#     401 AuthenticationRequiredRequesterPays（带了但没认证）
#   因此必须先具备：
#     ① 已开通计费的 Google Cloud 项目（拿到 PROJECT_ID）
#     ② 认证凭据（gcloud auth application-default login 或服务账号 JSON）
#     ③ 该项目对 broad-public-datasets 桶有 storage.objects.get 权限
#
# 用法：
#   # 1) 先认证（二选一）
#   #    gcloud auth application-default login
#   #    export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json
#   # 2) 跑本脚本，传你的计费项目 ID
#   bash scripts/download_funcotator_datasources.sh <YOUR_GCP_PROJECT_ID>
#
# 说明：Funcotator 在官方 Session 14 里是「可选」项（数据源约 6 GB，全量可达数十 GB）。
#       本项目主线用 VEP；本脚本备好以便日后需要 GATK 原生注释时直接使用。

set -euo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
OUT_DIR=$WS/rawData/funcotator_dataSources
LOG=$WS/logs/funcotator_download.log

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

mkdir -p "$(dirname "$OUT_DIR")" "$WS/logs"

echo "[$(date)] === Funcotator data sources download ===" | tee "$LOG"
echo "Output: $OUT_DIR" | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 0) 前置检查：项目 ID + 凭据 ────────────────────────────────────
GCP_PROJECT="${1:-}"

if [ -z "$GCP_PROJECT" ]; then
  echo "ERROR: no GCP project id given." | tee -a "$LOG"
  echo "  Usage: bash scripts/download_funcotator_datasources.sh <YOUR_GCP_PROJECT_ID>" | tee -a "$LOG"
  echo "  Reason: broad-public-datasets is a Requester Pays bucket — a billing project is mandatory." | tee -a "$LOG"
  exit 1
fi
echo "GCP project (billing for requester pays): $GCP_PROJECT" | tee -a "$LOG"

# 凭据检查（ADC 文件 / 服务账号 JSON 任一存在即可）
ADC_FILE="$HOME/.config/gcloud/application_default_credentials.json"
HAS_CRED=0
if [ -f "$ADC_FILE" ]; then
  echo "Credentials: application default credentials found" | tee -a "$LOG"
  HAS_CRED=1
elif [ -n "${GOOGLE_APPLICATION_CREDENTIALS:-}" ] && [ -f "${GOOGLE_APPLICATION_CREDENTIALS}" ]; then
  echo "Credentials: service account JSON ($GOOGLE_APPLICATION_CREDENTIALS)" | tee -a "$LOG"
  HAS_CRED=1
fi

if [ "$HAS_CRED" -eq 0 ]; then
  echo "ERROR: no Google credentials found." | tee -a "$LOG"
  echo "  Authenticate first, one of:" | tee -a "$LOG"
  echo "    gcloud auth application-default login" | tee -a "$LOG"
  echo "    export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json" | tee -a "$LOG"
  echo "  (Without credentials the bucket returns 401 AuthenticationRequiredRequesterPays.)" | tee -a "$LOG"
  exit 1
fi
echo | tee -a "$LOG"

# ── 1) 断点重入：已有数据源目录则跳过 ──────────────────────────────
if [ -d "$OUT_DIR" ] && [ -n "$(ls -A "$OUT_DIR" 2>/dev/null)" ]; then
  echo "[$(date)] Data sources already present, skip download:" | tee -a "$LOG"
  du -sh "$OUT_DIR" | tee -a "$LOG"
  ls "$OUT_DIR" | head -5 | tee -a "$LOG"
  exit 0
fi

# ── 2) 下载 Germline HG38 数据源 ───────────────────────────────────
# --germline            下载 germline（种系）预打包数据源
# --hg38                使用 HG38 数据源桶（本教程参考为 GRCh38）
# --validate-integrity  下载后用 sha256 校验完整性
# --extract-after-download  下载后自动解压到同级目录
# --gcs-project-for-requester-pays  访问 requester pays 桶所需的计费项目
echo "[$(date)] Downloading germline HG38 data sources (several GB)..." | tee -a "$LOG"
echo "  This may take a long time depending on bandwidth." | tee -a "$LOG"

gatk FuncotatorDataSourceDownloader \
  --germline \
  --hg38 \
  --validate-integrity \
  --extract-after-download \
  --gcs-project-for-requester-pays "$GCP_PROJECT" \
  -O "$OUT_DIR" \
  2>>"$LOG"

echo "[$(date)] Download finished." | tee -a "$LOG"
echo | tee -a "$LOG"

# ── 3) 结果校验 ────────────────────────────────────────────────────
if [ -d "$OUT_DIR" ] && [ -n "$(ls -A "$OUT_DIR" 2>/dev/null)" ]; then
  echo "Data sources installed:" | tee -a "$LOG"
  du -sh "$OUT_DIR" | tee -a "$LOG"
  ls -lh "$OUT_DIR" | head -15 | tee -a "$LOG"
  echo | tee -a "$LOG"
  echo "Next: annotate with Funcotator (tutorial's actual command), e.g." | tee -a "$LOG"
  echo "  gatk Funcotator -R rawData/ref/Homo_sapiens_assembly38.fasta \\" | tee -a "$LOG"
  echo "    -V 07.var.filtered/NA12878_chr20_filtered.vcf.gz \\" | tee -a "$LOG"
  echo "    -O 08.annotation/NA12878_chr20_Funcotator.vcf.gz \\" | tee -a "$LOG"
  echo "    --output-file-format VCF \\" | tee -a "$LOG"
  echo "    --data-sources-path $OUT_DIR \\" | tee -a "$LOG"
  echo "    --ref-version hg38 -L chr20" | tee -a "$LOG"
  echo "  NOTE: tutorial's comment mentions --variant-type SNP_INDEL, but that flag is" | tee -a "$LOG"
  echo "        NOT part of the tutorial's actual command and GATK 4.6 rejects it" | tee -a "$LOG"
  echo "        (\"variant-type is not a recognized option\"). Funcotator annotates all" | tee -a "$LOG"
  echo "        variants present in the input VCF, so no type filter is needed." | tee -a "$LOG"
else
  echo "ERROR: download appears to have failed — directory missing or empty." | tee -a "$LOG"
  echo "  Check log: $LOG" | tee -a "$LOG"
  echo "  Common causes: missing billing project, no credentials, project lacks bucket permission." | tee -a "$LOG"
  exit 1
fi

# ── 备选路径（若将来找到非 requester-pays 的镜像）────────────────────
# Funcotator 支持用自备的 tar.gz 绕过 GCS：
#   gatk FuncotatorDataSourceDownloader \
#     --testing-override-path-for-datasources /path/to/funcotator_dataSources.tar.gz \
#     --validate-integrity --extract-after-download -O $OUT_DIR
# 若该 tar.gz 的 sha256 与官方不符，可再加 --testing-override-path-for-datasources-sha256 指向自备校验文件。