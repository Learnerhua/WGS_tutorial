---
name: snakemake-sandbox-pipeline
description: Wrap existing pipeline scripts in a sandboxed Snakemake DAG.
version: 1.0.0
author: <你的用户名>, Hermes Agent
license: MIT
platforms: [linux, macos]
metadata:
  hermes:
    tags: [snakemake, pipeline, sandbox, reproducibility, bioinformatics]
    related_skills: [snakemake-pipeline-debugging, bash-pipeline-scheduler]
---

# Snakemake 沙箱化包裹既有流水线

把一条**已跑通的顺序脚本流水线**（`NN.run_xxx.sh` × N）改造成 Snakemake 编排，
且**所有产物写入独立沙箱目录**，绝不覆盖原结果。

核心原则：**包裹已验证的脚本，绝不重写命令**。

## When to Use

- 已有 N 个顺序执行的流水线脚本，想获得断点续跑 / 并行 / 一条命令重建
- 改造**不能破坏**原产物（原目录要留着对比，或含手工产物无法自动重建）
- 需要证明"重跑结果与原结果一致"（可复现性要求）
- **不适用**：从零写新流程（直接写细粒度 rule 更清晰）

## Prerequisites

- Snakemake ≥ 7；`graphviz`（`dot`，生成 DAG 图用）
- 步骤脚本自身的输入输出路径**可参数化**（见 Procedure 1）
- 明确区分"只读输入"（原始数据/参考/数据库）与"可写产物"

## How to Run

```
terminal(command="bash scripts/NN.run_snakemake.sh")         # 路径校验 + dry-run
terminal(command="bash scripts/NN.run_snakemake.sh --run")   # 实际执行
terminal(command="bash scripts/NN.run_snakemake.sh --dag")   # 生成 DAG 图
terminal(command="bash scripts/NN.verify.sh full")           # 核验 + 与原结果对比
```

## Procedure

**1. 给脚本加输出根重定向（向后兼容）**
完成标准：**不设环境变量时行为与改造前完全一致**。

每个脚本在 `WS=...` 后加一行，并把 `$WS/<步骤>` 改为 `$ROOT/<步骤>`
（`$WS/rawData`、`$WS/scripts`、`$WS/logs` 等只读/日志路径**保持不动**）：

```bash
ROOT="${WGS_OUT_ROOT:-$WS}"      # 默认原地；设环境变量则重定向到沙箱
```

**2. 写 config.yaml（四组）**
完成标准：**路径与参数零硬编码在 Snakefile 中**。

| 组 | 内容 |
|---|---|
| 1 样本与参考 | 样本名、染色体、参考 FASTA、输入 FASTQ |
| 2 工具绝对路径 | 各环境可执行文件**绝对路径**（绕开 `conda activate`）|
| 3 已下载资源 | 参考数据库 / cache / BED（**假设已存在，不触网**）|
| 4 步骤参数 | 线程数、各步阈值 |

另加 `out_root: "<沙箱目录>"` 作为沙箱根。

**3. 写 Snakefile：一个 rule 包一个脚本**
完成标准：`rule all` 覆盖全部终产物，dry-run 无 `MissingInputException`。

```python
configfile: "scripts/NN.config.yaml"
O     = config["out_root"]
O_ABS = os.path.abspath(O)                 # ★ 必须绝对路径
def S(*p): return os.path.join(O, *p)

rule step5:
    input:  S("04.xxx/out.bam"), ref=REF
    output: S("05.yyy/out.vcf.gz"), S("05.yyy/out.vcf.gz.tbi")
    log:    "logs/NN.step5.log"
    threads: config["threads"]
    resources: mem_mb=32*1024, runtime=3600
    shell:
        "set +eo pipefail; export WGS_OUT_ROOT=" + O_ABS + "; "
        "bash scripts/05.run_yyy.sh 2>> {log}"
```

**rule 粒度**：一个 rule 包一个**已跑通的脚本**。刻意重写细粒度命令会引入参数漂移。

**3. 让样本名可配置（否则流程锁死单样本）**
完成标准：换样本只改 config + 输入文件，脚本不用动。

在 config 里定义（本项目实例，长短前缀因历史命名并存）：
```yaml
sample:       "NA12878_chr20"    # 主前缀（多数产物）
sample_short: "NA12878"          # 短前缀（少数产物）
sample_label: "NA12878 chr20"    # 图标题等展示文本
```

每个脚本加（**默认值 = 改造前的值 → 行为不变**）：
```bash
SAMPLE="${SAMPLE:-NA12878_chr20}"
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"
```

Snakefile 统一 export 给所有 rule：
```python
EXPORT = ("set +eo pipefail; export WGS_OUT_ROOT=" + O_ABS
          + " SAMPLE=" + SMP + " SAMPLE_SHORT=" + SMPS
          + ' SAMPLE_LABEL="' + SML + '"; ')
shell: EXPORT + "bash scripts/NN.run_xxx.sh 2>> {log}"
```

Snakefile 自身路径也用 f-string：`S(f"07.var.filtered/{SMP}_filtered.vcf.gz")`。

**关键点**：`<<'PYEOF'` 引号 heredoc **不展开** shell 变量 → 内嵌 Python 块必须
改走环境变量（`SMP="$SAMPLE" python3 - ... <<'PYEOF'` + `os.environ['SMP']`）。

**5. 写启动脚本**
完成标准：**默认 dry-run，加 `--run` 才执行**。

必须含**路径校验**（逐条校验 config 里的路径存在，缺失即退出），
并支持 `--out DIR` / `WGS_OUT_ROOT` 覆盖沙箱根，优先级 CLI > 环境变量 > config。

**6. 验证一致性**
完成标准：与原结果的**逐条内容 md5 相同**（不只是计数相同）。

顺序：dry-run 确认 DAG 完整 → 真跑 → 对比计数 / 摘要指标 / **逐条 md5** / 图片字节大小 / 泄漏检测。

## Pitfalls

| 症状 | 根因 | 修法 |
|---|---|---|
| `Specified output directory does not exist` | 沙箱根传了**相对路径**，脚本内部 `cd` 后失效 | 传 `os.path.abspath(O)` |
| 规则报 "exited with non-zero exit code" 但日志无错 | Snakemake 给 shell 块加 `-euo pipefail`；脚本内 `tee` broken pipe 返回非 0 | shell 开头加 `set +eo pipefail;` |
| **结果与原流程差百分之几** | **手写命令参数漂移** | **包裹脚本，不重写命令** |
| GATK `NumberFormatException: "55.74"` | 按表达式字面量推断类型：`QD < 2` 被当整数 | 表达式带小数点 `QD < 2.0`（FILTER 名可仍为 `QD2`）|
| `MissingOutputException`（`.bam.bai`）| GATK `--CREATE_INDEX true` 产 `X.bai` 而非 `X.bam.bai` | 补 `samtools index X.bam`，两种命名都留 |
| 某规则**静默不执行** | 产物既不在 `rule all` 也无下游依赖 → 孤岛规则不被调度 | 把产物加进 `rule all` |
| 补了新依赖却不重跑 | Snakemake 只验**产物存在**，**不检查 input** | 删掉该规则产物强制重跑 |
| 后台进程启动即死 | `pkill -f snakemake` **匹配到自己的命令行** | kill 与启动**分两次调用** |
| md5 对比"假通过" | 查询命令出错返回空 → 两边都是空串哈希 | 确保查询非空；`bcftools query` 样本字段写 `[%GT]` |
| **流程锁死单样本** | 脚本内 `SAMPLE=xxx` 硬编码、Snakefile 路径写死 | 三变量参数化（sample/sample_short/sample_label）+ EXPORT 统一 export；**引号 heredoc 内的 Python 块须用 `os.environ` 取值** |
| 沙箱泄漏写到原目录 | Snakefile 的 output 漏加沙箱前缀 | `find . -newermt <今天> -not -path "./<沙箱>/*"` 定期扫 |
| **根目录散落中间文件** | Snakemake 的 CWD = 项目根；脚本内相对路径（或早期 output 未加前缀）把中间产物写到根目录 | `find . -maxdepth 1 -type f` 扫；确认无引用后删。**根目录只应留文档** |

## Verification

1. dry-run 输出 `Nothing to be done` → 全部产物就绪
2. 路径校验全部通过（无缺失）
3. **逐条内容 md5 相同**（非仅计数相同）
4. 泄漏检测：原目录无文件被改动（文档 / 已还原文件需白名单排除）
5. 手工产物（人工截图等）**不在任何 rule 内** → 天然不受影响
6. 根目录零散落文件（`find . -maxdepth 1 -type f` 仅返回文档）
7. 换样本验证：改 config 的 sample → 用 `--config sample=TEST_chr20` dry-run，路径与环境变量应全部改名、无残留旧名
8. **边界确认**：明确记录流程**不能**做什么（本项目：锁 chr20、truth 基准绑定 HG001）——能力边界写进文档，别让下游误用
