#!/usr/bin/env snakemake
# -*- mode: python -*-
# 13.snakefile — WGS Snakemake 流水线（Session 20）
#
# ★ 沙箱设计 ★
#   所有流水线产物写入 config["out_root"]（默认 13.snakemake/），绝不覆盖原有 01–12 步结果。
#   • 只读输入：rawData/（FASTQ 入口 / 参考 / 资源 / 数据库 / VEP cache）
#   • 写入沙箱：{out_root}/01.qc … {out_root}/12.figures
#   • 被包裹的步骤脚本通过 WGS_OUT_ROOT 环境变量重定向输出
#   • 命令行可用 --out DIR 或 WGS_OUT_ROOT 环境变量覆盖
#
# 用法：
#   bash scripts/13.run_snakemake.sh              # 路径校验 + dry-run
#   bash scripts/13.run_snakemake.sh --run        # 实际执行（写沙箱）
#   bash scripts/13.run_snakemake.sh --dag        # 生成 DAG 图

import os

configfile: "scripts/13.config.yaml"

# ─────────────────────────────────────────────────
# config 简写
# ─────────────────────────────────────────────────
WS      = config["project_root"]
O       = config["out_root"]          # ★ 产物沙箱根

SAMPLE  = config["sample"]
CHROM   = config["chrom"]
REF     = os.path.join(WS, config["ref"])            # 只读
R1      = os.path.join(WS, config["fastq_r1"])       # 只读入口
R2      = os.path.join(WS, config["fastq_r2"])       # 只读入口

BWA       = config["bwa"]
SAMTOOLS  = config["samtools"]
GATK      = config["gatk"]
BCFTOOLS  = config["bcftools"]
DELLY     = config["delly"]
SURVIVOR  = config["survivor"]

DNASEQ    = os.path.join(WS, config["mills_indels"])   # 只读
KNOWN_IND = os.path.join(WS, config["known_indels"])   # 只读
EXCL      = os.path.join(WS, config["excl_regions"])   # 只读

SMP  = config["sample"]            # 主前缀（01–09 产物）
SMPS = config["sample_short"]      # 短前缀（10–12 产物）
SML  = config["sample_label"]      # 展示文本（图标题）

MANTA_DIR = config["manta_dir"]
MANTA_BIN = os.path.join(MANTA_DIR, "bin")


O_ABS = os.path.abspath(O)   # 绝对路径（供被包裹脚本用，避免 cd 后相对路径失效）

# 传给被包裹脚本的环境变量（含样本名 → 换样本只需改 config.yaml）
EXPORT = ("set +eo pipefail; export WGS_OUT_ROOT=" + O_ABS
          + " SAMPLE=" + SMP + " SAMPLE_SHORT=" + SMPS
          + ' SAMPLE_LABEL="' + SML + '"'
          + " RG_FLOWCELL=" + config["rg_flowcell"]
          + " RG_PU=" + config["rg_pu"] + "; ")


def ffmt(v):
    """浮点格式：整数补 .0（GATK VariantFiltration 按字面量推断字段类型，必须带小数点）"""
    try:
        f = float(v)
        return f"{f:.1f}" if f == int(f) else str(f)
    except (TypeError, ValueError):
        return str(v)


def fmt(v):
    """数值格式化：整数去掉 .0（对齐原脚本的 FILTER 命名，如 -8 而非 -8.0）"""
    try:
        f = float(v)
        return str(int(f)) if f == int(f) else str(f)
    except (TypeError, ValueError):
        return str(v)


def S(*parts):
    """沙箱路径：{out_root}/<step>/..."""
    return os.path.join(O, *parts)


def WGS(*parts):
    """只读资源路径（相对项目根）"""
    return os.path.join(WS, *parts)


# ─────────────────────────────────────────────────
# 末态
# ─────────────────────────────────────────────────
FIG_NAMES = [
    f"{SMP}_variant_density",         # 长前缀
    f"{SMP}_circos",                  # 长前缀
    f"{SMPS}_F1_benchmark",           # 短前缀
    f"{SMPS}_SV_length_distribution", # 短前缀
    f"{SMPS}_annotation_summary",     # 短前缀
    f"{SMPS}_variant_qc_dashboard",   # 短前缀
]
FIG_R_NAMES = [f"{SMP}_consequence_summary"]

rule all:
    input:
        # 4.1 FastQC + MultiQC
        S(f"01.qc/fastqc/{SMP}_R1_fastqc.html"),
        S(f"01.qc/fastqc/{SMP}_R2_fastqc.html"),
        S("01.qc/multiqc/multiqc_report.html"),
        # 4.2 fastp 修剪
        S(f"02.trimmed/{SMP}_R1_trimmed.fastq.gz"),
        S(f"02.trimmed/{SMP}_R2_trimmed.fastq.gz"),
        S(f"02.trimmed/multiqc/{SMPS}_pre_post_trimming.html"),
        # 4.4 比对后 QC
        S(f"04.post_alignment_qc/{SMP}.flagstat.txt"),
        S("04.post_alignment_qc/chr20.interval_list"),
        # 4.5 BAM
        S(f"05.bqsr/{SMP}_BQSR.bam"),
        S(f"05.bqsr/{SMP}_BQSR.bam.bai"),
        # 4.6 GVCF
        S(f"06.gatk/{SMP}.g.vcf.gz"),
        S(f"06.gatk/{SMP}.g.vcf.gz.tbi"),
        # 4.7 过滤
        S(f"07.var.filtered/{SMP}_filtered.vcf.gz"),
        S(f"07.var.filtered/{SMP}_filtered.vcf.gz.tbi"),
        # 4.8 VEP
        S(f"08.annotation/{SMP}_VEP.vcf.gz"),
        S(f"08.annotation/{SMP}_VEP.vcf.gz.tbi"),
        # 4.9 数据库
        S("09.databases/gnomad_freq_tier.tsv"),
        S("09.databases/clinvar_classification.tsv"),
        S(f"09.databases/{SMP}_43candidates_clinically_annotated.tsv"),
        # 4.10 benchmark
        S(f"09.benchmark/hapy_{SMP}.summary.csv"),
        # 4.11 CNV
        S(f"10.cnv/results/{SMPS}_cnv_final.tsv"),
        # 4.12 SV
        S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz"),
        S(f"11.sv/merged/{SMPS}_sv_clinical_review.tsv"),
        # 4.13 可视化
        expand(S("12.figures/{n}.png"), n=FIG_NAMES),
        expand(S("12.figures/{n}.pdf"), n=FIG_NAMES),
        expand(S("12.figures/{n}.png"), n=FIG_R_NAMES),
        expand(S("12.figures/{n}.pdf"), n=FIG_R_NAMES),
        S(f"12.figures/{SMPS}_variant_report.html"),
        S(f"12.figures/{SMPS}_sv_breakpoint_report.html"),


# ═════════════════════════════════════════════════
# 4.1 FastQC + MultiQC（原始 FASTQ）
# ═════════════════════════════════════════════════
rule fastqc:
    input:
        r1=R1, r2=R2,
    output:
        h1=S(f"01.qc/fastqc/{SMP}_R1_fastqc.html"),
        h2=S(f"01.qc/fastqc/{SMP}_R2_fastqc.html"),
        z1=S(f"01.qc/fastqc/{SMP}_R1_fastqc.zip"),
        z2=S(f"01.qc/fastqc/{SMP}_R2_fastqc.zip"),
        mq=S("01.qc/multiqc/multiqc_report.html"),
    log: "logs/13.fastqc.log"
    threads: 4
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/01.run_fastqc.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.2 fastp 修剪
# ═════════════════════════════════════════════════
rule fastp:
    input:
        r1=R1, r2=R2,
    output:
        t1=S(f"02.trimmed/{SMP}_R1_trimmed.fastq.gz"),
        t2=S(f"02.trimmed/{SMP}_R2_trimmed.fastq.gz"),
        js=S(f"02.trimmed/fastp/{SMP}_fastp.json"),
        ht=S(f"02.trimmed/fastp/{SMP}_fastp.html"),
        mq=S(f"02.trimmed/multiqc/{SMPS}_pre_post_trimming.html"),
    log: "logs/13.fastp.log"
    threads: config["fastp_threads"]
    resources: mem_mb=8 * 1024, runtime=900
    shell:
        EXPORT + "bash scripts/02.run_fastp.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.3 BWA-MEM2 比对（输入 = fastp 修剪后 FASTQ）
# ═════════════════════════════════════════════════
rule bwa_mem2:
    """4.3 BWA-MEM2 比对（包裹原脚本：含正确 @RG + @PG 修复）"""
    input:
        r1=S(f"02.trimmed/{SMP}_R1_trimmed.fastq.gz"),
        r2=S(f"02.trimmed/{SMP}_R2_trimmed.fastq.gz"),
        ref=REF,
    output:
        bam=S(f"03.aligned/{SMP}_sorted.bam"),
        bai=S(f"03.aligned/{SMP}_sorted.bam.bai"),
    log: "logs/13.bwa_mem2.log"
    threads: config["bwa_threads"]
    resources: mem_mb=32 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/03.run_bwa_mem2.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.4 比对后 QC
# ═════════════════════════════════════════════════
rule post_alignment_qc:
    """4.4 比对后 QC（包裹原脚本）"""
    input:
        S(f"03.aligned/{SMP}_sorted.bam"),
    output:
        flagstat=S(f"04.post_alignment_qc/{SMP}.flagstat.txt"),
        ilist=S("04.post_alignment_qc/chr20.interval_list"),
    log: "logs/13.post_alignment_qc.log"
    threads: 4
    resources: mem_mb=8 * 1024, runtime=900
    shell:
        EXPORT + "bash scripts/04.run_post_alignment_qc.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.5 MarkDup + BQSR
# ═════════════════════════════════════════════════
rule markdup_bqsr:
    """4.5 MarkDup + BQSR（包裹原脚本，保证参数逐字一致）"""
    input:
        S(f"03.aligned/{SMP}_sorted.bam"),
        ref=REF,
    output:
        md_bam=S(f"05.bqsr/{SMP}_markdup_sorted.bam"),
        md_bai=S(f"05.bqsr/{SMP}_markdup_sorted.bam.bai"),
        metrics=S(f"05.bqsr/{SMP}_duplication_metrics.txt"),
        recal=S(f"05.bqsr/{SMP}_recal_data.table"),
        bam=S(f"05.bqsr/{SMP}_BQSR.bam"),
        bai=S(f"05.bqsr/{SMP}_BQSR.bam.bai"),
    log: "logs/13.markdup_bqsr.log"
    threads: config["bqsr_threads"]
    resources: mem_mb=config["bqsr_mem"] * 1024, runtime=3600
    shell:
        EXPORT + "bash scripts/05.run_markdup_bqsr.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.6 HaplotypeCaller（包裹原脚本）
# ═════════════════════════════════════════════════
rule haplotypecaller:
    input:
        bam=S(f"05.bqsr/{SMP}_BQSR.bam"),
        ref=REF,
    output:
        gvcf=S(f"06.gatk/{SMP}.g.vcf.gz"),
        idx=S(f"06.gatk/{SMP}.g.vcf.gz.tbi"),
    log: "logs/13.haplotypecaller.log"
    threads: config["hc_threads"]
    resources: mem_mb=config["hc_mem"] * 1024, runtime=2400
    shell:
        EXPORT + "bash scripts/06.run_haplotypecaller.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.7 基因分型 + 过滤 + 合并（包裹原脚本）
# ═════════════════════════════════════════════════
rule variant_filtering:
    input:
        gvcf=S(f"06.gatk/{SMP}.g.vcf.gz"),
        ref=REF,
    output:
        vcf=S(f"07.var.filtered/{SMP}_filtered.vcf.gz"),
        idx=S(f"07.var.filtered/{SMP}_filtered.vcf.gz.tbi"),
    log: "logs/13.variant_filtering.log"
    threads: 4
    resources: mem_mb=16 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/07.run_variant_filtering.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.8 VEP 注释（包裹脚本）
# ═════════════════════════════════════════════════
rule vep_annotate:
    input:
        vcf=S(f"07.var.filtered/{SMP}_filtered.vcf.gz"),
    output:
        vcf=S(f"08.annotation/{SMP}_VEP.vcf.gz"),
        idx=S(f"08.annotation/{SMP}_VEP.vcf.gz.tbi"),
    log: "logs/13.vep.log"
    threads: config["vep_fork"]
    resources: mem_mb=8 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/08.run_vep_annotation.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.9 数据库交叉
# ═════════════════════════════════════════════════
rule database_crossref:
    input:
        vcf=S(f"08.annotation/{SMP}_VEP.vcf.gz"),
    output:
        gnomad=S("09.databases/gnomad_freq_tier.tsv"),
        clinvar=S("09.databases/clinvar_classification.tsv"),
        candidates=S(f"09.databases/{SMP}_43candidates_clinically_annotated.tsv"),
    log: "logs/13.database_crossref.log"
    threads: 4
    resources: mem_mb=8 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/09.run_database_crossref.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.10 hap.py 基准
# ═════════════════════════════════════════════════
rule hapy_benchmark:
    input:
        vcf=S(f"07.var.filtered/{SMP}_filtered.vcf.gz"),
    output:
        summary=S(f"09.benchmark/hapy_{SMP}.summary.csv"),
    log: "logs/13.hapy.log"
    threads: 16
    resources: mem_mb=16 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/09.run_happy_benchmark.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.11 CNV
# ═════════════════════════════════════════════════
rule cnvkit_batch:
    input:
        bam=S(f"05.bqsr/{SMP}_BQSR.bam"),
    output:
        cnr=S(f"10.cnv/results/{SMPS}.cnr"),
        cns=S(f"10.cnv/results/{SMPS}.cns"),
        called=S(f"10.cnv/results/{SMPS}_called.cns"),
    log: "logs/13.cnvkit.log"
    threads: config["cnv_threads"]
    resources: mem_mb=8 * 1024, runtime=1200
    shell:
        EXPORT + "bash scripts/10.run_cnvkit.sh 2>> {log}"

rule cnvkit_annotate:
    input:
        called=S(f"10.cnv/results/{SMPS}_called.cns"),
    output:
        final=S(f"10.cnv/results/{SMPS}_cnv_final.tsv"),
        annotated=S(f"10.cnv/results/{SMPS}_cnv_annotated.txt"),
    log: "logs/13.cnv_annotate.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/10.run_cnv_annotation.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.12 SV
# ═════════════════════════════════════════════════
rule manta:
    input:
        bam=S(f"05.bqsr/{SMP}_BQSR.bam"),
    output:
        vcf=S(f"11.sv/manta/{SMPS}_manta_pass.vcf.gz"),
    log: "logs/13.manta.log"
    threads: config["sv_threads"]
    resources: mem_mb=8 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/11.run_sv.sh manta 2>> {log}"

rule delly:
    input:
        bam=S(f"05.bqsr/{SMP}_BQSR.bam"),
    output:
        vcf=S(f"11.sv/delly/{SMPS}_delly_pass.vcf.gz"),
    log: "logs/13.delly.log"
    threads: config["sv_threads"]
    resources: mem_mb=8 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/11.run_sv.sh delly 2>> {log}"

rule survivor_merge:
    input:
        manta=S(f"11.sv/manta/{SMPS}_manta_pass.vcf.gz"),
        delly=S(f"11.sv/delly/{SMPS}_delly_pass.vcf.gz"),
    output:
        vcf=S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz"),
        idx=S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz.tbi"),
    log: "logs/13.survivor.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/11.run_sv.sh merge 2>> {log}"

rule sv_gnomad_annotate:
    """4.12 gnomAD-SV 人群频率注释（11.run_sv.sh annotate，产出 svdb TSV）"""
    input:
        vcf=S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz"),
    output:
        tsv=S(f"11.sv/merged/{SMP}_SVs_svdb_annotated.tsv"),
    log: "logs/13.sv_gnomad_annotate.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=1800
    shell:
        EXPORT + "bash scripts/11.run_sv.sh annotate 2>> {log}"

rule sv_annotate:
    input:
        vcf=S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz"),
        svdb=S(f"11.sv/merged/{SMP}_SVs_svdb_annotated.tsv"),
    output:
        review=S(f"11.sv/merged/{SMPS}_sv_clinical_review.tsv"),
        gene=S(f"11.sv/merged/{SMPS}_sv_gene_annotated.tsv"),
        bed=S(f"11.sv/merged/{SMPS}_sv_calls.bed"),
    log: "logs/13.sv_annotate.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/11.run_sv_annotation.sh 2>> {log}"

# ═════════════════════════════════════════════════
# 4.13 可视化
# ═════════════════════════════════════════════════
rule plot_matplotlib:
    input:
        vcf=S(f"07.var.filtered/{SMP}_filtered.vcf.gz"),
        cnr=S(f"10.cnv/results/{SMPS}.cnr"),
        sv=S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz"),
        summary=S(f"09.benchmark/hapy_{SMP}.summary.csv"),
        gnomad=S("09.databases/gnomad_freq_tier.tsv"),
        clinvar=S("09.databases/clinvar_classification.tsv"),
        candidates=S(f"09.databases/{SMP}_43candidates_clinically_annotated.tsv"),
    output:
        expand(S("12.figures/{n}.png"), n=FIG_NAMES),
        expand(S("12.figures/{n}.pdf"), n=FIG_NAMES),
    log: "logs/13.matplotlib.log"
    threads: 4
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/12.run_matplotlib_plots.sh 2>> {log}"

rule plot_ggplot2:
    input:
        vcf=S(f"08.annotation/{SMP}_VEP.vcf.gz"),
    output:
        png=S(f"12.figures/{SMP}_consequence_summary.png"),
        pdf=S(f"12.figures/{SMP}_consequence_summary.pdf"),
    log: "logs/13.ggplot2.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/12.run_ggplot2_plots.sh 2>> {log}"

rule igv_reports_snv:
    input:
        vcf=S(f"07.var.filtered/{SMP}_filtered.vcf.gz"),
    output:
        S(f"12.figures/{SMPS}_variant_report.html"),
    log: "logs/13.igv_reports_snv.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/12.run_igv_reports.sh snv 2>> {log}"

rule igv_reports_sv:
    input:
        vcf=S(f"11.sv/merged/{SMP}_SVs_merged.vcf.gz"),
    output:
        S(f"12.figures/{SMPS}_sv_breakpoint_report.html"),
    log: "logs/13.igv_reports_sv.log"
    threads: 2
    resources: mem_mb=4 * 1024, runtime=600
    shell:
        EXPORT + "bash scripts/12.run_igv_reports.sh sv 2>> {log}"


rule clean:
    shell:
        "echo '清理沙箱请手动执行: rm -rf " + O + "'"
