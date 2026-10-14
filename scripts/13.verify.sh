#!/usr/bin/env bash
# 13.verify.sh — Session 20 沙箱流水线核验（完整性 + 与原目录对比 + 内容级）
#
# 用法：bash scripts/13.verify.sh [quick|full]
#   quick（默认）：产物计数 + 摘要指标对比（快，秒级）
#   full         ：再加内容级 md5 逐条对比（慢，分钟级）
export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH
cd /Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial || exit 1

MODE="${1:-quick}"
O=13.snakemake     # 沙箱
P=.                # 原目录（基准）

ok()   { printf "  %-34s %-14s %-14s %s\n" "$1" "$2" "$3" "$([ "$2" = "$3" ] && echo ✅ || echo ❌)"; }
oks()  { printf "  %-34s %-14s %-14s %s\n" "$1" "-" "$2" "$([ "$2" = same ] && echo ✅ || echo ❌)"; }
sep()  { printf '─%.0s' $(seq 1 78); echo; }

echo "════════════════════════════════════════════════════════════════════════"
echo " Session 20 沙箱核验  模式=$MODE   $(date '+%F %T')"
echo "════════════════════════════════════════════════════════════════════════"
echo

# ── 一、产物完整性 ──────────────────────────────────────────────
echo "【一】产物完整性"
sep
miss=0
while IFS= read -r f; do
  [ -z "$f" ] && continue
  if [ -f "$f" ]; then
    printf "  ✅ %-64s %s\n" "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else
    printf "  ❌ 缺失 %s\n" "$f"; miss=$((miss+1))
  fi
done <<'EOF'
13.snakemake/01.qc/fastqc/NA12878_chr20_R1_fastqc.html
13.snakemake/01.qc/fastqc/NA12878_chr20_R2_fastqc.html
13.snakemake/01.qc/multiqc/multiqc_report.html
13.snakemake/02.trimmed/NA12878_chr20_R1_trimmed.fastq.gz
13.snakemake/02.trimmed/NA12878_chr20_R2_trimmed.fastq.gz
13.snakemake/02.trimmed/multiqc/NA12878_pre_post_trimming.html
13.snakemake/03.aligned/NA12878_chr20_sorted.bam
13.snakemake/04.post_alignment_qc/NA12878_chr20.flagstat.txt
13.snakemake/04.post_alignment_qc/chr20.interval_list
13.snakemake/05.bqsr/NA12878_chr20_BQSR.bam
13.snakemake/05.bqsr/NA12878_chr20_recal_data.table
13.snakemake/06.gatk/NA12878_chr20.g.vcf.gz
13.snakemake/06.gatk/NA12878_chr20.g.vcf.gz.tbi
13.snakemake/07.var.filtered/NA12878_chr20_filtered.vcf.gz
13.snakemake/07.var.filtered/NA12878_chr20_filtered.vcf.gz.tbi
13.snakemake/08.annotation/NA12878_chr20_VEP.vcf.gz
13.snakemake/08.annotation/NA12878_chr20_VEP.vcf.gz.tbi
13.snakemake/09.databases/gnomad_freq_tier.tsv
13.snakemake/09.databases/clinvar_classification.tsv
13.snakemake/09.databases/NA12878_chr20_43candidates_clinically_annotated.tsv
13.snakemake/09.benchmark/hapy_NA12878_chr20.summary.csv
13.snakemake/10.cnv/results/NA12878_cnv_final.tsv
13.snakemake/11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz
13.snakemake/11.sv/merged/NA12878_sv_clinical_review.tsv
13.snakemake/12.figures/NA12878_chr20_variant_density.png
13.snakemake/12.figures/NA12878_chr20_variant_density.pdf
13.snakemake/12.figures/NA12878_chr20_circos.png
13.snakemake/12.figures/NA12878_chr20_circos.pdf
13.snakemake/12.figures/NA12878_F1_benchmark.png
13.snakemake/12.figures/NA12878_F1_benchmark.pdf
13.snakemake/12.figures/NA12878_SV_length_distribution.png
13.snakemake/12.figures/NA12878_SV_length_distribution.pdf
13.snakemake/12.figures/NA12878_annotation_summary.png
13.snakemake/12.figures/NA12878_annotation_summary.pdf
13.snakemake/12.figures/NA12878_variant_qc_dashboard.png
13.snakemake/12.figures/NA12878_variant_qc_dashboard.pdf
13.snakemake/12.figures/NA12878_chr20_consequence_summary.png
13.snakemake/12.figures/NA12878_chr20_consequence_summary.pdf
13.snakemake/12.figures/NA12878_variant_report.html
13.snakemake/12.figures/NA12878_sv_breakpoint_report.html
EOF
echo
echo "  缺失总数: $miss"
echo

# ── 二、与原目录指标对比 ────────────────────────────────────────
echo "【二】与原目录指标对比"
sep
printf "  %-34s %-14s %-14s %s\n" "指标" "原目录" "沙箱" "一致"

a=$(grep "mapped (" "$P"/04.post_alignment_qc/NA12878_chr20.flagstat.txt 2>/dev/null | head -1 | awk '{print $1}')
b=$(grep "mapped (" "$O"/04.post_alignment_qc/NA12878_chr20.flagstat.txt 2>/dev/null | head -1 | awk '{print $1}')
ok "比对读段数" "$a" "$b"

V1="$P/07.var.filtered/NA12878_chr20_filtered.vcf.gz"; V2="$O/07.var.filtered/NA12878_chr20_filtered.vcf.gz"
n1=$(bcftools view -H "$V1" 2>/dev/null | wc -l); n2=$(bcftools view -H "$V2" 2>/dev/null | wc -l)
ok "过滤后变异总数" "$n1" "$n2"
s1=$(bcftools view -H -v snps "$V1" 2>/dev/null | wc -l); s2=$(bcftools view -H -v snps "$V2" 2>/dev/null | wc -l)
ok "SNP 数" "$s1" "$s2"
i1=$(bcftools view -H -v indels "$V1" 2>/dev/null | wc -l); i2=$(bcftools view -H -v indels "$V2" 2>/dev/null | wc -l)
ok "INDEL 数" "$i1" "$i2"

calc_ts() { bcftools view -H -v snps "$1" 2>/dev/null | awk -F'\t' '{r=$4;a=substr($5,1,1); if((r=="A"&&a=="G")||(r=="G"&&a=="A")||(r=="C"&&a=="T")||(r=="T"&&a=="C")) t++; else v++} END{printf "%.3f",(v>0?t/v:0)}'; }
t1=$(calc_ts "$V1"); t2=$(calc_ts "$V2")
ok "Ts/Tv" "$t1" "$t2"

h1=$(bcftools query -f '[%GT]\n' "$V1" 2>/dev/null | grep -cE '^(0[/|]1)$')
h2=$(bcftools query -f '[%GT]\n' "$V2" 2>/dev/null | grep -cE '^(0[/|]1)$')
ok "杂合 het 数" "$h1" "$h2"

ok "hap.py summary" "$(md5sum "$P"/09.benchmark/hapy_NA12878_chr20.summary.csv 2>/dev/null|cut -c1-8)" "$(md5sum "$O"/09.benchmark/hapy_NA12878_chr20.summary.csv 2>/dev/null|cut -c1-8)"

for t in manta delly; do
  x=$(bcftools view -H "$P"/11.sv/$t/*pass*.vcf.gz 2>/dev/null | wc -l)
  y=$(bcftools view -H "$O"/11.sv/$t/*pass*.vcf.gz 2>/dev/null | wc -l)
  ok "$t PASS" "$x" "$y"
done
x=$(bcftools view -H "$P"/11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz 2>/dev/null | wc -l)
y=$(bcftools view -H "$O"/11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz 2>/dev/null | wc -l)
ok "SV 共识" "$x" "$y"
echo

# ── 三、内容级核验 ──────────────────────────────────────────────
if [ "$MODE" = "full" ]; then
  echo "【三】内容级核验（逐条 md5）"
  sep
  c1=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$V1" 2>/dev/null | md5sum | cut -c1-8)
  c2=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$V2" 2>/dev/null | md5sum | cut -c1-8)
  ok "变异逐条内容 md5" "$c1" "$c2"

  v1=$(bcftools query -f '%CHROM:%POS:%REF:%ALT:%INFO/CSQ\n' "$P"/08.annotation/NA12878_chr20_VEP.vcf.gz 2>/dev/null | md5sum | cut -c1-8)
  v2=$(bcftools query -f '%CHROM:%POS:%REF:%ALT:%INFO/CSQ\n' "$O"/08.annotation/NA12878_chr20_VEP.vcf.gz 2>/dev/null | md5sum | cut -c1-8)
  ok "VEP CSQ 内容 md5" "$v1" "$v2"

  sv1=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SVTYPE\n' "$P"/11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz 2>/dev/null | md5sum | cut -c1-8)
  sv2=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SVTYPE\n' "$O"/11.sv/merged/NA12878_chr20_SVs_merged.vcf.gz 2>/dev/null | md5sum | cut -c1-8)
  ok "SV 共识内容 md5" "$sv1" "$sv2"

  for f in 09.databases/gnomad_freq_tier.tsv 09.databases/clinvar_classification.tsv \
           09.databases/NA12878_chr20_43candidates_clinically_annotated.tsv \
           10.cnv/results/NA12878_cnv_final.tsv \
           11.sv/merged/NA12878_sv_clinical_review.tsv \
           09.benchmark/hapy_NA12878_chr20.summary.csv; do
    r=$(diff -q "$P/$f" "$O/$f" >/dev/null 2>&1 && echo same || echo diff)
    oks "$(basename "$f")" "$r"
  done

  # PNG：真字节级（md5）；PDF：内嵌 CreationDate 每次不同 → 去时间戳后比内容
  psame=0; ptot=0; dsame=0; dtot=0
  for n in NA12878_chr20_variant_density NA12878_chr20_circos NA12878_F1_benchmark \
           NA12878_SV_length_distribution NA12878_annotation_summary \
           NA12878_variant_qc_dashboard NA12878_chr20_consequence_summary; do
    ptot=$((ptot+1))
    [ "$(md5sum "$P/12.figures/$n.png" 2>/dev/null|cut -d' ' -f1)" = "$(md5sum "$O/12.figures/$n.png" 2>/dev/null|cut -d' ' -f1)" ] && psame=$((psame+1))
    dtot=$((dtot+1))
    a=$(strings "$P/12.figures/$n.pdf" 2>/dev/null | grep -vE "CreationDate|ModDate" | md5sum | cut -d' ' -f1)
    b=$(strings "$O/12.figures/$n.pdf" 2>/dev/null | grep -vE "CreationDate|ModDate" | md5sum | cut -d' ' -f1)
    [ -n "$a" ] && [ "$a" = "$b" ] && dsame=$((dsame+1))
  done
  printf "  %-34s %-14s %-14s %s\n" "图 PNG 字节级一致（md5）" "-" "$psame/$ptot" "$([ "$psame" = "$ptot" ] && echo ✅ || echo ❌)"
  printf "  %-34s %-14s %-14s %s\n" "图 PDF 内容一致（去时间戳）" "-" "$dsame/$dtot" "$([ "$dsame" = "$dtot" ] && echo ✅ || echo ❌)"
  echo
fi

# ── 四、沙箱泄漏检测 ────────────────────────────────────────────
echo "【四】沙箱泄漏检测（原目录中今天被改动的文件）"
sep
# 白名单：文档（人工编辑）+ 05.bqsr recal 表（曾误写、已还原，内容已正确）
leak=$(find . -maxdepth 3 -type f -newermt "$(date +%Y-%m-%d) 00:00" \
  -not -path "./13.snakemake/*" -not -path "./scripts/*" -not -path "./logs/*" \
  -not -path "./backups/*" -not -path "./.snakemake/*" \
  -not -path "./WGS_tutorial_git/*" -not -path "./WGS_tutorial_git" \
  -not -name "*.md" -not -name "*.tgz" 2>/dev/null | grep -v "05.bqsr/NA12878_chr20_recal_data.table" | head -20)
if [ -z "$leak" ]; then
  echo "  ✅ 无泄漏（文档与已还原文件已白名单）"
else
  echo "$leak" | sed 's/^/  ⚠️  /'
  echo "  ⚠️  共 $(echo "$leak" | wc -l) 个文件 —— 检查是否有脚本越界写入"
fi
echo
echo "【五】根目录残留检测（散落的中间文件）"
sep
# 允许出现在根目录的：文档 + 已知目录
STRAY=$(find . -maxdepth 1 -type f \
  ! -name "WGS_tutorial.md" ! -name "DEVELOPMENT.md" ! -name "SKILL*.md" \
  ! -name ".gitignore" 2>/dev/null)
if [ -z "$STRAY" ]; then
  echo "  ✅ 根目录无散落文件（仅 3 个文档）"
else
  echo "$STRAY" | sed 's/^/  ⚠️  /'
  echo "  ⚠️  共 $(echo "$STRAY" | wc -l) 个 —— 通常是脚本以相对路径写出（工作目录 = 项目根）"
  echo "      确认无引用后可删：rm -f <文件>"
fi
echo
echo "  ℹ️  .snakemake/ 为 Snakemake 自带元数据目录（CWD 下），可安全删除，重跑会重建"
echo
echo "════════════════════════════════════════════════════════════════════════"
echo " 汇总：产物缺失 $miss 个 | 泄漏 $([ -z "$leak" ] && echo 0 || echo "$(echo "$leak"|wc -l)") 个 | 根目录残留 $([ -z "$STRAY" ] && echo 0 || echo "$(echo "$STRAY"|wc -l)") 个"
echo "════════════════════════════════════════════════════════════════════════"
