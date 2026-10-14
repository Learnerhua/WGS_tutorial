#!/usr/bin/env bash
# Session 19（第 12 步）主要结果图：每张同时输出 PNG + PDF
#
# 只画「主要分析结果」（覆盖 4.10 / 4.12 / 4.13 三步的核心产出）：
#   图 1  4.10  变异密度 + Ts/Tv 双面板          （SNV/INDEL 检出质量）
#   图 2  4.10  QC 仪表板 VAF / GQ / 深度        （基因型质量）
#   图 3  4.10  F1 基准柱状图（hap.py）          （准确性）
#   图 4  4.12  chr20 环形 CNV log2 总览         （CNV 结果）
#   图 5  4.13  SURVIVOR 共识 SV 长度 + 类型分布 （SV 结果）
#
# 用法：bash scripts/12.run_matplotlib_plots.sh [density|qc|f1|circos|svlen|all]
#   已存在的图（PNG + PDF 均在场）会跳过

set -uo pipefail

WS=/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial
# 输出根：默认原地（$WS）；设为 13.snakemake 时全部产物重定向到沙箱
ROOT="${WGS_OUT_ROOT:-$WS}"

# ── 样本名（可用环境变量覆盖 → 换样本复用本流程）──────────────────
SAMPLE="${SAMPLE:-NA12878_chr20}"                                  # 主前缀（01–09 产物）
SAMPLE_SHORT="${SAMPLE_SHORT:-${SAMPLE%%_chr*}}"                   # 短前缀（10–12 产物）
SAMPLE_LABEL="${SAMPLE_LABEL:-$(printf '%s' "$SAMPLE" | tr '_' ' ')}"  # 展示文本（图标题）

export PATH=/oldhome/ouyjh/miniforge3/envs/WGS/bin:$PATH

OUT="$ROOT/12.figures"
LOG="$WS/logs/12.run_matplotlib_plots.log"
ONLY="${1:-all}"
mkdir -p "$OUT" "$WS/logs"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
want() { [ "$ONLY" = "all" ] || [ "$ONLY" = "$1" ]; }

echo "=============================================="
echo " Session 19: matplotlib main-result figures"
echo "=============================================="
echo " Out : $OUT"
echo " Run : $ONLY"
echo

command -v python3 >/dev/null || { echo "ERROR: python3" | tee -a "$LOG"; exit 1; }
log "matplotlib $(python3 -c 'import matplotlib;print(matplotlib.__version__)')"
echo

# ══════════════════════════════════════════════════════
# 图 1  4.10 变异密度 + Ts/Tv
# ══════════════════════════════════════════════════════
if want density; then
PNG="$OUT/${SAMPLE}_variant_density.png"; PDF="$OUT/${SAMPLE}_variant_density.pdf"
if [ -f "$PNG" ] && [ -f "$PDF" ]; then
  log "SKIP  图1 变异密度"
else
log "[1/5] 4.10 变异密度（500 kb 窗口 + Ts/Tv）..."
SMP="$SAMPLE" SMPS="$SAMPLE_SHORT" SML="$SAMPLE_LABEL" python3 - "$PNG" "$PDF" <<'PYEOF' 2>&1 | grep -vE "UserWarning|Glyph" | tee -a "$LOG"
import sys, gzip
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt

png_path, pdf_path = sys.argv[1], sys.argv[2]
import os as _os                       # 样本名由 bash 侧 export 传入
SMP  = _os.environ.get('SMP',  'NA12878_chr20')
SMPS = _os.environ.get('SMPS', 'NA12878')
SML  = _os.environ.get('SML',  'NA12878 chr20')
WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'
vcf_path = f'{WS}/07.var.filtered/{SMP}_filtered.vcf.gz'
WINDOW = 500_000
TRANSITIONS = {("A","G"),("G","A"),("C","T"),("T","C")}

counts = defaultdict(lambda: defaultdict(int))
tstv   = defaultdict(lambda: [0, 0])
opener = gzip.open if vcf_path.endswith('.gz') else open
with opener(vcf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 8 or f[6] != 'PASS': continue
        pos, ref, alt = int(f[1]), f[3], f[4].split(',')[0]
        w = pos // WINDOW
        counts[w]['SNV' if (len(ref)==1 and len(alt)==1) else 'INDEL'] += 1
        if len(ref)==1 and len(alt)==1:
            pair = tuple(sorted([ref.upper(), alt.upper()]))
            tstv[w][0 if pair in TRANSITIONS else 1] += 1

windows = sorted(counts.keys())
xs   = [w * WINDOW / 1e6 for w in windows]
snv  = [counts[w]['SNV']   for w in windows]
ind  = [counts[w]['INDEL'] for w in windows]
tv_w = [w for w in windows if w in tstv]
tv_v = [tstv[w][0] / (tstv[w][1] if tstv[w][1] else 1) for w in tv_w]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 7), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})
ax1.bar(xs, snv, color="#2E86C1", alpha=0.85, label="SNV",   width=0.48)
ax1.bar(xs, ind, color="#D97706", alpha=0.85, label="INDEL", width=0.48, bottom=snv)
ax1.axvspan(26, 30, alpha=0.12, color="#B91C1C", label="Centromere-proximal")
ax1.set_ylabel("Variants per 500 kb window")
ax1.set_title(f"Fig 1 | {SML} Germline Variant Density (PASS)", fontweight="bold")
ax1.legend(fontsize=9, loc="upper right")
ax1.set_ylim(bottom=0)

ax2.bar([w * WINDOW / 1e6 for w in tv_w], tv_v, color="#5B21B6", alpha=0.7, width=0.48)
ax2.axhline(2.1, color="#0E7C5B", ls="--", lw=1.2, label="Expected Ts/Tv = 2.1")
ax2.set_xlabel("chr20 position (Mb)")
ax2.set_ylabel("Ts/Tv")
ax2.legend(fontsize=8)
ax2.set_ylim(0, 4)

plt.tight_layout()
plt.savefig(png_path, dpi=150, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
ts = sum(tstv[w][0] for w in tstv); tv = sum(tstv[w][1] for w in tstv)
print(f"  → {png_path}")
print(f"  → {pdf_path}")
print(f"  窗口 {len(windows)}｜SNV {sum(snv)}｜INDEL {sum(ind)}｜Ts/Tv={ts/tv:.3f}")
PYEOF
fi
echo
fi

# ══════════════════════════════════════════════════════
# 图 2  4.10 QC 仪表板
# ══════════════════════════════════════════════════════
if want qc; then
PNG="$OUT/${SAMPLE_SHORT}_variant_qc_dashboard.png"; PDF="$OUT/${SAMPLE_SHORT}_variant_qc_dashboard.pdf"
if [ -f "$PNG" ] && [ -f "$PDF" ]; then
  log "SKIP  图2 QC 仪表板"
else
log "[2/5] 4.10 QC 仪表板（VAF + GQ + 深度）..."
SMP="$SAMPLE" SMPS="$SAMPLE_SHORT" SML="$SAMPLE_LABEL" python3 - "$PNG" "$PDF" <<'PYEOF' 2>&1 | grep -vE "UserWarning|Glyph" | tee -a "$LOG"
import sys, gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

png_path, pdf_path = sys.argv[1], sys.argv[2]
import os as _os                       # 样本名由 bash 侧 export 传入
SMP  = _os.environ.get('SMP',  'NA12878_chr20')
SMPS = _os.environ.get('SMPS', 'NA12878')
SML  = _os.environ.get('SML',  'NA12878 chr20')
WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'
vcf_path = f'{WS}/07.var.filtered/{SMP}_filtered.vcf.gz'

vafs, gqs, depths = [], [], []
opener = gzip.open if vcf_path.endswith('.gz') else open
with opener(vcf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 10 or f[6] != 'PASS': continue
        fmt = dict(zip(f[8].split(':'), f[9].split(':')))
        try:
            ad = [int(x) for x in fmt.get('AD', '0,0').split(',')]
            if sum(ad) > 0: vafs.append(ad[1] / sum(ad))
        except: pass
        try: gqs.append(int(fmt.get('GQ', '.')))
        except: pass
        try: depths.append(int(fmt.get('DP', '.')))
        except: pass

fig = plt.figure(figsize=(15, 5))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.32)

ax1 = fig.add_subplot(gs[0])
ax1.hist(vafs, bins=50, color="#2E86C1", alpha=0.8, edgecolor='white')
ax1.axvline(0.5, color="#B91C1C", ls='--', lw=1.5, label='Expected het 0.5')
ax1.axvline(1.0, color="#0E7C5B", ls='--', lw=1.5, label='Expected hom 1.0')
ax1.set_xlabel('Variant Allele Fraction'); ax1.set_ylabel('Count')
ax1.set_title(f'VAF Distribution\n(n={len(vafs)}, PASS)', fontweight='bold')
ax1.legend(fontsize=8)

ax2 = fig.add_subplot(gs[1])
ax2.hist(gqs, bins=50, color="#0E7C5B", alpha=0.8, edgecolor='white')
ax2.axvline(20, color="#B91C1C", ls='--', lw=1.5, label='GQ=20 threshold')
ax2.set_xlabel('Genotype Quality (GQ)'); ax2.set_ylabel('Count')
ax2.set_title('Genotype Quality Distribution', fontweight='bold')
ax2.legend(fontsize=8)

ax3 = fig.add_subplot(gs[2])
ax3.hist([min(d, 80) for d in depths], bins=50, color="#5B21B6", alpha=0.8, edgecolor='white')
ax3.axvline(np.median(depths), color="#D97706", ls='--', lw=1.5,
            label=f'Median {np.median(depths):.0f}x')
ax3.set_xlabel('Read Depth (capped at 80x)'); ax3.set_ylabel('Count')
ax3.set_title('Read Depth at PASS Variants', fontweight='bold')
ax3.legend(fontsize=8)

plt.suptitle(f'Fig 2 | {SML} Variant QC Dashboard', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(png_path, dpi=150, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
het = sum(1 for v in vafs if 0.2 <= v <= 0.8)
print(f'  → {png_path}')
print(f'  → {pdf_path}')
print(f'  n={len(vafs)}｜het 占比 {het/max(len(vafs),1)*100:.1f}%｜深度中位 {np.median(depths):.0f}×')
PYEOF
fi
echo
fi

# ══════════════════════════════════════════════════════
# 图 3  4.10 F1 基准
# ══════════════════════════════════════════════════════
if want f1; then
PNG="$OUT/${SAMPLE_SHORT}_F1_benchmark.png"; PDF="$OUT/${SAMPLE_SHORT}_F1_benchmark.pdf"
if [ -f "$PNG" ] && [ -f "$PDF" ]; then
  log "SKIP  图3 F1 基准"
else
log "[3/5] 4.10 F1 基准（hap.py）..."
SMP="$SAMPLE" SMPS="$SAMPLE_SHORT" SML="$SAMPLE_LABEL" python3 - "$PNG" "$PDF" <<'PYEOF' 2>&1 | grep -vE "UserWarning|Glyph" | tee -a "$LOG"
import sys, csv
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

png_path, pdf_path = sys.argv[1], sys.argv[2]
import os as _os                       # 样本名由 bash 侧 export 传入
SMP  = _os.environ.get('SMP',  'NA12878_chr20')
SMPS = _os.environ.get('SMPS', 'NA12878')
SML  = _os.environ.get('SML',  'NA12878 chr20')
WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'
summary = f'{WS}/09.benchmark/hapy_{SMP}.summary.csv'

rows = []
with open(summary) as fh:
    for r in csv.DictReader(fh):
        rows.append({
            'Type': r['Type'], 'Filter': r['Filter'],
            'Recall':    float(r['METRIC.Recall']) * 100,
            'Precision': float(r['METRIC.Precision']) * 100,
            'F1':        float(r['METRIC.F1_Score']) * 100,
        })

categories = ['SNP_PASS', 'SNP_ALL', 'INDEL_PASS', 'INDEL_ALL']
data = {c: [0, 0, 0] for c in categories}
for r in rows:
    k = f"{r['Type']}_{r['Filter']}"
    if k in data: data[k] = [r['Recall'], r['Precision'], r['F1']]

x = np.arange(len(categories)); width = 0.25
fig, ax = plt.subplots(figsize=(11, 6))
b1 = ax.bar(x - width, [data[c][0] for c in categories], width, label='Recall', color="#2E86C1", alpha=0.85)
b2 = ax.bar(x,         [data[c][1] for c in categories], width, label='Precision', color="#0E7C5B", alpha=0.85)
b3 = ax.bar(x + width, [data[c][2] for c in categories], width, label='F1', color="#D97706", alpha=0.85)
for bars in [b1, b2, b3]:
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.15,
                f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=8)

ax.set_xticks(x); ax.set_xticklabels(categories)
ax.set_ylabel('Percentage (%)'); ax.set_ylim(95, 100.2)
ax.set_title(f'Fig 3 | {SML} F1 Benchmark (hap.py)\nSNP / INDEL x ALL / PASS', fontweight='bold')
ax.axhline(99, color="#B91C1C", ls='--', lw=1, alpha=0.6)
ax.text(2.5, 99.05, '99% target', ha='center', fontsize=8, color="#B91C1C")
ax.legend(loc='lower right'); ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(png_path, dpi=150, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
print(f'  → {png_path}')
print(f'  → {pdf_path}')
for r in rows: print(f"  {r['Type']:5s} {r['Filter']:4s}  R={r['Recall']:.2f}  P={r['Precision']:.2f}  F1={r['F1']:.2f}")
PYEOF
fi
echo
fi

# ══════════════════════════════════════════════════════
# 图 4  4.12 chr20 环形 CNV 总览
# ══════════════════════════════════════════════════════
if want circos; then
PNG="$OUT/${SAMPLE}_circos.png"; PDF="$OUT/${SAMPLE}_circos.pdf"
if [ -f "$PNG" ] && [ -f "$PDF" ]; then
  log "SKIP  图4 环形 CNV"
else
log "[4/5] 4.12 chr20 环形 CNV log2 总览..."
SMP="$SAMPLE" SMPS="$SAMPLE_SHORT" SML="$SAMPLE_LABEL" python3 - "$PNG" "$PDF" <<'PYEOF' 2>&1 | grep -vE "UserWarning|Glyph" | tee -a "$LOG"
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import pandas as pd

png_path, pdf_path = sys.argv[1], sys.argv[2]
import os as _os                       # 样本名由 bash 侧 export 传入
SMP  = _os.environ.get('SMP',  'NA12878_chr20')
SMPS = _os.environ.get('SMPS', 'NA12878')
SML  = _os.environ.get('SML',  'NA12878 chr20')
WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'
cnr_path = f'{WS}/10.cnv/results/{SMPS}.cnr'
CHR20_LEN = 64_444_167

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
theta = np.linspace(0, 2*np.pi, 1000)
ax.plot(theta, [1.0]*1000, color="#BDD9F2", lw=20, solid_capstyle="round")

n_bins = n_del = n_dup = 0
if os.path.exists(cnr_path):
    cnr = pd.read_csv(cnr_path, sep="\t")
    c20 = cnr[cnr["chromosome"] == "chr20"].copy()
    c20["theta"] = c20["start"] / CHR20_LEN * 2 * np.pi
    colors = c20["log2"].apply(lambda v: "#B91C1C" if v < -0.3 else "#0E7C5B" if v > 0.3 else "#2E86C1")
    ax.scatter(c20["theta"], 1.0 + c20["log2"] * 0.15, c=colors, s=2, alpha=0.5, zorder=3)
    n_bins = len(c20)
    n_del = int((c20["log2"] < -0.3).sum()); n_dup = int((c20["log2"] > 0.3).sum())

cen_s = 26e6 / CHR20_LEN * 2 * np.pi
cen_e = 30e6 / CHR20_LEN * 2 * np.pi
ax.plot([cen_s, cen_e], [1.0, 1.0], lw=25, color="#D97706", solid_capstyle="round", zorder=4)
ax.set_ylim(0, 1.5); ax.set_yticks([]); ax.set_xticks([])
ax.set_title(f"Fig 4 | {SML} CNV log2 Overview\n"
             "(red=DEL, green=DUP, blue=neutral, orange=centromere)", fontsize=11, pad=20)

plt.tight_layout()
plt.savefig(png_path, dpi=150, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
print(f"  → {png_path}")
print(f"  → {pdf_path}")
print(f"  bins={n_bins}｜DEL(log2<-0.3)={n_del}｜DUP(log2>0.3)={n_dup}")
PYEOF
fi
echo
fi

# ══════════════════════════════════════════════════════
# 图 5  4.13 SV 长度 + 类型分布
# ══════════════════════════════════════════════════════
if want svlen; then
PNG="$OUT/${SAMPLE_SHORT}_SV_length_distribution.png"; PDF="$OUT/${SAMPLE_SHORT}_SV_length_distribution.pdf"
if [ -f "$PNG" ] && [ -f "$PDF" ]; then
  log "SKIP  图5 SV 长度分布"
else
log "[5/5] 4.13 SV 长度 + 类型分布..."
SMP="$SAMPLE" SMPS="$SAMPLE_SHORT" SML="$SAMPLE_LABEL" python3 - "$PNG" "$PDF" <<'PYEOF' 2>&1 | grep -vE "UserWarning|Glyph" | tee -a "$LOG"
import sys, gzip
from collections import Counter, defaultdict
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

png_path, pdf_path = sys.argv[1], sys.argv[2]
import os as _os                       # 样本名由 bash 侧 export 传入
SMP  = _os.environ.get('SMP',  'NA12878_chr20')
SMPS = _os.environ.get('SMPS', 'NA12878')
SML  = _os.environ.get('SML',  'NA12878 chr20')
WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'
vcf_path = f'{WS}/11.sv/merged/{SMP}_SVs_merged.vcf.gz'

svs = []
opener = gzip.open if vcf_path.endswith('.gz') else open
with opener(vcf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        p = line.rstrip().split('\t')
        if len(p) < 8: continue
        info = dict(kv.split('=') for kv in p[7].split(';') if '=' in kv)
        svs.append((info.get('SVTYPE', '?'), abs(int(info.get('SVLEN', 0)))))

type_colors = {'DEL': '#B91C1C', 'DUP': '#0E7C5B', 'INS': '#2E86C1',
               'INV': '#D97706', 'BND': '#5B21B6'}
labels = ['50-100bp', '100-500bp', '0.5-1kb', '1-5kb', '5-10kb', '>10kb']
edges  = [50, 100, 500, 1000, 5000, 10000, 10**9]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax1 = axes[0]
by_type = defaultdict(list)
for t, l in svs: by_type[t].append(l)
for t, lens in sorted(by_type.items()):
    ax1.hist(lens, bins=np.logspace(1, 6, 30), alpha=0.7,
             label=f'{t} (n={len(lens)})', color=type_colors.get(t, '#999999'))
ax1.set_xscale('log')
ax1.set_xlabel('SV length (bp, log scale)'); ax1.set_ylabel('Count')
ax1.set_title('SURVIVOR consensus SV length distribution', fontweight='bold')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2 = axes[1]
bin_data = defaultdict(lambda: defaultdict(int))
for t, l in svs:
    for i in range(len(edges)-1):
        if edges[i] < l <= edges[i+1]:
            bin_data[labels[i]][t] += 1; break
bottom = np.zeros(len(labels))
for t in sorted(set(x for x, _ in svs)):
    vals = [bin_data[lbl].get(t, 0) for lbl in labels]
    ax2.bar(labels, vals, bottom=bottom, label=t,
            color=type_colors.get(t, '#999999'), edgecolor='white')
    bottom += np.array(vals)
ax2.set_xlabel('Length bin'); ax2.set_ylabel('Count')
ax2.set_title(f'Length bin x Type (n={len(svs)} SVs)', fontweight='bold')
ax2.legend(); ax2.tick_params(axis='x', rotation=20)
ax2.set_ylim(0, bottom.max() * 1.15)   # 留出标签空间
for i, total in enumerate(bottom):
    if total > 0:
        ax2.text(i, total + 0.4, str(int(total)), ha='center', va='bottom',
                 fontsize=9, fontweight='bold')

plt.suptitle(f'Fig 5 | {SML} Structural Variants (Manta + DELLY consensus)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(png_path, dpi=150, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
print(f'  → {png_path}')
print(f'  → {pdf_path}')
for t, c in Counter(x for x, _ in svs).most_common(): print(f'  {t}: {c}')
PYEOF
fi
echo
fi

# ══════════════════════════════════════════════════════
# 图 6  4.9–4.10 变异注释结果总览（4 面板）
# ══════════════════════════════════════════════════════
if want anno; then
PNG="$OUT/${SAMPLE_SHORT}_annotation_summary.png"; PDF="$OUT/${SAMPLE_SHORT}_annotation_summary.pdf"
if [ -f "$PNG" ] && [ -f "$PDF" ]; then
  log "SKIP  图6 注释总览"
else
log "[6/6] 4.9–4.10 变异注释结果总览（gnomAD/ClinVar/SIFT-PolyPhen/LOEUF）..."
SMP="$SAMPLE" SMPS="$SAMPLE_SHORT" SML="$SAMPLE_LABEL" python3 - "$PNG" "$PDF" <<'PYEOF' 2>&1 | grep -vE "UserWarning|Glyph" | tee -a "$LOG"
import sys, gzip, os
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec

png_path, pdf_path = sys.argv[1], sys.argv[2]
import os as _os                       # 样本名由 bash 侧 export 传入
SMP  = _os.environ.get('SMP',  'NA12878_chr20')
SMPS = _os.environ.get('SMPS', 'NA12878')
SML  = _os.environ.get('SML',  'NA12878 chr20')
WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'

fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.30)

# ── 面板 A：gnomAD 频率分层（4.10 步产物）──
axA = fig.add_subplot(gs[0, 0])
tier_p = f'{WS}/09.databases/gnomad_freq_tier.tsv'
tiers = []
if os.path.exists(tier_p):
    with open(tier_p) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                tiers.append((p[0], int(p[1]), float(p[2].rstrip('%'))))
if tiers:
    tier_en = {'novel (AF=0/缺失)': 'novel (AF=0/missing)',
           'ultra-rare (AF<0.001)': 'ultra-rare (AF<0.001)',
           'rare (0.001<=AF<0.01)': 'rare (0.001<=AF<0.01)',
           'common (AF>=0.01)': 'common (AF>=0.01)'}
    labels = [tier_en.get(t[0], t[0]) for t in tiers]
    vals   = [t[1] for t in tiers]
    pcts   = [t[2] for t in tiers]
    colors = ['#D97706', '#B91C1C', '#5B21B6', '#2E86C1']
    bars = axA.barh(labels, vals, color=colors[:len(labels)], alpha=0.85)
    axA.set_xscale('log')
    for bar, v, pc in zip(bars, vals, pcts):
        axA.text(bar.get_width() * 1.15, bar.get_y() + bar.get_height()/2,
                 f'{v:,} ({pc:.1f}%)', va='center', fontsize=9)
    axA.set_xlabel('Variant count (log scale)')
    axA.set_title('A | gnomAD Population Frequency Tiers\n(4.10 database crossref)', fontweight='bold')
    axA.set_xlim(right=max(vals) * 5)
    axA.grid(True, alpha=0.3, axis='x')

# ── 面板 B：ClinVar 致病性分类（4.10 步产物）──
axB = fig.add_subplot(gs[0, 1])
cv_p = f'{WS}/09.databases/clinvar_classification.tsv'
cvs = []
if os.path.exists(cv_p):
    with open(cv_p) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                cvs.append((p[0].replace('_', ' '), int(p[1])))
if cvs:
    cvs.sort(key=lambda x: x[1])
    labels = [c[0] for c in cvs]; vals = [c[1] for c in cvs]
    cmap = {'Benign': '#0E7C5B', 'Benign/Likely benign': '#2E86C1',
            'Likely benign': '#5B21B6', 'association': '#D97706',
            'Uncertain significance': '#B91C1C',
            'Conflicting classifications of pathogenicity': '#D97706',
            'confers sensitivity': '#B91C1C'}
    colors = [cmap.get(l, '#999999') for l in labels]
    bars = axB.barh(labels, vals, color=colors, alpha=0.85)
    for bar, v in zip(bars, vals):
        axB.text(bar.get_width() + max(vals)*0.015, bar.get_y() + bar.get_height()/2,
                 str(v), va='center', fontsize=9)
    axB.set_xlabel('Variant count')
    axB.set_title('B | ClinVar Clinical Significance\n(no Pathogenic; Benign = 93%)', fontweight='bold')
    axB.set_xlim(right=max(vals) * 1.18)
    axB.grid(True, alpha=0.3, axis='x')

# ── 面板 C：SIFT vs PolyPhen 预测（VEP CSQ 30/31 字段）──
axC = fig.add_subplot(gs[1, 0])
vcf = f'{WS}/08.annotation/{SMP}_VEP.vcf.gz'
sift_vals, pp_vals = [], []
opener = gzip.open if vcf.endswith('.gz') else open
with opener(vcf, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 8: continue
        csq = dict(kv.split('=', 1) for kv in p[7].split(';') if kv.startswith('CSQ='))
        if 'CSQ' not in csq: continue
        best = None
        for seg in csq['CSQ'].split(','):
            f = seg.split('|')
            if len(f) < 32: continue
            s, pp = f[29], f[30]
            # VEP 格式：deleterious(0.02) / probably_damaging(0.99) → 取括号内数值
            def _num(x):
                if not x: return None
                x = x.strip()
                if '(' in x and x.endswith(')'):
                    x = x[x.index('(')+1:-1]
                try:
                    return float(x)
                except ValueError:
                    return None
            sv, pv = _num(s), _num(pp)
            try:
                if sv is None or pv is None: continue
                if best is None: best = (sv, pv)
            except (ValueError, TypeError): continue
        if best: sift_vals.append(best[0]); pp_vals.append(best[1])

if sift_vals:
    axC.scatter(sift_vals, pp_vals, s=14, alpha=0.5, color='#5B21B6')
    axC.axvline(0.05, color='#B91C1C', ls='--', lw=1, label='SIFT 0.05 (deleterious)')
    axC.axhline(0.909, color='#0E7C5B', ls='--', lw=1, label='PolyPhen 0.909 (deleterious)')
    axC.set_xlabel('SIFT score (lower = more deleterious)')
    axC.set_ylabel('PolyPhen score (higher = more deleterious)')
    axC.set_title(f'C | SIFT vs PolyPhen\n(n={len(sift_vals)} variants with predictions)',
                  fontweight='bold')
    axC.legend(fontsize=8); axC.grid(True, alpha=0.3)
    harmful = sum(1 for s, pp in zip(sift_vals, pp_vals) if s < 0.05 and pp > 0.909)
    axC.text(0.02, 0.97, f'Concordant deleterious: {harmful}', transform=axC.transAxes,
             va='top', fontsize=9, color='#B91C1C', fontweight='bold')

# ── 面板 D：43 候选的 LOEUF 分布（基因约束）──
axD = fig.add_subplot(gs[1, 1])
cand_p = f'{WS}/09.databases/{SMP}_43candidates_clinically_annotated.tsv'
loeuf = []
if os.path.exists(cand_p):
    with open(cand_p) as fh:
        hdr = next(fh, '').rstrip('\n').split('\t')
        try:
            gi = hdr.index('gnomAD_LOEUF'); pi = hdr.index('gnomAD_pLI')
        except ValueError:
            gi = pi = None
        if gi is not None:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) > gi and p[gi]:
                    try: loeuf.append((p[0], float(p[gi])))
                    except ValueError: pass
if loeuf:
    loeuf.sort(key=lambda x: x[1])
    names = [x[0][:14] for x in loeuf]; vals = [x[1] for x in loeuf]
    colors = ['#B91C1C' if v < 0.35 else '#D97706' if v < 0.7 else '#0E7C5B' for v in vals]
    bars = axD.barh(names, vals, color=colors, alpha=0.85)
    axD.axvline(0.35, color='#B91C1C', ls='--', lw=1.2, label='LOEUF 0.35 (strongly constrained)')
    axD.set_xlabel('LOEUF (lower = more constrained)')
    axD.set_title(f'D | LOEUF Constraint of 43 Candidate Genes\n(n={len(loeuf)} genes with LOEUF)',
                  fontweight='bold')
    axD.legend(fontsize=8, loc='lower right')
    axD.grid(True, alpha=0.3, axis='x')
    axD.tick_params(axis='y', labelsize=7)

plt.suptitle(f'Fig 6 | {SML} Variant Annotation Summary (4.9 VEP/ClinVar + 4.10 crossref)',
             fontsize=14, fontweight='bold', y=0.985)
plt.savefig(png_path, dpi=150, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
print(f'  → {png_path}')
print(f'  → {pdf_path}')
print(f'  gnomAD 分层: {[(t[0], t[1]) for t in tiers]}')
print(f'  ClinVar 分类: {cvs}')
print(f'  SIFT/PolyPhen: n={len(sift_vals)}')
print(f'  LOEUF: n={len(loeuf)}')
PYEOF
fi
echo
fi

log "=== 产出（PNG + PDF 成对）==="
for f in "$OUT"/*.png "$OUT"/*.pdf; do
  [ -f "$f" ] && printf "  %-52s %9s B\n" "$(basename "$f")" "$(stat -c%s "$f")" | tee -a "$LOG"
done
echo
log "Done."