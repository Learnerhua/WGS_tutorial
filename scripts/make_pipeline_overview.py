#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_pipeline_overview.py — 生成 README 用的「流程总览」矢量图（中 / 英双语）。

用法：
    python3 scripts/make_pipeline_overview.py            # 写入仓库根目录
    python3 scripts/make_pipeline_overview.py --out DIR  # 指定输出目录

产物：
    pipeline_overview_zh.svg   简体中文版
    pipeline_overview_en.svg   English version

设计要点：
    • 白底浅色主题，颜色按分析阶段编码，底部附颜色图例；
    • 三个并行分支面板：拷贝数变异 CNV / 生殖系 SNV-INDEL / 结构变异 SV；
    • 分支汇聚到「可视化」与「Snakemake 编排」；
    • 纯手写 SVG，无第三方依赖，可在 GitHub README 中直接引用。
"""

import argparse
import math
import os

# ─────────────────────────────────────────────────────────────
# 画布与栅格
# ─────────────────────────────────────────────────────────────
W, H = 1480, 1720

TITLE_Y = 50

NODE_H = 76
NODE_GAP = 36
NODE_W = 380
PANEL_W = 420
PANEL_Y = 704

MID_CX = 740          # 画布中线，也是主脊（预处理 + SNV/INDEL）的轴
LEFT_CX, RIGHT_CX = 250, 1230
LEFT_PX, MID_PX, RIGHT_PX = 40, 530, 1020   # 面板左上角 x

# 主脊（预处理）各节点顶边 y
Y_INPUT = 110
Y_FASTP = Y_INPUT + NODE_H + 42          # 228
Y_BWA = Y_FASTP + NODE_H + 42            # 346
Y_POSTQC = Y_BWA + NODE_H + 42           # 464
Y_MARKDUP = Y_POSTQC + NODE_H + 42       # 582

BUS_Y = 682                              # 分支总线
NODE_Y0 = 760                            # 各分支首个节点顶边

Y_VIZ = 1390
Y_SNAKE = 1530
Y_LEGEND = 1660

ARROW = "#94A3B8"
LINE = "#CBD5E1"

# ─────────────────────────────────────────────────────────────
# 阶段配色（浅色主题）
# ─────────────────────────────────────────────────────────────
PHASE = {
    "slate":  {"border": "#94A3B8", "accent": "#475569", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "blue":   {"border": "#3B82F6", "accent": "#2563EB", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "cyan":   {"border": "#06B6D4", "accent": "#0891B2", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "green":  {"border": "#10B981", "accent": "#047857", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "indigo": {"border": "#6366F1", "accent": "#4F46E5", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "amber":  {"border": "#F59E0B", "accent": "#B45309", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "purple": {"border": "#A855F7", "accent": "#7E22CE", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "rose":   {"border": "#F43F5E", "accent": "#BE123C", "node": "#FFFFFF", "title": "#0F172A", "sub": "#64748B"},
    "dark":   {"border": "#0F172A", "accent": "#38BDF8", "node": "#1E293B", "title": "#F8FAFC", "sub": "#CBD5E1"},
}

PANEL_FILL = {
    "amber":  ("#FFFBEB", "#FDE68A"),
    "green":  ("#F0FDF4", "#BBF7D0"),
    "purple": ("#FAF5FF", "#E9D5FF"),
}

FONT = ("system-ui,-apple-system,'PingFang SC','Hiragino Sans GB','Microsoft YaHei',"
        "'Noto Sans CJK SC','Source Han Sans SC','Segoe UI',Roboto,Helvetica,Arial,sans-serif")


# ─────────────────────────────────────────────────────────────
# 基础绘图工具
# ─────────────────────────────────────────────────────────────
def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=15, weight="400", fill="#0F172A", anchor="middle", spacing=None):
    sp = f' letter-spacing="{spacing}"' if spacing else ""
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{sp}>{esc(s)}</text>')


def rect(x, y, w, h, fill="#FFFFFF", stroke="none", sw=1.5, rx=12, dash=None, filt=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    f = f' filter="url(#{filt})"' if filt else ""
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{d}{f}/>')


def rounded_path(pts, r=14):
    """折线 → 圆角路径（用于汇流总线）。"""
    if len(pts) < 2:
        return ""
    d = f"M {pts[0][0]},{pts[0][1]}"
    for i in range(1, len(pts) - 1):
        px, py = pts[i - 1]
        cx, cy = pts[i]
        nx, ny = pts[i + 1]

        def unit(ax, ay):
            L = math.hypot(ax, ay) or 1.0
            return ax / L, ay / L

        ux1, uy1 = unit(px - cx, py - cy)
        ux2, uy2 = unit(nx - cx, ny - cy)
        rr = min(r, math.hypot(px - cx, py - cy) / 2, math.hypot(nx - cx, ny - cy) / 2)
        d += (f" L {cx + ux1 * rr},{cy + uy1 * rr}"
              f" Q {cx},{cy} {cx + ux2 * rr},{cy + uy2 * rr}")
    d += f" L {pts[-1][0]},{pts[-1][1]}"
    return d


def edge(pts, dash=None, arrow=True, color=None):
    color = color or ARROW
    d = rounded_path(pts)
    da = f' stroke-dasharray="{dash}"' if dash else ""
    ar = ' marker-end="url(#ah)"' if arrow else ""
    return (f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2" '
            f'stroke-linecap="round"{da}{ar}/>')


def node(x, y, w, h, phase, title, sub, dash=None):
    c = PHASE[phase]
    out = [rect(x, y, w, h, fill=c["node"], stroke=c["border"], sw=1.6, rx=14,
                dash=dash, filt="sh")]
    # 左侧强调条
    out.append(f'<rect x="{x + 16}" y="{y + 18}" width="5" height="{h - 36}" rx="2.5" '
               f'fill="{c["accent"]}"/>')
    cy = y + h / 2
    out.append(text(x + w / 2, cy - 4, title, size=18, weight="700", fill=c["title"]))
    out.append(text(x + w / 2, cy + 20, sub, size=13, weight="400", fill=c["sub"]))
    return "".join(out)


def panel(x, y, w, h, phase, label):
    fill, stroke = PANEL_FILL[phase]
    accent = PHASE[phase]["accent"]
    out = [rect(x, y, w, h, fill=fill, stroke=stroke, sw=1.5, rx=18)]
    out.append(f'<circle cx="{x + 30}" cy="{y + 32}" r="6" fill="{accent}"/>')
    out.append(text(x + 46, y + 38, label, size=16, weight="700", fill=accent, anchor="start",
                    spacing="0.4"))
    return "".join(out)


# ─────────────────────────────────────────────────────────────
# 中英文文案
# ─────────────────────────────────────────────────────────────
CONTENT = {
    "zh": {
        "title": "NA12878 · chr20 生殖系 WGS 全流程",
        "out": "pipeline_overview_zh.svg",
        "input": ("原始数据 rawData/", "FASTQ · chr20 · ~30× 双端测序"),
        "fastqc": ("4.1 FastQC + MultiQC", "原始测序质量报告"),
        "ref": ("参考基因组 GRCh38", "Homo_sapiens_assembly38"),
        "fastp": ("4.2 fastp 修剪", "接头 / 低质量碱基"),
        "bwa": ("4.3 BWA-MEM2 比对", "→ GRCh38 坐标排序 BAM"),
        "postqc": ("4.4 比对后 QC", "mosdepth / Picard / MultiQC"),
        "markdup": ("4.5 MarkDuplicates + BQSR", "→ 分析就绪 BAM"),
        "cnv_label": "拷贝数变异 · CNV",
        "cnv": [("4.12 CNVkit", "拷贝数检测 → 分段"),
                ("4.12 CNV 注释", "基因 / 区间映射")],
        "snv_label": "生殖系 SNV / INDEL",
        "snv": [("4.6 HaplotypeCaller", "→ GVCF"),
                ("4.8 基因分型 + 硬过滤", "SNV / INDEL → 合并 VCF"),
                ("4.9 VEP 注释", "功能后果预测"),
                ("4.10 数据库交叉查询", "ClinVar / gnomAD / OMIM / HPO"),
                ("4.10 hap.py F1 基准", "vs GIAB 真集")],
        "snv_phases": ["green", "green", "indigo", "indigo", "green"],
        "sv_label": "结构变异 · SV",
        "sv": [("4.13 Manta + DELLY", "PASS 过滤"),
               ("4.13 SURVIVOR 合并", "共识 SV 集合"),
               ("4.13 gnomAD-SV 注释", "人群频率"),
               ("4.13 SV 基因注释", "临床 review")],
        "viz": ("4.14 可视化", "matplotlib / ggplot2 / igv-reports / IGV"),
        "snake": ("4.15 Snakemake 编排", "23 条规则 · 沙箱可复现"),
        "legend": [("输入", "slate"), ("质控", "blue"), ("比对", "cyan"),
                   ("SNV/INDEL", "green"), ("注释", "indigo"), ("CNV", "amber"),
                   ("SV", "purple"), ("可视化", "rose"), ("编排", "dark")],
    },
    "en": {
        "title": "NA12878 · chr20 Germline WGS Pipeline",
        "out": "pipeline_overview_en.svg",
        "input": ("Raw data  rawData/", "FASTQ · chr20 · ~30× paired-end"),
        "fastqc": ("4.1 FastQC + MultiQC", "Raw sequencing QC report"),
        "ref": ("Reference  GRCh38", "Homo_sapiens_assembly38"),
        "fastp": ("4.2 fastp trimming", "Adapters / low-quality bases"),
        "bwa": ("4.3 BWA-MEM2 alignment", "→ GRCh38 coordinate-sorted BAM"),
        "postqc": ("4.4 Post-alignment QC", "mosdepth / Picard / MultiQC"),
        "markdup": ("4.5 MarkDuplicates + BQSR", "→ Analysis-ready BAM"),
        "cnv_label": "Copy-number · CNV",
        "cnv": [("4.12 CNVkit", "Copy-number calling → segments"),
                ("4.12 CNV annotation", "Gene / interval mapping")],
        "snv_label": "Germline SNV / INDEL",
        "snv": [("4.6 HaplotypeCaller", "→ GVCF"),
                ("4.8 Genotyping + hard filter", "SNV / INDEL → merged VCF"),
                ("4.9 VEP annotation", "Functional consequence"),
                ("4.10 Database cross-reference", "ClinVar / gnomAD / OMIM / HPO"),
                ("4.10 hap.py F1 benchmark", "vs GIAB truth set")],
        "snv_phases": ["green", "green", "indigo", "indigo", "green"],
        "sv_label": "Structural · SV",
        "sv": [("4.13 Manta + DELLY", "PASS filtering"),
               ("4.13 SURVIVOR merge", "Consensus SV set"),
               ("4.13 gnomAD-SV annotation", "Population frequency"),
               ("4.13 SV gene annotation", "Clinical review")],
        "viz": ("4.14 Visualization", "matplotlib / ggplot2 / igv-reports / IGV"),
        "snake": ("4.15 Snakemake orchestration", "23 rules · sandboxed & reproducible"),
        "legend": [("Input", "slate"), ("QC", "blue"), ("Alignment", "cyan"),
                   ("SNV/INDEL", "green"), ("Annotation", "indigo"), ("CNV", "amber"),
                   ("SV", "purple"), ("Visualization", "rose"), ("Orchestration", "dark")],
    },
}


# ─────────────────────────────────────────────────────────────
# 组装 SVG
# ─────────────────────────────────────────────────────────────
def build(lang):
    C = CONTENT[lang]
    o = []
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
             f'width="{W}" height="{H}" role="img" aria-label="{esc(C["title"])}">')

    # defs
    o.append('<defs>'
             '<filter id="sh" x="-25%" y="-25%" width="150%" height="150%">'
             '<feDropShadow dx="0" dy="1.5" stdDeviation="2.6" flood-color="#334155" flood-opacity="0.10"/>'
             '</filter>'
             '<marker id="ah" viewBox="0 0 12 12" refX="10" refY="6" markerWidth="12" '
             'markerHeight="12" markerUnits="userSpaceOnUse" orient="auto">'
             f'<path d="M1,1 L11,6 L1,11 z" fill="{ARROW}"/></marker>'
             '</defs>')

    # 背景
    o.append(rect(0, 0, W, H, fill="#FFFFFF", rx=0))

    # 标题
    o.append(text(MID_CX, TITLE_Y, C["title"], size=25, weight="700", fill="#0F172A",
                  spacing="0.4"))

    # ── 连线（先画，压在节点下方）──────────────────────────
    # 主脊
    o.append(edge([(MID_CX, Y_INPUT + NODE_H), (MID_CX, Y_FASTP)]))
    o.append(edge([(MID_CX, Y_FASTP + NODE_H), (MID_CX, Y_BWA)]))
    o.append(edge([(MID_CX, Y_BWA + NODE_H), (MID_CX, Y_POSTQC)]))
    o.append(edge([(MID_CX, Y_POSTQC + NODE_H), (MID_CX, Y_MARKDUP)]))
    # 输入 → FastQC（侧向 QC，虚线）
    o.append(edge([(MID_PX + NODE_W + 20, Y_INPUT + NODE_H / 2),
                   (RIGHT_PX, Y_INPUT + NODE_H / 2)], dash="7 5"))
    # 参考 → BWA-MEM2（虚线）
    o.append(edge([(LEFT_PX + PANEL_W, Y_BWA + NODE_H / 2),
                   (MID_PX + 20, Y_BWA + NODE_H / 2)], dash="7 5"))
    # 总线分流
    o.append(edge([(MID_CX, Y_MARKDUP + NODE_H), (MID_CX, BUS_Y),
                   (LEFT_CX, BUS_Y), (LEFT_CX, PANEL_Y)]))
    o.append(edge([(MID_CX, BUS_Y), (RIGHT_CX, BUS_Y), (RIGHT_CX, PANEL_Y)]))
    o.append(edge([(MID_CX, BUS_Y), (MID_CX, PANEL_Y)]))
    # 汇流到可视化
    o.append(edge([(MID_CX, PANEL_Y + 600), (MID_CX, Y_VIZ)]))
    o.append(edge([(LEFT_CX, PANEL_Y + 264), (LEFT_CX, Y_VIZ - 30),
                   (490, Y_VIZ - 30), (490, Y_VIZ)]))
    o.append(edge([(RIGHT_CX, PANEL_Y + 488), (RIGHT_CX, Y_VIZ - 30),
                   (990, Y_VIZ - 30), (990, Y_VIZ)]))
    # 可视化 → 编排
    o.append(edge([(MID_CX, Y_VIZ + 88), (MID_CX, Y_SNAKE)]))

    # ── 分支面板 ─────────────────────────────────────────
    o.append(panel(LEFT_PX, PANEL_Y, PANEL_W, 264, "amber", C["cnv_label"]))
    o.append(panel(MID_PX, PANEL_Y, PANEL_W, 600, "green", C["snv_label"]))
    o.append(panel(RIGHT_PX, PANEL_Y, PANEL_W, 488, "purple", C["sv_label"]))

    # ── 节点 ─────────────────────────────────────────────
    o.append(node(MID_PX + 20, Y_INPUT, NODE_W, NODE_H, "slate", *C["input"]))
    o.append(node(RIGHT_PX, Y_INPUT, NODE_W, NODE_H, "blue", *C["fastqc"], dash="7 5"))
    o.append(node(LEFT_PX, Y_BWA, NODE_W, NODE_H, "slate", *C["ref"], dash="7 5"))
    o.append(node(MID_PX + 20, Y_FASTP, NODE_W, NODE_H, "blue", *C["fastp"]))
    o.append(node(MID_PX + 20, Y_BWA, NODE_W, NODE_H, "cyan", *C["bwa"]))
    o.append(node(MID_PX + 20, Y_POSTQC, NODE_W, NODE_H, "cyan", *C["postqc"]))
    o.append(node(MID_PX + 20, Y_MARKDUP, NODE_W, NODE_H, "cyan", *C["markdup"]))

    for i, (t, s) in enumerate(C["cnv"]):
        o.append(node(LEFT_PX + 20, NODE_Y0 + i * (NODE_H + NODE_GAP), NODE_W, NODE_H,
                      "amber", t, s))
    for i, (t, s) in enumerate(C["snv"]):
        o.append(node(MID_PX + 20, NODE_Y0 + i * (NODE_H + NODE_GAP), NODE_W, NODE_H,
                      C["snv_phases"][i], t, s))
    for i, (t, s) in enumerate(C["sv"]):
        o.append(node(RIGHT_PX + 20, NODE_Y0 + i * (NODE_H + NODE_GAP), NODE_W, NODE_H,
                      "purple", t, s))

    o.append(node(290, Y_VIZ, 900, 88, "rose", *C["viz"]))
    o.append(node(290, Y_SNAKE, 900, 88, "dark", *C["snake"]))

    # ── 图例 ─────────────────────────────────────────────
    items = C["legend"]
    step = 132
    x0 = MID_CX - step * len(items) / 2
    for i, (lab, ph) in enumerate(items):
        c = PHASE[ph]
        bx = x0 + i * step
        o.append(rect(bx, Y_LEGEND - 11, 22, 22, fill=c["node"], stroke=c["border"],
                      sw=1.6, rx=6))
        o.append(f'<rect x="{bx + 6}" y="{Y_LEGEND - 5}" width="10" height="10" rx="3" '
                 f'fill="{c["accent"]}"/>')
        o.append(text(bx + 30, Y_LEGEND + 5, lab, size=14, weight="500", fill="#475569",
                      anchor="start"))

    o.append("</svg>")
    return "\n".join(o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    help="输出目录（默认仓库根目录）")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for lang in ("zh", "en"):
        p = os.path.join(a.out, CONTENT[lang]["out"])
        with open(p, "w", encoding="utf-8") as f:
            f.write(build(lang))
        print(f"✔ {p}  ({os.path.getsize(p)} bytes)")


if __name__ == "__main__":
    main()
