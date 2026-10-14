#!/usr/bin/env python3
"""重画 Picard insert size 直方图。

Picard 默认 X 轴 = 数据实际范围（~1.1 Mb），导致主分布（200-800 bp）的柱状
条宽度在 1.1 Mb 的 X 轴上只有 ~0.07%，肉眼几乎不可见。

修复：X 轴限制在 0-1000 bp，丢弃 >1000 bp 的离群配对（来自跨染色体配对），
主分布就会清晰显示。
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PICARD_DIR = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial/04.post_alignment_qc/picard'
METRICS = f'{PICARD_DIR}/NA12878_chr20_insert_size_metrics.txt'
OUT_PDF = f'{PICARD_DIR}/NA12878_chr20_insert_size_histogram_fixed.pdf'
OUT_PNG = f'{PICARD_DIR}/NA12878_chr20_insert_size_histogram_fixed.png'

MAX_X = 1000   # X 轴上限；>1000 的配对视为跨染色体离群


def parse_picard_histogram(path):
    """读取直方图段：FR / RF / TANDEM 三组各自的 insert_size 分布。"""
    hist = {'FR': [], 'RF': [], 'TANDEM': []}
    in_hist = False
    with open(path) as f:
        for line in f:
            if line.startswith('## HISTOGRAM'):
                in_hist = True
                continue
            if in_hist:
                if line.startswith('##') or line.startswith('insert_size'):
                    continue
                if not line.strip():
                    break
                parts = line.rstrip('\n').split('\t')
                try:
                    size = int(parts[0])
                    fr = int(parts[1])
                    rf = int(parts[2])
                    tandem = int(parts[3])
                    hist['FR'].append((size, fr))
                    hist['RF'].append((size, rf))
                    hist['TANDEM'].append((size, tandem))
                except (ValueError, IndexError):
                    continue
    return hist


def main():
    hist = parse_picard_histogram(METRICS)

    # 拆分为 numpy 数组（按方向）
    arrays = {}
    for direction, pairs in hist.items():
        if not pairs:
            arrays[direction] = (np.array([]), np.array([]))
            continue
        sizes = np.array([p[0] for p in pairs])
        counts = np.array([p[1] for p in pairs])
        # 只保留 ≤ MAX_X 的部分（FR 主分布）
        mask = sizes <= MAX_X
        arrays[direction] = (sizes[mask], counts[mask])

    # 画图（FR 主分布 + RF/TANDEM 离群对比）
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 左图：FR 线性坐标（主分布细节）
    ax = axes[0]
    sizes, counts = arrays['FR']
    ax.bar(sizes, counts, width=1.0, color='steelblue', edgecolor='none', alpha=0.85)
    ax.axvline(397, color='red', linestyle='--', linewidth=1, label=f'MODE = 397 bp')
    ax.axvline(434, color='orange', linestyle='--', linewidth=1, label=f'MEDIAN = 434 bp')
    ax.set_xlabel('Insert size (bp)', fontsize=11)
    ax.set_ylabel('Number of read pairs', fontsize=11)
    ax.set_title(f'FR read pairs (main distribution)\nX-axis limited to 0-{MAX_X} bp', fontsize=12)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xlim(0, MAX_X)
    # Y 轴顶值给文字标头一点空间
    ymax = counts.max() * 1.15
    ax.set_ylim(0, ymax)

    # 右图：FR + RF + TANDEM log 坐标（全分布对比）
    ax = axes[1]
    for direction, color in [('FR', 'steelblue'), ('RF', 'crimson'), ('TANDEM', 'darkgreen')]:
        sizes, counts = arrays[direction]
        if len(sizes) == 0:
            continue
        # 真实总数（包括 >1000 bp 的离群配对，label 显示真实数量）
        total = sum(p[1] for p in hist[direction])
        # 用 bin=10 聚合显示（避免每个 bin 画一根线）
        bin_size = 10
        bins = np.arange(0, MAX_X + bin_size, bin_size)
        binned, _ = np.histogram(sizes, bins=bins, weights=counts)
        ax.plot(bins[:-1] + bin_size / 2, binned, color=color, linewidth=1.5,
                label=f'{direction} (n={total:,} total; showing <= {MAX_X} bp)')
    ax.set_xlabel('Insert size (bp)', fontsize=11)
    ax.set_ylabel('Number of read pairs (log scale)', fontsize=11)
    ax.set_yscale('log')
    ax.set_title(f'Three pair orientations (FR/RF/TANDEM)\nlog scale to show all three', fontsize=12)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(0, MAX_X)

    plt.suptitle('NA12878 chr20 — Insert size distribution (re-drawn from Picard histogram)',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_PDF, bbox_inches='tight')
    plt.savefig(OUT_PNG, bbox_inches='tight', dpi=150)
    plt.close()

    # 简要统计
    fr_sizes, fr_counts = arrays['FR']
    orig_total = sum(p[1] for p in hist['FR'])
    in_range = sum(p[1] for p in hist['FR'] if p[0] <= MAX_X)
    truncated = orig_total - in_range

    print(f"FR 总 read 对数: {orig_total:,}")
    print(f"  <= {MAX_X} bp: {in_range:,} ({in_range/orig_total*100:.2f}%)")
    print(f"  > {MAX_X} bp (truncated): {truncated:,} ({truncated/orig_total*100:.2f}%)")
    print(f"\n写入: {OUT_PDF}")
    print(f"写入: {OUT_PNG}")
    print(f"\nFR 主分布 (<={MAX_X} bp) 概览:")
    print(f"  峰位: {fr_sizes[fr_counts.argmax()]}")
    print(f"  最大计数值: {fr_counts.max():,}")
    print(f"  median (434 bp) 与 mode (397 bp) 见 matplotlib 虚线")


if __name__ == '__main__':
    main()
