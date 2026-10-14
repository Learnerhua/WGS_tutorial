#!/usr/bin/env python3
"""教程结构自检：验证 WGS_tutorial.md 结构规范。

用法：
    python3 scripts/check_tutorial_structure.py

检查项：
  1. 第 4 章每节（4.N）恰好 5 个三级标题，编号连续（4.N.1 – 4.N.5）
  2. 五段命名依次为 目的 / 基本概念 / 方法 / 结果 / 总结
  3. 全文档：所有父标题先于其子标题出现（防父子错位）
  4. 全文档：同一父下的子编号连续（无缺号）
  5. 全文档：无重复编号
  6. 不存在五级标题（##### ）
  7. 无教程口吻残留

退出码：0 = 全部通过，1 = 有问题
"""

import re
import sys
import pathlib
from collections import defaultdict

DOC = pathlib.Path(__file__).resolve().parent.parent / 'WGS_tutorial.md'

EXPECTED = ['目的', '基本概念', '方法', '结果', '总结']
SECTIONS = ['4.1', '4.2', '4.3', '4.4', '4.5', '4.6', '4.7', '4.8', '4.9', '4.10', '4.11', '4.12', '4.13', '4.14', '4.15']
BAD_TONE = ['我们将在', '让你能', '你将在']


def main() -> int:
    text = DOC.read_text(encoding='utf-8')
    lines = text.splitlines()
    problems = []

    # 收集所有编号标题：(行号, 层级, 编号, 标题)
    heads = []
    for idx, line in enumerate(lines, 1):
        m = re.match(r'^(#{2,6}) ((\d+)\.\d+(?:\.\d+)*)(.*)$', line.rstrip())
        if m:
            heads.append((idx, len(m.group(1)), m.group(2), m.group(4).strip()))

    # ── 1 & 2: 第 4 章每节 5 个三级标题 + 命名 ──
    for sec in SECTIONS:
        items = [(n, t) for _, lvl, n, t in heads
                 if lvl == 4 and n.rsplit('.', 1)[0] == sec]
        nums = [n for n, _ in items]
        want = [f'{sec}.{i}' for i in range(1, 6)]
        if nums != want:
            problems.append(f'{sec}: 三级标题编号异常 → {nums}（期望 {want}）')
            continue
        for i, (_, title) in enumerate(items, 1):
            if not title.startswith(EXPECTED[i - 1]):
                problems.append(f'{sec}.{i}: 命名应为「{EXPECTED[i-1]}…」，实际「{title}」')

    # ── 3: 父标题先于子标题 ──
    seen = set()
    for idx, lvl, num, _ in heads:
        parts = num.split('.')
        for k in range(2, len(parts) + 1):
            parent = '.'.join(parts[:k - 1])
            if parent.count('.') >= 1 and parent not in seen:
                problems.append(f'行 {idx}: 「{num}」的父标题「{parent}」尚未出现（父子错位）')
        seen.add(num)

    # ── 4: 同父下子编号连续 ──
    groups = defaultdict(list)
    for idx, lvl, num, _ in heads:
        parts = num.split('.')
        if len(parts) >= 3:
            groups['.'.join(parts[:-1])].append(int(parts[-1]))
    for parent, ks in sorted(groups.items()):
        ks_sorted = sorted(ks)
        if ks_sorted != list(range(1, len(ks_sorted) + 1)):
            problems.append(f'{parent}.x: 子编号不连续 → {ks_sorted}')

    # ── 5: 无重复编号 ──
    cnt = defaultdict(int)
    for _, _, num, _ in heads:
        cnt[num] += 1
    for num, c in cnt.items():
        if c > 1:
            problems.append(f'编号重复：{num}（{c} 次）')

    # ── 6: 无五级标题 ──
    h5 = [i for i, line in enumerate(lines, 1) if line.startswith('##### ')]
    if h5:
        problems.append(f'存在五级标题（##### ）：行 {h5[:10]}')

    # ── 7: 教程口吻 ──
    for bad in BAD_TONE:
        for idx, line in enumerate(lines, 1):
            if bad in line:
                problems.append(f'行 {idx}: 教程口吻残留「{bad}」')

    # ── 输出 ──
    if problems:
        print(f'❌ 发现 {len(problems)} 个问题：')
        for p in problems:
            print(f'   - {p}')
        return 1

    print('✅ 教程结构自检全部通过：')
    print(f'   - 第 4 章 {len(SECTIONS)} 节，每节 5 个三级标题（目的/基本概念/方法/结果/总结）')
    print('   - 全文档无重复编号、无缺号、父标题均先于子标题')
    print('   - 无五级标题、无教程口吻残留')
    print(f'   - 总行数: {len(lines)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
