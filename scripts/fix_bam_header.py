#!/usr/bin/env python3
"""修复 BWA-MEM2 写入 BAM header 的 @PG 行未转义 tab 问题。

背景：
  BWA-MEM2 把 `-R` 传入的 Read Group 字符串原样写进 @PG 行的 CL 字段，
  其中的 tab 未转义。samtools 只报警告（跳过），但 GATK/Picard 用的 htsjdk
  严格解析会直接抛异常：
      SAMFormatException: Problem parsing @PG key:value pair
      ID:xxx clashes with ID:bwa-mem2
  导致所有 GATK 步骤（MarkDuplicates/BQSR/HaplotypeCaller）失败。

修复：
  把每条 @PG 行「CL:」字段内的 tab 替换成空格，再用 samtools reheader 写回。
  只改 header，不触碰比对数据。

用法：
  python3 fix_bam_header.py <in.bam> <out.bam>
"""
import subprocess
import sys
import tempfile
import os

SAMTOOLS = 'samtools'   # 需在 PATH 中（或用绝对路径）


def export_header(bam: str) -> str:
    """导出 header（--no-PG 避免追加本次命令的 @PG 行）。"""
    r = subprocess.run([SAMTOOLS, 'view', '-H', '--no-PG', bam],
                       capture_output=True, text=True, check=True)
    return r.stdout


def sanitize_pg_lines(header: str) -> tuple[str, int]:
    """把 @PG 行 CL: 字段内的 tab 换成空格。返回 (新header, 修复行数)。"""
    out, n = [], 0
    for line in header.split('\n'):
        if line.startswith('@PG') and '\tCL:' in line:
            i = line.find('\tCL:')
            head, cl = line[:i], line[i + 4:]
            if '\t' in cl:
                cl = cl.replace('\t', ' ')
                n += 1
            line = head + '\tCL:' + cl
        out.append(line)
    return '\n'.join(out), n


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    in_bam, out_bam = sys.argv[1], sys.argv[2]

    header = export_header(in_bam)
    fixed, n = sanitize_pg_lines(header)
    print(f"@PG 行修复数: {n}")

    with tempfile.NamedTemporaryFile('w', suffix='.sam', delete=False) as tf:
        tf.write(fixed)
        hdr_path = tf.name
    try:
        with open(out_bam, 'wb') as fo:
            subprocess.run([SAMTOOLS, 'reheader', hdr_path, in_bam],
                           stdout=fo, check=True)
    finally:
        os.unlink(hdr_path)
    print(f"已写出: {out_bam}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
