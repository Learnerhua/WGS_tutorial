#!/usr/bin/env python3
"""13.validate_paths.py — 校验 13.config.yaml 里所有路径"""
import os, sys, yaml

WS = '/Nas3Data/users/ouyjh/WorkSpace/WGS_tutorial'
os.chdir(WS)

CFG = os.path.join(WS, 'scripts/13.config.yaml')
with open(CFG) as fh:
    cfg = yaml.safe_load(fh)

G1 = [('ref', cfg['ref']),
      ('ref_fai', cfg['ref_fai'])]

G2 = [('bwa',          cfg['bwa']),
      ('samtools',     cfg['samtools']),
      ('gatk',         cfg['gatk']),
      ('bcftools',     cfg['bcftools']),
      ('delly',        cfg['delly']),
      ('survivor',     cfg['survivor']),
      ('cnvkit',       cfg['cnvkit']),
      ('vep_perl',     cfg['vep_perl']),
      ('vep_bin',      cfg['vep_bin']),
      ('igv_reports',  cfg['igv_reports']),
      ('manta/configManta.py', os.path.join(cfg['manta_dir'], 'bin', 'configManta.py'))]

G3 = [('mills_indels', cfg['mills_indels']),
      ('known_indels', cfg['known_indels']),
      ('dgv_chr20',    cfg['dgv_chr20']),
      ('gnomad_sv',    cfg['gnomad_sv']),
      ('vep_cache',    cfg['vep_cache']),
      ('excl_regions', cfg['excl_regions']),
      ('genes_bed',    cfg['genes_bed'])]

RED  = '\033[0;31m'; GRN = '\033[0;32m'; YEL = '\033[1;33m'; NC = '\033[0m'

def fmt_sz(n):
    for u in ('B','K','M','G','T'):
        if n < 1024: return f'{n:.1f}{u}'
        n /= 1024
    return f'{n:.1f}P'

def show(label, path, kind='file'):
    full = path if os.path.isabs(path) else os.path.join(WS, path)
    if kind == 'dir':
        ok = os.path.isdir(full)
    elif kind == 'exec':
        ok = os.path.isfile(full)
    else:
        ok = os.path.exists(full)
    if ok:
        if kind in ('file', 'exec') and os.path.isfile(full):
            sz = os.path.getsize(full)
            print(f'  {GRN}✅{NC} {label:30s} {path}  ({fmt_sz(sz)})')
        else:
            print(f'  {GRN}✅{NC} {label:30s} {path}')
        return True
    else:
        print(f'  {RED}❌{NC} {label:30s} {path}  {YEL}(不存在){NC}')
        return False

miss = 0
print('  --- 组 1：样本与参考 ---')
for lbl, p in G1:
    if not show(lbl, p, 'file'): miss += 1

print('  --- 组 2：工具绝对路径 ---')
for lbl, p in G2:
    if not show(lbl, p, 'exec'): miss += 1

print('  --- 组 3：已下载资源 ---')
for lbl, p in G3:
    if not show(lbl, p, 'file'): miss += 1

# vep_cache 是 dir
if not show('vep_cache (dir)',  cfg['vep_cache'], 'dir'): miss += 1
# manta_dir 是 dir
if not show('manta_dir (dir)',   cfg['manta_dir'], 'dir'): miss += 1

print()
if miss > 0:
    print(f'{RED}❌ 共 {miss} 处缺失或不可用{NC}')
    print('  请在 config.yaml 中修正路径，或下载缺失的资源')
    sys.exit(1)
print(f'{GRN}✅ 全部 {len(G1)+len(G2)+len(G3)} 处路径就绪{NC}')
