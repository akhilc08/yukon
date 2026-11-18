#!/usr/bin/env python3
"""hits.py: scan rankg outputs (*/rin_0?.out) for curves that beat a cell target, print mkcand commands for new ones.
Target of cell (tors, r) = min(frozen baseline from records.json, board holder H if held (HELD below)).
Record line: p q H rlo rhi npts nsec ms tors [iso list [[H,[tors]],...]] (torsion printed like [4] or [2, 2])."""
import glob, json, os, re, sys
REC = json.load(open('/work/ec-rank/ec_verify/records.json'))['cells']
HELD = {'3/9': 9.512795, '3/11': 14.74411, '3/12': 17.864161, '4/9': 12.166195, '2x2/9': 9.460889,
        '2x2/12': 18.004452, '2x2/13': 20.25267, '2x2/14': 24.139479}
def tkey(t):
    t = [int(x) for x in re.findall(r'\d+', t)]
    return 'x'.join(map(str, t))
def target(cell):
    b = REC.get(cell, {}).get('record', {}).get('naive_height')
    if b is None: return None
    return min(b, HELD.get(cell, 1e9))
seen = set(open('hits.seen').read().split()) if os.path.exists('hits.seen') else set()
out = []
for f in glob.glob('*/rin_0?.out'):
    job = f.split('/')[0]
    for line in open(f):
        m = re.match(r'(-?\d+) (\d+) (\S+) (\d+) (\d+) \d+ \d+ \S+ (\[[\d, ]+\])\s*(.*)$', line.strip())
        if not m: continue
        p, q, H, rlo = m.group(1), m.group(2), float(m.group(3)), int(m.group(4))
        if rlo < 9: continue
        cands = [(H, tkey(m.group(6)), 0)]
        for h, t in re.findall(r'\[(\d+\.\d+), (\[[\d, ]+\])\]', m.group(7)):
            k = tkey(t); cands.append((float(h), k, 4 if k == '4' else 22))
        for h, k, iso in cands:
            if iso and k == tkey(m.group(6)): iso = 0
            cell = f'{k}/{rlo}'; tg = target(cell)
            if tg is None or h >= tg - 1e-4: continue  # records.json heights are rounded to 1e-6; skip ties (the record itself)
            key = f'{job}:{p}:{q}:{k}'
            if key in seen: continue
            out.append((cell, h, tg, job, p, q, iso, key))
for cell, h, tg, job, p, q, iso, key in sorted(out, key=lambda x: x[1] - x[2]):
    name = f"{cell.replace('/', '-')}-{job}-{p.replace('-', 'm')}-{q}"
    print(f'{key}\t{cell}\t{h:.4f}\t{tg:.4f}\t{job}\t{p}\t{q}\t{iso}\t{name}')
