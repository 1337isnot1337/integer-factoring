#!/usr/bin/env python3
"""Normalised first hit z = j*dens_N grouped by floor(log2 p), plus hit clumping.
Usage: python3 small_p.py data/scan_alltargets.jsonl data/scan_tail.jsonl"""
import json, math, statistics as st, sys
from collections import defaultdict
b = defaultdict(list); adj = defaultdict(lambda: [0, 0])
for path in sys.argv[1:]:
    for line in open(path):
        r = json.loads(line)
        d = sum(len(t['hits']) for t in r['targets']) / sum(t['scanned'] for t in r['targets'])
        for t in r['targets']:
            k = min(int(math.log2(t['p'])), 11)
            if t['hits']:
                b[k].append(t['hits'][0] * d)
            if len(t['hits']) >= 2:
                adj[k][0] += 1; adj[k][1] += (t['hits'][1] - t['hits'][0] == 1)
print('| floor(log2 p) | pairs | mean z | P(z>6) | P(2nd hit is next step) |')
print('|---|---|---|---|---|')
for k in sorted(b):
    v = b[k]
    print(f'| {k}{"+" if k == 11 else ""} | {len(v)} | {st.mean(v):.3f} | {sum(x > 6 for x in v)/len(v):.4f} | {adj[k][1]/adj[k][0]:.3f} |')
