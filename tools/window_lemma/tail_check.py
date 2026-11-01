#!/usr/bin/env python3
"""Tail test of the geometric model: z = j * dens_N over every (N, target) pair.
Usage: python3 tail_check.py data/scan_alltargets.jsonl data/scan_tail.jsonl"""
import json, math, sys
z = []
for path in sys.argv[1:]:
    for line in open(path):
        r = json.loads(line)
        hits = sum(len(t['hits']) for t in r['targets']); steps = sum(t['scanned'] for t in r['targets'])
        d = hits / steps
        z += [t['hits'][0] * d for t in r['targets'] if t['hits']]
n = len(z); z.sort()
print(f'pairs={n} max z={z[-1]:.2f} expected max ~ ln(n)+0.577 = {math.log(n)+0.5772:.2f}')
for t in (1, 2, 4, 6, 8, 10):
    obs = sum(1 for v in z if v > t) / n
    print(f'P(z>{t}) observed {obs:.2e}  Exp(1) model {math.exp(-t):.2e}')
