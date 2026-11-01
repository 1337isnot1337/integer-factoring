#!/usr/bin/env python3
"""Exact paired-frame, local inverse, and independent word-kernel controls."""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
from random import Random

from bounded_spin_crossing import paths
from paired_path_crossing import MULTIPLY as M,power,splitting_blocks,path_entry,polynomial_entry
from check_fixed_field_descent import Field

ROOT = Path(__file__).resolve().parents[1]


def dot(a,b):
    result = 0
    for x,y in zip(a,b):
        result ^= M[x][y]
    return result


def word_polynomial(four_paths,q):
    n = len(four_paths[0])
    words = [w for w in range(1 << n) if w.bit_count() == n//2]
    frame = [[path_entry(p,w,q,dual=i >= 2) for w in words]
             for i,p in enumerate(four_paths)]
    groups = defaultdict(list)
    for i,x in enumerate(words):
        for j,y in enumerate(words):
            groups[x & y,x ^ y].append((i,j))
    result = [0]*(n+1)
    for (_,difference),pairs in groups.items():
        s = difference.bit_count()
        for i,j in pairs:
            primal = M[frame[0][i]][frame[1][j]]
            if not primal:
                continue
            for v,u in pairs:
                scalar = M[primal][M[frame[2][u]][frame[3][v]]]
                if not scalar:
                    continue
                b = (words[i] ^ words[u]).bit_count()
                for k in range(s-b+1):
                    if comb(s-b,k)&1:
                        result[b+k] ^= scalar
    return result


def check():
    models,blocks,frame_pairs,entries = [],0,0,0
    for q in range(2,16):
        try:
            splits = splitting_blocks(q)
        except ValueError:
            continue
        models.append(q)
        for primal,inverse in splits.values():
            for i,j in product(range(2),repeat=2):
                assert dot(primal[i],[inverse[k][j] for k in range(2)]) == int(i == j)
            blocks += 1
        for n in (2,4,6,8):
            ps = paths(n)
            words = [w for w in range(1 << n) if w.bit_count() == n//2]
            f = [[path_entry(p,w,q) for w in words] for p in ps]
            d = [[path_entry(p,w,q,True) for w in words] for p in ps]
            for i,j in product(range(len(ps)),repeat=2):
                assert dot(f[i],d[j]) == int(i == j)
                frame_pairs += 1
        random = Random(q)
        for n in (2,4,6,8):
            ps = paths(n)
            labels = [[random.randrange(len(ps)) for _ in range(4)] for _ in range(3)]
            labels += [[0,len(ps)-1,len(ps)-1,0],[0,0,0,0]]
            for ids in labels:
                four = tuple(ps[i] for i in ids)
                actual = polynomial_entry(four,q)
                assert actual == word_polynomial(four,q)
                assert not any(actual[1::2])
                total = 0
                for value in actual:
                    total ^= value
                assert total == int(ids[0] == ids[3] and ids[1] == ids[2])
                entries += 1
    assert len(models) == 12
    assert power(8,5) == 1 and all(power(8,k) != 1 for k in range(1,5))
    certificates = []
    rank = Field(0b10011).rank
    for n in (8,10,12):
        ps = paths(n)
        random = Random(n)
        rows,labels = [],[]
        for i in range(150):
            ids = [random.randrange(len(ps)) for _ in range(4)]
            if i % 3 == 0:
                ids[2:] = ids[1::-1]
            elif i % 3 == 1:
                ids[2:] = ids[:2]
            value = polynomial_entry(tuple(ps[j] for j in ids),8)
            assert not any(value[1::2])
            rows.append(value[::2])
            labels.append(ids)
            if rank(rows) == n//2+1:
                break
        assert rank(rows) == n//2+1
        certificates.append({'length':n,'entry_labels':labels,'span_dimension':n//2+1})
    files = ('tools/paired_path_crossing.py','tools/check_paired_path_crossing.py',
             'tools/check_fourier_multiplication_gates.py','tools/bounded_spin_crossing.py',
             'tools/check_fixed_field_descent.py')
    return {'status':'pass','field':'F2[a]/(a^4+a+1)','parameters':models,
            'local_inverse_blocks':blocks,'biorthogonal_frame_pairs':frame_pairs,
            'independent_word_polynomials':entries,
            'q8_coefficient_span_certificates':certificates,
            'scope':'Uniform paired path crossing construction; no favorable uniform rank decomposition.',
            'source_sha256':{f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
