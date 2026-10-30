#!/usr/bin/env python3
"""Recompute exact F2 coefficient detectors for 7500 twisted-pairing catalogue choices."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import time


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    d = 5
    assert 24**2 < 2*17**2 and 97**17 < 10**17*5**24
    assert 97*13 > 125*10
    assert 10**2 > 2*7**2 and 11**7 > 5**10
    identity = [[int(i == j) for j in range(d)] for i in range(d)]
    shift = [[int(j == (i+1) % d) for j in range(d)] for i in range(d)]
    diagonal = [[int(i == j and i < 2) for j in range(d)] for i in range(d)]
    alternating = [[int((i,j) in [(0,1),(1,0),(2,3),(3,2)])
                    for j in range(d)] for i in range(d)]
    matrices = [identity, shift,
                [[a^b for a,b in zip(row_a,row_b)] for row_a,row_b in zip(identity,shift)],
                diagonal, alternating]
    names = ['I','P','I+P','D','A']
    pairs = list(product(range(d), repeat=2))

    def cuts(entry):
        return ([sum(entry(x,y,u,v) << (d*u+v) for u,v in pairs) for x,y in pairs],
                [sum(entry(x,y,u,v) << (d*y+v) for y,v in pairs) for x,u in pairs])

    crossing = cuts(lambda x,y,u,v: int(x == v and y == u))
    rows = [cuts(lambda x,y,u,v,a=a,b=b: matrices[a][x][y]*matrices[b][u][v])
            for a,b in pairs]
    columns = [cuts(lambda x,y,u,v,a=a,b=b: matrices[a][x][u]*matrices[b][y][v])
               for a,b in pairs]
    records = []
    rank_stream_hash = hashlib.sha256()
    witness = None
    for row_a,row_b in combinations(range(25),2):
        for column in range(25):
            bases = [crossing,rows[row_a],rows[row_b],columns[column]]
            cut_ranks = []
            for coefficients in range(1,16):
                ranks = []
                for cut in (0,1):
                    matrix_rows = [0]*25
                    for index in range(4):
                        if (coefficients >> index) & 1:
                            matrix_rows = [a^b for a,b in zip(matrix_rows,bases[index][cut])]
                    ranks.append(rank(matrix_rows))
                assert ranks[0] > 0 and ranks[1] > 0  # Basis independence over F2.
                cut_ranks.append(ranks)
                if coefficients & 1:
                    assert ranks[0]*ranks[1] >= 121
                    if ranks[0]*ranks[1] == 121 and witness is None:
                        witness = {'indices':[row_a,row_b,column],
                                   'coefficients':[(coefficients >> j) & 1 for j in range(4)],
                                   'cut_ranks':ranks}
            minimum = min(a*b for coefficient,(a,b) in enumerate(cut_ranks,1) if coefficient & 1)
            records.append([row_a,row_b,column,minimum])
            rank_stream_hash.update(json.dumps([row_a,row_b,column,cut_ranks],
                                               separators=(',',':')).encode())
    assert len(records) == 7500
    assert len({tuple(record[:3]) for record in records}) == 7500
    assert witness is not None
    assert min(record[3] for record in records) == 121
    result = {
        'field_size':2, 'alphabet':5, 'catalogue':names,
        'matrix_index_order':'lexicographic (A,B), each index0..4',
        'records_format':['first_R_index','second_R_index','C_index','minimum_E_nonzero_rank_product'],
        'catalogue_records':records, 'distinct_catalogue_choices':len(records),
        'independent_bases':len(records), 'nonzero_tensor_evaluations':15*len(records),
        'cut_rank_evaluations':30*len(records),
        'minimum_product_distribution':dict(sorted(Counter(record[3] for record in records).items())),
        'minimum_observed_E_rank_product':121, 'attaining_witness':witness,
        'uniform_coefficient_detector_weight':11,
        'rank_table_stream_sha256':rank_stream_hash.hexdigest(),
        'universal_lower_base':'125/13', 'consistent_target_upper_base':'97/10',
        'scope':'The specified finite F2 catalogue only; individual-word observable spectral transfers, every finite auxiliary dimension.',
    }
    if args.output:
        args.output.write_text(json.dumps(result,separators=(',',':'),sort_keys=True)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key != 'catalogue_records'},indent=2))
    print('seconds:',time.monotonic()-started)


if __name__ == '__main__':
    main()
