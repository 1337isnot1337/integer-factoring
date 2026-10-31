#!/usr/bin/env python3
"""Actual finite tensor-power jet ranks, beyond the sufficient Smith bound.

This searches one retained witness per cubic Smith signature. That is not
a complete search of coefficient tensors with those signatures.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from math import log2
from pathlib import Path

from check_multipoint_hermite_crossings import (
    flatten,interpolation_weights,local_polynomial,rank,tensor_power,
)
from search_multipoint_rank_one_seeds import KEYS,SWAP,pure_tensors,representatives

ROOT = Path(__file__).resolve().parents[1]


def schedules(copies):
    total = 3*copies+1
    result = []
    for first in range(total+1):
        for second in range(total-first+1):
            counts = (first,second,total-first-second)
            weights = interpolation_weights(total-1,[0,2,None],counts,3)
            selected,offset = [],0
            for place,count in enumerate(counts):
                selected.extend((place,k) for k in range(count) if weights[offset+k])
                offset += count
            result.append((counts,weights,selected))
    return result


def analyze(seed,copies,plans):
    table = tensor_power(seed,2,copies,3)
    degree = 3*copies
    locals_ = [{key:local_polynomial(value,point,degree,3)
                for key,value in table.items()} for point in (0,2,None)]
    ranks = [[[rank(flatten(local,2**copies,cut,k),3) for cut in (0,1)]
              for k in range(degree+1)] for local in locals_]
    best = None
    for counts,weights,selected in plans:
        maximum = max(ranks[p][k][0]*ranks[p][k][1] for p,k in selected)
        width = max(sum((r-1).bit_length() for r in ranks[p][k])
                    for p,k in selected if all(ranks[p][k]))
        pieces = sum(bool(ranks[p][k][0]) for p,k in selected)
        value = maximum,width,pieces,counts
        if best is None or value < best[0]:
            best = value,weights
    (maximum,width,pieces,counts),weights = best
    for labels in product(range(2**copies),repeat=4):
        x,y,u,v = labels
        expected = int(x == v and y == u)
        assert sum(table.get(labels,())) % 3 == expected
        reconstructed,offset = 0,0
        for local,count in zip(locals_,counts):
            values = local.get(labels,[0]*(degree+1))
            reconstructed += sum(weights[offset+k]*values[k] for k in range(count))
            offset += count
        assert reconstructed % 3 == expected
    return {'copies':copies,'complete_jet_cut_ranks':ranks,
            'best_rank_product':maximum,'rounded_width_sum_at_best_product':width,
            'nonzero_pieces':pieces,'multiplicities':counts,'weights':weights,
            'finite_normalized_log_rank_product':log2(maximum)/(4*copies),
            'complete_reconstructed_entries':2**(4*copies)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--copies',type=int,nargs='+',default=[1,2])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    assert all(1 <= copies <= 4 for copies in args.copies)
    catalogue = json.loads((ROOT/'evidence/multipoint-pure-cubics.json').read_text())
    tensors,reps = pure_tensors(),representatives()
    plans = {copies:schedules(copies) for copies in args.copies}
    records = []
    for number,record in enumerate(catalogue['catalogue']):
        ai,bi,ci = record['indices']
        a,b,c = reps[ai],tensors[bi],tensors[ci]
        seed = {key:(x,(y-e-z)%3,(-x-y-e)%3,z)
                for key,x,y,z,e in zip(KEYS,a,b,c,SWAP)}
        records.append({'signature_index':number,'indices':record['indices'],
                        'results':[analyze(seed,m,plans[m]) for m in args.copies]})
    result = {'scope': 'Actual finite jets of one witness per retained cubic Smith signature; no asymptotic bound or complete coefficient search.',
              'witnesses':len(records),'copies':args.copies,'records':records}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for copies in args.copies:
        histogram = Counter(record['results'][args.copies.index(copies)]['best_rank_product']
                            for record in records)
        print(json.dumps({'copies':copies,'best_rank_product_histogram':dict(sorted(histogram.items())),
                          'best_finite_normalized_log_rank_product':log2(min(histogram))/(4*copies)}))


if __name__ == '__main__':
    main()
