#!/usr/bin/env python3
"""Finite Hermite schedules for the explicit matching-value curves."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import (
    tensor_power,local_polynomial,flatten,rank,interpolation_weights,
)
from matching_curve import matching_curve


def compositions(total,count):
    if count == 1:
        yield (total,)
        return
    for first in range(total+1):
        for rest in compositions(total-first,count-1):
            yield (first,)+rest


def prefix_schedule(ranks,total):
    states = {0:(0,())}
    for row in ranks:
        prefix = [0]
        for a,b in row:
            prefix.append(max(prefix[-1],a*b))
        updated = {}
        for used,(cost,counts) in states.items():
            for h in range(total-used+1):
                candidate = max(cost,prefix[h]),counts+(h,)
                if used+h not in updated or candidate < updated[used+h]:
                    updated[used+h] = candidate
        states = updated
    return states[total][1]


def analyze(seed,prime,copies):
    degree = prime*copies
    points = [a for a in range(prime) if a != 1]+[None]
    table = tensor_power(seed,2,copies,prime)
    locals_ = [{key:local_polynomial(value,point,degree,prime)
                for key,value in table.items()} for point in points]
    ranks = [[[rank(flatten(local,2**copies,cut,k),prime) for cut in (0,1)]
              for k in range(degree+1)] for local in locals_]
    plans = compositions(degree+1,len(points)) if prime == 5 else [prefix_schedule(ranks,degree+1)]
    best = None
    schedules = 0
    for counts in plans:
        weights = interpolation_weights(degree,points,counts,prime)
        selected,offset = [],0
        for j,count in enumerate(counts):
            selected.extend((j,k) for k in range(count) if weights[offset+k])
            offset += count
        maximum = max(ranks[j][k][0]*ranks[j][k][1] for j,k in selected)
        pieces = sum(bool(ranks[j][k][0]) for j,k in selected)
        candidate = maximum,pieces,counts
        schedules += 1
        if best is None or candidate < best[0]:
            best = candidate,weights
    (maximum,pieces,counts),weights = best
    for labels in product(range(2**copies),repeat=4):
        x,y,u,v = labels
        expected = int(x == v and y == u)
        value,offset = 0,0
        for local,count in zip(locals_,counts):
            poly = local.get(labels,[0]*(degree+1))
            value += sum(weights[offset+k]*poly[k] for k in range(count))
            offset += count
        assert value % prime == expected
    return {'copies':copies,'maximum_rank_product':maximum,'nonzero_pieces':pieces,
            'multiplicities':counts,'weights':weights,'points':points,
            'schedules_checked':schedules,'complete_jet_cut_ranks':ranks,
            'reconstructed_entries':2**(4*copies),
            'schedule_scope':'All minimal Hermite allocations' if prime == 5 else
                             'One prefix-minimax dynamic-programming schedule; not exhaustive after zero interpolation weights'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime',type=int,choices=(5,11),default=5)
    parser.add_argument('--copies',type=int,nargs='+',default=[1,2,3])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    assert all(1 <= m <= 3 for m in args.copies)
    row_points = (0,2,3) if args.prime == 5 else (0,2,3,4,5,6)
    seed,_ = matching_curve(args.prime,row_points)
    records = []
    for copies in args.copies:
        record = analyze(seed,args.prime,copies)
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k not in
                          ('complete_jet_cut_ranks','weights','points')},sort_keys=True),flush=True)
    result = {'prime':args.prime,'degree':args.prime,'records':records,
              'scope':'Finite exact rank and reconstruction data; no uniform improvement.',
              'source_sha256':{name:sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                               for name in ('matching_curve.py','check_matching_curve_jets.py')}}
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n')


if __name__ == '__main__':
    main()
