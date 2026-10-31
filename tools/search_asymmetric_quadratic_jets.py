#!/usr/bin/env python3
"""Actual jets of the retained unequal-determinant quadratic witnesses."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import log2
from pathlib import Path

from check_joint_determinant_barrier import determinant
from check_multipoint_hermite_crossings import flatten
from search_multipoint_actual_jets import analyze,schedules
from search_multipoint_rank_one_seeds import seed_from_values

ROOT = Path(__file__).resolve().parents[1]


def determinants(seed):
    result = []
    for cut in (0,1):
        coefficients = [flatten(seed,2,cut,k) for k in range(3)]
        matrix = [[tuple(coefficients[k][i][j] for k in range(3))
                   for j in range(4)] for i in range(4)]
        result.append(determinant(matrix,3))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--copies',type=int,nargs='+',default=[1,2,3])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    assert all(1 <= m <= 4 for m in args.copies)
    catalogue = ROOT/'evidence/multipoint-rank-two-seeds.json'
    saved = json.loads(catalogue.read_text())
    plans = {m:schedules(m,2) for m in args.copies}
    records,excluded = [],[]
    for index,record in enumerate(saved['catalogue']):
        seed = seed_from_values(*record['values'])
        f,g = determinants(seed)
        if f == g:
            excluded.append(index)
            continue
        records.append({'signature_index':index,
                        'results':[analyze(seed,m,plans[m],2) for m in args.copies]})
    assert len(records) == 195 and len(excluded) == 126
    result = {'scope':'One retained witness per unequal-determinant quadratic profile signature; finite actual ranks and complete Hermite reconstruction, not a uniform rate or complete coefficient search.',
              'input_catalogue_sha256':sha256(catalogue.read_bytes()).hexdigest(),
              'excluded_equal_determinant_signatures':excluded,
              'tested_witnesses':len(records),'copies':args.copies,'records':records,
              'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('search_asymmetric_quadratic_jets.py',
                                            'search_multipoint_actual_jets.py')}}
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    for m in args.copies:
        histogram = Counter(record['results'][args.copies.index(m)]['best_rank_product']
                            for record in records)
        print(json.dumps({'copies':m,'best_rank_product_histogram':dict(sorted(histogram.items())),
                          'best_finite_normalized_log_rank_product':log2(min(histogram))/(4*m)}))


if __name__ == '__main__':
    main()
