#!/usr/bin/env python3
"""Complete small block-braid coefficient tensors and both actual cut profiles."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from bounded_spin_crossing import paths
from check_fixed_field_descent import Field
from quantum_block_braid import polynomial_entry

ROOT = Path(__file__).resolve().parents[1]


def check():
    rank = Field(0b10011).rank
    records = []
    for n,expected_dimension in ((4,5),(6,33)):
        ps = paths(n)
        d = len(ps)
        tensor = {}
        for ids in product(range(d),repeat=4):
            values = polynomial_entry(tuple(ps[i] for i in ids))
            assert values[0] == int(ids[0] == ids[2] and ids[1] == ids[3])
            total = 0
            for value in values:
                total ^= value
            assert total == int(ids[0] == ids[3] and ids[1] == ids[2])
            tensor[ids] = values
        coefficients = [[tensor[ids][j] for ids in product(range(d),repeat=4)]
                        for j in range(n*n+2)]
        assert rank(coefficients) == expected_dimension
        pairs = []
        for j in range(n*n+2):
            row = [[tensor[p,q,r,s][j] for r,s in product(range(d),repeat=2)]
                   for p,q in product(range(d),repeat=2)]
            col = [[tensor[p,q,r,s][j] for q,s in product(range(d),repeat=2)]
                   for p,r in product(range(d),repeat=2)]
            pairs.append([rank(row),rank(col)])
        records.append({'length':n,'labels':d,'complete_tensor_entries':d**4,
                        'coefficient_span_dimension':expected_dimension,'rank_pairs':pairs})
    files = ('tools/check_quantum_block_braid_profile.py','tools/quantum_block_braid.py',
             'tools/paired_path_crossing.py','tools/bounded_spin_crossing.py',
             'tools/bounded_spin_support.py','tools/check_fourier_multiplication_gates.py',
             'tools/check_fixed_field_descent.py')
    return {'status':'pass','field':'F16','q':8,'profiles':records,
            'scope':'Complete finite coefficient profiles, without a rank-sum optimum or a uniform favorable rank bound.',
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
