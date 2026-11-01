#!/usr/bin/env python3
"""Full independent CSS contraction and exact Frobenius equivalence control."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import numpy as np

from check_cyclic_pauli_crossings import independent_columns
from collective_css_pauli import code, compress
from cyclic_pauli_crossings import construct, monomial
from native_matrix_smith_f4 import MatrixOracle
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply
from windowed_f4_smith import pack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    begun = time.monotonic()
    planes, metadata = compress()
    cell, cell_metadata = construct(1)
    assert independent_columns(cell, cell_metadata, range(9)) == 81
    cell = [[pack((int(cell[0][i, j]), int(cell[1][i, j])))
             for j in range(9)] for i in range(9)]
    # Rebuild the split from literal F3 orbits; do not reuse code()'s words.
    cosets, occupied = [], set()
    for s in product(range(3), repeat=4):
        if (s[0]+s[1]+2*s[2]+2*s[3]) % 3 or s in occupied:
            continue
        orbit = [tuple((si+t) % 3 for si in s) for t in range(3)]
        cosets.append(orbit)
        occupied.update(orbit)
    assert len(cosets) == 9 and len(occupied) == 27
    words, embedding = code()
    assert cosets == [[tuple(v) for v in orbit] for orbit in words.tolist()]
    checked, summands, nonzero_summands = 0, 0, 0
    for x, y, v, u in product(range(9), repeat=4):
        result = 0
        for wx, wy, wv, wu in product(cosets[x], cosets[y], cosets[v], cosets[u]):
            summands += 1
            values = [cell[3*wx[i]+wy[i]][3*wv[i]+wu[i]] for i in range(4)]
            if not all(values):
                continue
            value = 1
            for a in values:
                value = polynomial_multiply(value, a)
            result ^= value
            nonzero_summands += 1
        assert result == pack((int(planes[0][9*x+y, 9*v+u]),
                               int(planes[1][9*x+y, 9*v+u])))
        checked += 1
    # Check every single-site Weyl compression using independent field sums.
    detections = 0
    for site, shift, phase in product(range(4), range(3), range(3)):
        if shift == phase == 0:
            continue
        a, b = [0]*4, [0]*4
        a[site], b[site] = shift, phase
        permutation, factors = monomial(tuple(a), tuple(b))
        action = np.zeros((81, 81), dtype=np.uint8)
        action[permutation, np.arange(81)] = factors
        for i, j in product(range(9), repeat=2):
            result = 0
            for r in np.flatnonzero(embedding[:, i]):
                for c in np.flatnonzero(embedding[:, j]):
                    result ^= int(action[r, c])
            assert result == 0
        detections += 1
    # T_code(z;x,y,v,u)=sigma(T_cell(z^2;Mx,My,Mv,Mu)^tensor2).
    # Here sigma squares F4 coefficients and M(i,j)=(i,i+j) over F3.
    permutation = [3*(i//3)+(i//3+i % 3) % 3 for i in range(9)]
    equivalences = 0
    for x, y, v, u in product(range(9), repeat=4):
        mx, my, mv, mu = [permutation[i] for i in (x, y, v, u)]
        value = polynomial_multiply(cell[3*(mx//3)+my//3][3*(mv//3)+mu//3],
                                    cell[3*(mx % 3)+my % 3][3*(mv % 3)+mu % 3])
        expected = 0
        for j in range(7):
            c = (value >> (2*j)) & 3
            expected |= (0, 1, 3, 2)[c] << (4*j)
        assert expected == pack((int(planes[0][9*x+y, 9*v+u]),
                                 int(planes[1][9*x+y, 9*v+u])))
        equivalences += 1
    # Independent untruncated profiles, including genuinely unequal heights.
    records = []
    rng = np.random.default_rng(42220261101)
    for infinity in (False, True):
        p = ([a.reshape(9,9,9,9).transpose(0,3,2,1).reshape(81,81)
              for a in planes] if infinity else planes)
        matrix = [[(int(p[0][i,j]), int(p[1][i,j])) for j in range(81)] for i in range(81)]
        packed = [[pack(value) for value in row] for row in matrix]
        native = MatrixOracle(matrix, 81)
        for weighted in (False, True):
            row = rng.integers(-3, 4, 81).tolist() if weighted else [0]*81
            col = rng.integers(-3, 4, 81).tolist() if weighted else [0]*81
            exact, _ = WeightedSmithOracle(packed, infinity=infinity).profile(row, col, 1)
            got = native.profile(row, col, infinity=infinity)['orders']
            assert list(map(int, got)) == list(exact)
            records.append({'infinity': infinity, 'row_weights': row, 'column_weights': col,
                            'orders': exact, 'histogram': dict(sorted(Counter(exact).items()))})
    names = ['collective_css_pauli.py', 'check_collective_css_pauli.py',
             'cyclic_pauli_crossings.py', 'weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'metadata': metadata,
              'independent_cell_entries': 81, 'complete_contracted_entries': checked,
              'literal_summands_examined': summands, 'nonzero_summands': nonzero_summands,
              'single_site_nonidentity_Weyl_annihilations': detections,
              'full_polynomial_equivalence_entries': equivalences,
              'equivalence_permutation': permutation, 'untruncated_orders': 324,
              'profiles': records, 'seconds': time.monotonic()-begun,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Exact finite polynomial identity over F4 and independent complete orders. '
                       'This code is a Frobenius/relabelled two-copy cell and gives no new support family '
                       'even at arbitrary physical heights. No stronger factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'complete_contracted_entries',
                     'full_polynomial_equivalence_entries', 'untruncated_orders', 'seconds')}))


if __name__ == '__main__':
    main()
