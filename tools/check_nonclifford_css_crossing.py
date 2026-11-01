#!/usr/bin/env python3
"""Independent complete polynomials and untruncated orders for the fixed CSS twist."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import time

import numpy as np

from native_matrix_smith_f4 import MatrixOracle
from nonclifford_css_crossing import construct
from physical_split_crossings import F4_PRODUCT
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply, scalar_multiply
from windowed_f4_smith import pack

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    begun = time.monotonic()
    planes,metadata,original,changed = construct()
    O = metadata['rotation']
    original = [[pack((int(original[0][i,j]),int(original[1][i,j]))) for j in range(9)] for i in range(9)]
    changed_packed = [[pack((int(changed[0][i,j]),int(changed[1][i,j]))) for j in range(9)] for i in range(9)]
    K = [[F4_PRODUCT[O[i//3][j//3]][O[i % 3][j % 3]] for j in range(9)] for i in range(9)]
    for i,j in product(range(9),repeat=2):
        value = 0
        for a,b in product(range(9),repeat=2):
            value ^= scalar_multiply(original[a][b],F4_PRODUCT[K[a][i]][K[b][j]])
        assert value == changed_packed[i][j]
    # Literal F3 orbits and packed scalar polynomial arithmetic, separate from
    # the vectorized field-plane constructor.
    cosets,occupied = [],set()
    for word in product(range(3),repeat=4):
        if (word[0]+word[1]+2*word[2]+2*word[3]) % 3 or word in occupied:
            continue
        orbit = [tuple((s+t) % 3 for s in word) for t in range(3)]
        cosets.append(orbit)
        occupied.update(orbit)
    assert len(cosets) == 9
    checked,summands = 0,0
    for x,y,v,u in product(range(9),repeat=4):
        result = 0
        for wx,wy,wv,wu in product(cosets[x],cosets[y],cosets[v],cosets[u]):
            summands += 1
            values = [(changed_packed if site == 0 else original)[3*wx[site]+wy[site]][3*wv[site]+wu[site]]
                      for site in range(4)]
            if all(values):
                value = 1
                for a in values:
                    value = polynomial_multiply(value,a)
                result ^= value
        assert result == pack((int(planes[0][9*x+y,9*v+u]),int(planes[1][9*x+y,9*v+u])))
        checked += 1
    profiles = []
    randomizer = random.Random(20261101)
    for infinity in (False,True):
        p = ([a.reshape(9,9,9,9).transpose(0,3,2,1).reshape(81,81) for a in planes]
             if infinity else planes)
        matrix = [[(int(p[0][i,j]),int(p[1][i,j])) for j in range(81)] for i in range(81)]
        packed = [[pack(value) for value in row] for row in matrix]
        row = [randomizer.randrange(-2,3) for _ in range(81)]
        column = [randomizer.randrange(-2,3) for _ in range(81)]
        exact,_ = WeightedSmithOracle(packed,infinity=infinity).profile(row,column,1)
        actual = MatrixOracle(matrix,81).profile(row,column,infinity=infinity)['orders']
        assert list(map(int,actual)) == list(exact)
        profiles.append(dict(infinity=infinity,row_weights=row,column_weights=column,orders=exact))
    names = ['tools/check_nonclifford_css_crossing.py','tools/nonclifford_css_crossing.py',
             'tools/collective_css_pauli.py','tools/cyclic_pauli_crossings.py',
             'tools/weighted_smith_oracle_f4.py','tools/native_matrix_smith_f4.py',
             'tools/native_matrix_smith_f4.cpp','tools/native_rectangular_smith_f4.cpp']
    result = dict(scope='Exact fixed-code polynomial crossing and independent complete weighted orders; no stronger bound.',
                  metadata=metadata,complete_contracted_entries=checked,literal_coset_summands=summands,
                  independent_twisted_cell_entries=81,untruncated_orders=162,profiles=profiles,
                  source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    result = json.loads(json.dumps(result))
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/nonclifford-css-controls.json').read_text())
    print(json.dumps({k:result[k] for k in ('scope','complete_contracted_entries','literal_coset_summands','untruncated_orders')}))
    print(json.dumps({'seconds':time.monotonic()-begun,'structural_gate':metadata['structural_gate']}))


if __name__ == '__main__':
    main()
