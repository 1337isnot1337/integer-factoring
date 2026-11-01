#!/usr/bin/env python3
"""Primary-character code certificate for all unprojected fused tensor powers."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_fused_scalar_projection import coefficient_matrices
from f4_lattice_linear import f16_multiply

ROOT = Path(__file__).resolve().parents[1]
MUL = [[f16_multiply(a, b) for b in range(16)] for a in range(16)]
INV = [0]+[next(b for b in range(1, 16) if MUL[a][b] == 1) for a in range(1, 16)]
COMPONENTS = (
    ((0, 0, 1), ((0, 0, 1, 0, 0, 0, 0, 0, 0),)),
    ((0, 12, 13), ((0, 1, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0))),
    ((1, 0, 0), ((0, 0, 0, 0, 0, 0, 1, 0, 0),)),
    ((7, 5, 3), ((1, 0, 0, 0, 14, 0, 0, 0, 1),)),
    ((11, 0, 10), ((14, 0, 0, 0, 1, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0, 0, 0, 1))),
    ((12, 13, 0), ((0, 0, 0, 1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 1, 0))),
)


def eliminate(matrix):
    matrix = [list(row) for row in matrix]
    pivots = []
    for c in range(len(matrix[0])):
        p = next((i for i in range(len(pivots), len(matrix)) if matrix[i][c]), None)
        if p is None:
            continue
        r = len(pivots)
        matrix[r], matrix[p] = matrix[p], matrix[r]
        inverse = INV[matrix[r][c]]
        matrix[r] = [MUL[inverse][v] for v in matrix[r]]
        for i in range(len(matrix)):
            if i != r and matrix[i][c]:
                scalar = matrix[i][c]
                matrix[i] = [a ^ MUL[scalar][b] for a, b in zip(matrix[i], matrix[r])]
        pivots.append(c)
        if len(pivots) == len(matrix):
            break
    return matrix, pivots


def multiply(first, second):
    result = [[0]*len(second[0]) for _ in first]
    for r, row in enumerate(first):
        for k, scalar in enumerate(row):
            for c, value in enumerate(second[k]):
                result[r][c] ^= MUL[scalar][value]
    return result


def certificate():
    _, operators = coefficient_matrices()
    vectors = [v for _, basis in COMPONENTS for v in basis]
    change = [list(row) for row in zip(*vectors)]
    identity = [[int(i == j) for j in range(9)] for i in range(9)]
    reduced, pivots = eliminate([a+b for a, b in zip(change, identity)])
    assert pivots == list(range(9))
    inverse = [row[9:] for row in reduced]
    assert multiply(change, inverse) == identity
    transformed = [multiply(multiply(inverse, matrix), change) for matrix in operators]
    widths = [len(basis) for _, basis in COMPONENTS]
    ownership = [i for i, width in enumerate(widths) for _ in range(width)]
    assert all(matrix[r][c] == 0 for matrix in transformed
               for r, c in product(range(9), repeat=2) if ownership[r] != ownership[c])
    offset, nilpotent_ranks = 0, []
    for (character, _), width in zip(COMPONENTS, widths):
        radical = []
        for scalar, matrix in zip(character, transformed):
            radical.append([[matrix[offset+r][offset+c] ^ (scalar if r == c else 0)
                             for c in range(width)] for r in range(width)])
        zero = [[0]*width for _ in range(width)]
        assert all(multiply(a, b) == zero for a, b in product(radical, repeat=2))
        nilpotent_ranks.append([len(eliminate(matrix)[1]) for matrix in radical])
        offset += width
    characters = [list(character) for character, _ in COMPONENTS]
    # Every possible support of total weight < 3 has an injective complement.
    supports = []
    for mask in range(1 << 6):
        if sum(widths[i] for i in range(6) if mask >> i & 1) >= 3:
            continue
        complement = [row for i, row in enumerate(characters) if not (mask >> i & 1)]
        rank = len(eliminate(complement)[1])
        assert rank == 3
        supports.append(dict(support=[i for i in range(6) if mask >> i & 1],
                             complement_rank=rank))
    assert len(supports) == 10
    witness = [1, 15, 9]
    evaluations = []
    for row in characters:
        value = 0
        for a, b in zip(row, witness):
            value ^= MUL[a][b]
        evaluations.append(value)
    assert sum(width for width, value in zip(widths, evaluations) if value) == 3
    names = ('check_fused_primary_code.py', 'check_fused_scalar_projection.py',
             'f4_lattice_linear.py')
    return dict(physical_alphabet=3, primary_characters=characters, primary_dimensions=widths,
                primary_basis=vectors, square_zero_radical_pair_products=54,
                primary_radical_ranks=nilpotent_ranks, weighted_code_distance=3,
                low_weight_support_certificates=supports,
                minimum_word=dict(coefficients=witness, evaluations=evaluations),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in names},
                scope='Every positive tensor length and every field extension of the fixed F16 encoding, for the full unprojected coefficient tensor-product span.',
                conclusion='Every nonzero piece has opposite-cut rank at least 3^m; every crossing decomposition in that span has rank-sum cost at least 3^(3m/2), giving rho at least 3/4.',
                limit='Entrywise scalar projections, physical type restrictions and other crossings are not excluded.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/fused-primary-code.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'primary_basis')}))


if __name__ == '__main__':
    main()
