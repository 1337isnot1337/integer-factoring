#!/usr/bin/env python3
"""Complete short-word coefficient recombinations for the Steiner rank barrier."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_steiner_cubic_crossings import curve, table
from weighted_smith_oracle import polynomial_multiply


def two_word_coefficients(matrix):
    n = len(matrix)
    rows = [[(j, value) for j, value in enumerate(row) if value] for row in matrix]
    result = [[0]*(n*n) for _ in range(7)]
    for r, s in product(range(n), repeat=2):
        for (j, first), (k, second) in product(rows[r], rows[s]):
            polynomial = polynomial_multiply(first, second)
            bit = 1 << (n*j+k)
            while polynomial:
                term = polynomial & -polynomial
                result[term.bit_length()-1][n*r+s] ^= bit
                polynomial ^= term
    return result


def check():
    q, t = 7, 7
    matrix = curve(table(q, 'boolean'))
    column = [[matrix[q*x+y][q*u+v] for y, v in product(range(q), repeat=2)]
              for x, u in product(range(q), repeat=2)]
    records = []
    for length in (1, 2):
        if length == 1:
            sides = [[[sum(((value >> k) & 1) << j for j, value in enumerate(row)) for row in side]
                      for k in range(4)] for side in (matrix, column)]
        else:
            sides = [two_word_coefficients(side) for side in (matrix, column)]
        dimension, size = 3*length+1, len(sides[0][0])
        detected, minimum_row, minimum_column, minimum_product = 0, size+1, size+1, size*size+1
        ranks = []
        for code in range(1, 1 << dimension):
            pair = []
            for side in sides:
                rows = [0]*size
                for k in range(dimension):
                    if code >> k & 1:
                        rows = [a ^ b for a, b in zip(rows, side[k])]
                rank = binary_rank(rows)
                if code == 1 << (dimension-1):
                    assert binary_rank(rows, low=True) == rank
                pair.append(rank)
            ranks.append((code, *pair))
            if code >> (dimension-1) & 1:
                detected += 1
                r, c = pair
                assert r >= (2*t)**length and c >= (4*t)**length
                minimum_row, minimum_column = min(minimum_row, r), min(minimum_column, c)
                minimum_product = min(minimum_product, r*c)
        assert minimum_row == (2*t)**length and minimum_column == (4*t)**length
        assert minimum_product == (8*t*t)**length
        # Record the complete univariate coefficient span's threshold as a
        # separate finite observation; the theorem covers multivariate words.
        pivots = {}
        def reduce(code):
            while code:
                bit = code.bit_length()-1
                if bit not in pivots:
                    break
                code ^= pivots[bit]
            return code
        for code, r, c in sorted(ranks, key=lambda a: a[1]*a[2]):
            residual = reduce(code)
            if residual:
                pivots[residual.bit_length()-1] = residual
            if reduce((1 << dimension)-1) == 0:
                threshold = r*c
                break
        records.append(dict(word_length=length, physical_alphabet=q**length,
                            complete_nonzero_binary_recombinations=(1 << dimension)-1,
                            detected_recombinations=detected,
                            minimum_detected_row_rank=minimum_row,
                            minimum_detected_column_rank=minimum_column,
                            minimum_detected_product=minimum_product,
                            exact_univariate_spanning_threshold_squared=threshold))
    assert 8*t*t > q**3
    return dict(records=records, exact_seven_letter_comparison={'8T_squared': 8*t*t, 'q_cubed': q**3},
                scope='Finite complete binary univariate controls; the written block proof covers arbitrary multivariate recombinations and characteristic-two extensions.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
