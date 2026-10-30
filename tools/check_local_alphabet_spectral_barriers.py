#!/usr/bin/env python3
"""Exact controls for local-alphabet observable-transfer bounds."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import time

from check_binary_two_piece_f4_boundary import Field, rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    rank_pairs = 0
    for d in range(1, 13):
        n = d*d
        for r, c in product(range(1, n+1), repeat=2):
            assert (n+1)**2*r*c-n*(r+c)**2 == (n*r-c)*(n*c-r) >= 0
            rank_pairs += 1
    direct_checks = []
    for cardinality in (2, 4):
        field = Field(cardinality)
        for d in range(2, 8):
            pairs = list(product(range(d), repeat=2))
            q, n, delta = d*(d-1)//2, d*d, d % 2
            minimum_squared_cost = n*n
            for b, c in product(range(cardinality), repeat=2):
                def tensor(x, y, u, v):
                    value = int(x == v and y == u)
                    value = field.add(value, field.mul(b, int(x == y and u == v)))
                    return field.add(value, field.mul(c, int(x == u and y == v)))
                row = [[tensor(x, y, u, v) for u, v in pairs] for x, y in pairs]
                column = [[tensor(x, y, u, v) for y, v in pairs] for x, u in pairs]
                rr, cr = rank(row, field), rank(column, field)
                loss = delta and field.add(field.add(1, b), c) == 0
                expected_r = q+int(b != 0) if c == 1 else n-int(loss)
                expected_c = q+int(c != 0) if b == 1 else n-int(loss)
                assert (rr, cr) == (expected_r, expected_c)
                assert rr*cr >= (q+1)**2
                if b == c == 1:
                    assert rr == cr == q+1
                minimum_squared_cost = min(minimum_squared_cost, rr*cr)
            assert minimum_squared_cost == (q+1)**2
            direct_checks.append({'field_size': cardinality, 'alphabet': d,
                                  'coefficient_pairs': cardinality**2,
                                  'minimum_squared_cost': minimum_squared_cost})
    assert 256 > 25*8  # (16/5)^2 > 2^3
    assert 729 > 25*27  # (27/5)^2 > 3^3
    assert 2**29 > 17**7
    assert 11**7 > 5**10
    for d in range(6, 101):
        assert (d-1)**2 >= 4*d
    result = {
        'positive_rank_pairs_checked': rank_pairs,
        'direct_flattening_checks': direct_checks,
        'direct_coefficient_cases_total': sum(item['coefficient_pairs'] for item in direct_checks),
        'small_alphabet_rational_bounds': {'2': '16/5', '3': '27/5', '4': '128/17'},
        'exact_integer_comparisons': {'2^29': 2**29, '17^7': 17**7,
                                      '11^7': 11**7, '5^10': 5**10},
        'scope': 'Bounds the observable word-rank spectral method, not arbitrary grouped-piece true ranks.',
        'seconds': time.monotonic()-started,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
