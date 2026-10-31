#!/usr/bin/env python3
"""Complete F3 quadratic crossings with rank-one values at two places.

The endpoint is Swap at one.  Both cut ranks at zero and minus one are
at most one, so each nonzero value is a pure four-factor tensor.  The
Swap stabilizer reduces the first value to five nonzero orbits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import flatten, local_polynomial


KEYS = tuple(product(range(2), repeat=4))
SWAP = tuple(int(x == v and y == u) for x, y, u, v in KEYS)


def convolution(a, b, precision, prime=3):
    out = [0] * precision
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:precision-i]):
                out[i+j] = (out[i+j] + x*y) % prime
    return out


def smith_orders(matrix, degree, prime=3):
    """Minimum-order pivots in F_p[[t]], with determinant degree guard."""
    n = len(matrix)
    precision = n * degree + 1
    work = [[list(a) + [0] * (precision-len(a)) for a in row] for row in matrix]
    orders = []
    for k in range(n):
        candidates = [(next((j for j, a in enumerate(work[r][c]) if a), precision), r, c)
                      for r in range(k, n) for c in range(k, n)]
        order, r, c = min(candidates)
        assert order < precision, 'endpoint proves generic invertibility'
        orders.append(order)
        work[k], work[r] = work[r], work[k]
        for row in work:
            row[k], row[c] = row[c], row[k]
        unit = work[k][k][order:] + [0] * order
        inverse = [pow(unit[0], -1, prime)] + [0] * (precision-1)
        for i in range(1, precision):
            inverse[i] = -inverse[0] * sum(unit[j]*inverse[i-j]
                                          for j in range(1, i+1)) % prime
        for j in range(k+1, n):
            work[k][j] = convolution(work[k][j], inverse, precision, prime)
        for i in range(k+1, n):
            factor = work[i][k][order:] + [0] * order
            for j in range(k+1, n):
                correction = convolution(factor, work[k][j], precision, prime)
                work[i][j] = [(a-b) % prime for a, b in zip(work[i][j], correction)]
    assert orders == sorted(orders)
    assert sum(orders) <= n*degree
    return tuple(orders)


def seed_from_values(first, second):
    # Lagrange interpolation at 0,1,2 in F3.
    return {key: (a, (b-e) % 3, (-a-b-e) % 3)
            for key, a, b, e in zip(KEYS, first, second, SWAP)}


def profiles(seed):
    out = []
    for point in (0, 2, None):
        local = {key: local_polynomial(value, point, 2, 3)
                 for key, value in seed.items()}
        cuts = []
        for cut in (0, 1):
            coefficients = [flatten(local, 2, cut, k) for k in range(3)]
            matrix = [[tuple(coefficients[k][i][j] for k in range(3))
                       for j in range(4)] for i in range(4)]
            cuts.append(smith_orders(matrix, 2))
        out.append(tuple(cuts))
    return tuple(out)


def pure_tensors():
    directions = ((1, 0), (0, 1), (1, 1), (1, 2))
    answer = {(0,) * 16}
    for vectors in product(directions, repeat=4):
        for scalar in (1, 2):
            answer.add(tuple(scalar * vectors[0][x] * vectors[1][y]
                             * vectors[2][u] * vectors[3][v] % 3
                             for x, y, u, v in KEYS))
    assert len(answer) == 513
    return sorted(answer)


def representatives():
    answer = [(0,) * 16]
    # Normalize x,y factors to e0. The paired u,v factor is either
    # e0 (nonzero pairing) or e1 (zero pairing). Only when both
    # pairings are nonzero is the total scalar an invariant.
    for u, v in product(range(2), repeat=2):
        for scalar in ((1, 2) if (u, v) == (0, 0) else (1,)):
            answer.append(tuple(scalar * int(key == (0, 0, u, v)) for key in KEYS))
    assert len(answer) == 6
    return answer


def uniform_prefix_threshold(pair, rank_product=8):
    options = []
    for k, l in product(range(1, 5), repeat=2):
        if k*l >= rank_product:
            threshold = max(Fraction(sum(pair[0][:k]), k),
                            Fraction(sum(pair[1][:l]), l))
            options.append((threshold, k, l))
    return min(options)


def run():
    tensors, reps = pure_tensors(), representatives()
    signatures = Counter()
    witnesses = {}
    for first in reps:
        for second in tensors:
            seed = seed_from_values(first, second)
            signature = profiles(seed)
            signatures[signature] += 1
            witnesses.setdefault(signature, (first, second))
    pairs = {pair for sig in signatures for pair in sig}
    thresholds = {pair: uniform_prefix_threshold(pair) for pair in pairs}
    failures = [sig for sig in signatures if sum(thresholds[pair][0] for pair in sig) > 2]
    return {'scope': 'All F3 degree-at-most-two crossings with Swap endpoint and both cut ranks at most one at 0 and 2.',
            'representatives': len(reps), 'second_values': len(tensors),
            'profiled_seeds': sum(signatures.values()), 'signatures': len(signatures),
            'local_profile_pairs': len(pairs),
            'uniform_prefix_rho_floor': '3/4' if not failures else None,
            'prefix_certificate_failures': len(failures),
            'catalogue': [{'profiles': sig, 'count': count,
                           'values': witnesses[sig],
                           'prefix_certificates': [(str(thresholds[pair][0]),
                                                    *thresholds[pair][1:]) for pair in sig]}
                          for sig,count in sorted(signatures.items())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key,value in result.items() if key != 'catalogue'}, indent=2))


if __name__ == '__main__':
    main()
