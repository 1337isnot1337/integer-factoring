#!/usr/bin/env python3
"""Verify local tensor-power profiles and the bounded-weight permutation lemma."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random

from hamming_pfaffian import labels, profile
from hamming_tensor_power import local_terms, matchings, permuted_label, tensor_profile
from modular_pfaffian import polynomial_product

ROOT = Path(__file__).resolve().parents[1]


def determinant(matrix):
    result = 0
    for permutation in permutations(range(len(matrix))):
        term = 1
        for i, j in enumerate(permutation):
            term = polynomial_product(term, matrix[i][j])
        result ^= term
    return result


def polynomial_rank(matrix):
    for size in range(4, 0, -1):
        for rows in combinations(range(4), size):
            for columns in combinations(range(4), size):
                if determinant([[matrix[i][j] for j in columns] for i in rows]):
                    return size
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(1729)
    matching_count = sum(1 for _ in matchings(tuple(range(8))))
    assert matching_count == 60 and len(local_terms()) == 23
    local_data = []
    for local, scalar in local_terms():
        row = [list(local[4*i:4*i+4]) for i in range(4)]
        column = [[local[(2*x+y)*4+2*u+v] for y, v in product(range(2), repeat=2)]
                  for x, u in product(range(2), repeat=2)]
        local_data.append({'low_scalar': scalar, 'matrix': row,
                           'generic_ranks': [polynomial_rank(row), polynomial_rank(column)]})
    rank_counts = Counter(tuple(item['generic_ranks']) for item in local_data)
    assert rank_counts == {(4, 1): 2, (1, 4): 2, (3, 3): 18, (4, 4): 1}
    checked = 0
    for width in range(1, 7):
        domain = labels(width)
        d = len(domain)
        indices = list(product(range(d), repeat=4)) if width <= 2 else [
            tuple(rng.randrange(d) for _ in range(4)) for _ in range(300)]
        for weights in [(1,)*(width+2), (1,)+tuple(1+rng.randrange(3) for _ in range(width+1))]:
            for ix in indices:
                x, y, u, v = (domain[a] for a in ix)
                assert tensor_profile(width, ix, weights) == profile(x, v, u, y, weights)
                checked += 1
    # Exchange the leading bit with a free bit: this exercises representative
    # complementation, which a permutation of just the free bits would miss.
    permutation_checks = 0
    for width in range(1, 5):
        d = 1 << width
        weights = (1,)+tuple(1+j % 3 for j in range(width+1))
        for permutation in permutations(range(width+1)):
            relabel = [permuted_label(width, a, permutation) for a in range(d)]
            assert sorted(relabel) == list(range(d))
            changed_weights = (1,)+tuple(weights[j+1] for j in permutation)
            indices = list(product(range(d), repeat=4)) if width <= 2 else [
                tuple(rng.randrange(d) for _ in range(4)) for _ in range(24)]
            for ix in indices:
                changed = tuple(relabel[a] for a in ix)
                assert tensor_profile(width, ix, weights) == tensor_profile(
                    width, changed, changed_weights)
                permutation_checks += 1
    paths = [Path(__file__), ROOT/'tools/hamming_tensor_power.py']
    result = {'status': 'exact local formula and permutation reduction; no improved exponent',
              'allowed_matchings': matching_count, 'local_terms': local_data,
              'profile_comparisons': checked, 'permutation_comparisons': permutation_checks,
              'sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('local_terms', 'sha256')}))


if __name__ == '__main__':
    main()
