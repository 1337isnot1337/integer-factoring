#!/usr/bin/env python3
"""Exact simplex cancellations and full coherent matching-profile controls."""

from itertools import product
from collections import Counter
import json

from matching_coherent_projections import canonical_matching_split, canonical_projection
from perfect_matching_crossing import dense_seed, binary_rank, toeplitz_profile
from check_perfect_matching_crossings import literal_cut, reverse


def literal_kronecker_contraction(source, left, right):
    """Contract each coefficient through literal row/column Kronecker maps."""
    q, n = len(left), len(left[0])
    row_maps = [sum((a[i] & b[j]) << (n * i + j)
                    for i in range(n) for j in range(n))
                for a in left for b in left]
    col_maps = [sum((a[i] & b[j]) << (n * i + j)
                    for i in range(n) for j in range(n))
                for a in right for b in right]
    result = [0] * q ** 4
    for power in range(4):
        source_rows = [0] * n ** 2
        for index, (value, high) in enumerate(source['values']):
            assert high == 0
            if not (value >> power & 1):
                continue
            x, y, u, v = (index // n ** k % n for k in range(4))
            source_rows[n * x + y] ^= 1 << (n * u + v)
        for r, row in enumerate(row_maps):
            combined = 0
            while row:
                bit = row & -row
                row ^= bit
                combined ^= source_rows[bit.bit_length() - 1]
            x, y = divmod(r, q)
            for c, column in enumerate(col_maps):
                u, v = divmod(c, q)
                bit = (combined & column).bit_count() & 1
                result[x + q * y + q ** 2 * u + q ** 3 * v] ^= bit << power
    return result


def verify():
    source = dense_seed(3)
    left, right = canonical_matching_split(True)
    patterns = [sum(right[j][i] << j for j in range(4)) for i in range(15)]
    assert sorted(patterns) == list(range(1, 16))
    triple_checks = 0
    for a, b, c in product(range(5), repeat=3):
        total = sum(right[a][i] * right[b][i] * right[c][i] for i in range(15)) & 1
        assert total == int(a == b == c == 4)
        triple_checks += 1
    results, coefficient_entries, profile_ranks = [], 0, 0
    for include_trivial in (False, True):
        seed = canonical_projection(include_trivial)
        q = seed['dimensions'][0]
        direct = literal_kronecker_contraction(source, seed['left'], seed['right'])
        assert direct == [value[0] for value in seed['values']]
        coefficient_entries += 4 * q ** 4
        degree = seed['degree']
        assert degree == (3 if include_trivial else 1)
        assert all(not (value >> 2 & 1) for value in direct)
        if include_trivial:
            assert [i for i, value in enumerate(direct) if value >> 3 & 1] == [q ** 4 - 1]
        ranks = []
        profiles = []
        matrices = [literal_cut(seed, side) for side in (False, True)]
        for power in range(degree + 1):
            ranks.append([power] + [binary_rank([
                sum(((value >> power) & 1) << j for j, value in enumerate(row))
                for row in matrix]) for matrix in matrices])
        assert ranks == ([[0, 8, 8], [1, 20, 20], [2, 0, 0], [3, 1, 1]]
                         if include_trivial else [[0, 8, 8], [1, 12, 12]])
        for matrix in matrices:
            zero, transcript = toeplitz_profile(matrix, 6)
            profile_ranks += len(transcript) - 1
            reciprocal = [[reverse(value, degree) for value in row] for row in matrix]
            infinity_cleared, transcript = toeplitz_profile(reciprocal, 6)
            profile_ranks += len(transcript) - 1
            infinity = {order - degree: count for order, count in infinity_cleared.items()}
            assert zero == ({0: 8, 1: 16, 3: 1} if include_trivial else {0: 8, 1: 8})
            assert infinity == ({-3: 1, -1: 20, 1: 4} if include_trivial else {-1: 12, 1: 4})
            profiles.append({'zero': zero, 'infinity': infinity})
        results.append({'alphabet': q, 'degree': degree, 'coefficient_rank_pairs': ranks,
                        'complete_profiles': profiles,
                        'left': seed['left'], 'right': seed['right']})
    return {'status': 'exact coherent matching splits, simplex cancellation, and profile controls passed',
            'simplex_column_patterns': patterns, 'triple_product_controls': triple_checks,
            'independent_coefficient_contractions': coefficient_entries,
            'independent_block_toeplitz_ranks': profile_ranks,
            'cases': results, 'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
