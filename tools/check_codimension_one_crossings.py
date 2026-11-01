#!/usr/bin/env python3
"""Exhaust every constant rank-seven retraction of three binary sites."""

import argparse
from collections import Counter, deque
from hashlib import sha256
import json
from math import comb, log
from pathlib import Path

from check_block_smith_crossings import reverse_polynomial, toeplitz_rank
from linear_physical_restriction import (
    linear_restriction_matrices, linear_restriction_profiles,
)


def permute_vector(value, permutation):
    return sum(((value >> index) & 1) << permutation[index] for index in range(8))


def symmetry_generators():
    """Local GL2 swaps/transvections and adjacent site permutations."""
    result = []
    for site in range(3):
        stride = 1 << site
        lower = sum(1 << index for index in range(8) if not index & stride)
        upper = 255 ^ lower
        swap = tuple(((value & lower) << stride) | ((value & upper) >> stride)
                     for value in range(256))
        forward = tuple(value ^ ((value & upper) >> stride) for value in range(256))
        dual = tuple(value ^ ((value & lower) << stride) for value in range(256))
        result.append((swap, swap))
        result.append((forward, dual))
    for first, second in ((0, 1), (1, 2)):
        permutation = []
        for index in range(8):
            if ((index >> first) ^ (index >> second)) & 1:
                index ^= (1 << first) | (1 << second)
            permutation.append(index)
        table = tuple(permute_vector(value, permutation) for value in range(256))
        result.append((table, table))
    return result


def pair_retraction(vector, functional):
    """Image ker(f), kernel span(v), with f(v)=1."""
    if not 1 <= vector < 256 or not 1 <= functional < 256:
        raise ValueError('The vector and functional must be nonzero eight-bit masks.')
    if (vector & functional).bit_count() % 2 != 1:
        raise ValueError('The functional must take value one on the kernel vector.')
    pivot = (functional & -functional).bit_length()-1
    free = [index for index in range(8) if index != pivot]
    embedding = [(1 << index) | (((functional >> index) & 1) << pivot) for index in free]
    retraction = [(1 << index) ^ (((vector >> index) & 1)*functional) for index in free]
    for index in range(8):
        image = 0
        for column, covector in zip(embedding, retraction):
            if covector >> index & 1:
                image ^= column
        assert image == (1 << index) ^ (((functional >> index) & 1)*vector)
    return embedding, retraction


def quadratic_invariant(value):
    return sum(((value >> index) & 1)*((value >> (7-index)) & 1)
               for index in range(4)) % 2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/codimension-one-crossings.json')
    args = parser.parse_args()
    generators = symmetry_generators()
    identity = [1 << index for index in range(8)]
    ambient = linear_restriction_matrices(3, identity, identity)
    for forward, dual in generators:
        changed = linear_restriction_matrices(
            3, [forward[1 << index] for index in range(8)],
            [dual[1 << index] for index in range(8)],
        )
        assert changed['row'] == ambient['row'] and changed['column'] == ambient['column']
    row_profiles = {
        0: (1, 6, 12, 13, 11, 6),
        1: (1, 6, 15, 7, 14, 6),
    }
    column_profile = (1, 6, 15, 18, 9)
    visited, records, examples = set(), [], {}
    type_pairs, type_orbits = Counter(), Counter()
    action_checks = 0
    for vector in range(1, 256):
        for functional in range(1, 256):
            if (vector & functional).bit_count() % 2 != 1:
                continue
            start = 256*vector+functional
            if start in visited:
                continue
            visited.add(start)
            queue = deque([start])
            orbit_size = 0
            kind = quadratic_invariant(vector)*quadratic_invariant(functional)
            while queue:
                point = queue.popleft()
                current_vector, current_functional = divmod(point, 256)
                orbit_size += 1
                assert quadratic_invariant(current_vector)*quadratic_invariant(current_functional) == kind
                for forward, dual in generators:
                    next_vector, next_functional = forward[current_vector], dual[current_functional]
                    assert (next_vector & next_functional).bit_count() % 2 == 1
                    neighbor = 256*next_vector+next_functional
                    action_checks += 1
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            embedding, retraction = pair_retraction(vector, functional)
            profiles = linear_restriction_profiles(3, embedding, retraction)
            expected_row = tuple(order for order, count in enumerate(row_profiles[kind]) for _ in range(count))
            expected_column = tuple(order for order, count in enumerate(column_profile) for _ in range(count))
            assert profiles['row_infinity_orders'] == expected_row
            assert profiles['column_zero_orders'] == expected_column
            examples.setdefault(kind, (embedding, retraction, profiles))
            type_pairs[kind] += orbit_size
            type_orbits[kind] += 1
            records.append({'representative': [vector, functional], 'orbit_size': orbit_size,
                            'quadratic_invariant_product': kind,
                            'row_infinity_counts': row_profiles[kind],
                            'column_zero_counts': column_profile,
                            'row_precision': profiles['row_precision'],
                            'column_precision': profiles['column_precision']})
    assert len(visited) == 32640 and len(records) == 60
    assert dict(type_pairs) == {0: 25920, 1: 6720}
    assert dict(type_orbits) == {0: 51, 1: 9}
    toeplitz_controls = 0
    for embedding, retraction, profiles in examples.values():
        matrices = linear_restriction_matrices(3, embedding, retraction)
        reversed_row = [[reverse_polynomial(value, 3) for value in line] for line in matrices['row']]
        for matrix, orders in ((reversed_row, profiles['row_infinity_orders']),
                               (matrices['column'], profiles['column_zero_orders'])):
            for precision in range(1, max(orders)+2):
                assert toeplitz_rank(matrix, precision) == sum(max(0, precision-order) for order in orders)
                toeplitz_controls += 1
    # P_B-P_A has positive real pressure, although the coefficient field is F2.
    assert [right-left for left, right in zip(row_profiles[0], row_profiles[1])] == [0, 0, 3, -6, 3, 0]
    polynomial = [0]*10
    for i, left in enumerate(row_profiles[0]):
        for j, right in enumerate(column_profile):
            polynomial[i+j] += left*right
    polynomial[3] -= 686
    assert polynomial == [1, 12, 63, -493, 386, 537, 543, 405, 207, 54]
    degree = 48
    bernstein = [sum(coefficient*comb(degree-j, k-j)
                     for j, coefficient in enumerate(polynomial) if j <= k)
                 for k in range(degree+1)]
    assert min(bernstein) > 0
    root = Path(__file__).resolve().parents[1]
    report = {
        'status': 'pass',
        'scope': 'All constant rank-seven retractions of the three-site binary elementary tensor, without additional gauges.',
        'projectors_checked': len(visited), 'symmetry_orbits': len(records),
        'local_symmetry_generators_checked_on_full_tensor': len(generators),
        'projector_symmetry_actions_checked': action_checks,
        'independent_toeplitz_rank_controls': toeplitz_controls,
        'profile_types': 2, 'type_projector_counts': dict(type_pairs),
        'type_orbit_counts': dict(type_orbits),
        'uniform_pressure_lower_bound': 686,
        'rho_lower_bound_exact': 'log(686)/(4 log(7))',
        'rho_lower_bound_display_only': log(686)/(4*log(7)),
        'bernstein_polynomial': polynomial,
        'bernstein_degree': degree, 'scaled_bernstein_coefficients': bernstein,
        'floating_arithmetic_used_for_certification': False,
        'factoring_exponent_improved': False,
        'records': records,
        'engine_sha256': sha256((root/'tools/linear_physical_restriction.py').read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root/output
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in (
        'status', 'projectors_checked', 'symmetry_orbits', 'profile_types',
        'uniform_pressure_lower_bound', 'rho_lower_bound_display_only',
    )}))


if __name__ == '__main__':
    main()
