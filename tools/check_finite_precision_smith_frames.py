#!/usr/bin/env python3
"""Independent exact controls for finite regular-frame representatives."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_block_smith_crossings import toeplitz_rank
from check_full_alphabet_two_cut_lattices import independent_profiles
from finite_precision_smith_frames import regular_smith_frame, build_lattice
from full_alphabet_two_cut_lattices import CompressedFullAlphabetFiber
from search_full_alphabet_two_cut_lattices import root_for_profiles
from search_full_alphabet_reference_walk import paired_energies
from weyl_lattice_filter_synthesis import matrix_multiply, regular_smith_basis


def constant_rank(matrix):
    pivots = {}
    for line in matrix:
        row = sum((value & 1) << j for j, value in enumerate(line))
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--search-records', type=Path, nargs='*', default=[])
    args = parser.parse_args()
    rng = random.Random(311849)
    cases = literal_products = toeplitz_checks = exact_comparisons = 0
    for size in range(2, 8):
        for trial in range(6):
            orders = sorted(rng.randrange(10) for _ in range(size))
            left = [[int(i == j) for j in range(size)] for i in range(size)]
            right = [[int(i == j) for j in range(size)] for i in range(size)]
            for i in range(size):
                for j in range(i):
                    left[i][j] = rng.randrange(1 << 6)
                    right[j][i] = rng.randrange(1 << 6)
            diagonal = [[(1 << orders[i]) if i == j else 0 for j in range(size)]
                        for i in range(size)]
            source = matrix_multiply(matrix_multiply(left, diagonal), right)
            rng.shuffle(source)
            for allowance in (1, 4, 9):
                precision = max(orders)+allowance
                actual, a, b, reduced = regular_smith_frame(source, precision)
                assert list(actual) == orders
                assert constant_rank(a) == constant_rank(b) == size
                product = matrix_multiply(matrix_multiply(a, source), b)
                mask = (1 << precision)-1
                assert [[value & mask for value in row] for row in product] == reduced
                literal_products += 1
                for length in sorted({1, orders[size//2]+1, max(orders)+1}):
                    assert toeplitz_rank(product, length) == sum(max(0, length-v) for v in orders)
                    toeplitz_checks += 1
                cases += 1
            if size <= 3:
                assert list(regular_smith_basis(source)[0]) == orders
                exact_comparisons += 1
    rejected = 0
    for source, precision in [([[1, 0], [0, 1 << 5]], 5), ([[1, 1], [1, 1]], 8)]:
        try:
            regular_smith_frame(source, precision)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('An unresolved final channel was accepted.')
    root = Path(__file__).resolve().parents[1]
    reference_rows = []
    for index in (0, 1):
        path = root/'evidence'/f'affine-degenerate-witness-{index}-integer-chamber.json'
        reference = json.loads(path.read_text())
        values = tuple((value[0], 0) for value in reference['values'])
        lattice = build_lattice(values, reference['q'], reference['degree'],
                                reference['heights_infinity_then_zero'], raise_cap=12)
        expected = tuple(map(tuple, reference['profiles_zero_then_infinity'][::-1]))
        assert lattice.profiles(values) == expected
        assert lattice.profiles(reference['values']) == expected
        independent, counts = independent_profiles(lattice, values)
        assert independent == expected
        toeplitz_checks += counts
        floors = tuple(max(0, value-4) for value in lattice.orders)
        fiber = CompressedFullAlphabetFiber(lattice, floors)
        for targets in (lattice.orders, floors):
            solution = fiber.solve(targets)
            assert solution is not None
            sample = fiber.sample(solution, rng)
            lattice.check_witness(sample, targets)
        cost = root_for_profiles(expected, *paired_energies(reference['heights_infinity_then_zero'], reference['q']))
        assert abs(cost['rho']-reference['gibbs_root_search']['rho']) < 1e-9
        reference_rows.append({'reference': path.name, 'complete_profiles': expected,
                               'binary_and_F4_profiles_equal': True, 'replayed_rho': cost['rho'],
                               'coefficient_witnesses': 2, 'fiber_variables': fiber.variables})
    searches = []
    for path in args.search_records:
        record = json.loads(path.read_text())
        assert record['completed_requested_trials']
        for name, digest in record['source_sha256'].items():
            assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest
        reference = record['reference']
        lattice = build_lattice(record['binary_reference_values'], record['physical_alphabet'],
                                record['degree'], record['heights_infinity_then_zero'],
                                record['arguments']['upper_window'])
        for witness in record['retained_witnesses']:
            actual = lattice.check_witness(witness['values'], witness['targets'])
            assert actual == tuple(map(tuple, witness['profiles']))
            cost = root_for_profiles(actual, *record['physical_paired_energies'])
            assert abs(cost['rho']-witness['rho']) < 1e-9
        searches.append({'record': path.name, 'completed_trials': record['completed_trials'],
                         'distinct_profile_classes': record['distinct_profile_classes'],
                         'retained_literal_replays': len(record['retained_witnesses'])})
    result = {'status': 'Finite-frame precision, actual unit representatives, and literal coefficient witnesses passed.',
              'known_complete_profile_cases': cases, 'full_untruncated_product_checks': literal_products,
              'independent_toeplitz_ranks': toeplitz_checks,
              'untruncated_fraction_free_comparisons': exact_comparisons,
              'rejected_insufficient_precision_or_rank': rejected,
              'rounded_gauge_references': reference_rows, 'search_replays': searches,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'finite_precision_smith_frames.py',
                                             'full_alphabet_two_cut_lattices.py')}}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
