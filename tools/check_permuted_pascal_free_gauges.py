#!/usr/bin/env python3
"""Independent complete Fock controls outside the explicit Pascal gauge cone.

The pure Python polynomial elimination is independent of the native weighted
pivot engine used by the search. All heights and resulting orders here are
integers. These controls validate the evaluator, not a universal pressure floor.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_weighted_pascal_flags import (
    exact_weighted_profile, original_fock_profile, subset_sums,
)
from permuted_pascal_kernels import permuted_pascal_kernel
from tensor_pascal_box_crossings import exterior_tensor, matrix_product
from weighted_pascal_flags import (
    change_flag, clifford_data, inverse_kernel, natural_heights,
)


def check_case(images, perturbation, shear, native=False):
    seed = permuted_pascal_kernel(images)
    n = seed['n']
    if shear:
        flag = [[int(i == j) for j in range(n)] for i in range(n)]
        flag[0][-1] ^= 1
        seed = change_flag(seed, flag)
    assert matrix_product(seed['matrix'], inverse_kernel(seed)) == [
        [int(i == j) for j in range(2 * n)] for i in range(2 * n)]
    tensor = exterior_tensor(seed)
    cuts = []
    for place, inherited in enumerate(natural_heights(seed)):
        heights = [value + perturbation * (((5 * i + 3 * place) % 7) - 3)
                   for i, value in enumerate(inherited)]
        values, shifts, rows, columns = clifford_data(seed, place)
        row = [heights[i] for i in rows]
        column = [heights[i] for i in columns]
        row += [-value for value in row]
        column += [-value for value in column]
        orders = exact_weighted_profile(values, shifts, row, column, bool(place))
        assert orders == [-value for value in orders[::-1]]
        positive = orders[2 * n:]
        twice_mean = (2 * seed['outer'] * n if place == 0 else 0) + sum(heights)
        assert (twice_mean - sum(positive)) % 2 == 0
        minimum = (twice_mean - sum(positive)) // 2
        reconstructed = [minimum + value for value in subset_sums(positive)]
        original = original_fock_profile(seed, tensor, heights, place)
        assert reconstructed == original
        record = {'place': place, 'heights': heights, 'clifford_orders': orders,
                  'complete_fock_orders_compared': len(original),
                  'minimum_order': original[0], 'maximum_order': original[-1]}
        if native:
            from native_matrix_smith_f4 import MatrixOracle
            oracle = MatrixOracle([[(entry, 0) for entry in line] for line in values],
                                  4 * n, shifts, 1024)
            result = oracle.profile(row, column, 1., infinity=bool(place))
            assert result['orders'] == orders
            record['native_orders_compared'] = len(orders)
        cuts.append(record)
    return {'images': images, 'zero_energies': seed['zero_energies'],
            'perturbation': perturbation, 'constant_flag_shear': shear, 'cuts': cuts}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', action='store_true',
                        help='Also compare the exact integer orders with the search engine.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    samples = [([0, 1], 1, False), ([1, 0], 2, False),
               ([1, 0], 1, True), ([0, 1, 2, 3], 1, False),
               ([0, 1, 3, 2], 1, False), ([0, 3, 1, 2], 2, False),
               ([2, 0, 3, 1], 1, True)]
    cases = [check_case(images, perturbation, shear, args.native)
             for images, perturbation, shear in samples]
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/check_permuted_pascal_free_gauges.py',
               'tools/permuted_pascal_kernels.py', 'tools/weighted_pascal_flags.py',
               'tools/check_weighted_pascal_flags.py', 'tools/tensor_pascal_box_crossings.py',
               'tools/cube_symmetric_exterior_blocks.py']
    if args.native:
        sources += ['tools/search_weighted_pascal_flags.py', 'tools/native_matrix_smith_f4.py',
                    'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp',
                    'tools/native_weighted_smith_f4.py']
    report = {'status': 'pass', 'scope': 'Complete finite original-polynomial controls for '
              'arbitrary exponent permutations, independent physical gauges and one flag shear. '
              'No pressure improvement or universal exclusion is claimed.',
              'native_comparison_requested': args.native, 'cases': cases,
              'complete_fock_orders_compared': sum(cut['complete_fock_orders_compared']
                                                  for case in cases for cut in case['cuts']),
              'source_sha256': {path: sha256((root / path).read_bytes()).hexdigest()
                                for path in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'native_comparison_requested', 'complete_fock_orders_compared')}))


if __name__ == '__main__':
    main()
