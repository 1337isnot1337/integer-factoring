#!/usr/bin/env python3
"""Exact dressed input repairs and common-flag two-sided crossing controls."""
import argparse
import json
from pathlib import Path
import random

from check_rectangular_f4_repair import build, matrix_multiply, times_binary, binary_power
from check_rectangular_refined_repairs import height, single_shear, postcorrect_seed, precorrect_seed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = build((0, 1), 2, 1, 2)
    counts = {'strict_height': 0, 'equal_height': 0}
    cases = []
    for row in range(12):
        for column in range(1, 12):
            if row == column or height(row) > height(column):
                continue
            counts['strict_height' if height(row) < height(column) else 'equal_height'] += 1
            for coefficient in (1, 2, 3):
                matrix, wx, wy = single_shear(row, column, coefficient)
                seed = precorrect_seed(base, matrix, wx, wy)
                cases.append({'row': row, 'column': column, 'coefficient': coefficient,
                              'denominator_degree': seed['denominator_degree'],
                              'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                              'two_copy_type_checks': seed['two_copy_type_checks'],
                              'coefficient_sha256': seed['coefficient_sha256']})
    source = build((0, 1, 2), 4, 1, 2, check_types=False)
    background = []
    for row in range(1, 12):
        if height(row) != height(0):
            continue
        for coefficient in (1, 2, 3):
            matrix, wx, wy = single_shear(row, 0, coefficient)
            seed = precorrect_seed(base, matrix, wx, wy, keep_background=False)
            assert not seed['input_background_fixed']
            background.append({'row': row, 'column': 0, 'coefficient': coefficient,
                               'denominator_degree': seed['denominator_degree'],
                               'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                               'two_copy_type_checks': seed['two_copy_type_checks'],
                               'coefficient_sha256': seed['coefficient_sha256']})
    rng = random.Random(970212)
    joint = []
    for case in range(3):
        wx, wy = [rng.randrange(-50, 51) for _ in range(4)], [rng.randrange(-50, 51) for _ in range(3)]
        keys = [(height(i), wx[i % 4]+wy[i // 4]) for i in range(12)]
        V = [[int(i == j) if i == j else rng.randrange(4) if i and keys[i] < keys[j] else 0
              for j in range(12)] for i in range(12)]
        W = [[int(i == j) if i == j else rng.randrange(4) if j and keys[i] < keys[j] else 0
              for j in range(12)] for i in range(12)]
        seed = precorrect_seed(postcorrect_seed(source, V, wx, wy, check_types=False), W, wx, wy)
        # The all-low core block, including every bulk assignment, is fixed.
        A, B, bulk = source['q_x'], source['q_u'], source['bulk_size']
        checks = 0
        for x in range(bulk):
            for y in range(bulk):
                for u in range(bulk):
                    for v in range(bulk):
                        index = x+A*(y+B*(u+B*v))
                        assert times_binary(seed['values'][index], binary_power(3, source['denominator_degree'])) == \
                               times_binary(source['values'][index], binary_power(3, seed['denominator_degree']))
                        checks += 1
        joint.append({'postcorrection': V, 'input_correction': W, 'secondary_x': wx, 'secondary_y': wy,
                      'denominator_degree': seed['denominator_degree'],
                      'ramified_endpoint_checks': seed['ramified_endpoint_checks'],
                      'two_copy_type_checks': seed['two_copy_type_checks'],
                      'fixed_background_entries_checked': checks, 'coefficient_sha256': seed['coefficient_sha256']})
    V, vx, vy = single_shear(3, 8)
    W, wx, wy = single_shear(8, 3)
    try:
        precorrect_seed(postcorrect_seed(base, V, vx, vy, check_types=False), W, wx, wy)
    except ValueError:
        rejected = True
    else:
        raise AssertionError('Incompatible secondary flags were accepted.')
    invalid_endpoint = matrix_multiply(matrix_multiply(V, base['ramified_endpoint_matrix']), W)
    bad_diagonal = [i for i in range(12) if invalid_endpoint[i][i] != 1]
    assert bad_diagonal == [3]
    result = {'status': 'exact input and common-flag two-sided polynomial repairs; no pressure claim',
              'single_input_directions_with_column_zero_fixed': counts,
              'single_input_cases': cases, 'joint_flag_cases': joint,
              'background_changing_input_cases': background,
              'incompatible_flag_control': {'rejected': rejected, 'output_direction': [3, 8],
                                            'input_direction': [8, 3], 'failed_endpoint_diagonal': bad_diagonal}}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
