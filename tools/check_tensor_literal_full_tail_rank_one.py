#!/usr/bin/env python3
"""Nonhomogeneous full-tail cores versus the original raw Laurent tensor."""
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import check_tensor_hot_compounds_f16 as controls
import tensor_hot_compounds_f16 as field
from check_tensor_literal_hot_pfaffians_f16 import reference_pfaffian
from native_matrix_smith_f4 import MatrixOracle
from tensor_literal_full_tail_rank_one import core_profile, maps
from tensor_literal_rank_one_cell_search import all_cases


@lru_cache(None)
def raw_coefficients(modes):
    original_regular = reference_pfaffian(modes)

    @lru_cache(None)
    def coefficient(mask):
        value = original_regular(mask)
        if not value:
            return 0, 0, 0, 0
        row_height = sum(2*(i % modes)+1 for i in range(2*modes) if mask >> i & 1)
        degree = value.bit_length()-1
        reversed_value = int(bin(value)[2:][::-1], 2)
        return value, -row_height, reversed_value, row_height-degree
    return coefficient


def reference_entry(code, x, y, u, v, cut, clearing):
    p, left, modes, heights = (code[k] for k in ('embedding', 'retraction', 'modes', 'heights'))
    raw = raw_coefficients(modes)
    planes = [0]*4
    for tx, ty, tu, tv in product(left[x].items(), left[y].items(), p[u].items(), p[v].items()):
        terms = (tx, ty, tu, tv)
        mask = sum(w << (leg*modes) for leg, ((w, _), _) in enumerate(terms))
        scalar, exponent = 1, 0
        for (_, e), c in terms:
            scalar = field.mul(scalar, c)
            exponent += e
        value, offset, reversed_value, reversed_offset = raw(mask)
        if not value:
            continue
        if cut == 'zero':
            power = offset+exponent+heights[x]+heights[y]+clearing
        else:
            value = reversed_value
            power = reversed_offset-exponent+heights[u]+heights[v]+clearing
        assert power >= 0
        for bit in range(4):
            if scalar >> bit & 1:
                planes[bit] ^= value << power
    return tuple(planes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-full-tail-rank-one-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    records, entries, order_checks, toeplitz_checks, endpoint_checks = [], 0, 0, 0, 0
    for case in all_cases(4):
        p, left = maps(case)
        words = [tail << 1 for tail in range(8)]
        code = {'modes': 4, 'q': 8, 'heights': [field.height(w) for w in words],
                'embedding': [{(w | bit, 0): c for bit, c in enumerate(p) if c} for w in words],
                'retraction': [{(w | bit, 0): c for bit, c in enumerate(left) if c} for w in words]}
        cut_records = []
        for cut in ('zero', 'infinity'):
            clearing = 2*int((left if cut == 'zero' else p)[1] != 0)
            matrix = [[(0, 0, 0, 0) for _ in range(64)] for _ in range(64)]
            for x, y, u, v in product(range(8), repeat=4):
                value = reference_entry(code, x, y, u, v, cut, clearing)
                assert value[2:] == (0, 0)
                assert sum((plane.bit_count() % 2) << bit for bit, plane in enumerate(value)) == int(x == v and y == u)
                if cut == 'zero':
                    matrix[8*x+y][8*v+u] = value
                else:
                    matrix[8*x+u][8*v+y] = value
                entries += 1
                endpoint_checks += 1
            converted = [[value[:2] for value in row] for row in matrix]
            precision = 128
            while True:
                try:
                    expected = MatrixOracle(converted, 64, precision=precision).profile()['orders']
                    break
                except ArithmeticError:
                    precision *= 2
                    assert precision <= 8192
            assert all(value == round(value) for value in expected)
            expected = sorted(map(round, expected))
            rank_checks = controls.toeplitz_controls(matrix, expected)
            actual = core_profile(4, case, cut)
            assert actual['clearing_power'] == clearing
            assert actual['cleared_orders'] == expected
            assert actual['orders'] == [value-clearing for value in expected]
            order_checks += 64
            toeplitz_checks += len(rank_checks)
            cut_records.append({'active_profile': actual, 'literal_cleared_orders': expected,
                                'literal_smith_precision': precision, 'toeplitz_checks': rank_checks})
        records.append({'case': case, 'embedding': p, 'retraction': left,
                        'assigned_low_virtual_height': 0, 'cuts': cut_records})
    names = [Path(__file__).name, 'tensor_literal_full_tail_rank_one.py',
             'tensor_literal_full_tail_core.py', 'tensor_literal_rank_one_cell_search.py',
             'check_tensor_literal_hot_pfaffians_f16.py', 'check_tensor_hot_compounds_f16.py',
             'reflection_dual_rail_gaussian.py', 'tensor_hot_compounds_f16.py',
             'tensor_orthogonal_scaling.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_weighted_smith_f4.py',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'cases': records, 'frobenius_representatives': len(records),
              'normalized_f4_pairs_covered': 20, 'raw_laurent_literal_entries': entries,
              'endpoint_checks': endpoint_checks, 'complete_smith_order_comparisons': order_checks,
              'independent_toeplitz_rank_checks': toeplitz_checks,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'seconds': time.monotonic()-started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('cases', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
