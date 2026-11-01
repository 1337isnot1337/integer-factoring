#!/usr/bin/env python3
"""Search triple-single output corrections outside the cubic syndrome theorem."""

import argparse
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from math import isqrt, sqrt
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local, minimum_decomposition
from check_double_endpoint_retracts import independent_coefficients


def run(pairs, orientations=8):
    assert pairs >= 3 and pairs % 2
    local = evaluate_local((2, 16, 1, 1, 1, 1))
    width = pairs * pairs
    degrees = [p + 4 * q for q in range(4) for p in range(4 - q)]

    def monomial(i, j):
        return tuple(1 if h == i else 3 if h == j else 0 for h in range(pairs))

    inputs = [[monomial(i, j) for j in range(pairs) if i != j] for i in range(pairs)]
    outputs = [monomial(i, (i + 1) % pairs) for i in range(pairs)]

    @lru_cache(None)
    def pattern_value(pattern):
        value = 1
        for index in pattern:
            value = multiply(value, local[index])
            if not value:
                break
        return value

    @lru_cache(None)
    def kernel(first, second, y, v):
        syndrome = sum(((a.bit_count() + b.bit_count()) % 2) << i
                       for i, (a, b) in enumerate(zip(first, second)))
        if syndrome ^ (1 << y) ^ (1 << v):
            return 0
        value = 0
        for yy, vv in product(inputs[y], inputs[v]):
            pattern = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + d
                                   for a, b, c, d in zip(first, yy, second, vv)))
            value ^= pattern_value(pattern)
        return value

    baseline = [[kernel(outputs[x], outputs[u], y, v)
                 for v in range(pairs) for u in range(pairs)]
                for x in range(pairs) for y in range(pairs)]

    def profile(tensor):
        rows = [[0] * width for _ in degrees]
        columns = [[0] * width for _ in degrees]
        for x, y, v, u in product(range(pairs), repeat=4):
            value = tensor[x * pairs + y][v * pairs + u]
            assert (value & 1) == int(x == v and y == u)
            assert not (value >> 3) & 1 and not (value >> 12) & 1
            for i, degree in enumerate(degrees):
                if value >> degree & 1:
                    rows[i][x * pairs + y] ^= 1 << (v * pairs + u)
                    columns[i][x * pairs + u] ^= 1 << (v * pairs + y)
        chosen = independent_coefficients([[matrix] for matrix in rows])
        assert chosen[0] == 0
        row, column = [rows[i] for i in chosen], [columns[i] for i in chosen]
        ranks = [[0, 0]]
        for mask in range(1, 1 << len(chosen)):
            matrices = []
            for basis in (row, column):
                matrix = [0] * width
                for i, source in enumerate(basis):
                    if mask >> i & 1:
                        matrix = [a ^ b for a, b in zip(matrix, source)]
                matrices.append(rank_two(matrix))
            ranks.append(matrices)
        cost, witness = minimum_decomposition([sqrt(r * c) for r, c in ranks], target=1)
        lower, _ = minimum_decomposition([isqrt(r * c) for r, c in ranks], target=1)
        return {'basis_degrees': [degrees[i] for i in chosen], 'rank_table': ranks,
                'cost': cost, 'integer_lower': lower, 'witness': witness,
                'exclusion': lower ** 70 > pairs ** 99}

    base_profile = profile(baseline)
    best = dict(base_profile, corrected_label=None, triple=None)
    records = []
    for i in range(pairs):
        for triple, orientation in product(combinations(range(pairs), 3), range(orientations)):
            extra = tuple(1 + ((orientation >> triple.index(j)) & 1) if j in triple else 0
                          for j in range(pairs))
            tensor = [row.copy() for row in baseline]
            for y, v, u in product(range(pairs), repeat=3):
                tensor[i * pairs + y][v * pairs + u] ^= kernel(extra, outputs[u], y, v)
            for x, y, v in product(range(pairs), repeat=3):
                tensor[x * pairs + y][v * pairs + i] ^= kernel(outputs[x], extra, y, v)
            for y, v in product(range(pairs), repeat=2):
                tensor[i * pairs + y][v * pairs + i] ^= kernel(extra, extra, y, v)
            changed = sum(a != b for row, old in zip(tensor, baseline) for a, b in zip(row, old))
            result = dict(profile(tensor) if changed else base_profile,
                          corrected_label=i, triple=list(triple), orientation=orientation,
                          changed_polynomial_entries=changed)
            result['rank_table_sha256'] = hashlib.sha256(
                json.dumps(result['rank_table'], separators=(',', ':')).encode()).hexdigest()
            if result['cost'] < best['cost']:
                best = result
            records.append({key: value for key, value in result.items() if key != 'rank_table'})
    untouched = pairs - 1
    lower_square_numerator = untouched ** 2 * (untouched - 1) ** 2
    return {'pairs': pairs, 'labels': pairs, 'length': 3,
            'untouched_output_labels': untouched,
            'all_field_rank_cost_lower_square_numerator': lower_square_numerator,
            'all_field_rank_cost_lower_square_denominator': 2,
            'all_single_output_changes_excluded_by_integer_gate':
                lower_square_numerator ** 35 > 2 ** 35 * pairs ** 99,
            'complete_single_corrections': len(records), 'orientations_per_triple': orientations,
            'all_recorded_spaces_excluded': all(row['exclusion'] for row in records),
            'baseline': base_profile, 'best': best, 'records': records,
            'distinct_local_patterns': pattern_value.cache_info().currsize,
            'contracted_physical_entries': kernel.cache_info().currsize}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', type=int, default=7)
    parser.add_argument('--orientations', type=int, default=8, choices=(1, 8))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.pairs, args.orientations)
    result['scope'] = ('Every single triple-mode correction in the recorded orientations to one cyclic coordinate output; '
                       'not all output maps or combinations of corrections.')
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('records', 'baseline', 'best')}))
    print(json.dumps({key: value for key, value in result['best'].items() if key != 'rank_table'}))


if __name__ == '__main__':
    main()
