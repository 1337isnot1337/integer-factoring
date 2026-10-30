#!/usr/bin/env python3
"""Exhaust the P10 fan's five lowest-y coefficient matrices by syndrome."""

import argparse
import hashlib
import json
from math import isqrt, log2
from pathlib import Path
import random
import time

from check_modular_crossing import rank_two
from check_ward_fan_leading_block import (
    FORMS, LENGTH, PAIRS, ROOT, RestrictedOracle, build_code, high_pivot_rank,
    independent_original_entry,
)
from explore_complement_edge_pencils import evaluate_local


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    code = build_code()
    groups = {tau: [(a, b) for a in range(len(code)) for b in range(len(code))
                    if code[a]['type_syndrome'] ^ code[b]['type_syndrome'] == tau]
              for tau in range(32)}
    groups = {tau: pairs for tau, pairs in groups.items() if pairs}
    assert len(groups) == 31
    assert sum(len(pairs) ** 2 for pairs in groups.values()) == 1917184
    lookups = {tau: {pair: i for i, pair in enumerate(pairs)} for tau, pairs in groups.items()}
    # Coefficient matrices remain split into blocks of size at most592.
    columns = {tau: [[0] * len(pairs) for _ in range(5)] for tau, pairs in groups.items()}
    oracle = RestrictedOracle(code)
    rng = random.Random(981274)
    original = evaluate_local((3, 1 ^ (1 << 16), 1, 1, 1, 1))
    truncation = sum(1 << (a + 16 * b) for a in range(8) for b in range(4))
    original = [value & truncation for value in original]
    reports, controls = [], []
    row_totals = [0] * 5
    nonzero_polynomials = 0
    for tau, pairs in groups.items():
        matrices = [[] for _ in range(5)]
        digest = hashlib.sha256()
        first_nonzero = None
        sample = {(rng.randrange(len(pairs)), rng.randrange(len(pairs))) for _ in range(2)}
        for i, (x, y) in enumerate(pairs):
            rows = [0] * 5
            for j, (v, u) in enumerate(pairs):
                value = oracle.polynomial(x, y, u, v)
                digest.update(bytes([value]))
                if value:
                    nonzero_polynomials += 1
                    if first_nonzero is None:
                        first_nonzero = [x, y, u, v, value]
                    column_tau = code[x]['type_syndrome'] ^ code[u]['type_syndrome']
                    a = lookups[column_tau][x, u]
                    b = lookups[column_tau][v, y]
                    for degree in range(3, 8):
                        if value >> degree & 1:
                            rows[degree - 3] |= 1 << j
                            columns[column_tau][degree - 3][a] |= 1 << b
                if (i, j) in sample:
                    actual = independent_original_entry(code, x, y, u, v, original)
                    assert actual == value << 48
                    controls.append({'type_pair_syndrome': tau, 'labels': [x, y, u, v],
                                     'univariate_quotient': value})
            oracle.logical_cache.clear()
            for matrix, row in zip(matrices, rows):
                matrix.append(row)
        ranks = [rank_two(matrix) for matrix in matrices]
        assert ranks == [high_pivot_rank(matrix) for matrix in matrices]
        row_totals = [a + b for a, b in zip(row_totals, ranks)]
        reports.append({'type_pair_syndrome': tau, 'dimension': len(pairs),
                        'row_coefficient_ranks_a3_through_a7': ranks,
                        'univariate_quotient_block_sha256': digest.hexdigest(),
                        'first_nonzero_entry': first_nonzero,
                        'row_coefficient_rows': [[hex(row) for row in matrix] for matrix in matrices]})
        print(json.dumps({'type_pair_syndrome': tau, 'dimension': len(pairs), 'ranks': ranks,
                          'physical_queries': oracle.counts['physical_query_requests'],
                          'distinct_physical_polynomials': len(oracle.physical_cache),
                          'elapsed_seconds': round(time.monotonic() - started, 3)}), flush=True)
    column_totals = [0] * 5
    for report in reports:
        matrices = columns[report['type_pair_syndrome']]
        ranks = [rank_two(matrix) for matrix in matrices]
        assert ranks == [high_pivot_rank(matrix) for matrix in matrices]
        report['column_coefficient_ranks_a3_through_a7'] = ranks
        report['column_coefficient_rows'] = [[hex(row) for row in matrix] for matrix in matrices]
        column_totals = [a + b for a, b in zip(column_totals, ranks)]
    old_path = ROOT / 'evidence/ward-fan-first-profile.json'
    old = next(record for record in json.loads(old_path.read_text()) if record['orientation_selected'])
    assert old['labels'] == len(code)
    updated = []
    for a, b, row, column in old['profile']:
        if b == 3 and 3 <= a <= 7:
            row, column = min(row, row_totals[a - 3]), min(column, column_totals[a - 3])
        updated.append([a, b, row, column])
    total = 0
    for _, _, row, column in updated:
        product_value = row * column
        root = isqrt(product_value)
        total += root + int(root * root < product_value)
    names = (Path(__file__).name, 'check_ward_fan_leading_block.py',
             'ward_fan_entry_queries.py', 'check_ward_fan_entry_queries.py',
             'check_ward_fan_family.py', 'check_ward_stable_operators.py',
             'complement_shell_tensor.py', 'explore_complement_edge_pencils.py',
             'pfaffian_particle_hole.py')
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names},
              'pairs': PAIRS, 'length': LENGTH, 'forms': FORMS, 'labels': len(code),
              'code_basis_sha256': hashlib.sha256(json.dumps(code, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
              'original_profile_sha256': hashlib.sha256(old_path.read_bytes()).hexdigest(),
              'blocks': reports, 'row_coefficient_ranks_a3_through_a7': row_totals,
              'column_coefficient_ranks_a3_through_a7': column_totals,
              'coefficient_entries_exhausted': 1917184,
              'nonzero_quotient_polynomials': nonzero_polynomials,
              'entire_y3_boundary_zero': nonzero_polynomials == 0,
              'independent_original_kernel_controls': controls,
              'query_counts': dict(oracle.counts),
              'complete_upper_profile': {'previous_integer_rank_sum_upper_bound': old['integer_rank_sum_upper_bound'],
                                         'integer_rank_sum_upper_bound': total,
                                         'rho_of_upper_bound': log2(total) / (2 * log2(len(code))),
                                         'profile': updated},
              'scope': 'Exact complete y^3 boundary for the specified finite fan code and output map; no improved factoring exponent asserted.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({'row_ranks': row_totals, 'column_ranks': column_totals,
                      'nonzero_polynomials': nonzero_polynomials, 'integer_rank_sum_upper_bound': total,
                      'elapsed_seconds': round(time.monotonic() - started, 3)}), flush=True)


if __name__ == '__main__':
    main()
