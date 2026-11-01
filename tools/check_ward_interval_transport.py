#!/usr/bin/env python3
"""Active-support transport and optimal fixed-order binary form bases."""

import argparse
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import random

from check_modular_crossing import rank_two
from check_ward_form_operator_transport import insertion_transitions
from check_ward_stable_operators import bits, closed_operator_interface


def minimal_span_basis(rows, pairs):
    rows = list(rows)
    forms = len(rows)
    change = [1 << i for i in range(forms)]
    pivot = 0
    for column in range(pairs):
        found = next((i for i in range(pivot, forms) if rows[i] >> column & 1), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        change[pivot], change[found] = change[found], change[pivot]
        for i in range(pivot + 1, forms):
            if rows[i] >> column & 1:
                rows[i] ^= rows[pivot]
                change[i] ^= change[pivot]
        pivot += 1
        if pivot == forms:
            break
    assert pivot == forms
    reductions = 0
    while True:
        endings = {}
        for i, value in enumerate(rows):
            last = value.bit_length() - 1
            if last in endings:
                earlier = endings[last]
                first = (rows[earlier] & -rows[earlier]).bit_length() - 1
                rows[earlier] ^= value
                change[earlier] ^= change[i]
                assert (rows[earlier] & -rows[earlier]).bit_length() - 1 == first
                assert rows[earlier].bit_length() - 1 < last
                reductions += 1
                break
            endings[last] = i
        else:
            break
    assert rank_two(change) == forms
    return rows, change, reductions


def interval_profile(rows, pairs):
    endpoints = [((row & -row).bit_length() - 1, row.bit_length() - 1) for row in rows]
    return [sum(first <= i <= last for first, last in endpoints) for i in range(pairs)]


def optimal_basis_controls():
    rng = random.Random(674319)
    checked = exhaustive = 0
    examples = []
    candidates = []
    for forms, maximum in ((1, 5), (2, 5), (3, 3)):
        for pairs in range(forms, maximum + 1):
            for columns in product(range(1 << forms), repeat=pairs):
                if rank_two(columns) == forms:
                    rows = [sum(((column >> a) & 1) << i for i, column in enumerate(columns)) for a in range(forms)]
                    candidates.append((rows, pairs, True))
    for _ in range(500):
        forms = rng.randrange(1, 6)
        pairs = rng.randrange(forms, 11)
        rows = [rng.randrange(1 << pairs) for _ in range(forms)]
        if rank_two(rows) == forms:
            candidates.append((rows, pairs, False))
    changes = {k: [rowset for rowset in product(range(1, 1 << k), repeat=k) if rank_two(rowset) == k]
               for k in range(1, 4)}
    for rows, pairs, complete in candidates:
        optimized, transform, reductions = minimal_span_basis(rows, pairs)
        forms = len(rows)
        for row, coefficient in zip(optimized, transform):
            expected = 0
            for i in bits(coefficient):
                expected ^= rows[i]
            assert row == expected
        lower = [rank_two([row & ((1 << (i + 1)) - 1) for row in rows]) +
                 rank_two([row & ~((1 << i) - 1) for row in rows]) - forms for i in range(pairs)]
        assert interval_profile(optimized, pairs) == lower
        if complete:
            for basis in changes[forms]:
                altered = []
                for coordinate in basis:
                    value = 0
                    for i in bits(coordinate):
                        value ^= rows[i]
                    altered.append(value)
                assert all(a >= b for a, b in zip(interval_profile(altered, pairs), lower))
                exhaustive += 1
        checked += 1
        if rows == [3, 12, 48]:
            examples.append({'rows': rows, 'optimized_rows': optimized, 'width': max(lower)})
    columns = (3, 3, 5, 5, 7, 7)
    rows = [sum(((column >> a) & 1) << i for i, column in enumerate(columns)) for a in range(3)]
    optimized, change, reductions = minimal_span_basis(rows, len(columns))
    return {'full_rank_form_matrices': checked, 'all_invertible_row_changes_checked': exhaustive,
            'example': {'columns': list(columns), 'original_rows': rows, 'optimized_rows': optimized,
                        'row_change': change, 'original_width': max(interval_profile(rows, len(columns))),
                        'optimal_width': max(interval_profile(optimized, len(columns))), 'trailing_reductions': reductions}}


def transport_controls():
    base = closed_operator_interface()
    inserted = insertion_transitions()
    ordinary = base['transport_transitions']
    incidence = inserted['physical_insertion_transitions_incidence']
    creation = inserted['physical_insertion_transitions_creation']
    basis = base['homogeneous_basis_matrices']
    matrices = [[[[(value >> position) & 1 for value in row] for row in table]
                 for table in (ordinary, incidence, creation)] for position in range(256)]
    rng = random.Random(880213)
    g_grades = ((1, 0), (0, 1), (2, 1), (1, 2))
    f_grades = ((1, 0), (0, 1), (0, 0), (1, 1))
    pair_grades = [(a + c, b + d) for a, b in f_grades for c, d in f_grades]
    def act(matrix, vector):
        answer = []
        for row in matrix:
            value = 0
            for bit, coefficient in zip(row, vector):
                if bit:
                    value ^= coefficient
            answer.append(value)
        return answer
    def acted(r, o, vector, scalar):
        first, second = act(r, vector), act(o, vector)
        return [a ^ (b if scalar else b << 1) for a, b in zip(first, second)]
    def full_ring(columns, positions, initial):
        forms = max(columns).bit_length()
        values = [initial] + [[0] * 8 for _ in range((1 << forms) - 1)]
        for column, position in zip(columns, positions):
            b, r, o = matrices[position]
            following = [act(b, vector) for vector in values]
            for mask, vector in enumerate(values):
                addition = acted(r, o, vector, False)
                for form in bits(column):
                    bit = 1 << form
                    if not mask & bit:
                        following[mask | bit] = [a ^ c for a, c in zip(following[mask | bit], addition)]
            values = following
        return values[-1][0]
    def active_ring(columns, positions, initial, scalar=False):
        forms = max(columns).bit_length()
        supports = [[i for i, column in enumerate(columns) if column >> a & 1] for a in range(forms)]
        if any(not support for support in supports):
            return 0, 0
        first, last = [s[0] for s in supports], [s[-1] for s in supports]
        active, values, peak = [], [initial], 0
        for site, (column, position) in enumerate(zip(columns, positions)):
            enlarged = sorted(active + [a for a in range(forms) if first[a] == site])
            lifted = [[0] * 8 for _ in range(1 << len(enlarged))]
            for mask, vector in enumerate(values):
                new_mask = sum(((mask >> j) & 1) << enlarged.index(a) for j, a in enumerate(active))
                lifted[new_mask] = vector
            active, values = enlarged, lifted
            peak = max(peak, len(active))
            b, r, o = matrices[position]
            following = [act(b, vector) for vector in values]
            for mask, vector in enumerate(values):
                addition = acted(r, o, vector, scalar)
                for a in bits(column):
                    bit = 1 << active.index(a)
                    if not mask & bit:
                        following[mask | bit] = [x ^ y for x, y in zip(following[mask | bit], addition)]
            remaining = [a for a in active if last[a] > site]
            required = sum(1 << j for j, a in enumerate(active) if last[a] == site)
            values = [[0] * 8 for _ in range(1 << len(remaining))]
            for mask, vector in enumerate(following):
                if mask & required == required:
                    new_mask = sum(((mask >> active.index(a)) & 1) << j for j, a in enumerate(remaining))
                    values[new_mask] = [x ^ y for x, y in zip(values[new_mask], vector)]
            active = remaining
        assert not active
        return values[0][0], peak
    cases = []
    for columns in ((1, 3, 2), (1, 3, 2, 6, 4), (1, 3, 2, 6, 4, 12, 8),
                    (1, 2, 4, 3, 6, 4), (3,), (1, 2, 7), (1, 0, 3, 0, 2), (3, 3, 5, 5, 7, 7)):
        pairs, forms = len(columns), max(columns).bit_length()
        choices = [set(selected) for selected in combinations(range(pairs), forms)
                   if rank_two([columns[i] for i in selected]) == forms]
        original_rows = [sum(((column >> a) & 1) << i for i, column in enumerate(columns)) for a in range(forms)]
        optimized_columns = None
        if rank_two(columns) == forms:
            optimized, _, _ = minimal_span_basis(original_rows, pairs)
            optimized_columns = [sum(((row >> i) & 1) << a for a, row in enumerate(optimized)) for i in range(pairs)]
        count = nonzero = grade_checks = 0
        attempts = 0
        width = max(interval_profile(original_rows, pairs))
        while count < 100:
            attempts += 1
            assert attempts < 100000
            selected = rng.choice(choices) if choices else set(range(min(forms, pairs)))
            positions, channel = [0] * pairs, 0
            for site in reversed(range(pairs)):
                tables = (incidence, creation) if site in selected else (ordinary,)
                options = [(b, table[channel][b]) for table in tables for b in range(8) if table[channel][b]]
                if not options:
                    break
                channel, value = rng.choice(options)
                positions[site] = rng.choice(list(bits(value)))
            else:
                boundary = rng.choice(list(bits(basis[channel])))
                initial = [(value >> boundary) & 1 for value in basis]
                expected = full_ring(columns, positions, initial)
                actual, peak = active_ring(columns, positions, initial)
                scalar, scalar_peak = active_ring(columns, positions, initial, True)
                assert actual == expected
                assert scalar == (expected.bit_count() & 1)
                assert peak == scalar_peak == width
                if optimized_columns is not None:
                    changed, _ = active_ring(optimized_columns, positions, initial)
                    assert changed == expected
                if expected:
                    assert expected.bit_count() == 1
                    nonzero += 1
                    degree = expected.bit_length() - 1
                    row, col = boundary // 4, boundary % 4
                    first = g_grades[row][0] - g_grades[col][0]
                    second = g_grades[row][1] - g_grades[col][1]
                    for position in positions:
                        row, col = position // 16, position % 16
                        first += pair_grades[row][0] - pair_grades[col][0]
                        second += pair_grades[row][1] - pair_grades[col][1]
                    assert first == forms + degree
                    assert abs(second) <= forms - degree and (second - forms + degree) % 2 == 0
                    grade_checks += 1
                count += 1
        cases.append({'columns': list(columns), 'forms': forms, 'form_rank': rank_two(columns),
                      'active_interval_width': width, 'operator_channel_bound': 8 * 2 ** width,
                      'full_ring_vs_active_memory_entries': count, 'scalar_vs_polynomial_entries': count,
                      'nonzero_entries': nonzero, 'grade_controls': grade_checks,
                      'optimized_columns': optimized_columns,
                      'basis_changed_entry_controls': count if optimized_columns is not None else 0})
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'check_ward_stable_operators.py',
                                            'check_ward_form_operator_transport.py')},
              'fixed_order_basis_optimization': optimal_basis_controls(),
              'interval_transport': transport_controls()}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
