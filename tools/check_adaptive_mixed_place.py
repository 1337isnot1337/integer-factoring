#!/usr/bin/env python3
"""Independent exact controls for adaptive opposite Hasse-jet flags.

Finite identities test the implementation; the research note supplies the
uniform proof. The rational interval example deliberately uses nonnested
comparison allocations. No improved factoring exponent is claimed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import lcm
from pathlib import Path

import adaptive_mixed_place as implementation
from check_mixed_place_hermite import (
    apply_functional, combine, jet_vector, nullspace,
)
from check_multipoint_hermite_crossings import (
    affine_seed, flatten, local_polynomial, local_profiles,
    mixing_seed, multiply, rank, solve, tensor_power,
)


def convolution_counts(orders, copies, cutoff):
    """Independent ordinary convolution, with no packed arithmetic."""
    exact = [1] + [0] * cutoff
    for _ in range(copies):
        updated = [0] * (cutoff + 1)
        for total, count in enumerate(exact):
            for order in orders:
                if total + order <= cutoff:
                    updated[total + order] += count
        exact = updated
    result, cumulative = [], 0
    for count in exact:
        cumulative += count
        result.append(cumulative)
    return result


def count_controls():
    profiles = [[0], [1], [0, 0], [0, 2], [0, 0, 1, 3],
                [2, 2, 7], [0, 1, 1, 4], [3, 5]]
    word_cases, words, convolution_cases, coefficients = 0, 0, 0, 0
    for orders in profiles:
        for copies in range(1, 6):
            cutoff = 17
            actual = implementation.channel_counts(orders, copies, cutoff)
            counts = Counter(map(sum, product(orders, repeat=copies)))
            brute = [sum(value for order, value in counts.items() if order <= j)
                     for j in range(cutoff + 1)]
            assert actual == brute == convolution_counts(orders, copies, cutoff)
            # A smaller truncation must agree with the same initial segment.
            for short in (0, 1, 5, 9):
                assert implementation.channel_counts(orders, copies, short) == brute[:short + 1]
            word_cases += 1
            words += len(orders) ** copies
            coefficients += cutoff + 1
    for orders in ([0, 0, 1, 3], [2, 2, 7], [0] * 3 + [1] * 3 + [2] * 3,
                   [0, 1, 1, 1, 1, 2, 4, 4, 4]):
        for copies in (6, 17, 31, 64):
            cutoff = 2 * copies + 3
            assert implementation.channel_counts(orders, copies, cutoff) == convolution_counts(
                orders, copies, cutoff)
            convolution_cases += 1
            coefficients += cutoff + 1
    return {'word_cases': word_cases, 'enumerated_words': words,
            'larger_convolution_cases': convolution_cases,
            'cumulative_coefficients_checked': coefficients}


def compositions(total, places):
    if places == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, places - 1):
            yield (first,) + tail


def greedy_controls():
    profile_sets = [
        [[0, 1]], [[0, 0], [1, 2]], [[0, 3], [0, 1]],
        [[2, 2], [1, 3]], [[0, 1], [0, 1], [0, 1]],
        [[0, 0, 2], [1, 1, 1], [0, 3, 3]],
        [[2, 2, 2], [0, 1, 2], [3, 3, 4]],
    ]
    flags, prefixes, allocations = 0, 0, 0
    for profiles in profile_sets:
        for copies in (1, 2, 4):
            for degree in range(8):
                flag = implementation.greedy_flag(profiles, copies, degree)
                independent = [convolution_counts(p, copies, degree) for p in profiles]
                costs = [[(j + 1) * count for j, count in enumerate(counts)]
                         for counts in independent]
                seen = [0] * len(profiles)
                cumulative = 0
                for size, jet in enumerate(flag, 1):
                    place, order = jet['place_index'], jet['order']
                    assert order == seen[place]
                    seen[place] += 1
                    assert jet['cost'] == costs[place][order]
                    cumulative += jet['cost']
                    assert jet['prefix_cost'] == min(cumulative, len(profiles[0]) ** copies)
                    candidate_costs = []
                    for counts in compositions(size, len(profiles)):
                        candidate_costs.append(sum(sum(c[:r]) for c, r in zip(costs, counts)))
                        allocations += 1
                    minimum = min(candidate_costs)
                    assert cumulative == minimum
                    assert jet['prefix_cost'] == min(minimum, len(profiles[0]) ** copies)
                    prefixes += 1
                flags += 1
    return {'flags': flags, 'prefixes': prefixes,
            'initial_segment_allocations': allocations}


def adaptive_basis(row_profiles, column_profiles, copies, degree,
                   row_points, column_points, prime):
    """Build the actual functional intersection, independently of core code."""
    assert set(row_points).isdisjoint(column_points)
    assert 1 not in row_points and 1 not in column_points
    row_flag = implementation.greedy_flag(row_profiles, copies, degree)
    column_flag = implementation.greedy_flag(column_profiles, copies, degree)
    for flag, points in ((row_flag, row_points), (column_flag, column_points)):
        for jet in flag:
            jet['point'] = points[jet['place_index']]
            jet['vector'] = jet_vector(degree, jet['point'], jet['order'], prime)
    rows = [jet['vector'] for jet in row_flag]
    columns = [jet['vector'] for jet in column_flag]
    assert rank(rows, prime) == rank(columns, prime) == degree + 1
    dimensions = 0
    for r in range(degree + 2):
        for c in range(degree + 2):
            actual = r + c - rank(rows[:r] + columns[:c], prime)
            assert actual == max(0, r + c - degree - 1)
            dimensions += 1
    basis = []
    for k in range(degree + 1):
        selected_rows, selected_columns = rows[:k + 1], columns[:degree - k + 1]
        matrix = [[v[j] for v in selected_rows] + [-v[j] % prime for v in selected_columns]
                  for j in range(degree + 1)]
        kernel = nullspace(matrix, prime)
        assert len(kernel) == 1
        left, right = kernel[0][:k + 1], kernel[0][k + 1:]
        vector = combine(left, selected_rows, prime)
        assert vector == combine(right, selected_columns, prime)
        inverse = pow(next(value for value in vector if value), -1, prime)
        vector = [value * inverse % prime for value in vector]
        left = [value * inverse % prime for value in left]
        right = [value * inverse % prime for value in right]
        assert left[-1] and right[-1]
        basis.append({'vector': vector, 'row_coefficients': left,
                      'column_coefficients': right})
    vectors = [entry['vector'] for entry in basis]
    assert rank(vectors, prime) == degree + 1
    weights = solve(list(map(list, zip(*vectors))), [1] * (degree + 1), prime)
    assert combine(weights, vectors, prime) == [1] * (degree + 1)
    return row_flag, column_flag, basis, weights, dimensions


def reconstruction_control(name, seed, alphabet, degree, copies, prime,
                           row_points, column_points):
    full_degree = degree * copies
    all_points = row_points + column_points
    profiles = {point: local_profiles(seed, alphabet, degree, point, prime)
                for point in all_points}
    row_profiles = [profiles[point][0] for point in row_points]
    column_profiles = [profiles[point][1] for point in column_points]
    rows, columns, basis, weights, dimensions = adaptive_basis(
        row_profiles, column_profiles, copies, full_degree, row_points, column_points, prime)
    table = tensor_power(seed, alphabet, copies, prime)
    locals_by_point = {
        point: {key: local_polynomial(value, point, full_degree, prime)
                for key, value in table.items()} for point in all_points}
    dimension, pieces, records, representations = alphabet ** copies, [], [], 0
    for k, (entry, weight) in enumerate(zip(basis, weights)):
        piece = {}
        for key, polynomial in table.items():
            value = weight * apply_functional(entry['vector'], polynomial, prime) % prime
            left = sum(c * locals_by_point[jet['point']][key][jet['order']]
                       for c, jet in zip(entry['row_coefficients'], rows)) % prime
            right = sum(c * locals_by_point[jet['point']][key][jet['order']]
                        for c, jet in zip(entry['column_coefficients'], columns)) % prime
            assert value == weight * left % prime == weight * right % prime
            representations += 2
            if value:
                piece[key] = [value]
        ranks = [rank(flatten(piece, dimension, cut), prime) for cut in (0, 1)]
        bounds = [rows[k]['prefix_cost'], columns[full_degree - k]['prefix_cost']]
        assert all(actual <= bound for actual, bound in zip(ranks, bounds))
        pieces.append(piece)
        records.append({'index': k, 'weight': weight, 'ranks': ranks,
                        'adaptive_prefix_bounds': bounds})
    for key in product(range(dimension), repeat=4):
        x, y, u, v = key
        expected = int(x == v and y == u)
        assert sum(table.get(key, ())) % prime == expected
        assert sum(piece.get(key, [0])[0] for piece in pieces) % prime == expected
    return {'name': name, 'prime': prime, 'alphabet': alphabet,
            'seed_degree_bound': degree, 'copies': copies,
            'row_profiles': [{'point': p, 'orders': profiles[p][0]} for p in row_points],
            'column_profiles': [{'point': p, 'orders': profiles[p][1]} for p in column_points],
            'row_flag': rows, 'column_flag': columns, 'pieces': records,
            'complete_endpoint_entries': dimension ** 4,
            'representation_entries': representations,
            'flag_dimension_checks': dimensions, 'exact_rank_checks': 2 * len(pieces)}


def scalar_controls():
    original = [[0, 0, 0, 1, 1, 1, 2, 2, 2]] * 2
    records = []
    for shifts in ([1, 0], [0, 2], [1, 2]):
        shifted = [[order + shift for order in profile]
                   for profile, shift in zip(original, shifts)]
        for copies in (1, 2, 5, 13):
            degree = 2 * copies
            added_free = copies * sum(shifts)
            before = implementation.greedy_flag(original, copies, degree)
            after = implementation.greedy_flag(shifted, copies, degree + added_free)
            assert all(jet['cost'] == 0 for jet in after[:added_free])
            polynomial_factor = 1 + copies * max(shifts)
            for i, jet in enumerate(before):
                old, new = jet['prefix_cost'], after[i + added_free]['prefix_cost']
                assert old <= new <= polynomial_factor * old
            records.append({'copies': copies, 'local_content_orders': list(shifts),
                            'added_zero_jets': added_free,
                            'prefix_comparisons': len(before),
                            'cost_factor_bound': polynomial_factor})
    return records


def rational_certificate():
    """An independent all-local-pairs checker, not the proposal verifier."""
    row = [[0] * 3 + [1] * 3 + [2] * 3] * 2
    column = [[0, 1, 1, 1, 1, 2, 4, 4, 4]]
    certificate = {
        'alphabet': 3, 'degree': 2, 'rho': '99/100',
        'row_profiles': row, 'column_profiles': column,
        'intervals': [
            {'left': '0', 'right': '1/2',
             'row': [(0, '1/2', '2')], 'column': [(0, '2', '1')]},
            {'left': '1/2', 'right': '1',
             'row': [(1, '1', '1')], 'column': [(0, '3/2', '5/4')]},
            {'left': '1', 'right': '3/2',
             'row': [(0, '3/4', '4/3'), (1, '3/4', '4/3')],
             'column': [(0, '1', '3/2')]},
            {'left': '3/2', 'right': '2',
             'row': [(0, '1', '1'), (1, '1', '1')],
             'column': [(0, '1/2', '2')]},
        ],
    }

    def verify(data):
        q, degree, rho = data['alphabet'], data['degree'], Fraction(data['rho'])
        assert q >= 2 and Fraction(1, 2) <= rho < 1
        previous, comparisons, records = Fraction(0), 0, []
        for interval in data['intervals']:
            left, right = Fraction(interval['left']), Fraction(interval['right'])
            assert left == previous and left < right <= degree
            previous = right
            parsed = {}
            exponent = (4 * rho).denominator
            for side, target in (('row', right), ('column', degree - left)):
                parsed[side] = [(i, Fraction(a), Fraction(s)) for i, a, s in interval[side]]
                assert sum(a for _, a, _ in parsed[side]) >= target > 0
                assert len(set(i for i, _, _ in parsed[side])) == len(parsed[side])
                for i, amount, scale in parsed[side]:
                    assert 0 <= i < len(data[side + '_profiles']) and amount > 0 and scale >= 1
                    exponent = lcm(exponent, amount.denominator)
            interval_records = []
            for i, amount, scale in parsed['row']:
                u = sum((scale ** -order for order in data['row_profiles'][i]), Fraction())
                for j, other_amount, other_scale in parsed['column']:
                    v = sum((other_scale ** -order for order in data['column_profiles'][j]), Fraction())
                    bound = ((u * v) ** exponent * scale ** int(exponent * amount)
                             * other_scale ** int(exponent * other_amount))
                    target = q ** int(4 * rho * exponent)
                    assert bound < target
                    interval_records.append({
                        'row_place': i, 'column_place': j, 'common_power': exponent,
                        'left_numerator': str(bound.numerator),
                        'right_after_clearing_denominator': str(target * bound.denominator)})
                    comparisons += 1
            records.append(interval_records)
        assert previous == degree
        return comparisons, records

    comparisons, powered = verify(certificate)
    # Place zero's allocation drops from 1/2 to zero at the next checkpoint.
    assert certificate['intervals'][0]['row'][0][0] == 0
    assert all(i != 0 for i, _, _ in certificate['intervals'][1]['row'])
    invalid = deepcopy(certificate)
    invalid['intervals'][1]['column'][0] = (0, '1', '5/4')
    try:
        verify(invalid)
    except AssertionError:
        rejected = True
    else:
        raise AssertionError('an uncovered complementary budget was accepted')
    return {'certificate': certificate, 'exact_local_pair_comparisons': comparisons,
            'powered_inequalities': powered, 'nonnested_allocation_accepted': True,
            'underfunded_complementary_budget_rejected': rejected}


def audit_exported_certificate(path, profile_path):
    """Check an exported witness against its exact source and all intervals.

    This does not call the generator's verifier or any numerical routine.
    Finite-field Smith-profile verification belongs to the cited source
    checker; here every selected profile must match that source verbatim.
    """
    certificate = json.loads(path.read_text())
    source_bytes = profile_path.read_bytes()
    assert certificate['source_profiles_sha256'] == sha256(source_bytes).hexdigest()
    cases = json.loads(source_bytes)['cases']
    seed = certificate['seed']
    source = next(case for case in cases if all(case[key] == value for key, value in seed.items()))
    assert certificate['alphabet'] == source['alphabet']
    assert certificate['degree'] == source['degree']
    q, degree = certificate['alphabet'], certificate['degree']
    rho = Fraction(*certificate['rho'])
    assert q >= 2 and degree >= 1 and Fraction(1, 2) <= rho < 1
    by_point = {json.dumps(profile['point']): profile for profile in source['profiles']}
    points = {}
    for side, cut in (('row', 0), ('column', 1)):
        profiles = certificate[side + '_profiles']
        points[side] = [json.dumps(profile['point']) for profile in profiles]
        assert profiles and len(set(points[side])) == len(profiles)
        for profile in profiles:
            assert profile['orders'] == by_point[json.dumps(profile['point'])]['cuts'][cut]
            assert len(profile['orders']) == q * q
            assert all(isinstance(order, int) and order >= 0 for order in profile['orders'])
    assert set(points['row']).isdisjoint(points['column'])
    previous, bounds, pairs, largest_bits = Fraction(), 0, 0, 0
    nonnested = False
    previous_amounts = None
    for interval in certificate['intervals']:
        left, right = Fraction(*interval['left']), Fraction(*interval['right'])
        assert left == previous and left < right <= degree
        previous = right
        exponent = (4 * rho).denominator
        parsed, current_amounts = {}, {}
        for side, required in (('row', right), ('column', degree - left)):
            profiles = certificate[side + '_profiles']
            allocations = interval[side]
            parsed[side] = []
            seen = set()
            amounts = [Fraction()] * len(profiles)
            for item in allocations:
                index = item['profile']
                amount, scale = Fraction(*item['amount']), Fraction(*item['scale'])
                assert isinstance(index, int) and 0 <= index < len(profiles) and index not in seen
                assert amount >= 0 and scale >= 1
                seen.add(index)
                amounts[index] = amount
                exponent = lcm(exponent, amount.denominator)
                parsed[side].append((profiles[index]['orders'], amount, scale))
            assert sum(amounts) >= required > 0
            current_amounts[side] = amounts
        if previous_amounts is not None:
            # Row budgets grow with the interval, but their individual place
            # allocations are free to decrease. Columns have the reverse order.
            nonnested |= any(a < b for a, b in zip(current_amounts['row'], previous_amounts['row']))
            nonnested |= any(a > b for a, b in zip(current_amounts['column'], previous_amounts['column']))
        previous_amounts = current_amounts
        maxima = []
        for side in ('row', 'column'):
            local_bounds = []
            for orders, amount, scale in parsed[side]:
                partition = sum((scale ** -order for order in orders), Fraction())
                value = partition ** exponent * scale ** int(exponent * amount)
                local_bounds.append(value)
                largest_bits = max(largest_bits, value.numerator.bit_length(), value.denominator.bit_length())
                bounds += 1
            maxima.append(max(local_bounds))
        target_power = 4 * rho * exponent
        assert target_power.denominator == 1
        assert maxima[0] * maxima[1] < q ** int(target_power)
        pairs += len(parsed['row']) * len(parsed['column'])
    assert previous == degree
    return {'status': 'pass', 'certificate_sha256': sha256(path.read_bytes()).hexdigest(),
            'source_profiles_sha256': sha256(source_bytes).hexdigest(),
            'seed': seed, 'rho': str(rho), 'intervals': len(certificate['intervals']),
            'powered_local_bounds': bounds, 'implied_local_pair_comparisons': pairs,
            'largest_integer_bits': largest_bits,
            'comparison_allocations_are_nonnested': nonnested,
            'scope': 'Exact interval and source-binding audit; the source checker establishes the local profiles.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--certificate', type=Path,
                        help='also independently audit an exported rational interval certificate')
    parser.add_argument('--profiles', type=Path,
                        default=Path(__file__).resolve().parent.parent / 'evidence/fused-closed-profiles.json')
    args = parser.parse_args()
    counts, greedy = count_controls(), greedy_controls()
    seed = mixing_seed()
    with_content = {key: multiply(value, [2, 2], 3) for key, value in seed.items()}
    affine_content = {key: multiply(value, [2, 4], 5)
                      for key, value in affine_seed(5).items()}
    cases = [
        reconstruction_control('binary affine', affine_seed(2), 2, 1, 3, 2, [None], [0]),
        reconstruction_control('ternary', seed, 3, 2, 1, 3, [2, None], [0]),
        reconstruction_control('ternary tensor square', seed, 3, 2, 2, 3, [2, None], [0]),
        reconstruction_control('ternary scalar content', with_content, 3, 3, 2, 3, [2, None], [0]),
        reconstruction_control('affine scalar content', affine_content, 2, 2, 2, 5, [0, 2], [None, 3]),
        reconstruction_control('surplus degree at infinity', affine_seed(5), 2, 3, 2, 5, [0, 2], [None, 3]),
    ]
    content_case = cases[3]
    assert content_case['row_profiles'][0]['orders'] == [v + 1 for v in cases[2]['row_profiles'][0]['orders']]
    assert [jet['cost'] for jet in content_case['row_flag'][:2]] == [0, 0]
    assert sum(piece['ranks'] != [0, 0] for piece in content_case['pieces'][:2]) == 0
    assert all(jet['cost'] == 0 for jet in cases[5]['column_flag'][:4])
    totals = {key: sum(case[key] for case in cases) for key in
              ('complete_endpoint_entries', 'representation_entries',
               'flag_dimension_checks', 'exact_rank_checks')}
    result = {'schema_version': 1, 'status': 'pass',
              'claim': 'Exact adaptive flags, original-tensor crossing controls, content accounting, and a uniform rational interval example; factoring score unchanged.',
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in
                                (Path(__file__), Path(implementation.__file__))},
              'count_controls': counts, 'greedy_optimality_controls': greedy,
              'reconstruction_cases': cases, 'scalar_content_controls': scalar_controls(),
              'rational_interval_control': rational_certificate(), 'totals': totals}
    if args.certificate:
        result['exported_certificate_audit'] = audit_exported_certificate(args.certificate, args.profiles)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'totals': totals,
                      'count_controls': counts, 'greedy_controls': greedy}, sort_keys=True))


if __name__ == '__main__':
    main()
