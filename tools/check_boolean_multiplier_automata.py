#!/usr/bin/env python3
"""Exact independent controls for the unbounded multiplier reduction."""

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import prod
from pathlib import Path
import random

import numpy as np

from boolean_multiplier_pressure import MultiplierEnvelope
from boolean_multiplier_profiles import (
    code_group_signatures_for_types, complete_multiplier_automaton,
    finite_degree_representatives, homogeneous_layer_patterns, layer_pattern_excess,
    minimal_divergence_multiplier, multiplier_profile, multiply, order,
    rational_target, series_prefix, signature_patterns,
)
from check_symmetric_binary_filter_synthesis import full_transform, xor
from search_global_code_multipliers import CODE_MATRICES


def coefficient_dictionary(polynomial, spacing=1):
    return {bit * spacing: 1 for bit in range(polynomial.bit_length())
            if polynomial & (1 << bit)}


def finite_covers():
    comparisons = represented = 0
    for boolean_degree in range(1, 5):
        for maximum_degree in range(1, 10):
            exhaustive = set()
            for polynomial in range(1, 1 << (maximum_degree + 1)):
                if polynomial.bit_count() % 2:
                    continue
                profile = multiplier_profile(polynomial, boolean_degree)
                exhaustive.add((profile['degree'], profile['excesses']))
                represented += 1
            compressed = set()
            for polynomial in finite_degree_representatives(boolean_degree, maximum_degree):
                profile = multiplier_profile(polynomial, boolean_degree)
                compressed.add((profile['degree'], profile['excesses']))
            assert exhaustive == compressed
            comparisons += 1
    return comparisons, represented


def route_integral_multipliers():
    counts = Counter()
    for boolean_degree in range(1, 6):
        automaton = complete_multiplier_automaton(boolean_degree)
        tails = {tail['pattern']: tail for tail in automaton['singleton_tails']}
        patterns = tuple(range(1, 1 << (boolean_degree + 1)))
        targets = {pattern: series_prefix(pattern, 16) for pattern in patterns}
        for polynomial in range(1, 1 << 13):
            if polynomial.bit_count() % 2:
                continue
            counts['integral_multiplier_routes'] += 1
            original = multiplier_profile(polynomial, boolean_degree)
            active, depth, prefix = patterns, 0, 0
            while len(active) > 1:
                coefficient = (polynomial >> depth) & 1
                prefix |= coefficient << depth
                active = tuple(pattern for pattern in active
                               if ((targets[pattern] >> depth) & 1) == coefficient)
                depth += 1
                if not active:
                    candidate = prefix
                    if candidate.bit_count() % 2:
                        candidate ^= 1 << depth
                    if not candidate:
                        assert not any(original['excesses'])
                        counts['baseline_routes'] += 1
                    else:
                        reduced = multiplier_profile(candidate, boolean_degree)
                        assert reduced['excesses'] == original['excesses']
                        assert reduced['degree'] <= original['degree']
                        counts['finite_branch_routes'] += 1
                    break
            if not active:
                continue
            pattern = active[0]
            degree, denominator = rational_target(pattern)
            divergence = order((1 << degree) ^ multiply(polynomial, denominator))
            assert divergence >= depth
            candidate = minimal_divergence_multiplier(pattern, divergence)
            assert candidate
            reduced = multiplier_profile(candidate, boolean_degree)
            assert reduced['excesses'] == original['excesses']
            assert reduced['degree'] <= original['degree']
            offset = reduced['degree'] - divergence
            if offset >= 0:
                assert divergence <= original['degree']
                counts['continuous_dominance_routes'] += 1
            else:
                tail = tails[pattern]
                if divergence < tail['cycle_start']:
                    row = next(row for row in tail['startup'] if row['divergence'] == divergence)
                else:
                    base = tail['cycle_start'] + (divergence - tail['cycle_start']) % tail['period']
                    row = next(row for row in tail['cycle'] if row['divergence'] == base)
                assert row['degree_offset'] == offset
                counts['negative_gap_routes'] += 1
    return counts


def periodic_tails():
    count = 0
    manifests = []
    for boolean_degree in range(1, 7):
        automaton = complete_multiplier_automaton(boolean_degree)
        for tail in automaton['singleton_tails']:
            for part in ('startup', 'cycle'):
                for row in tail['full_tail'][part]:
                    for repeat in ((0, 1, 2) if part == 'cycle' else (0,)):
                        divergence = row['divergence'] + repeat * tail['period']
                        candidate = minimal_divergence_multiplier(tail['pattern'], divergence)
                        if not candidate:
                            assert row['degree_offset'] is None
                        else:
                            assert candidate.bit_count() % 2 == 0
                            assert candidate.bit_length() - 1 - divergence == row['degree_offset']
                            degree, denominator = rational_target(tail['pattern'])
                            assert order((1 << degree) ^ multiply(candidate, denominator)) == divergence
                        count += 1
        manifests.append({
            'maximum_boolean_degree': boolean_degree,
            'continuous_pieces': len(automaton['continuous_profiles']),
            'negative_gap_periodic_rows': sum(len(tail['cycle']) for tail in automaton['singleton_tails']),
            'maximum_period': max(tail['period'] for tail in automaton['singleton_tails']),
        })
    return count, manifests


def literal_colored_cutoffs(generator):
    states = witnesses = 0
    for modes in (2, 3, 4):
        for _ in range(4):
            truth = tuple(generator.randrange(2) for _ in range(1 << modes))
            if len(set(truth)) == 1:
                truth = (truth[0] ^ 1,) + truth[1:]
            layers = homogeneous_layer_patterns(truth, modes)
            boolean_degree = max(pattern.bit_length() - 1 for pattern in layers)
            automaton = complete_multiplier_automaton(boolean_degree)
            for row in generator.sample(list(automaton['continuous_profiles']),
                                        min(12, len(automaton['continuous_profiles']))):
                for fractional in (Fraction(0), Fraction(1, 3), Fraction(1)):
                    cutoff = row['lower_degree'] + fractional
                    if row['upper_degree'] is not None:
                        assert cutoff <= row['upper_degree']
                    if len(row['active_patterns']) == 1:
                        prefix = series_prefix(row['active_patterns'][0], cutoff.numerator // cutoff.denominator + 1)
                    else:
                        prefix = row['prefix']
                    spacing = 12
                    exponent = int(spacing * cutoff)
                    multiplier = coefficient_dictionary(prefix, spacing)
                    multiplier[exponent] = multiplier.get(exponent, 0) ^ 2
                    f = []
                    for value in truth:
                        polynomial = {0: 1}
                        if value:
                            for power, coefficient in multiplier.items():
                                polynomial[power] = polynomial.get(power, 0) ^ coefficient
                                if not polynomial[power]:
                                    del polynomial[power]
                        assert xor(polynomial.values()) in (1, 2, 3)
                        f.append(polynomial)
                    beta = full_transform(tuple(f), modes, spacing)
                    for state, pattern in enumerate(layers):
                        if pattern in row['active_patterns']:
                            expected = spacing * state.bit_count() + exponent - spacing * (pattern.bit_length() - 1)
                        else:
                            expected = spacing * (state.bit_count() - row['outside_exponents'][pattern])
                        assert min(beta[state]) == expected
                        assert max(f[state]) == exponent * truth[state]
                        assert xor(beta[state].values()) == xor(f[state].values())
                        states += 1
                    witnesses += 1
    return witnesses, states


def code_layer_controls():
    controls = 0
    for rank, columns in CODE_MATRICES:
        types = tuple(sorted(set(columns)))
        signatures = code_group_signatures_for_types(rank, types)
        syndromes, parity_states = [], []
        for state in range(1 << len(columns)):
            syndrome = parity_state = 0
            for bit, column in enumerate(columns):
                if state & (1 << bit):
                    syndrome ^= column
                    parity_state ^= 1 << types.index(column)
            syndromes.append(syndrome)
            parity_states.append(parity_state)
        assert set(syndromes) == set(range(1 << rank))
        for syndrome in range(1 << rank):
            truth = tuple(int(value == syndrome) for value in syndromes)
            layers = homogeneous_layer_patterns(truth, len(columns))
            grouped = signature_patterns(signatures, rank, syndrome)
            assert all(pattern == grouped[parity_state] for pattern, parity_state in zip(layers, parity_states))
            controls += len(layers)
    return controls


def exact_heterogeneous_distributions(generator):
    controls = 0
    denominator = 101
    for matrix_index in (0, 1, 4, 8, 10, 12, 13, 14):
        rank, columns = CODE_MATRICES[matrix_index]
        types = tuple(sorted(set(columns)))
        signatures = code_group_signatures_for_types(rank, types)
        numerators = tuple(generator.randrange(1, denominator) for _ in columns)
        groups = [[index for index, column in enumerate(columns) if column == group] for group in types]
        compressed_weights = []
        group_weights = []
        for group in groups:
            total = denominator ** len(group)
            signed = prod(denominator - 2 * numerators[index] for index in group)
            assert (total - signed) % 2 == 0
            group_weights.append(((total + signed) // 2, (total - signed) // 2))
        for state in range(1 << len(types)):
            compressed_weights.append(prod(pair[(state >> bit) & 1] for bit, pair in enumerate(group_weights)))
        assert sum(compressed_weights) == denominator ** len(columns)
        actual_syndromes = []
        physical_weights = []
        for state in range(1 << len(columns)):
            syndrome = 0
            for bit, column in enumerate(columns):
                if state & (1 << bit):
                    syndrome ^= column
            actual_syndromes.append(syndrome)
            physical_weights.append(prod(numerators[bit] if state & (1 << bit)
                                         else denominator - numerators[bit] for bit in range(len(columns))))
        for syndrome in range(1 << rank):
            truth = tuple(int(value == syndrome) for value in actual_syndromes)
            layers = homogeneous_layer_patterns(truth, len(columns))
            direct = Counter()
            for pattern, weight in zip(layers, physical_weights):
                direct[pattern] += weight
            grouped = Counter()
            for pattern, weight in zip(signature_patterns(signatures, rank, syndrome), compressed_weights):
                grouped[pattern] += weight
            assert direct == grouped
            controls += 1
    return controls


def numerical_envelope_controls(generator):
    cases = 0
    for boolean_degree in range(1, 6):
        envelope = MultiplierEnvelope(boolean_degree)
        for _ in range(8):
            q = np.array([generator.random() + .03 for _ in envelope.pattern_degrees])
            q /= q.sum()
            p, tau = generator.uniform(.03, .97), generator.uniform(.02, 3)
            value, selected, divergence = envelope.evaluate([q], [p], tau)
            finite_minimum = 0.0
            for multiplier in range(1, 1 << 10):
                if multiplier.bit_count() % 2:
                    continue
                profile = multiplier_profile(multiplier, boolean_degree)
                cost = float(np.log(q @ np.exp(-2 * tau * np.array(profile['excesses'])))
                             + np.log(1 - p + p * np.exp(2 * tau * profile['degree'])))
                finite_minimum = min(finite_minimum, cost)
            assert value[0] <= finite_minimum + 1e-10
            if selected[0] >= 0:
                row = envelope.describe(int(selected[0]), float(divergence[0]))
                if row['kind'] in ('startup', 'cycle'):
                    multiplier = minimal_divergence_multiplier(row['pattern'], int(divergence[0]))
                    profile = multiplier_profile(multiplier, boolean_degree)
                    cost = float(np.log(q @ np.exp(-2 * tau * np.array(profile['excesses'])))
                                 + np.log(1 - p + p * np.exp(2 * tau * profile['degree'])))
                    assert abs(cost - value[0]) < 1e-9
            cases += 1
        envelope.exponential_arrays.cache_clear()
    return cases


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root / 'evidence' / 'boolean-multiplier-automata.json')
    args = parser.parse_args()
    generator = random.Random(801113)
    covers, represented = finite_covers()
    routes = route_integral_multipliers()
    tail_controls, manifests = periodic_tails()
    witnesses, states = literal_colored_cutoffs(generator)
    group_controls = code_layer_controls()
    distributions = exact_heterogeneous_distributions(generator)
    numerical = numerical_envelope_controls(generator)
    result = {
        'status': 'Exact finite covers, unbounded-tail routes, literal F4 cutoffs, and code distributions passed.',
        'finite_cover_comparisons': covers, 'exhaustively_represented_multipliers': represented,
        **routes, 'periodic_tail_controls': tail_controls,
        'literal_colored_witnesses': witnesses, 'literal_colored_state_controls': states,
        'independent_code_homogeneous_layer_controls': group_controls,
        'exact_heterogeneous_pattern_distributions': distributions,
        'floating_envelope_sanity_cases': numerical, 'automaton_sizes': manifests,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
