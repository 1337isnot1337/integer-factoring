#!/usr/bin/env python3
"""Original-state and exhaustive endpoint controls for joint syndrome filters."""

import argparse
from fractions import Fraction
import json
from math import ceil, log2
from pathlib import Path
import random

from check_symmetric_binary_filter_synthesis import full_transform, xor
from joint_syndrome_filter_synthesis import (
    add_equation, basis_release_targets, coefficient_tables, greedy_synthesize,
    synthesize, synthesize_compressed,
)
from weyl_lattice_filter_synthesis import ones_endpoint_possible


def brute_feasible(rank, types, spacing, degrees, targets):
    variables, masks, equations, _ = coefficient_tables(rank, types, spacing, degrees)
    rows = tuple(row for target, values in zip(targets, equations) for row in values[:target])
    endpoint_support = 0
    for candidate in range(1 << len(variables)):
        if any((candidate & row).bit_count() & 1 for row in rows):
            continue
        for syndrome, mask in enumerate(masks):
            if (candidate & mask).bit_count() & 1:
                endpoint_support |= 1 << syndrome
    return endpoint_support == (1 << (1 << rank)) - 1


def exhaustive_controls(generator):
    count = feasible = 0
    for rank, types in ((1, (1,)), (2, (1, 2, 3))):
        for spacing in (1, 3):
            for _ in range(40):
                degrees = tuple(generator.randrange(3) for _ in range(1 << rank))
                maximum = max(degrees) + spacing * rank
                if count % 4 == 0:
                    targets = (spacing * rank,) * (1 << len(types))
                else:
                    targets = tuple(generator.randrange(maximum + 2) for _ in range(1 << len(types)))
                expected = brute_feasible(rank, types, spacing, degrees, targets)
                actual = synthesize(rank, types, spacing, degrees, targets)
                assert bool(actual) == expected
                feasible += int(expected)
                count += 1
    return count, feasible


def endpoint_color_control():
    rank, types, spacing = 2, (1, 2, 3), 1
    degrees, targets = (0, 1, 2, 0), (2, 2, 2, 3, 2, 2, 2, 3)
    witness = synthesize(rank, types, spacing, degrees, targets)
    assert witness
    _, masks, equations, _ = coefficient_tables(rank, types, spacing, degrees)
    pivots, remainders = {}, masks
    for target, rows in zip(targets, equations):
        for row in rows[:target]:
            remainders = add_equation(row, pivots, remainders)
            assert remainders
    assert not ones_endpoint_possible(remainders)
    assert all(witness['endpoint_values'])
    return {key: witness[key] for key in ('degree_bounds', 'order_bounds', 'f', 'delta',
                                         'actual_degrees', 'actual_orders', 'endpoint_values')}


def state_data(columns, types):
    syndromes, parity_states = [], []
    for state in range(1 << len(columns)):
        syndrome = parity = 0
        for bit, column in enumerate(columns):
            if state & (1 << bit):
                syndrome ^= column
                parity ^= 1 << types.index(column)
        syndromes.append(syndrome)
        parity_states.append(parity)
    return tuple(syndromes), tuple(parity_states)


def original_polynomial_controls(generator):
    state_count = witness_count = rational_cost_count = 0
    sources = ((2, (1, 2, 3, 1)), (2, (1, 2, 2, 3, 3)),
               (3, (1, 2, 4, 3, 5)), (3, (1, 2, 4, 7, 1)))
    for rank, columns in sources:
        types = tuple(sorted(set(columns)))
        syndromes, parity_states = state_data(columns, types)
        for spacing in (1, 3):
            for query in range(8):
                degrees = tuple(generator.randrange(0, 3 * spacing + 1) for _ in range(1 << rank))
                probabilities = tuple(generator.random() + .01 for _ in range(1 << len(types)))
                tau = generator.uniform(.1, 5)
                witness = greedy_synthesize(rank, types, spacing, degrees, probabilities, tau,
                                            seed=830719 + witness_count, priority_noise=query % 3)
                replay = synthesize(rank, types, spacing, degrees, witness['order_bounds'])
                assert replay
                for key in ('f', 'delta', 'actual_orders', 'actual_degrees', 'endpoint_values'):
                    assert replay[key] == witness[key]
                f = tuple(witness['f'][syndrome] for syndrome in syndromes)
                beta = full_transform(f, len(columns), spacing)
                for state, parity in enumerate(parity_states):
                    shift = spacing * (state.bit_count() - rank)
                    expected = {degree + shift: coefficient
                                for degree, coefficient in witness['delta'][parity].items()}
                    assert all(degree >= 0 for degree in expected)
                    assert beta[state] == expected
                    assert xor(beta[state].values()) == witness['endpoint_values'][syndromes[state]]
                    state_count += 1
                if query < 3:
                    n = 1 << ceil(log2(len(columns)))
                    energies = tuple(range(2 * n - 1, 2 * n - 2 * len(columns), -2))
                    parameter = Fraction(49, 50)
                    for charge in (-2, 0, 2):
                        p_weights = [parameter ** (spacing * sum(energy for bit, energy in enumerate(energies)
                                                                  if state & (1 << bit)))
                                     for state in range(1 << len(columns))]
                        q_weights = [parameter ** (spacing * sum(energy + charge for bit, energy in enumerate(energies)
                                                                  if state & (1 << bit)))
                                     for state in range(1 << len(columns))]
                        original_zero = sum(weight * parameter ** (2 * n * (min(beta[state]) - spacing * state.bit_count()))
                                            for state, weight in enumerate(q_weights)) / sum(q_weights)
                        original_infinity = sum(weight * parameter ** (-2 * n * max(f[state]))
                                                for state, weight in enumerate(p_weights)) / sum(p_weights)
                        q_groups = [Fraction(0)] * (1 << len(types))
                        p_syndromes = [Fraction(0)] * (1 << rank)
                        for state in range(1 << len(columns)):
                            q_groups[parity_states[state]] += q_weights[state]
                            p_syndromes[syndromes[state]] += p_weights[state]
                        reduced_zero = sum(weight * parameter ** (2 * n * (witness['actual_orders'][state] - spacing * rank))
                                           for state, weight in enumerate(q_groups)) / sum(q_groups)
                        reduced_infinity = sum(weight * parameter ** (-2 * n * witness['actual_degrees'][syndrome])
                                               for syndrome, weight in enumerate(p_syndromes)) / sum(p_syndromes)
                        assert original_zero == reduced_zero and original_infinity == reduced_infinity
                        rational_cost_count += 1
                witness_count += 1
    return witness_count, state_count, rational_cost_count


def residue_compression_controls(generator):
    comparisons = 0
    for rank, types in ((2, (1, 2, 3)), (3, tuple(range(1, 8)))):
        for spacing in (1, 3, 5, 17):
            for query in range(16):
                degrees = tuple(generator.randrange(2 * spacing + 1) for _ in range(1 << rank))
                targets = ((spacing * rank,) * (1 << len(types)) if query % 4 == 0 else
                           tuple(generator.randrange((rank + 2) * spacing + 1)
                                 for _ in range(1 << len(types))))
                full = synthesize(rank, types, spacing, degrees, targets)
                compressed = synthesize_compressed(rank, types, spacing, degrees, targets)
                assert bool(full) == bool(compressed)
                if full:
                    assert full['constraint_rank'] == compressed['constraint_rank']
                    assert full['endpoint_image_dimension'] == compressed['endpoint_image_dimension']
                comparisons += 1
    large = []
    base_degrees, base_targets = (0, 1, 2, 0), (2, 2, 2, 3, 2, 2, 2, 3)
    for spacing in (127, 257, 65537):
        witness = synthesize_compressed(2, (1, 2, 3), spacing,
                                       tuple(spacing * value for value in base_degrees),
                                       tuple(spacing * value for value in base_targets))
        assert witness
        assert witness['actual_degrees'] == tuple(spacing * value for value in base_degrees)
        assert witness['actual_orders'] == tuple(spacing * value for value in base_targets)
        syndromes, parity_states = state_data((1, 2, 3, 1), (1, 2, 3))
        full = tuple(witness['f'][syndrome] for syndrome in syndromes)
        beta = full_transform(full, 4, spacing)
        for state, parity in enumerate(parity_states):
            expected = {degree + spacing * (state.bit_count() - 2): coefficient
                        for degree, coefficient in witness['delta'][parity].items()}
            assert beta[state] == expected
        large.append({key: witness[key] for key in (
            'spacing', 'residue_pattern_count', 'coefficient_variables',
            'maximum_local_coefficient_variables', 'endpoint_values')})
    return comparisons, large


def target_release_control():
    rank, types, degrees = 3, tuple(range(1, 8)), (4, 0, 0, 0, 0, 0, 0, 0)
    free_mask = (1 << 0) | (1 << 1) | (1 << 3)
    probabilities = tuple(Fraction(1, 8) if state & ~free_mask == 0 else Fraction(0)
                          for state in range(128))
    source = greedy_synthesize(rank, types, 1, degrees, probabilities, .1)
    releases = basis_release_targets(rank, types, 1, degrees, source['order_bounds'])
    parameter = Fraction(9, 10)  # exp(-2*tau) for this exact cost comparison.

    def cost(witness):
        return (sum(probability * parameter ** (value - rank)
                    for probability, value in zip(probabilities, witness['actual_orders']))
                * sum(parameter ** (-value) for value in witness['actual_degrees']) / 8)

    original = cost(source)
    assert original > 1
    best = None
    for release in releases:
        assert synthesize(rank, types, 1, degrees, release['order_bounds'], construct=False)
        for state, probability in enumerate(probabilities):
            if not probability:
                continue
            targets = list(release['order_bounds'])
            targets[state] = source['order_bounds'][state] + 1
            witness = greedy_synthesize(rank, types, 1, degrees, probabilities, .1,
                                         initial_targets=targets)
            if witness and (best is None or cost(witness) < best[0]):
                best = cost(witness), release['released_pivot'], state, witness
    assert best and best[0] < 1
    return {
        'scope': 'Uniform-bias algorithmic control, not a physical factoring candidate.',
        'feasible_basis_releases': len(releases),
        'original_relative_cost': [original.numerator, original.denominator],
        'improved_relative_cost': [best[0].numerator, best[0].denominator],
        'released_pivot': best[1], 'forced_output_state': best[2],
        'original_f': source['f'], 'improved_f': best[3]['f'],
    }


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root / 'evidence' / 'joint-syndrome-filter-synthesis.json')
    args = parser.parse_args()
    generator = random.Random(830717)
    decisions, feasible = exhaustive_controls(generator)
    endpoint_control = endpoint_color_control()
    witnesses, states, rational_costs = original_polynomial_controls(generator)
    compressed, large_spacings = residue_compression_controls(generator)
    release_control = target_release_control()
    result = {
        'status': 'Joint coefficient feasibility, original-state transforms, endpoints, and full costs passed.',
        'exhaustive_binary_kernel_endpoint_decisions': decisions,
        'feasible_exhaustive_cases': feasible,
        'all_ones_endpoint_restriction_negative_control': endpoint_control,
        'canonical_literal_witness_replays': witnesses,
        'independent_original_polynomial_state_controls': states,
        'exact_rational_full_cost_controls': rational_costs,
        'independent_residue_compression_comparisons': compressed,
        'large_spacing_literal_controls': large_spacings,
        'target_walk_negative_control_and_basis_release': release_control,
        'scope': 'Exact feasibility and witness interface; the weighted target walk is a search heuristic.',
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
