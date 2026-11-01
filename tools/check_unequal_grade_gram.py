#!/usr/bin/env python3
"""Independent checks of unequal-grade complete folded-Gram channels."""

from collections import defaultdict
from itertools import combinations, product
import argparse
import hashlib
import json
from pathlib import Path
import random

from check_complement_shell_crossings import shell
from check_folded_gram_channels import incoming_generators
from check_homogeneous_central_profiles import python_profile
from check_homogeneous_central_symbols import binomial_two
from check_modular_crossing import column_bits, rank_two
from complement_shell_tensor import folded_profile, pair_occupations
from folded_gram_aggregate import complete_gram_channels as equal_grade_channels
from folded_gram_ranks import binomial_parity, tensor_power_blocks
from unequal_grade_gram import (
    complete_gram_channels, component_block_rank, weighted_exterior_blocks,
)


def verify_weighted_blocks():
    """Enumerate residue occupancies, independently of Gaussian characters."""
    checked = odd = 0
    for modes in range(1, 7):
        expected = defaultdict(lambda: defaultdict(int))
        for occupations in product(range(3), repeat=modes):
            total = sum(occupations)
            singles = occupations.count(1)
            doubles = occupations.count(2)
            energy = sum((2 * index + 1) * value
                         for index, value in enumerate(occupations))
            for highest, twist, multiplicity in tensor_power_blocks(singles):
                expected[total, energy][highest, twist + doubles] += multiplicity
        for total in range(2 * modes + 1):
            for energy in range(2 * modes ** 2 + 1):
                actual = {(m, d): n for m, d, n in weighted_exterior_blocks(
                    modes, total, energy)}
                assert actual == dict(expected[total, energy]), (modes, total, energy)
                checked += 1
                odd += total & 1
    return {"all_weighted_multiplicities": checked, "odd_exterior_grades": odd}


def verify_noncentral_ranks():
    rng = random.Random(202610301723)
    checked = noncentral = 0
    while checked < 240:
        length = rng.randrange(1, 9)
        states = rng.randrange(1, 9)
        delta = states - length
        creators = rng.randrange(max(0, -delta), length + 1)
        double = rng.randrange(min(creators, creators + delta) + 1)
        p = creators - double
        singles = 2 * p + delta
        if singles > 9:
            continue
        shift = rng.randrange(-length, length + 1)
        masks = [sum(1 << bit for bit in subset)
                 for subset in combinations(range(singles), p)]
        order = length + shift - creators - double
        rows = [sum(binomial_parity(length - double - (a & b).bit_count(), order) << j
                    for j, b in enumerate(masks)) for a in masks]
        expected = rank_two(rows)
        actual = sum(multiplicity * component_block_rank(
            length, states, creators, highest, shift)
            for highest, _, multiplicity in tensor_power_blocks(singles))
        assert actual == expected, (length, states, creators, double, shift)
        checked += 1
        noncentral += length != states
    return {"dense_Johnson_comparisons": checked, "noncentral_weight_slices": noncentral}


def merged_polynomial(first, second, length, states, reverse=False):
    delta = states - length
    value = 0
    for (double, single), rows in first.items():
        i = double.bit_count()
        creators = i + (single.bit_count() - delta) // 2
        for a in rows:
            for b in second.get((double, single), ()):
                # incoming_generators stores unchanged-only sets. Their
                # intersection exceeds the creator-only intersection by delta.
                overlap = (a & b).bit_count() - delta
                assert 0 <= overlap <= creators - i
                value ^= binomial_two(length - i - overlap) << (creators + i)
    if reverse:
        value = sum(((value >> degree) & 1) << (2 * length - degree)
                    for degree in range(2 * length + 1))
    return value


def unfold(value, modulus, length, parameter_energy, state_energy):
    result = 0
    offset = state_energy + parameter_energy - length * modulus
    while value:
        bit = value & -value
        degree = (bit.bit_length() - 1) * modulus + offset
        assert degree >= 0
        result |= 1 << degree
        value ^= bit
    return result


def verify_complete_entries():
    cases = []
    for modulus, length, states, j, energy in (
            (16, 2, 4, 16, 32), (16, 4, 2, 32, 16),
            (16, 3, 4, 23, 32), (16, 2, 3, 16, 23),
            (16, 3, 2, 23, 16), (16, 2, 2, 12, 20)):
        parameters, state_labels = shell(modulus, length, j), shell(modulus, states, energy)
        dimension = min(len(parameters), len(state_labels))
        parameters, state_labels = parameters[:dimension], state_labels[:dimension]
        partner_labels = [tuple(sorted(modulus - r for r in label)) for label in parameters]
        parameter_occupations = [pair_occupations(label, modulus) for label in parameters]
        state_occupations = [pair_occupations(label, modulus) for label in state_labels]
        parameter_masks = [sum(1 << ((r - 1) // 2) for r in label) for label in parameters]
        state_masks = [sum(1 << ((r - 1) // 2) for r in label) for label in state_labels]
        first = {(x, y): incoming_generators(parameter_masks[x], state_masks[y], modulus // 2, False)
                 for x, y in product(range(dimension), repeat=2)}
        second = {(x, y): incoming_generators(parameter_masks[x], state_masks[y], modulus // 2, True)
                  for x, y in product(range(dimension), repeat=2)}
        rows = [[0] * dimension ** 2 for _ in range(2 * length + 1)]
        native_checks = 0
        stride = max(1, dimension ** 4 // 32)
        for index, (x, y, u, v) in enumerate(product(range(dimension), repeat=4)):
            expected = folded_profile(parameter_occupations[x], state_occupations[y],
                                      state_occupations[u], parameter_occupations[v])
            row = merged_polynomial(first[x, y], first[v, u], length, states)
            column = merged_polynomial(second[x, u], second[v, y], length, states, True)
            assert row == column == expected, (modulus, length, states, j, energy, x, y, u, v)
            assert expected.bit_count() % 2 == int(x == v and y == u)
            if index % stride == 0:
                assert unfold(expected, modulus, length, j, energy) == python_profile(
                    parameters[x], state_labels[u], partner_labels[v], state_labels[y])
                native_checks += 1
            for degree in range(2 * length + 1):
                if expected >> degree & 1:
                    rows[degree][x * dimension + y] ^= 1 << (v * dimension + u)
        bounds = {v: (r, c) for v, r, c in complete_gram_channels(
            modulus, length, states, j, energy)}
        profile = []
        for degree, matrix in enumerate(rows):
            r, c = rank_two(matrix), rank_two(column_bits(matrix, dimension))
            upper_r, upper_c = bounds.get(degree - length, (0, 0))
            assert r <= upper_r and c <= upper_c
            if r or c:
                profile.append([degree - length, r, c])
        cases.append({
            "mode_cutoff": modulus, "operator_length": length, "state_length": states,
            "operator_energy": j, "state_energy": energy, "selected_labels": dimension,
            "complete_entry_polynomials_in_both_cuts": dimension ** 4,
            "independent_native_Pfaffian_comparisons": native_checks,
            "exact_nonzero_profile_in_modulus_units": profile,
        })
    return cases


def verify_specializations_and_symmetries():
    equal = []
    for args in ((16, 4, 32), (32, 3, 47), (32, 8, 128), (64, 16, 512)):
        modulus, length, energy = args
        assert complete_gram_channels(modulus, length, length, energy, energy) == equal_grade_channels(*args)
        equal.append(list(args))
    exchange = reflection = 0
    for modulus, length, states, j, energy in (
            (16, 2, 3, 12, 19), (16, 3, 4, 19, 28),
            (32, 4, 7, 52, 91), (64, 8, 12, 192, 288)):
        original = complete_gram_channels(modulus, length, states, j, energy)
        swapped = complete_gram_channels(modulus, states, length, energy, j)
        # The normalized polynomial can have different zero margins.
        # If either cut is zero, the complete coefficient tensor is zero.
        nonzero = lambda value: tuple(row for row in value if row[1] and row[2])
        assert nonzero(original) == nonzero(swapped)
        exchange += 1
        reflected = complete_gram_channels(modulus, length, states, j, states * modulus - energy)
        assert nonzero(original) == nonzero(tuple((-v, c, r) for v, r, c in reversed(reflected)))
        reflection += 1
    return {"equal_grade_profile_comparisons": equal,
            "operator_state_exchange_checks": exchange,
            "state_energy_reflection_checks": reflection}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "all-grade complete folded-Gram channels; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "module_sha256": hashlib.sha256(Path(__file__).with_name("unequal_grade_gram.py").read_bytes()).hexdigest(),
        "weighted_multiplicity_checks": verify_weighted_blocks(),
        "noncentral_block_rank_checks": verify_noncentral_ranks(),
        "complete_source_cases": verify_complete_entries(),
        "specializations_and_symmetries": verify_specializations_and_symmetries(),
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
