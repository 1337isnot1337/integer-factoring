#!/usr/bin/env python3
"""Exact independent controls for the moving-support eight-mode core."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_cancellation_module_automata import pmul
from coherent_full_tail_core import binary_rank
from coherent_moving_pair_core import core, images
from coherent_physical_tail_core import concentration_numerator


def dot(left, right):
    value = 0
    for a, b in zip(left, right):
        if a and b:
            value ^= pmul(a, b)
    return value


def shifted(value, exponent):
    if exponent >= 0:
        return value << exponent
    assert not value & ((1 << -exponent) - 1)
    return value >> -exponent


def double(matrix):
    size = len(matrix)
    return [row + [0] * size for row in matrix] + [[0] * size + row for row in matrix]


def concentration_controls():
    """Check the full-mode equations without inverting a rational matrix.

    The source of the full-mode vectors is the literal L/R shear formula.
    C*transported_creator=original_creator and transported_annihilator
    =original_annihilator*C are checked after clearing a local unit.
    The column equations use the separately J-weighted concentration.
    """
    counts = dict(weighted_polynomial_entries=0, constant_unit_ranks=0,
                  full_operator_entries=0, canonical_pairings=0)
    for M in (4, 8, 16, 32, 64, 128):
        numerator, denominator = concentration_numerator(M, 1)
        assert denominator & 1
        weighted = [[shifted(numerator[i][j], 2 * (j - i))
                     for j in range(M)] for i in range(M)]
        counts['weighted_polynomial_entries'] += M * M
        for matrix in (numerator, weighted):
            assert binary_rank([sum((value & 1) << j for j, value in enumerate(row))
                                for row in matrix]) == M
            counts['constant_unit_ranks'] += 1
        C, V = double(numerator), double(weighted)
        target, diagonal = images(M)
        retained = [0, 1, M - 2, M - 1]
        active = retained + [M + j for j in retained]
        spectators = [2 * j + 1 for j in range(2, M - 2)]
        spectators += [2 * M + 2 * j + 1 for j in range(2, M - 2)]
        assert sorted(diagonal + spectators) == list(range(1, 4 * M, 2))
        for local, a in enumerate((0, M - 1)):
            unit = lambda index: [int(j == index) for j in range(2 * M)]
            rx = (unit(a), unit(a))
            ry = (unit(M + a), unit(M + a))
            cv = (unit(a), unit(a))
            cu = (unit(M + a), unit(M + a))
            for b in range(M):
                Fba, Fab = 3 << (2 * (b & ~a)), 3 << (2 * (a & ~b))
                rx[0][M + b] ^= Fba
                ry[1][b] ^= Fab
                cv[1][M + b] ^= Fba
                cu[0][b] ^= Fab
            for group, original in enumerate((rx, ry, cv, cu)):
                change = C if group < 2 else V
                expected = [[0] * (2 * M) for _ in range(2)]
                for part in range(2):
                    for index, value in zip(active, target[group][local][part]):
                        expected[part][index] = value
                assert [dot(row, expected[0]) for row in change] \
                    == [pmul(denominator, value) for value in original[0]]
                assert [dot(original[1], column) for column in zip(*change)] \
                    == [pmul(denominator, value) for value in expected[1]]
                counts['full_operator_entries'] += 4 * M
        for side in (target[:2], target[2:]):
            operators = side[0] + side[1]
            for i, (creator, _) in enumerate(operators):
                for j, (_, annihilator) in enumerate(operators):
                    assert dot(creator, annihilator) == int(i == j)
                    counts['canonical_pairings'] += 1
    return counts


def literal_control():
    from check_block_smith_crossings import smith_profile
    from coherent_hot_compounds import profiles
    from eight_letter_reflection_retraction import ambient_normalized_tensor

    M = 4
    words = [word for word in range(16) if word & 9 != 9]
    reverse = lambda word: sum(((word >> i) & 1) << (3 - i) for i in range(4))
    selected = list(map(reverse, words))
    q = len(words)
    tensor = ambient_normalized_tensor()
    matrix = [[tensor[a | b << 4 | c << 8 | d << 12] for c in selected for d in selected]
              for a in selected for b in selected]
    checked = 0
    for x, y, u, v in product(range(q), repeat=4):
        assert matrix[q * x + y][q * u + v].bit_count() % 2 == int(x == v and y == u)
        checked += 1
    actual = []
    for parity in (0, 1):
        indices = [i for i in range(q * q)
                   if (words[i // q].bit_count() + words[i % q].bit_count()) % 2 == parity]
        block = [[matrix[i][j] for j in indices] for i in indices]
        degree = max(value.bit_length() - 1 for row in block for value in row)
        orders, _ = smith_profile(block, degree)
        actual.extend(orders)
    expected, _, _ = core(M)
    binary_code = [1 << word for word in words]
    compound = profiles(binary_code, binary_code, M)
    assert sorted(actual) == expected == compound['zero_orders_t'] == compound['infinity_orders_t']
    return dict(literal_endpoint_entries=checked, independent_python_smith_orders=len(actual),
                original_compound_orders=2 * q * q, active_core_orders=len(expected),
                physical_forbidden_pair=[0, M - 1])


def replay(path):
    evidence = json.loads(path.read_text())
    assert evidence['complete'] and not evidence['factoring_exponent_improved']
    for name, expected in evidence['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    count = 0
    first_differences = []
    for record in evidence['records']:
        M = record['M']
        orders, blocks, degree = core(M)
        assert orders == record['orders'] and blocks == record['blocks']
        assert degree == record['actual_coefficient_degree']
        assert record['physical_forbidden_pair'] == [0, M - 1]
        assert record['physical_alphabet_dimension'] == 3 * 2 ** (M - 2)
        assert record['physical_hot_heights'] == [0, 1, 2 * M - 1]
        assert record['physical_tail_heights'] == [2 * j + 1 for j in range(1, M - 1)]
        assert record['spectator_orders'] == [2 * j + 1 for j in range(2, M - 2)] \
            + [2 * M + 2 * j + 1 for j in range(2, M - 2)]
        assert record['coefficient_degree_bound'] == 28 * M
        assert record['raw_tensor_degree_w_bound'] == 8 * M * M
        gaussian = Counter({0: 1})
        _, diagonal = images(M)
        for shift in diagonal:
            for exponent, multiplicity in list(gaussian.items()):
                gaussian[exponent + shift] += multiplicity
        actual = Counter(orders)
        difference = {exponent: gaussian[exponent] - actual[exponent]
                      for exponent in gaussian.keys() | actual.keys()
                      if gaussian[exponent] != actual[exponent]}
        first = min(difference)
        assert first == 6 * M - 2 and difference[first] == 1
        first_differences.append(dict(M=M, degree=first, coefficient=1))
        count += len(orders)
    return dict(intervals=len(evidence['records']), recomputed_active_orders=count,
                finite_first_counting_differences=first_differences,
                finite_pressure_values_are_numerical=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='passed', factoring_exponent_improved=False,
                  concentration=concentration_controls(), literal_control=literal_control())
    if args.evidence:
        result['complete_profile_replay'] = replay(args.evidence)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
