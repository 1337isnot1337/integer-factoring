#!/usr/bin/env python3
"""Exact finite checks for the weighted divisor-prefix constructor.

This is not the large factoring algorithm or a modular clique evaluator.
It checks the pinned multiplication gates on small forced settings, and
independently compares cyclic-polynomial grouped tables with state sums.
"""

from __future__ import annotations

from itertools import combinations, product
import json

try:
    from half_clique_codec import EqualityNetwork, MessageEdge
except ModuleNotFoundError:
    from tools.half_clique_codec import EqualityNetwork, MessageEdge


def ring_add(a, b, characteristic):
    return tuple((x + y) % characteristic for x, y in zip(a, b))


def ring_mul(a, b, characteristic):
    size = len(a)
    assert len(b) == size
    result = [0] * size
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[(i + j) % size] += x * y
    return tuple(x % characteristic for x in result)


def monomial(exponent, size):
    return tuple(int(i == exponent % size) for i in range(size))


def weighted_state_sum(network, exponents, characteristic, size):
    result = [0] * size
    for states in product(*(range(len(row)) for row in network.weights)):
        if not all(network.weights[v][states[v]] for v in range(len(states))):
            continue
        if not all(edge.first_messages[states[edge.first]] ==
                   edge.second_messages[states[edge.second]] for edge in network.edges):
            continue
        exponent = sum(exponents[v][states[v]] for v in range(len(states)))
        result[exponent % size] += 1
    return tuple(x % characteristic for x in result)


def weighted_group_tables(network, exponents, groups, k, width, characteristic, size):
    edge_order = sorted(network.edges, key=lambda edge: (edge.first, edge.second))
    cross = {pair: tuple(edge for edge in edge_order
                        if tuple(sorted((groups[edge.first], groups[edge.second]))) == pair)
             for pair in combinations(range(k), 2)}
    assert all(sum(edge.width for edge in edges) <= width for edges in cross.values())
    result = []
    for vertex in range(k):
        members = [v for v, group in enumerate(groups) if group == vertex]
        internal = [edge for edge in edge_order
                    if groups[edge.first] == groups[edge.second] == vertex]
        row = [[0] * size for _ in range(1 << ((k - 1) * width))]
        for choices in product(*(range(len(network.weights[v])) for v in members)):
            state = dict(zip(members, choices))
            if not all(network.weights[v][state[v]] for v in members):
                continue
            if any(edge.first_messages[state[edge.first]] !=
                   edge.second_messages[state[edge.second]] for edge in internal):
                continue
            index = block_index = 0
            for other in range(k):
                if other == vertex:
                    continue
                label = offset = 0
                for edge in cross[min(vertex, other), max(vertex, other)]:
                    if groups[edge.first] == vertex:
                        message = edge.first_messages[state[edge.first]]
                    else:
                        message = edge.second_messages[state[edge.second]]
                    label |= message << offset
                    offset += edge.width
                index |= label << (block_index * width)
                block_index += 1
            exponent = sum(exponents[v][state[v]] for v in members) % size
            row[index][exponent] = (row[index][exponent] + 1) % characteristic
        result.append(tuple(tuple(entry) for entry in row))
    return tuple(result)


def weighted_clique_value(rows, width, characteristic, size):
    k = len(rows)
    pairs = tuple(combinations(range(k), 2))
    total = (0,) * size
    for values in product(range(1 << width), repeat=len(pairs)):
        edges = dict(zip(pairs, values))
        term = monomial(0, size)
        for vertex in range(k):
            labels = [edges[min(vertex, other), max(vertex, other)]
                      for other in range(k) if other != vertex]
            index = sum(label << (i * width) for i, label in enumerate(labels))
            term = ring_mul(term, rows[vertex][index], characteristic)
        total = ring_add(total, term, characteristic)
    return total


def params(s):
    return 1 << (s + 1), 1 << s, 1 << (2 * s + 1)


def input_digits(value, s):
    a, h, digit_bits = params(s)
    mask = (1 << digit_bits) - 1
    return tuple(tuple((value >> (digit_bits * (r + a * i))) & mask
                       for i in range(h)) for r in range(a))


def source_coefficients(left, right, s):
    """The three pinned faOut/midOut/ibOut formulas, modulo their actual M."""
    a, h, digit_bits = params(s)
    modulus = (1 << (4 * digit_bits)) + 1
    omega = 16
    omega_inverse = pow(omega, a * a - 1, modulus)
    scale_inverse = pow(a * a, -1, modulus)
    p_digits, q_digits = input_digits(left, s), input_digits(right, s)

    def forward(digits, r, k):
        return pow(omega, r * k, modulus) * sum(
            digits[r][i] * pow(omega, a * i * k, modulus) for i in range(h)) % modulus

    middle = []
    for k in range(a):
        yp = [forward(p_digits, r, k) for r in range(a)]
        yq = [forward(q_digits, r, k) for r in range(a)]
        products = []
        for k2 in range(a):
            hp = sum(yp[r] * pow(omega, a * r * k2, modulus) for r in range(a))
            hq = sum(yq[r] * pow(omega, a * r * k2, modulus) for r in range(a))
            products.append(hp * hq % modulus)
        middle.append([pow(omega_inverse, k * j, modulus) * sum(
            products[k2] * pow(omega_inverse, a * k2 * j, modulus) for k2 in range(a)
        ) % modulus for j in range(a)])
    return tuple(scale_inverse * sum(
        middle[k][j % a] * pow(omega_inverse, a * k * (j // a), modulus)
        for k in range(a)) % modulus for j in range(a * a))


def carry_accepts(coefficients, value, s):
    _, _, digit_bits = params(s)
    base = 1 << digit_bits
    carry = 0
    for j, coefficient in enumerate(coefficients):
        combined = carry + coefficient
        if combined % base != (value >> (digit_bits * j)) % base:
            return False
        carry = combined // base
    return carry == 0


def source_prefix_accepts(digits, k, x, s):
    a, h, digit_bits = params(s)
    if x >= 1 << (digit_bits * digit_bits - k):
        return False
    for r in range(a):
        for i in range(h):
            for t in range(digit_bits):
                position = digit_bits * (r + a * i) + t
                if position >= k and ((digits[r][i] >> t) & 1) != ((x >> (position - k)) & 1):
                    return False
    return True


def local_weight_exponent(digits, r, s, size):
    a, h, digit_bits = params(s)
    return sum(digits[r][i] * pow(2, digit_bits * (r + a * i), size)
               for i in range(h)) % size


def verify():
    source_cases = carry_cases = prefix_cases = grouped_cases = digit_weight_cases = 0
    accepted = {value: [] for value in range(1, 16)}
    for left, right in product(range(16), repeat=2):
        coefficients = source_coefficients(left, right, 0)
        pd, qd = [left & 3, left >> 2], [right & 3, right >> 2]
        expected = tuple(sum(pd[i] * qd[j - i] for i in range(2) if 0 <= j - i < 2)
                         for j in range(4))
        assert coefficients == expected
        source_cases += 1
        for value in range(1, 16):
            good = carry_accepts(coefficients, value, 0)
            assert good == (left * right == value)
            carry_cases += 1
            if good:
                accepted[value].append(input_digits(left, 0))

    # The next scale uses the composite modulus 2^32+1 and two digits per input gate.
    larger_inputs = (0, 1, 255, 256, 257, 65535, (1 << 32) - 1, (1 << 64) - 1)
    for left, right in product(larger_inputs, repeat=2):
        coefficients = source_coefficients(left, right, 1)
        pd = [(left >> (8 * i)) & 255 for i in range(8)]
        qd = [(right >> (8 * i)) & 255 for i in range(8)]
        expected = tuple(sum(pd[i] * qd[j - i] for i in range(8) if 0 <= j - i < 8)
                         for j in range(16))
        assert coefficients == expected
        source_cases += 1
    for left, size in product(larger_inputs, (1, 3, 5, 7, 11)):
        digits = input_digits(left, 1)
        assert sum(local_weight_exponent(digits, r, 1, size) for r in range(4)) % size == left % size
        digit_weight_cases += 1

    for value in range(1, 16):
        n = value.bit_length()
        for k in range(n + 1):
            for x in range(1 << (n - k)):
                for characteristic, size in product((2, 3, 5), (1, 2, 3, 5, 7)):
                    actual = (0,) * size
                    for digits in accepted[value]:
                        if source_prefix_accepts(digits, k, x, 0):
                            term = monomial(0, size)
                            for r in range(2):
                                term = ring_mul(term, monomial(
                                    local_weight_exponent(digits, r, 0, size), size), characteristic)
                            actual = ring_add(actual, term, characteristic)
                    expected = [0] * size
                    for divisor in range(1, value + 1):
                        if value % divisor == 0 and divisor >> k == x:
                            expected[divisor % size] += 1
                    assert actual == tuple(c % characteristic for c in expected)
                    prefix_cases += 1

    for bits in product((0, 1), repeat=8):
        network = EqualityNetwork(((1, 1), (1, 1), (1, 1)), (
            MessageEdge(0, 1, 1, bits[:2], bits[2:4]),
            MessageEdge(1, 2, 1, bits[4:6], bits[6:]),
        ))
        exponents = ((0, 1), (1, 2), (2, 0))
        for characteristic, size in product((2, 3), (3, 5)):
            expected = weighted_state_sum(network, exponents, characteristic, size)
            for groups in ((0, 1, 0), (1, 0, 1)):
                rows = weighted_group_tables(network, exponents, groups, 2, 3, characteristic, size)
                assert weighted_clique_value(rows, 3, characteristic, size) == expected
                assert all(entry == (0,) * size for row in rows for entry in row[4:])
                grouped_cases += 1

    network = EqualityNetwork(((1, 1), (1, 0), (1, 1)), (
        MessageEdge(0, 1, 1, (0, 1), (0, 1)),
        MessageEdge(0, 2, 1, (1, 0), (0, 1)),
        MessageEdge(1, 2, 1, (1, 0), (0, 1)),
    ))
    exponents = ((0, 1), (1, 2), (2, 0))
    rows = weighted_group_tables(network, exponents, (1, 0, 1), 3, 2, 2, 5)
    assert rows[2][0] == monomial(0, 5)
    assert all(entry == (0,) * 5 for entry in rows[2][1:])
    assert weighted_clique_value(rows, 2, 2, 5) == weighted_state_sum(network, exponents, 2, 5)
    grouped_cases += 1

    # A nonempty set can vanish after both ordinary counting and cyclic reduction.
    assert ring_add(monomial(3, 1), monomial(9, 1), 2) == (0,)
    assert ring_add(monomial(3, 3), monomial(9, 3), 2) == (0, 0, 0)
    assert any(ring_add(monomial(3, 5), monomial(9, 5), 2))

    return {
        "scope": "finite weighted constructor checks; no modular evaluator or new factoring exponent",
        "source_three_stage_convolutions": source_cases,
        "source_carry_acceptance_cases": carry_cases,
        "distributed_digit_weight_cases": digit_weight_cases,
        "weighted_divisor_prefix_cases": prefix_cases,
        "weighted_grouped_contractions": grouped_cases,
        "coefficient_characteristics": [2, 3, 5],
        "cyclic_dimensions": [1, 2, 3, 5, 7],
        "negative_control": "nontrivial divisors of 9 vanish mod 2 after exponents reduce mod 1 or 3, but not 5",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
