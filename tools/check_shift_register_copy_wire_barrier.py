#!/usr/bin/env python3
"""Exact controls for the copied-wire obstruction for all full-cycle registers."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path

from check_controlled_permutation_smith import binary_rank, phi_rows
from check_cyclic_algebra_crossings import jet_orders
from check_shift_register_crossings import (
    full_cycle, joined_cycle, literal_jets, maximal_registers,
)
from check_quadratic_register_crossings import register, truth_tables
from cyclic_algebra_crossing import construct


def colours(bits, first):
    return [int(i in first) for i in range(bits)]


def agreements(pattern, shift):
    return sum(value == pattern[(i+shift) % len(pattern)]
               for i, value in enumerate(pattern))


def colour_controls():
    count = 0
    for bits in range(4, 15, 2):
        for first in combinations(range(bits), bits//2):
            if 0 not in first:
                continue
            pattern = colours(bits, first)
            values = [agreements(pattern, j) for j in range(1, bits)]
            for length in range(1, bits):
                assert max(values[:length])*2*length >= bits*(length-1)
                correlations = [2*c-bits for c in values[:length]]
                expansion = bits*(length+1)+2*sum(
                    (length-j)*value for j, value in enumerate(correlations))
                direct = sum(sum(1 if pattern[(i+t) % bits] else -1
                                 for t in range(length+1))**2 for i in range(bits))
                assert expansion == direct >= 0
                count += 1
    return count


def digit(state, position, radix):
    return state//radix**position % radix


def encode_digits(state, first, second, radix):
    a = sum(digit(state, i, radix)*radix**j for j, i in enumerate(first))
    b = sum(digit(state, i, radix)*radix**j for j, i in enumerate(second))
    return a*radix**len(second)+b


def euler_register(radix, digits):
    """An explicit full shift-register cycle by an Euler tour of its word graph."""
    nodes = radix**(digits-1)
    outgoing, vertices, edges, reverse_tour = [0]*nodes, [0], [], []
    while vertices:
        vertex = vertices[-1]
        if outgoing[vertex] < radix:
            symbol = outgoing[vertex]
            outgoing[vertex] += 1
            edge = radix*vertex+symbol
            vertices.append(edge % nodes)
            edges.append(edge)
        else:
            vertices.pop()
            if edges:
                reverse_tour.append(edges.pop())
    tour = reverse_tour[::-1]
    n = radix**digits
    assert sorted(tour) == list(range(n))
    cycle = [0]*n
    for i, state in enumerate(tour):
        cycle[state] = tour[(i+1) % n]
        assert cycle[state]//radix == state % nodes
    assert full_cycle(cycle)
    return cycle


def control_register(transition, radix, digits, first, coefficient_control=False):
    second = tuple(i for i in range(digits) if i not in first)
    q, n = radix**(digits//2), radix**digits
    encoding = [encode_digits(s, first, second, radix) for s in range(n)]
    cycle = [0]*n
    for state, target in enumerate(transition):
        cycle[encoding[state]] = encoding[target]
    assert full_cycle(cycle)
    assert all(transition[s]//radix == s % (n//radix) for s in range(n))
    pattern, canonical_power, ranks = colours(digits, first), list(range(n)), [n]
    power_cases = triple_cases = coefficient_cases = 0
    for j in range(1, 2*digits):
        canonical_power = [transition[s] for s in canonical_power]
        power = [0]*n
        for s, target in enumerate(canonical_power):
            power[encoding[s]] = encoding[target]
        rows = phi_rows(power, q)
        if radix in (2, 4):
            rank = binary_rank(rows)
        else:
            from check_multipoint_hermite_crossings import rank as prime_rank
            rank = prime_rank([[(row >> i) & 1 for i in range(n)] for row in rows], 3)
        ranks.append(rank)
        if j < digits:
            copied = [i for i in range(j, digits) if pattern[i] == pattern[i-j]]
            assert len(copied) >= agreements(pattern, j)-j
            blocks = set()
            for s, target in enumerate(canonical_power):
                row_label = column_label = 0
                for k, i in enumerate(copied):
                    a, b = digit(target, i, radix), digit(s, i-j, radix)
                    assert a == b
                    row_label += (a if pattern[i] else b)*radix**k
                    column_label += (b if pattern[i] else a)*radix**k
                assert row_label == column_label
                blocks.add(row_label)
            assert len(blocks) == radix**len(copied)
            assert rank >= len(blocks)
            power_cases += 1
    for j in range(1, digits):
        for k in range(1, digits):
            assert ranks[j]*ranks[k]*ranks[j+k] >= n
            triple_cases += 1
    if coefficient_control:
        _, _, infinity = construct(q, 2, cycle)
        orders = jet_orders(infinity, 16)
        for j in range(1, digits):
            early = sum(order < j for order in orders)
            assert ranks[j] <= j*(j+1)//2 * early
            coefficient_cases += 1
    return power_cases, triple_cases, coefficient_cases


def register_controls():
    power_cases = triple_cases = coefficient_cases = 0
    categories = dict(linear_binary=0, nonlinear_binary=0, ternary=0, quaternary=0)
    def add(transition, radix, digits, first, category, coefficients=False):
        nonlocal power_cases, triple_cases, coefficient_cases
        a, b, c = control_register(transition, radix, digits, first, coefficients)
        power_cases += a
        triple_cases += b
        coefficient_cases += c
        categories[category] += 1
    for bits in (4, 6):
        for index, (_, permutation) in enumerate(maximal_registers(bits)):
            joined = list(permutation)
            joined[0], joined[1 << (bits-1)] = joined[1 << (bits-1)], joined[0]
            for first in combinations(range(bits), bits//2):
                if 0 in first:
                    add(joined, 2, bits, first, 'linear_binary',
                        bits == 4 or (index == 0 and first == (0, 1, 3)))
    selected = 0
    for code, table in enumerate(truth_tables()):
        if not code >> 6:
            continue
        transition = register(table)
        if not full_cycle(transition):
            continue
        for first in combinations(range(6), 3):
            if 0 in first:
                add(transition, 2, 6, first, 'nonlinear_binary')
        selected += 1
        if selected == 8:
            break
    for radix, category in ((3, 'ternary'), (4, 'quaternary')):
        transition = euler_register(radix, 4)
        for first in combinations(range(4), 2):
            if 0 in first:
                add(transition, radix, 4, first, category)
    return dict(copied_wire_power_cases=power_cases, orthogonal_triples=triple_cases,
                coefficient_channel_cases=coefficient_cases, controls_by_class=categories)


def counterexample():
    bits, taps, q = 12, 2100, 64
    n = 1 << bits
    first = (0, 1, 3, 6, 7, 10)
    second = tuple(i for i in range(bits) if i not in first)
    permutation = [((state << 1) & (n-1)) | ((state & taps).bit_count() & 1)
                   for state in range(n)]
    cycle = joined_cycle(permutation, first, second)
    assert full_cycle(cycle)
    r, k, _, _ = literal_jets(cycle, q)
    r2 = binary_rank(phi_rows([cycle[value] for value in cycle], q))
    assert (r, r2, k) == (17, 18, 17)
    assert r*r2 < n <= r*r*r2
    return dict(bits=bits, taps=taps, first=list(first), full_cycle=n,
                rank_one_power=r, rank_two_power=r2, order_one_pivots=k,
                false_two_rank_product=r*r2, valid_triple_product=r*r*r2)


def integer_pressure_controls():
    # The per-channel n/j lower bound, without logs, at rational test scales.
    from fractions import Fraction
    scalar_cases = 0
    for n in (4, 9, 16, 64, 81, 256):
        for j in range(1, min(n, 17)):
            for scale in map(Fraction, ('1', '65/64', '9/8', '3/2', '2', '8')):
                value = scale**(n-j)+(n-1)*scale**(-j)
                assert j*value >= n
                scalar_cases += 1
    for bits in range(576, 1296, 2):
        length = isqrt(bits//2)
        exponent = bits//2-(bits+2*length-1)//(2*length)-length
        power = 12*exponent-5*bits
        assert power >= 0 and (1 << power) >= length**36
    large = []
    for bits in (1296, 2048, 4096, 16384):
        length = isqrt(bits)
        exponent = bits//2-(bits+2*length-1)//(2*length)-length
        # [d^(m+E)/J^3]^12 >= [q^(17/6)]^12; d=2 is the worst case.
        assert (1 << (12*exponent-5*bits)) >= length**36
        large.append(dict(bits=bits, J=length, E=exponent,
                          rational_rho_lower_bound='17/24', integer_comparison=True))
    assert 2*17**2 > 24**2
    return dict(scalar_cases=scalar_cases,
                complete_finite_cutoff=dict(first_even_bits=576, last_even_bits=1294, comparisons=360),
                large_size_checks=large)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(colour_window_cases=colour_controls(), registers=register_controls(),
                  counterexample=counterexample(), pressure=integer_pressure_controls(),
                  asymptotic_bare_rho_lower_bound='3/4',
                  explicit_cutoff=dict(even_register_digits_at_least=576, rho_lower_bound='17/24'),
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Every prime-power digit full-cycle shift register, arbitrary nonlinear feedback, coordinate partitions; bare two-place Smith criterion.',
                  factoring_exponent='1/2; unchanged')
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(colour_cases=result['colour_window_cases'], registers=result['registers'],
                          asymptotic_bare_rho='3/4', factoring_exponent=result['factoring_exponent'])))


if __name__ == '__main__':
    main()
