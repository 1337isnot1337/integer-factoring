#!/usr/bin/env python3
"""Exact low-jet identities and the arbitrary-rank pole-two criterion.

The main finite assertion checks every entry in degrees zero through seven
of the four- and sixteen-letter reflection seeds.  The ordinary proof then
applies to every square-zero generator, without enumerating those matrices.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from square_zero_coupled_reflection_frames import matrix_terms, from_terms


def multiply(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        a <<= 1
        b >>= 1
    return value ^ 7 if value & 4 else value


def product(first, second):
    output = [[0] * len(second[0]) for _ in first]
    for i in range(len(first)):
        for j in range(len(second[0])):
            for k in range(len(second)):
                output[i][j] ^= multiply(first[i][k], second[k][j])
    return output


def coefficient(value, exponent):
    if exponent < 0:
        return 0
    a, b = value
    return ((a >> exponent) & 1) | (((b >> exponent) & 1) << 1)


def low_jets(seed):
    q = seed['q_x']

    def entry(x, y, u, v, power):
        return coefficient(seed['values'][x + q * (y + q * (u + q * v))], power)

    h = [[int(i == j and i.bit_count() % 2) for j in range(q)] for i in range(q)]
    j = [[entry(x, y, 0, 0, 4) for y in range(q)] for x in range(q)]
    k = [[entry(x, y, 0, 0, 6) for y in range(q)] for x in range(q)]
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    dx, du = int(x == y), int(u == v)
                    hxy, huv = h[x][y], h[u][v]
                    jxy, juv = j[x][y], j[u][v]
                    expect = {
                        0: dx * du,
                        2: (hxy if du else 0) ^ (huv if dx else 0),
                        4: ((jxy if du else 0) ^ (juv if dx else 0)
                            ^ multiply(hxy, huv) ^ (huv if dx else 0)),
                        6: ((k[x][y] if du else 0) ^ (k[u][v] if dx else 0)
                            ^ multiply(hxy, juv) ^ multiply(jxy, huv)
                            ^ multiply(hxy, huv) ^ (huv if dx else 0)),
                    }
                    for power in range(8):
                        assert entry(x, y, u, v, power) == expect.get(power, 0)
    return h, j, k


def paired_matching(q, parity_preserving):
    if parity_preserving:
        groups = [[i for i in range(q) if i.bit_count() % 2 == parity] for parity in (0, 1)]
        pairs = [group[i:i + 2] for group in groups for i in range(0, len(group), 2)]
    else:
        even = [i for i in range(q) if i.bit_count() % 2 == 0]
        odd = [i for i in range(q) if i.bit_count() % 2]
        pairs = list(zip(even, odd))
    generator = [[0] * q for _ in range(q)]
    for i, j in pairs:
        for row in (i, j):
            for col in (i, j):
                generator[row][col] = 1
    return generator


def controls(seed, h):
    q = seed['q_x']
    generators = [('parity-preserving perfect matching', paired_matching(q, True)),
                  ('opposite-parity perfect matching', paired_matching(q, False))]
    # Image span(e0,e1) is H-invariant, although this G does not commute H.
    noncommuting = [[0] * q for _ in range(q)]
    noncommuting[0][2] = 2
    noncommuting[1][3] = 3
    generators.append(('noncommuting rank-two map with parity-invariant image', noncommuting))
    results = []
    for name, generator in generators:
        assert product(generator, generator) == [[0] * q for _ in range(q)]
        obstruction = product(product(generator, h), generator)
        permitted = not any(value for row in obstruction for value in row)
        terms = matrix_terms(seed, generator)
        for index in range(q ** 4):
            actual = 0
            for order in range(5):
                actual ^= coefficient(terms[order].get(index, (0, 0)), -4 + 2 * order)
            x, y, u, v = (index // q ** leg % q for leg in range(4))
            assert actual == multiply(obstruction[x][y], obstruction[u][v])
        try:
            changed = from_terms(seed, terms, pole=2)
        except ArithmeticError:
            assert not permitted
            changed = None
        else:
            assert permitted
            assert changed['degree'] == seed['degree']
            assert changed['monomial_content'] == 0
        results.append({'name': name, 'generator': generator, 'GHG': obstruction,
                        'pole_two_polynomial': permitted,
                        'unchanged_degree': None if changed is None else changed['degree']})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    for interval in (4, 8):
        source = tensor(interval)
        bell, factors = bell_normalize(source)
        h, j, k = low_jets(bell)
        q = source['q_x']
        cases.append({'interval': interval, 'alphabet': q, 'degree': source['degree'],
                      'H': h, 'J': j, 'K': k, 'bell_factors': factors,
                      'low_coefficient_entries_checked': 8 * q ** 4,
                      'obstruction_entries_checked': 3 * q ** 4,
                      'selected_generator_controls': controls(bell, h),
                      'source_seed_sha256': sha256(json.dumps(source['values'], separators=(',', ':')).encode()).hexdigest()})
    names = (Path(__file__).name, 'square_zero_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py')
    report = {'status': 'passed', 'scope': 'all low-jet entries for two fixed seeds, yielding an all-rank algebraic criterion; no pressure claim',
              'cases': cases,
              'low_coefficient_entries_checked': sum(case['low_coefficient_entries_checked'] for case in cases),
              'obstruction_entries_checked': sum(case['obstruction_entries_checked'] for case in cases),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
