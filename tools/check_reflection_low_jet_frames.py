#!/usr/bin/env python3
"""Check q-by-q reflection jets against independent complete crossings."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, coupled_admissible_poles
from reflection_low_jet_frames import reflection_jets, matrix_coefficient, pole_conditions
from square_zero_coupled_reflection_frames import matrix_terms, matching_generator, product


def polynomial_product(a, b, cutoff):
    result = 0
    while b:
        bit = b & -b
        b ^= bit
        result ^= a << (bit.bit_length() - 1)
    return result & ((1 << cutoff) - 1)


def product_entry(first, second, i, j, cutoff):
    result = 0
    for k in range(len(first)):
        result ^= polynomial_product(first[i][k], second[k][j], cutoff)
    return result


def verify_normalization(jets):
    q, cutoff = jets['alphabet'], jets['cutoff']
    assert all(not value & sum(1 << exponent for exponent in range(1, cutoff, 2))
               for row in jets['B'] for value in row)
    parity = [[int(i == j and i.bit_count() % 2) for j in range(q)] for i in range(q)]
    h = product(product(jets['Q0_inverse'], parity), jets['Q0'])
    b2 = matrix_coefficient(jets['B'], 2)
    scalar = b2[0][0] ^ h[0][0]
    assert b2 == [[h[i][j] ^ (scalar if i == j else 0) for j in range(q)] for i in range(q)]
    return h, scalar


def selected_generators(q):
    entries = [(0, 1), (0, q // 2), (1, 2), (q - 1, 0),
               (0, q - 1), (q // 2 - 1, q // 2)]
    if q == 16:
        entries += [(6, 9), (9, 6)]
    cases = []
    for index, (i, j) in enumerate(dict.fromkeys(entries)):
        generator = [[0] * q for _ in range(q)]
        generator[i][j] = 1 + index % 3
        cases.append((f'E_{i},{j} with nonzero F4 scalar', generator))
    for preserve in (False, True):
        permutation = list(range(q))
        groups = ([list(range(q))] if not preserve
                  else [[i for i in range(q) if i.bit_count() % 2 == parity] for parity in (0, 1)])
        for group in groups:
            for i in range(0, len(group), 2):
                a, b = group[i:i + 2]
                permutation[a], permutation[b] = b, a
        cases.append(('parity-preserving matching' if preserve else 'adjacent-letter matching',
                      matching_generator(permutation)))
    noncommuting = [[0] * q for _ in range(q)]
    noncommuting[0][2], noncommuting[1][3] = 2, 3
    cases.append(('noncommuting rank-two map', noncommuting))
    return cases


def complete_control(interval):
    source = tensor(interval)
    bell, factors = bell_normalize(source)
    jets = reflection_jets(interval)
    q, cutoff = jets['alphabet'], jets['cutoff']
    assert jets['Q0'] == factors['first_factor']
    assert jets['Q0_inverse'] == factors['inverse_factor']
    mask = (1 << cutoff) - 1
    for index, (value, high) in enumerate(source['values']):
        assert high == 0
        x, y, u, v = (index // q ** leg % q for leg in range(4))
        assert value & mask == polynomial_product(jets['Q_H'][x][y], jets['R_G'][u][v], cutoff)
        assert bell['values'][index][0] & mask == polynomial_product(jets['B'][x][y], jets['inverse_B'][u][v], cutoff)
    for i in range(q):
        for j in range(q):
            assert product_entry(jets['Q_H'], jets['R_G'], i, j, cutoff) == int(i == j)
            assert product_entry(jets['B'], jets['inverse_B'], i, j, cutoff) == int(i == j)
    h, scalar = verify_normalization(jets)
    generators = []
    for name, generator in selected_generators(q):
        conditions = pole_conditions(jets, generator)
        predicted = [case['pole'] for case in conditions if case['polynomial']]
        terms = matrix_terms(bell, generator)
        literal = coupled_admissible_poles(terms, interval // 2)
        assert predicted == literal
        generators.append({'name': name, 'generator': generator, 'allowed_poles': predicted,
                           'rejected_poles': [{'pole': case['pole'],
                                               'negative_exponents': sorted(case['negative_coefficients'])}
                                              for case in conditions if not case['polynomial']]})
    return {'interval': interval, 'alphabet': q, 'cutoff': cutoff,
            'Q0': jets['Q0'], 'Q0_inverse': jets['Q0_inverse'],
            'B': jets['B'], 'inverse_B': jets['inverse_B'], 'H': h,
            'B2_minus_H_scalar': scalar,
            'original_and_bell_tensor_entries': 2 * q ** 4,
            'inverse_product_entries': 2 * q ** 2,
            'generator_controls': generators,
            'literal_pole_comparisons': len(generators) * (interval // 2)}


def larger_control():
    """A q=256 matrix-jet control without constructing a q**4 tensor."""
    began = time.monotonic()
    jets = reflection_jets(16)
    q, cutoff = jets['alphabet'], jets['cutoff']
    h, scalar = verify_normalization(jets)
    pairs = [(0, 0), (q - 1, q - 1), (0, q - 1), (q - 1, 0)]
    pairs += [(i, (73 * i + 19) % q) for i in range(0, q, 17)]
    for i, j in pairs:
        assert product_entry(jets['Q_H'], jets['R_G'], i, j, cutoff) == int(i == j)
        assert product_entry(jets['B'], jets['inverse_B'], i, j, cutoff) == int(i == j)
    return {'interval': 16, 'alphabet': q, 'cutoff': cutoff,
            'selected_minor_matrix_entries': q ** 2,
            'sampled_inverse_product_entries': 2 * len(pairs),
            'B2_minus_H_scalar': scalar, 'H_is_letter_parity':
                h == [[int(i == j and i.bit_count() % 2) for j in range(q)] for i in range(q)],
            'B_sha256': sha256(json.dumps(jets['B'], separators=(',', ':')).encode()).hexdigest(),
            'inverse_B_sha256': sha256(json.dumps(jets['inverse_B'], separators=(',', ':')).encode()).hexdigest(),
            'elapsed_seconds': time.monotonic() - began,
            'scope': 'physical matrix jets only; no full q**4 crossing or pressure is constructed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-interval16', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [complete_control(interval) for interval in (4, 8)]
    larger = larger_control() if args.include_interval16 else None
    names = (Path(__file__).name, 'reflection_low_jet_frames.py',
             'reflection_dual_rail_gaussian.py', 'reflection_laurent_physical_frames.py',
             'square_zero_coupled_reflection_frames.py')
    report = {'status': 'passed', 'scope': 'uniform criterion proof with independent complete finite tensor controls and optional larger physical-jet control',
              'cases': cases, 'larger_control': larger,
              'original_and_bell_tensor_entries': sum(case['original_and_bell_tensor_entries'] for case in cases),
              'inverse_product_entries': sum(case['inverse_product_entries'] for case in cases),
              'literal_pole_comparisons': sum(case['literal_pole_comparisons'] for case in cases),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
