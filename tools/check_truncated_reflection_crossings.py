#!/usr/bin/env python3
"""Independent complete coefficient checks of low-jet degree truncations."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, digest
from truncate_reflection_crossings import truncate_even_crossing, truncate_double_endpoint


def divide_one_plus_w_squared(value):
    quotient = 0
    while value.bit_length() > 2:
        power = value.bit_length() - 3
        quotient ^= 1 << power
        value ^= 5 << power
    assert value == 0
    return quotient


def minimum_trace_path(q, power):
    output = []
    for index in range(q ** 4):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        crossing = int(x == v and y == u)
        trace_carrier = int(x == v == 0 and y == u) ^ int(x == v and y == u == 0)
        output.append((crossing ^ trace_carrier ^ (trace_carrier << power), 0))
    return output


def paired_trace(seed, power):
    q = seed['q_x']
    expected = (1 | (1 << power), 0)
    hashes = []
    for first, second in ((0, 3), (1, 2)):
        other = [leg for leg in range(4) if leg not in (first, second)]
        matrix = []
        for i in range(q):
            for j in range(q):
                letters, total = [0] * 4, (0, 0)
                letters[other[0]], letters[other[1]] = i, j
                for letter in range(q):
                    letters[first] = letters[second] = letter
                    index = sum(value * q ** leg for leg, value in enumerate(letters))
                    total = add(total, seed['values'][index])
                assert total == (expected if i == j else (0, 0))
                matrix.append(total)
        hashes.append(sha256(json.dumps(matrix, separators=(',', ':')).encode()).hexdigest())
    return {'trace_degree': power, 'polynomial_entries': 2 * q * q, 'trace_sha256': hashes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records, dictionary_entries, quotient_entries = [], 0, 0
    for interval, bounds in ((4, tuple(range(2, 25, 2))), (8, (16, 32, 48, 64, 68, 72, 80, 88, 96))):
        base = tensor(interval)
        q = base['q_x']
        for bound in bounds:
            trace_degree = min(bound, interval ** 2)
            minimum = minimum_trace_path(q, trace_degree)
            changed = truncate_even_crossing(base, bound)
            general = truncate_double_endpoint(base, bound)
            assert changed['values'] == general['values']
            assert changed['values'][0] == (1 << min(bound, interval ** 2 // 2), 0)
            assert any((a | b) & 1 for a, b in changed['values'])
            for index, (pair, source) in enumerate(zip(changed['values'], base['values'])):
                x, y, u, v = ((index // q ** leg) % q for leg in range(4))
                for plane in range(2):
                    original = {k for k in range(source[plane].bit_length()) if source[plane] >> k & 1}
                    kept = {k for k in original if k < bound}
                    if (len(kept) % 2) ^ int(plane == 0 and x == v and y == u):
                        kept.add(bound)
                    assert pair[plane] == sum(1 << k for k in kept)
                dictionary_entries += 1
            # The quotient construction is independent of the parity shortcut.
            mask = (1 << (bound - 1)) - 1
            reference = []
            for source, carrier in zip(base['values'], minimum):
                delta = add(source, carrier)
                quotient = tuple(divide_one_plus_w_squared(value) & mask for value in delta)
                reference.append(tuple(part ^ value ^ (value << 2)
                                       for part, value in zip(carrier, quotient)))
                quotient_entries += 1
            assert reference == changed['values']
            trace = paired_trace(changed, trace_degree)
            records.append({'interval': interval, 'alphabet': q, 'requested_degree': bound,
                            'actual_degree': changed['degree'],
                            'unchanged_coefficients_through': bound - 1,
                            'tensor_sha256': digest(changed['values']),
                            'endpoint_entries': q ** 4, 'paired_trace': trace})
    names = (Path(__file__).name, 'truncate_reflection_crossings.py',
             'reflection_dual_rail_gaussian.py', 'reflection_laurent_physical_frames.py')
    report = {'status': 'passed', 'scope': 'Exact degree reduction with unchanged low jets and paired traces; no pressure or factoring certificate.',
              'seeds_checked': len(records), 'independent_dictionary_entries': dictionary_entries,
              'independent_quotient_entries': quotient_entries, 'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('records', 'source_sha256')}))


if __name__ == '__main__':
    main()
