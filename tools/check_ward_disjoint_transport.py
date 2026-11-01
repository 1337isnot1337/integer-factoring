#!/usr/bin/env python3
"""Exact form-ring versus block-flag transport for disjoint Ward supports."""

import argparse
import hashlib
import json
from pathlib import Path
import random

from check_ward_form_operator_transport import insertion_transitions
from check_ward_stable_operators import bits, closed_operator_interface


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = closed_operator_interface()
    inserted = insertion_transitions()
    ordinary = base['transport_transitions']
    incidence = inserted['physical_insertion_transitions_incidence']
    creation = inserted['physical_insertion_transitions_creation']
    basis = base['homogeneous_basis_matrices']
    rng = random.Random(617401)
    records = []
    g_grades = ((1, 0), (0, 1), (2, 1), (1, 2))
    f_grades = ((1, 0), (0, 1), (0, 0), (1, 1))
    pair_grades = [(a + c, b + d) for a, b in f_grades for c, d in f_grades]

    def matrix_at(table, position):
        return [[(value >> position) & 1 for value in row] for row in table]

    def act(matrix, vector):
        result = [0] * 8
        for a in range(8):
            for b in range(8):
                if matrix[a][b]:
                    result[a] ^= vector[b]
        return result

    def inserted_action(r, o, vector):
        first, second = act(r, vector), act(o, vector)
        return [a ^ (b << 1) for a, b in zip(first, second)]

    for sizes in ((1, 1), (2, 2), (2, 2, 2), (3, 2, 4), (1, 3, 2, 2)):
        supports = [a for a, size in enumerate(sizes) for _ in range(size)]
        count = nonzero = grade_checks = 0
        attempts = 0
        while count < 160:
            attempts += 1
            assert attempts < 100000
            selected = {sum(sizes[:a]) + rng.randrange(size) for a, size in enumerate(sizes)}
            positions = [0] * len(supports)
            channel = 0  # The final ancillary observable is the identity basis element.
            for site in reversed(range(len(supports))):
                if site in selected:
                    options = [(b, table[channel][b]) for table in (incidence, creation)
                               for b in range(8) if table[channel][b]]
                else:
                    options = [(b, ordinary[channel][b]) for b in range(8) if ordinary[channel][b]]
                if not options:
                    break
                channel, value = rng.choice(options)
                positions[site] = rng.choice(list(bits(value)))
            else:
                boundary = rng.choice(list(bits(basis[channel])))
                initial = [(value >> boundary) & 1 for value in basis]
                full = [initial] + [[0] * 8 for _ in range((1 << len(sizes)) - 1)]
                for site, position in enumerate(positions):
                    b = matrix_at(ordinary, position)
                    r = matrix_at(incidence, position)
                    o = matrix_at(creation, position)
                    following = [act(b, vector) for vector in full]
                    bit = 1 << supports[site]
                    for mask, vector in enumerate(full):
                        if not mask & bit:
                            inserted_value = inserted_action(r, o, vector)
                            following[mask | bit] = [a ^ c for a, c in zip(following[mask | bit], inserted_value)]
                    full = following
                expected = full[-1][0]
                current, offset = initial, 0
                for size in sizes:
                    zero, one = current, [0] * 8
                    for position in positions[offset:offset + size]:
                        b = matrix_at(ordinary, position)
                        r = matrix_at(incidence, position)
                        o = matrix_at(creation, position)
                        no_insertion = act(b, zero)
                        previous = act(b, one)
                        insertion = inserted_action(r, o, zero)
                        zero, one = no_insertion, [a ^ c for a, c in zip(previous, insertion)]
                    current = one
                    offset += size
                actual = current[0]
                assert actual == expected
                scalar, offset = initial, 0
                for size in sizes:
                    zero, one = scalar, [0] * 8
                    for position in positions[offset:offset + size]:
                        b = matrix_at(ordinary, position)
                        r = matrix_at(incidence, position)
                        o = matrix_at(creation, position)
                        d = [[x ^ y for x, y in zip(first, second)] for first, second in zip(r, o)]
                        next_zero = act(b, zero)
                        previous, selected_value = act(b, one), act(d, zero)
                        zero, one = next_zero, [x ^ y for x, y in zip(previous, selected_value)]
                    scalar = one
                    offset += size
                assert scalar[0] == (actual.bit_count() & 1)
                if actual:
                    assert actual.bit_count() == 1
                    nonzero += 1
                    row, col = boundary // 4, boundary % 4
                    first = g_grades[row][0] - g_grades[col][0]
                    second = g_grades[row][1] - g_grades[col][1]
                    for position in positions:
                        row, col = position // 16, position % 16
                        first += pair_grades[row][0] - pair_grades[col][0]
                        second += pair_grades[row][1] - pair_grades[col][1]
                    for degree in bits(actual):
                        assert first == len(sizes) + degree
                        assert abs(second) <= len(sizes) - degree
                        assert (second - len(sizes) + degree) % 2 == 0
                        grade_checks += 1
                count += 1
        assert nonzero > 0
        records.append({'support_sizes': list(sizes), 'formal_ring_vs_block_flag_entries': count,
                        'nonzero_entries': nonzero, 'exact_grade_controls': grade_checks,
                        'streamed_eta_one_vs_full_polynomial_entries': count,
                        'scalar_operator_channels': 16})
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'check_ward_stable_operators.py',
                                            'check_ward_form_operator_transport.py')},
              'disjoint_support_cases': records}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
