#!/usr/bin/env python3
"""Exact checks for logarithmic blocking and two-cut padded factors."""

import argparse
from itertools import product
import json
from math import isqrt
import random

import blocked_modular_factors as blocked
import linear_modular_factors as local


def word(value, width):
    return format(value, '0%db' % width) if width else ''


def verify():
    from check_modular_crossing import normalized_piece
    rng = random.Random(184729)
    lengths = set(range(1, 4097))
    for power in range(1, 257):
        for delta in (-2, -1, 0, 1, 2):
            if (1 << power) + delta > 0:
                lengths.add((1 << power) + delta)
    schedule_cases = 0
    for length in sorted(lengths):
        logarithm = (length + 1).bit_length()
        assert 1 << (logarithm - 1) < length + 2 <= 1 << logarithm
        root = isqrt(logarithm)
        root += root * root < logarithm
        assert (root - 1) ** 2 < logarithm <= root ** 2
        r = blocked.chosen_block_size(length)
        assert r == min(length, root * (logarithm + 1).bit_length())
        blocks = (length + r - 1) // r
        assert blocks * r <= 2 * length
        # h_t<=8t² gives at most32t²+1 local shifts. This conservative
        # integer bound tests the entropy accounting without computing h_t.
        entropy = blocks * (32 * r * r).bit_length()
        assert entropy * entropy * logarithm <= 16 ** 2 * length ** 2
        schedule_cases += 1

    # Two one-bit blocks exercise all tuple codes independently of the
    # asymptotic schedule. The construction/factor code is identical.
    h = local.label_weight(1)
    assert h == 8
    labels = tuple(tuple(local.coded_parts(local.partition_unrank(h, i))) for i in range(2))
    shifts = tuple(range(-2 * h, 2 * h + 1))
    tables = {s: normalized_piece(2, h, s, labels) for s in shifts}
    for x, y, v, u in product(range(2), repeat=4):
        total = 0
        for s in shifts:
            total ^= tables[s][2 * x + y][2 * v + u]
        assert total == int(x == v and y == u)

    code = blocked.piece_code(2, (-2 * h, -2 * h), block_size=1)
    coding_cases = 0
    for first, second in product(shifts, repeat=2):
        assert code == blocked.piece_code(2, (first, second), block_size=1)
        decoded = tuple(blocked.decode_piece(2, code, block_size=1))
        assert decoded == ((1, h, first), (1, h, second))
        code = blocked.next_piece(2, code, block_size=1)
        coding_cases += 1
    assert code is None

    # Check both local factorizations, including explicit identity caps and
    # zero padding, against the independent normal-order source matrices.
    local_factorizations = capped_checks = 0
    for shift in (-2 * h, -1, 0, 2 * h):
        for row_cut in (False, True):
            width, count = blocked.local_width(1, h, shift, row_cut)
            assert width <= 2
            for x, y, u, v in product(range(2), repeat=4):
                a, b, c, d = (x, y, v, u) if row_cut else (x, u, v, y)
                value = 0
                for channel in range(1 << width):
                    value ^= (blocked.local_factor(1, h, shift, row_cut, True, a, b, channel, count)
                              * blocked.local_factor(1, h, shift, row_cut, False, c, d, channel, count))
                assert value == tables[shift][2 * x + y][2 * v + u]
                local_factorizations += 1
                capped_checks += count > 4

    # Exact complete two-block crossing, using independent local matrices.
    tensor_entries = crossing_entries = 0
    for x, y, u, v in product(range(4), repeat=4):
        total = 0
        for s, t in product(shifts, repeat=2):
            entry = (tables[s][2 * (x >> 1) + (y >> 1)][2 * (v >> 1) + (u >> 1)]
                     * tables[t][2 * (x & 1) + (y & 1)][2 * (v & 1) + (u & 1)])
            total ^= entry
            tensor_entries += 1
        assert total == int(x == v and y == u)
        crossing_entries += 1

    # Query the production bit-string wrapper on both cuts. This catches
    # channel-block offsets that a plain tensor-product identity cannot.
    wrapped_factorizations = wrapped_piece_entries = 0
    for _ in range(20):
        ss = tuple(rng.choice((-16, -1, 0, 16)) for _ in range(2))
        code = blocked.piece_code(2, ss, block_size=1)
        x, y, u, v = (rng.randrange(4) for _ in range(4))
        expected = (tables[ss[0]][2 * (x >> 1) + (y >> 1)][2 * (v >> 1) + (u >> 1)]
                    * tables[ss[1]][2 * (x & 1) + (y & 1)][2 * (v & 1) + (u & 1)])
        strings = tuple(word(value, 2) for value in (x, y, u, v))
        assert blocked.piece_entry(2, code, *strings, block_size=1) == expected
        wrapped_piece_entries += 1
        rw, cw = blocked.widths(2, code, block_size=1)
        assert rw <= 4 and cw <= 4
        for row_cut, width in ((True, rw), (False, cw)):
            a, b, c, d = (x, y, v, u) if row_cut else (x, u, v, y)
            value = 0
            for channel in range(1 << width):
                bits = word(channel, width)
                value ^= (blocked.factor(2, code, row_cut, True, word(a, 2), word(b, 2), bits, block_size=1)
                          * blocked.factor(2, code, row_cut, False, word(c, 2), word(d, 2), bits, block_size=1))
            assert value == expected
            wrapped_factorizations += 1

    # The default schedule is used directly at small lengths, while the
    # preceding tests exercise multiple blocks without a large finite table.
    default_cases = 0
    for length in (1, 2, 3):
        layout = tuple(blocked.block_layout(length))
        assert sum(t for t, _, _ in layout) == length
        ss = tuple(-2 * weight for _, weight, _ in layout)
        code = blocked.piece_code(length, ss)
        assert tuple((t, h, s) for t, h, s in blocked.decode_piece(length, code)) == tuple(
            (t, h, -2 * h) for t, h, _ in layout)
        default_cases += 1

    # A valid padded-width channel need not be a live source channel.
    padding_found = False
    for shift in shifts:
        for row_cut in (False, True):
            width, count = blocked.local_width(1, h, shift, row_cut)
            if 0 < count < (1 << width) and count <= 4:
                assert blocked.local_factor(1, h, shift, row_cut, True, 0, 0, count, count) == 0
                padding_found = True
    assert padding_found
    malformed = '1' * len(blocked.piece_code(2, (0, 0), block_size=1))
    try:
        tuple(blocked.decode_piece(2, malformed, block_size=1))
    except ValueError:
        pass
    else:
        raise AssertionError('An out-of-range tuple digit was accepted')

    return {
        'status': 'pass',
        'scope': 'exact block construction and indexing; asymptotic rate follows from the accompanying proof',
        'integer_schedule_and_entropy_cases': schedule_cases,
        'piece_tuple_codes_and_odometer_steps': coding_cases,
        'independent_local_factorizations': local_factorizations,
        'identity_cap_factorizations': capped_checks,
        'explicit_two_block_piece_entries': tensor_entries,
        'explicit_two_block_crossing_entries': crossing_entries,
        'production_wrapped_factorizations': wrapped_factorizations,
        'production_wrapped_piece_entries': wrapped_piece_entries,
        'default_schedule_small_cases': default_cases,
        'negative_controls': ['padded factor channel', 'out-of-range tuple digit'],
        'source_revision': '16a5ce41304fcdaea4e0be7c132e1c361dc26da5',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    encoded = json.dumps(verify(), indent=2, sort_keys=True) + '\n'
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end='')
