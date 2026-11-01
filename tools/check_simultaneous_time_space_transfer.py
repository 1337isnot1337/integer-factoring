#!/usr/bin/env python3
"""Exact controls for type-blocked crossings and packed transfer bounds.

All matrices here are small independent test fixtures. The asymptotic
oracle uses the streamed algorithms in the accompanying proof, not
these materialized test tables.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from math import comb
from pathlib import Path
import random

from check_streamed_seed_transfer import (
    flatten_column, multinomial, rank_factor, type_unrank,
)
from check_subexponential_piece_transfer import ceil_fraction, schedule_logs
from check_bit_tagged_witness_recovery import PrimeExtension


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, length - 1):
            yield (first,) + rest


def make_seed(bits):
    """An exact crossing with duplicate pieces to exercise real caps."""
    d = 1 << bits
    first = [[0] * (d * d) for _ in range(d * d)]
    other = [[0] * (d * d) for _ in range(d * d)]
    for x, y, u, v in product(range(d), repeat=4):
        r, e = int(x == y and u == v), int(x == v and y == u)
        first[x * d + y][u * d + v] = r
        other[x * d + y][u * d + v] = r ^ e
    # R + 3(E+R) = E over F2. The repeated pieces are intentional.
    entries = (first, other, other, other)
    cuts = []
    for cut in (0, 1):
        factors = []
        for matrix in entries:
            matrix = matrix if cut == 0 else flatten_column(matrix, d)
            left, right = rank_factor(matrix)
            actual = len(right)
            padded = 1 << (actual - 1).bit_length()
            left = [row + [0] * (padded - actual) for row in left]
            right = right + [[0] * (d * d) for _ in range(padded - actual)]
            factors.append((left, right))
        cuts.append(tuple(factors))
    return {"bits": bits, "dimension": d, "entries": entries, "cuts": cuts}


def split_label(value, block_bits, blocks, tail_bits):
    tail = value & ((1 << tail_bits) - 1)
    value >>= tail_bits
    result = [0] * blocks
    for j in range(blocks - 1, -1, -1):
        value, result[j] = divmod(value, 1 << block_bits)
    assert value == 0
    return tuple(result), tail


def direct_piece(seed, counts, coordinates, tail_seed=None, tail_piece=0):
    """Direct coefficient sum, independent of channel unranking."""
    blocks, k = sum(counts), len(counts)
    tail_bits = 0 if tail_seed is None else tail_seed["bits"]
    digits = [split_label(a, seed["bits"], blocks, tail_bits) for a in coordinates]
    result = 0
    for word in product(range(k), repeat=blocks):
        if tuple(word.count(i) for i in range(k)) != counts:
            continue
        term = 1
        for j, piece in enumerate(word):
            x, y, u, v = (a[0][j] for a in digits)
            d = seed["dimension"]
            term &= seed["entries"][piece][x * d + y][u * d + v]
        result ^= term
    if tail_seed is not None:
        x, y, u, v = (a[1] for a in digits)
        d = tail_seed["dimension"]
        result &= tail_seed["entries"][tail_piece][x * d + y][u * d + v]
    return result


def fixture_factors(seed, counts, cut, tail_seed=None, tail_piece=0):
    """Materialize independent scalar queries for finite verification only."""
    blocks = sum(counts)
    tail_bits = 0 if tail_seed is None else tail_seed["bits"]
    label_bits = blocks * seed["bits"] + tail_bits
    d = 1 << label_bits
    original = seed["cuts"][cut]
    ranks = tuple(len(right) for _, right in original)
    tail_factors = None if tail_seed is None else tail_seed["cuts"][cut][tail_piece]
    tail_rank = 1 if tail_factors is None else len(tail_factors[1])
    raw = multinomial(counts) * tail_rank
    for n, rank in zip(counts, ranks):
        raw *= rank ** n
    padded = 1 << (raw - 1).bit_length()
    capped = padded > d * d
    channels = min(padded, d * d)
    left_rows, right_columns = [], []
    zero_padding_checks = unrankings = 0
    if capped:
        left_rows = [1 << pair for pair in range(d * d)]
        for second in range(d * d):
            c, e = divmod(second, d)
            column = 0
            for first in range(d * d):
                a, b = divmod(first, d)
                coords = (a, b, c, e) if cut == 0 else (a, c, b, e)
                column |= direct_piece(seed, counts, coords, tail_seed, tail_piece) << first
            right_columns.append(column)
    else:
        decoded = []
        seen = set()
        for channel in range(raw):
            full_index, tail_index = divmod(channel, tail_rank)
            word, indices = type_unrank(counts, ranks, full_index)
            assert (word, indices, tail_index) not in seen
            seen.add((word, indices, tail_index))
            decoded.append((word, indices, tail_index))
            unrankings += 1
        for first, second in product(range(d), repeat=2):
            aa, at = split_label(first, seed["bits"], blocks, tail_bits)
            bb, bt = split_label(second, seed["bits"], blocks, tail_bits)
            left = right = 0
            for channel, (word, indices, ti) in enumerate(decoded):
                lv = rv = 1
                for j, piece in enumerate(word):
                    ll, rr = original[piece]
                    pair = aa[j] * seed["dimension"] + bb[j]
                    lv &= ll[pair][indices[j]]
                    rv &= rr[indices[j]][pair]
                if tail_factors is not None:
                    ll, rr = tail_factors
                    pair = at * tail_seed["dimension"] + bt
                    lv &= ll[pair][ti]
                    rv &= rr[ti][pair]
                left |= lv << channel
                right |= rv << channel
            assert left >> raw == right >> raw == 0
            zero_padding_checks += 2 * (padded - raw)
            left_rows.append(left)
            right_columns.append(right)
    assert len(left_rows) == len(right_columns) == d * d
    assert all(value.bit_length() <= channels for value in left_rows + right_columns)
    return left_rows, right_columns, {
        "raw_channels": raw, "padded_channels": channels,
        "capped": capped, "unrankings": unrankings,
        "zero_padding_checks": zero_padding_checks,
    }


def algebra_checks():
    seeds = {bits: make_seed(bits) for bits in (1, 2)}
    entries = unrankings = padding = caps = crossing = pieces = 0
    reports = []
    for block_bits, blocks, tail_bits in ((1, 1, 0), (1, 2, 0), (1, 3, 0), (2, 1, 1)):
        seed = seeds[block_bits]
        tail = None if not tail_bits else seeds[tail_bits]
        d = 1 << (block_bits * blocks + tail_bits)
        total = [[0] * (d * d) for _ in range(d * d)]
        local_pieces = local_caps = 0
        for counts in compositions(blocks, len(seed["entries"])):
            for tail_piece in range(1 if tail is None else len(tail["entries"])):
                tables = []
                for cut in (0, 1):
                    left, right, report = fixture_factors(seed, counts, cut, tail, tail_piece)
                    tables.append((left, right))
                    unrankings += report["unrankings"]
                    padding += report["zero_padding_checks"]
                    caps += report["capped"]
                    local_caps += report["capped"]
                for x, y, u, v in product(range(d), repeat=4):
                    expected = direct_piece(seed, counts, (x, y, u, v), tail, tail_piece)
                    for cut, (left, right) in enumerate(tables):
                        a, b = (x * d + y, u * d + v) if cut == 0 else (x * d + u, y * d + v)
                        assert ((left[a] & right[b]).bit_count() & 1) == expected
                        entries += 1
                    total[x * d + y][u * d + v] ^= expected
                local_pieces += 1
                pieces += 1
        assert local_pieces == comb(blocks + 3, 3) * (1 if tail is None else 4)
        for x, y, u, v in product(range(d), repeat=4):
            assert total[x * d + y][u * d + v] == int(x == v and y == u)
            crossing += 1
        reports.append({"block_bits": block_bits, "blocks": blocks,
                        "tail_bits": tail_bits, "label_dimension": d,
                        "pieces": local_pieces, "capped_cuts": local_caps})
    assert caps and padding and unrankings
    return {"two_cut_entry_checks": entries, "channel_unrankings": unrankings,
            "zero_padding_checks": padding, "capped_cuts": caps,
            "crossing_entries": crossing, "type_pieces": pieces,
            "finite_families": reports}


def combinatorial_checks():
    rng = random.Random(20261101)
    cases = 0
    for blocks in range(1, 129):
        for k in range(1, 25):
            assert comb(blocks + k - 1, k - 1) <= (blocks + 1) ** k
            counts = [0] * k
            for _ in range(blocks):
                counts[rng.randrange(k)] += 1
            assert multinomial(counts) <= k ** blocks
            # A factorial counter's bit length is at most b ceil(log2(b+1)).
            fact = 1
            for value in range(1, blocks + 1):
                fact *= value
            assert fact.bit_length() <= 1 + blocks * blocks.bit_length()
            cases += 1
    # A finite exact sample of the endpoint's subpolynomial index entropy.
    endpoint = []
    for exponent in (8, 12, 16, 20, 24):
        label = 1 << exponent
        block = max(1, exponent // 4)
        blocks, tail = divmod(label, block)
        k = (block + 1) ** 2
        index_upper = k * (blocks + 1).bit_length() + (tail + 2).bit_length()
        endpoint.append({"label_bits": label, "block_bits": block,
                         "seed_pieces": k, "index_bits_upper": index_upper,
                         "linear_workspace_budget": label})
        assert index_upper < label
    return {"type_and_factorial_inequalities": cases, "endpoint_samples": endpoint}


def extraction_checks():
    """Stream field evaluations; compare with a separate word expansion."""
    rng = random.Random(20261102)
    cases = points = 0
    fields = {}
    for prime in (2, 3):
        for blocks in range(1, 5):
            degree = 1
            while prime ** degree - 1 <= blocks:
                degree += 1
            field = fields.setdefault((prime, degree), PrimeExtension(prime, degree))
            for k in (1, 2, 3):
                for _ in range(2):
                    forms = [[rng.randrange(prime) for _ in range(k)] for _ in range(blocks)]
                    counts = [0] * k
                    for _ in range(blocks):
                        counts[rng.randrange(k)] += 1
                    expected = 0
                    for word in product(range(k), repeat=blocks):
                        if [word.count(i) for i in range(k)] != counts:
                            continue
                        term = 1
                        for j, i in enumerate(word):
                            term = term * forms[j][i] % prime
                        expected = (expected + term) % prime
                    actual = 0
                    for vector in product(range(field.size), repeat=k):
                        value = 1
                        for form in forms:
                            linear = 0
                            for coefficient, coordinate in zip(form, vector):
                                linear = field.add(linear, field.multiply(coefficient, coordinate))
                            value = field.multiply(value, linear)
                        for count, coordinate in zip(counts, vector):
                            value = field.multiply(value, field.power(coordinate, field.size - 1 - count))
                        actual = field.add(actual, value)
                        points += 1
                    actual = field.multiply(pow(-1, k, prime), actual)
                    assert actual == expected
                    cases += 1
    # The strict degree cutoff is necessary: Z^3 aliases the constant
    # coefficient if F4 is incorrectly used with degree bound three.
    field = PrimeExtension(2, 2)
    aliased = 0
    for value in range(field.size):
        aliased = field.add(aliased, field.power(value, 6))
    assert aliased == 1  # The actual constant coefficient of Z^3 is zero.
    # In odd characteristic the (-1)^k normalization is necessary.
    field = PrimeExtension(3, 1)
    raw = 0
    for value in range(field.size):
        raw = field.add(raw, field.power(value, 2))
    assert raw == 2 and field.multiply(2, raw) == 1
    return {"coefficient_checks": cases, "field_vectors_evaluated": points,
            "primes": [2, 3], "negative_controls": ["degree aliasing", "odd-characteristic sign"]}


def packed_cost_checks():
    logs = schedule_logs(5)
    paths = steps = entropy = 0
    for level in range(1, 6):
        root_log = logs[level]
        n = 1 << (4 * root_log)
        cutoff = 1 << (3 * root_log)
        smallest_label = 1 << root_log
        slack = Fraction(1, level + 1)
        rounding = Fraction(1, 100 * (level + 1))
        for rho in (Fraction(1, 2), Fraction(2, 3), Fraction(3, 4), Fraction(7, 8)):
            envelope = Fraction(1)
            for j in range(level):
                envelope *= 1 + (3 * Fraction(1, 1 << (j + 5)) + slack + rounding) / rho
            assert envelope < 27
            for forced_stop in sorted({0, level, max(1, level // 2)}):
                j, volume = level, n
                volume_sum = query_bits = width_bits = 0
                index_bits = {1: 0, 2: 0}
                while j and volume > cutoff:
                    parent_size_squared = 1 << (2 * logs[j])
                    child_size = 1 << logs[j - 1]
                    label = volume // parent_size_squared
                    assert volume % parent_size_squared == 0 and label > smallest_label
                    label_log = label.bit_length() - 1
                    assert label_log >= root_log
                    assert label <= volume // child_size
                    volume_sum += volume
                    query_bits += volume // child_size
                    width_bits += parent_size_squared * (2 * label + 1).bit_length()
                    for delta in index_bits:
                        per_piece = (label + label_log ** delta - 1) // label_log ** delta
                        index_bits[delta] += parent_size_squared * per_piece
                    zeta = Fraction(1, 1 << (j + 4))
                    raw = (rho + 3 * zeta + slack) * volume
                    child_square = child_size * child_size
                    next_volume = child_square * ceil_fraction(raw / child_square)
                    assert next_volume <= (rho + 3 * zeta + slack + rounding) * volume
                    if j == forced_stop:
                        next_volume = child_square * max(1, cutoff // (2 * child_square))
                    volume, j = next_volume, j - 1
                    steps += 1
                assert volume_sum <= envelope * n / (1 - rho)
                assert query_bits <= envelope * n * rho ** level / ((1 << 17) - 1)
                assert width_bits * smallest_label <= volume_sum * (4 * root_log + 4)
                for delta, bits in index_bits.items():
                    assert bits <= volume_sum * (Fraction(1, root_log ** delta)
                                                 + Fraction(1, smallest_label))
                    entropy += 1
                assert volume <= envelope * n * rho ** level + cutoff
                # Hash descriptors fit even at rho=1/2 (alpha=1).
                q = (9 * n + 4 * root_log - 1) // (4 * root_log)
                assert q <= n * rho ** level
                paths += 1
    return {"packed_paths": paths, "rounded_grid_steps": steps,
            "compact_index_sum_checks": entropy,
            "rhos_tested": ["1/2", "2/3", "3/4", "7/8"],
            "scope": "Exact declared cost envelopes, including rho=1/2; no new crossing rate."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"status": "PASS", "algebra": algebra_checks(),
              "combinatorics": combinatorial_checks(),
              "coefficient_extraction": extraction_checks(),
              "packed_costs": packed_cost_checks(),
              "scope": "Conditional transfer and exact type blocking; no better contracting seed is supplied."}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"],
                      "two_cut_entry_checks": result["algebra"]["two_cut_entry_checks"],
                      "channel_unrankings": result["algebra"]["channel_unrankings"],
                      "capped_cuts": result["algebra"]["capped_cuts"],
                      "packed_paths": result["packed_costs"]["packed_paths"]}))


if __name__ == "__main__":
    main()
