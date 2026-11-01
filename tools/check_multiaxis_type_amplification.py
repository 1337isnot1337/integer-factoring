#!/usr/bin/env python3
"""Check actual common-cut tensors under fixed-axis type amplification.

These finite seeds test identities and indexing. None certifies a better
factoring exponent. The asymptotic transfer has a separate written proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement, product
import json
from math import comb, factorial
from pathlib import Path
import random


def digits(value: int, count: int, base: int) -> tuple[int, ...]:
    out = []
    for _ in range(count):
        out.append(value % base)
        value //= base
    assert value == 0
    return tuple(out)


def encode(values: tuple[int, ...], base: int) -> int:
    out = 0
    for value in reversed(values):
        assert 0 <= value < base
        out = out * base + value
    return out


def endpoint(letters: tuple[int, ...]) -> int:
    return int(all(letters[v] == letters[len(letters) - 1 - v]
                   for v in range(len(letters) // 2)))


def halves(d: int, axis: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left = tuple(v for v in range(1 << d) if not (v >> axis) & 1)
    right = tuple(v for v in range(1 << d) if (v >> axis) & 1)
    return left, right


def factor_binary(rows: list[int], columns: int) -> tuple[int, list[int], list[int]]:
    """Produce a literal row-basis factorization, without rank guesses."""
    basis: list[int] = []
    pivots: list[int] = []
    coordinates = []
    for row in rows:
        remaining = row
        code = 0
        for i, (pivot, vector) in enumerate(zip(pivots, basis)):
            if (remaining >> pivot) & 1:
                remaining ^= vector
                code ^= 1 << i
        if remaining:
            pivots.append((remaining & -remaining).bit_length() - 1)
            basis.append(remaining)
            code ^= 1 << (len(basis) - 1)
        coordinates.append(code)
    for row, code in zip(rows, coordinates):
        reconstructed = 0
        for i, vector in enumerate(basis):
            if (code >> i) & 1:
                reconstructed ^= vector
        assert reconstructed == row
    assert all(row < (1 << columns) for row in basis)
    rank = len(basis)
    width = (rank - 1).bit_length() if rank else 0
    capacity = 1 << width
    return width, coordinates, basis + [0] * (capacity - rank)


def make_seed(d: int, t: int, kind: str) -> dict:
    legs = 1 << d
    m = legs // 2
    base = 1 << t
    entries = base ** legs
    target = [endpoint(digits(i, legs, base)) for i in range(entries)]
    if kind == "dense_three":
        rng = random.Random(1067 + 31 * d + t)
        first = [rng.randrange(2) for _ in range(entries)]
        second = [rng.randrange(2) for _ in range(entries)]
        tensors = [first, second, [e ^ a ^ b for e, a, b in zip(target, first, second)]]
    elif kind == "supported_basis":
        tensors = []
        for index, value in enumerate(target):
            if value:
                tensor = [0] * entries
                tensor[index] = 1
                tensors.append(tensor)
    elif kind == "support_split":
        first = [0] * entries
        first[next(i for i, value in enumerate(target) if value)] = 1
        tensors = [first, [a ^ b for a, b in zip(target, first)]]
    else:
        raise ValueError(kind)
    assert [sum(t[i] for t in tensors) % 2 for i in range(entries)] == target
    factors = []
    width_sums = []
    side = base ** m
    for tensor in tensors:
        cuts = []
        for axis in range(d):
            left, right = halves(d, axis)
            rows = []
            for row in range(side):
                lv = digits(row, m, base)
                packed = 0
                for col in range(side):
                    cv = digits(col, m, base)
                    all_letters = [0] * legs
                    for vertex, value in zip(left, lv):
                        all_letters[vertex] = value
                    for vertex, value in zip(right, cv):
                        all_letters[vertex] = value
                    packed |= tensor[encode(tuple(all_letters), base)] << col
                rows.append(packed)
            width, lfactor, rfactor = factor_binary(rows, side)
            assert width <= m * t
            cuts.append((width, lfactor, rfactor))
        factors.append(cuts)
        width_sums.append(sum(cut[0] for cut in cuts))
    pieces = len(tensors)
    entropy = (pieces - 1).bit_length()
    quality = Fraction(max(width_sums) + d * entropy, d * m * t)
    assert quality >= Fraction(1, d)
    return {"d": d, "t": t, "kind": kind, "pieces": tensors, "factors": factors,
            "base": base, "legs": legs, "m": m, "quality": quality,
            "width_sums": width_sums}


def ceil_root(value: int, degree: int) -> int:
    low, high = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while low + 1 < high:
        mid = (low + high) // 2
        if mid ** degree >= value:
            high = mid
        else:
            low = mid
    assert (high - 1) ** degree < value <= high ** degree
    return high


def rank_sum_controls() -> dict:
    seeds = [make_seed(2, 1, "dense_three"), make_seed(3, 1, "support_split")]
    inequalities = 0
    summaries = []
    for seed in seeds:
        d = seed["d"]
        actual_ranks = [[sum(bool(row) for row in cut[2]) for cut in cuts]
                        for cuts in seed["factors"]]
        products = []
        for ranks in actual_ranks:
            term = 1
            for rank in ranks:
                assert rank > 0
                term *= rank
            products.append(term)
        for denominator in (1, 16, 256):
            upper_terms = [Fraction(ceil_root(term * denominator ** d, d), denominator)
                           for term in products]
            assert all(z ** d >= term for z, term in zip(upper_terms, products))
            upper = sum(upper_terms, Fraction(0))
            for blocks in (1, 2, 3, 4):
                for typ in combinations_with_replacement(range(len(products)), blocks):
                    counts = [typ.count(i) for i in range(len(products))]
                    multinomial = factorial(blocks)
                    for count in counts:
                        multinomial //= factorial(count)
                    cut_product = multinomial ** d
                    for piece in typ:
                        cut_product *= products[piece]
                    assert cut_product <= upper ** (d * blocks)
                    inequalities += 1
        if seed["kind"] == "support_split":
            assert actual_ranks == [[1, 1, 1], [15, 15, 15]]
            assert sum(ceil_root(term, d) for term in products) == 16
            assert seed["quality"] == Fraction(5, 4)
            summaries.append({"axes": d, "actual_cut_ranks": actual_ranks,
                              "exact_rank_sum": 16, "rank_sum_quality": "1",
                              "rounded_maximum_quality": str(seed["quality"])})
    return {"multinomial_rank_sum_inequalities": inequalities,
            "strict_quality_comparison": summaries}


def type_words(typ: tuple[int, ...]) -> list[tuple[int, ...]]:
    """Independent exhaustive word filter for these deliberately tiny blocks."""
    if not typ:
        return [()]
    return [word for word in product(range(max(typ) + 1), repeat=len(typ))
            if tuple(sorted(word)) == typ]


def literal_piece(seed: dict, words: list[tuple[int, ...]], tail: int,
                  letters: tuple[int, ...]) -> int:
    t = seed["t"]
    b = len(words[0])
    mask = (1 << t) - 1
    tail_letters = tuple(value >> (b * t) for value in letters)
    if not endpoint(tail_letters):
        return 0
    answer = 0
    for word in words:
        value = 1
        for block, piece in enumerate(word):
            block_letters = tuple((letter >> (block * t)) & mask for letter in letters)
            value &= seed["pieces"][piece][encode(block_letters, 1 << t)]
        answer ^= value
    return answer


def explicit_cut(seed: dict, words: list[tuple[int, ...]], tail: int,
                 letters: tuple[int, ...], axis: int) -> int:
    """Contract the explicit type-word channel, without literal_piece."""
    d, t, m = seed["d"], seed["t"], seed["m"]
    left, right = halves(d, axis)
    b = len(words[0])
    mask = (1 << t) - 1
    rows = [encode(tuple((letters[v] >> (block * t)) & mask for v in left), 1 << t)
            for block in range(b)]
    cols = [encode(tuple((letters[v] >> (block * t)) & mask for v in right), 1 << t)
            for block in range(b)]
    row_tail = encode(tuple(letters[v] >> (b * t) for v in left), 1 << tail)
    col_tail = tuple(letters[v] >> (b * t) for v in right)
    right_position = {v: i for i, v in enumerate(right)}
    answer = 0
    for word in words:
        local_cuts = [seed["factors"][piece][axis] for piece in word]
        for channels in product(*(range(1 << cut[0]) for cut in local_cuts)):
            lval = 1
            rval = 1
            for block, (channel, cut) in enumerate(zip(channels, local_cuts)):
                _, lfactor, rfactor = cut
                lval &= (lfactor[rows[block]] >> channel) & 1
                rval &= (rfactor[channel] >> cols[block]) & 1
            if not lval or not rval:
                continue
            # The full tail channel is indexed by its left-half labels.
            for tail_channel in range(1 << (m * tail)):
                if tail_channel != row_tail:
                    continue
                decoded = digits(tail_channel, m, 1 << tail)
                matches = all(decoded[i] == col_tail[right_position[(1 << d) - 1 - v]]
                              for i, v in enumerate(left))
                answer ^= int(matches)
    return answer


def support_of_piece(seed: dict, words: list[tuple[int, ...]], tail: int) -> list[tuple[int, ...]]:
    """Actual support controls for supported-basis seeds, built from their arrays."""
    if seed["kind"] != "supported_basis":
        return []
    t, m, legs = seed["t"], seed["m"], seed["legs"]
    points = [digits(next(i for i, v in enumerate(tensor) if v), legs, 1 << t)
              for tensor in seed["pieces"]]
    answer = []
    for word in words:
        for tail_half in product(range(1 << tail), repeat=m):
            letters = [0] * legs
            for block, piece in enumerate(word):
                for v, value in enumerate(points[piece]):
                    letters[v] |= value << (block * t)
            for v in range(m):
                letters[v] |= tail_half[v] << (len(word) * t)
                letters[legs - 1 - v] |= tail_half[v] << (len(word) * t)
            answer.append(tuple(letters))
    return answer


def check_case(d: int, t: int, kind: str, blocks: int, tail: int,
               exhaustive: bool) -> dict:
    seed = make_seed(d, t, kind)
    assert 0 <= tail < t or tail == 0
    pieces = len(seed["pieces"])
    lam = blocks * t + tail
    full_base = 1 << lam
    types = list(combinations_with_replacement(range(pieces), blocks))
    all_points = (list(product(range(full_base), repeat=seed["legs"])) if exhaustive else None)
    rng = random.Random(7019 + 23 * d + 11 * t + tail)
    checks = caps = channel_words = 0
    for typ in types:
        words = type_words(typ)
        counts = [typ.count(i) for i in range(pieces)]
        multinomial = factorial(blocks)
        for count in counts:
            multinomial //= factorial(count)
        assert len(words) == multinomial
        assert multinomial <= pieces ** blocks
        channel_words += len(words)
        raw_widths = []
        for axis in range(d):
            exponent = sum(seed["factors"][i][axis][0] for i in typ) + seed["m"] * tail
            cardinality = multinomial * (1 << exponent)
            width = (cardinality - 1).bit_length()
            raw_widths.append(width)
            if width > seed["m"] * lam:
                caps += 1
        bound = (blocks * (max(seed["width_sums"]) + d * (pieces - 1).bit_length())
                 + d * seed["m"] * tail + d)
        assert sum(raw_widths) <= bound
        points = all_points if all_points is not None else (
            support_of_piece(seed, words, tail)
            + [tuple(rng.randrange(full_base) for _ in range(seed["legs"])) for _ in range(12)])
        for letters in points:
            expected = literal_piece(seed, words, tail, letters)
            for axis in range(d):
                actual = explicit_cut(seed, words, tail, letters, axis)
                assert actual == expected, (d, t, kind, typ, letters, axis)
                checks += 1
                if raw_widths[axis] > seed["m"] * lam:
                    left, _ = halves(d, axis)
                    row = encode(tuple(letters[v] for v in left), full_base)
                    assert row < 1 << (seed["m"] * lam)
                    # Identity cap has just this nonzero left channel; its
                    # right factor is the complete piece entry.
                    assert int(row == row) * expected == actual
                    assert int((row ^ 1) == row) == 0
    endpoint_checks = 0
    if all_points is not None:
        word_lists = [type_words(typ) for typ in types]
        for letters in all_points:
            value = sum(literal_piece(seed, words, tail, letters) for words in word_lists) % 2
            assert value == endpoint(letters)
            endpoint_checks += 1
    else:
        # Sum the independently generated sparse supports. At most one
        # type contributes to each endpoint point in these basis controls.
        support: dict[tuple[int, ...], int] = {}
        for typ in types:
            for point in support_of_piece(seed, type_words(typ), tail):
                support[point] = support.get(point, 0) ^ 1
        expected_size = 1 << (seed["m"] * lam)
        assert len(support) == expected_size
        assert all(value == 1 and endpoint(point) for point, value in support.items())
        endpoint_checks = len(support)
    return {"axes": d, "seed_length": t, "seed_kind": kind,
            "seed_pieces": pieces, "seed_quality": str(seed["quality"]),
            "blocks": blocks, "tail": tail, "types": len(types),
            "type_words": channel_words, "cut_comparisons": checks,
            "capped_cuts": caps, "endpoint_checks": endpoint_checks}


def integer_controls() -> dict:
    selector_checks = 0
    for c1 in (1, 3, 17):
        bound = 64 * c1
        for power in (1, 3, 10, 63, 64, 127, 1024, 4096):
            n = (1 << power) - 1
            tmax = max(1, ((n + 2).bit_length() - 1) // bound)
            assert (1 << (64 * c1 * tmax)) <= (1 << (64 * c1)) * (n + 2)
            selector_checks += 1
    type_count_checks = 0
    for m in range(1, 17):
        for b in range(1, 17):
            assert comb(b + m - 1, m - 1) <= (b + 1) ** m
            type_count_checks += 1
    return {"selector_inequalities": selector_checks,
            "type_count_inequalities": type_count_checks}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [check_case(2, 1, "dense_three", 2, 0, True),
             check_case(3, 1, "supported_basis", 2, 0, False),
             check_case(2, 2, "supported_basis", 2, 1, False)]
    result = {"status": "exact finite amplification controls; no improved seed supplied",
              "complete_factoring_score": "1/2 unchanged", "cases": cases,
              "integer_controls": integer_controls(), "rank_sum_controls": rank_sum_controls()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"cases": len(cases),
                      "cut_comparisons": sum(c["cut_comparisons"] for c in cases),
                      "endpoint_checks": sum(c["endpoint_checks"] for c in cases),
                      "capped_cuts": sum(c["capped_cuts"] for c in cases),
                      "status": "passed; no stronger factoring exponent claimed"}, sort_keys=True))


if __name__ == "__main__":
    main()
