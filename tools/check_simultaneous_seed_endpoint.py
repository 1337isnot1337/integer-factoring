#!/usr/bin/env python3
"""Exact checks for the simultaneous finite-seed endpoint interface.

This is a finite algebra/indexing check, not a factoring implementation and
not evidence of any crossing parameter below the half-rate threshold.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from math import comb
from pathlib import Path


def ceil_log2(value: int) -> int:
    assert value >= 1
    return (value - 1).bit_length()


def tuples(length: int, alphabet: int):
    return itertools.product(range(alphabet), repeat=length)


def flat_index(label: tuple[int, ...], d: int) -> int:
    out = 0
    for value in label:
        out = out * d + value
    return out


def crossing(d: int) -> list[int]:
    return [int(x == v and y == u) for x, y, u, v in tuples(4, d)]


def flatten(tensor: list[int], d: int, cut: str) -> list[list[int]]:
    order = (0, 1, 2, 3) if cut == "row" else (0, 2, 1, 3)
    matrix = [[0] * (d * d) for _ in range(d * d)]
    for labels in tuples(4, d):
        i, j, k, ell = (labels[v] for v in order)
        matrix[d * i + j][d * k + ell] = tensor[flat_index(labels, d)]
    return matrix


def pivot_columns(matrix: list[list[int]], p: int) -> list[int]:
    work = [[x % p for x in row] for row in matrix]
    rank = 0
    columns: list[int] = []
    for col in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][col], -1, p)
        work[rank] = [(x * inverse) % p for x in work[rank]]
        for i in range(len(work)):
            if i != rank and work[i][col]:
                scale = work[i][col]
                work[i] = [(x - scale * y) % p for x, y in zip(work[i], work[rank])]
        columns.append(col)
        rank += 1
        if rank == len(work):
            break
    return columns


def inverse_matrix(matrix: list[list[int]], p: int) -> list[list[int]]:
    size = len(matrix)
    work = [row[:] + [int(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for col in range(size):
        pivot = next(i for i in range(col, size) if work[i][col] % p)
        work[col], work[pivot] = work[pivot], work[col]
        inverse = pow(work[col][col] % p, -1, p)
        work[col] = [(x * inverse) % p for x in work[col]]
        for i in range(size):
            if i != col:
                scale = work[i][col]
                work[i] = [(x - scale * y) % p for x, y in zip(work[i], work[col])]
    return [row[size:] for row in work]


def rank_factors(matrix: list[list[int]], p: int):
    columns = pivot_columns(matrix, p)
    assert columns
    left = [[row[j] % p for j in columns] for row in matrix]
    rows = pivot_columns([list(row) for row in zip(*left)], p)
    inv = inverse_matrix([left[i] for i in rows], p)
    right = [
        [sum(inv[i][j] * matrix[rows[j]][k] for j in range(len(rows))) % p
         for k in range(len(matrix[0]))]
        for i in range(len(columns))
    ]
    assert all(
        sum(left[i][k] * right[k][j] for k in range(len(columns))) % p == matrix[i][j] % p
        for i in range(len(matrix)) for j in range(len(matrix[0]))
    )
    return len(columns), left, right


def seed(a: int, p: int, pieces: list[list[int]]) -> dict:
    d = 1 << a
    assert all(len(tensor) == d**4 for tensor in pieces)
    assert all(sum(tensor[j] for tensor in pieces) % p == v for j, v in enumerate(crossing(d)))
    pieces = [tensor for tensor in pieces if any(x % p for x in tensor)]
    facts = [{cut: rank_factors(flatten(tensor, d, cut), p) for cut in ("row", "col")}
             for tensor in pieces]
    w = max(ceil_log2(f["row"][0]) + ceil_log2(f["col"][0]) for f in facts)
    quality = Fraction(w + 2 * ceil_log2(len(pieces)), 4 * a)
    assert quality >= Fraction(1, 2)
    return {"a": a, "p": p, "d": d, "pieces": pieces, "factors": facts, "quality": quality}


def two_piece_seed(p: int) -> dict:
    d = 2
    first = [int(x == y and u == v) for x, y, u, v in tuples(4, d)]
    second = [(x - y) % p for x, y in zip(crossing(d), first)]
    answer = seed(1, p, [first, second])
    assert [(f["row"][0], f["col"][0]) for f in answer["factors"]] == [(1, 4), (4, 1)]
    assert answer["quality"] == 1
    return answer


def tensor_seed(base: dict, power: int) -> dict:
    d, p = base["d"], base["p"]
    full_d = d**power
    output = []
    for word in tuples(power, len(base["pieces"])):
        tensor = []
        for labels in tuples(4, full_d):
            value = 1
            for j, piece in enumerate(word):
                local = tuple((x // d**(power - j - 1)) % d for x in labels)
                value = value * base["pieces"][piece][flat_index(local, d)] % p
            tensor.append(value)
        output.append(tensor)
    return seed(base["a"] * power, p, output)


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for value in range(total + 1):
        for rest in compositions(total - value, length - 1):
            yield (value,) + rest


def bounded_multinomial(counts: tuple[int, ...]) -> tuple[int, int]:
    remaining = sum(counts)
    output, max_numerator = 1, 1
    for count in counts:
        value = 1
        # Only the increasing half of the binomial row has bounded partials.
        for j in range(1, min(count, remaining - count) + 1):
            numerator = value * (remaining - j + 1)
            max_numerator = max(max_numerator, output * numerator)
            assert numerator % j == 0
            value = numerator // j
        output *= value
        remaining -= count
    return output, max_numerator


def unrank_word(counts: tuple[int, ...], ordinal: int) -> tuple[int, ...]:
    counts = list(counts)
    words = bounded_multinomial(tuple(counts))[0]
    assert 0 <= ordinal < words
    result = []
    while sum(counts):
        length = sum(counts)
        for i, count in enumerate(counts):
            numerator = words * count
            assert numerator % length == 0
            starts = numerator // length
            if ordinal < starts:
                result.append(i)
                counts[i] -= 1
                words = starts
                break
            ordinal -= starts
        else:
            raise AssertionError("Unranking lost the ordinal")
    return tuple(result)


def direct_channel_matrix(source: dict, width: int, counts: tuple[int, ...], cut: str):
    a, p, d = source["a"], source["p"], source["d"]
    k, b = divmod(width, a)
    full_d, tail_d = 1 << width, 1 << b
    word_count = bounded_multinomial(counts)[0]
    channels_per_word = tail_d**2
    for i, multiplicity in enumerate(counts):
        channels_per_word *= source["factors"][i][cut][0]**multiplicity
    cardinality = word_count * channels_per_word
    result = [[0] * (full_d**2) for _ in range(full_d**2)]
    channel_records = []
    for channel in range(cardinality):
        ordinal, residual = divmod(channel, channels_per_word)
        word = unrank_word(counts, ordinal)
        residual, tail_channel = divmod(residual, tail_d**2)
        locals_rev = []
        for letter in reversed(word):
            residual, local = divmod(residual, source["factors"][letter][cut][0])
            locals_rev.append(local)
        assert residual == 0
        channel_records.append((word, tuple(reversed(locals_rev)), tail_channel))
    assert len(set(channel_records)) == cardinality
    # Keep left/right lists only in the finite checker, not in the proposed algorithm.
    left = [[0] * cardinality for _ in range(full_d**2)]
    right = [[0] * full_d**2 for _ in range(cardinality)]
    tail = flatten(crossing(tail_d), tail_d, cut)
    for channel, (word, local_channels, tail_channel) in enumerate(channel_records):
        for pair in tuples(2, full_d):
            row = flat_index(pair, full_d)
            left_value = int(flat_index(tuple(x % tail_d for x in pair), tail_d) == tail_channel)
            right_value = tail[tail_channel][flat_index(tuple(x % tail_d for x in pair), tail_d)]
            for j, (letter, local) in enumerate(zip(word, local_channels)):
                local_pair = tuple((x >> (b + a * (k - j - 1))) % d for x in pair)
                local_row = flat_index(local_pair, d)
                _, seed_left, seed_right = source["factors"][letter][cut]
                left_value = left_value * seed_left[local_row][local] % p
                right_value = right_value * seed_right[local][local_row] % p
            left[row][channel] = left_value
            right[channel][row] = right_value
    for i in range(full_d**2):
        for j in range(full_d**2):
            result[i][j] = sum(left[i][ell] * right[ell][j] for ell in range(cardinality)) % p
    return result, cardinality


def literal_type_tensor(source: dict, width: int, counts: tuple[int, ...]) -> list[int]:
    a, p, d = source["a"], source["p"], source["d"]
    k, b = divmod(width, a)
    words = [word for word in tuples(k, len(source["pieces"]))
             if tuple(word.count(i) for i in range(len(source["pieces"]))) == counts]
    answer = []
    for labels in tuples(4, 1 << width):
        tail = tuple(x % (1 << b) for x in labels)
        if tail[0] != tail[3] or tail[1] != tail[2]:
            answer.append(0)
            continue
        total = 0
        for word in words:
            product = 1
            for j, letter in enumerate(word):
                local = tuple((x >> (b + a * (k - j - 1))) % d for x in labels)
                product = product * source["pieces"][letter][flat_index(local, d)] % p
            total += product
        answer.append(total % p)
    return answer


def run_checks() -> dict:
    tensor_entries = type_count = channel_count = word_checks = scalar_checks = 0
    seed_rows = []
    cases = []
    for p in (2, 3):
        base = two_piece_seed(p)
        doubled = tensor_seed(base, 2)
        seed_rows.append({"p": p, "base_ranks": [[f[cut][0] for cut in ("row", "col")] for f in base["factors"]],
                          "double_quality": str(doubled["quality"])})
        cases.extend((base, width) for width in (1, 2, 3))
        cases.append((doubled, 3))
    for source, width in cases:
        k, b = divmod(width, source["a"])
        summed = [0] * ((1 << width)**4)
        for counts in compositions(k, len(source["pieces"])):
            tensor = literal_type_tensor(source, width, counts)
            ranks = []
            for cut in ("row", "col"):
                result, cardinality = direct_channel_matrix(source, width, counts, cut)
                assert result == flatten(tensor, 1 << width, cut)
                ranks.append(cardinality)
                tensor_entries += len(tensor)
                channel_count += cardinality
            assert sum(ceil_log2(r) for r in ranks) <= 4 * source["quality"] * width + 4 * source["a"] + 2
            summed = [(x + y) % source["p"] for x, y in zip(summed, tensor)]
            type_count += 1
        assert summed == crossing(1 << width)

    for m in range(1, 6):
        for k in range(0, 9):
            for counts in compositions(k, m):
                value, numerator = bounded_multinomial(counts)
                independent = 1
                remaining = k
                for count in counts:
                    independent *= comb(remaining, count)
                    remaining -= count
                assert value == independent
                assert value <= m**k
                assert numerator <= max(1, value * max(1, k))
                scalar_checks += 1
                if m**k <= 256:
                    literal = [word for word in tuples(k, m)
                               if tuple(word.count(i) for i in range(m)) == counts]
                    assert [unrank_word(counts, i) for i in range(value)] == literal
                    word_checks += value

    # Test bounded delayed eligibility without treating abstract qualities as seeds.
    binary = two_piece_seed(2)
    ternary = two_piece_seed(3)
    records = [(4, binary), (17, ternary)]
    eligibility = []
    for budget in (0, 3, 4, 16, 17, 100):
        available = [(delay, item) for delay, item in records if delay <= budget]
        selected = min(available, key=lambda item: (item[1]["quality"], item[0])) if available else None
        eligibility.append({"budget": budget, "valid_records": len(available),
                            "selected_characteristic": selected[1]["p"] if selected else None})
    assert [row["valid_records"] for row in eligibility] == [0, 0, 1, 1, 2, 2]

    bad = [row[:] for row in binary["pieces"]]
    bad[0][0] ^= 1
    try:
        seed(1, 2, bad)
    except AssertionError:
        rejected_bad_crossing = True
    else:
        raise AssertionError("Invalid crossing was accepted")

    # Slow convergence example: exact rational recurrence only, no tensor assertion.
    slow = []
    for m in range(4, 41):
        beta = Fraction(1, 2) + Fraction(1, 8 * m)
        depth = m * m
        assert beta**10 * 2**9 < 1  # Below the 0.9-score contraction threshold.
        ratio = (2 * beta)**depth
        assert ratio >= 1 + Fraction(m, 4)  # Bernoulli lower bound.
        scalar_checks += 2
        if m in (4, 8, 16, 32, 40):
            slow.append({"m": m, "depth": depth, "beta": str(beta),
                         "ratio_lower_bound": str(1 + Fraction(m, 4))})

    path = Path(__file__)
    return {
        "claim": "Finite verification of direct type channels and bounded-space seed endpoint interfaces",
        "new_factoring_exponent": None,
        "coefficient_fields": [2, 3],
        "seeds": seed_rows,
        "amplification_cases": len(cases),
        "type_pieces": type_count,
        "flattening_entries_compared": tensor_entries,
        "explicit_channels_unranked": channel_count,
        "type_words_unranked": word_checks,
        "exact_scalar_checks": scalar_checks,
        "delayed_selection": eligibility,
        "negative_control_invalid_crossing_rejected": rejected_bad_crossing,
        "slow_scalar_convergence_only": slow,
        "script_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "scope": "Finite algebra/indexing evidence; no new asymptotic rank profile and no execution of the factoring machine.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
