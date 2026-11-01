#!/usr/bin/env python3
"""Exact ranks for a characteristic-two reciprocal pencil on energy shells.

This exploratory checker computes genuine restricted matrices.  It never
replaces their ranks by a graded-state dimension.  All arithmetic is over F2.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
from math import comb, log2, prod, sqrt
import json
import time


@lru_cache(None)
def layer_count(counts: tuple[int, ...], energy: int) -> int:
    if energy < 0:
        return 0
    if not counts:
        return int(energy == 0)
    weight = 1 << (len(counts) - 1)
    return sum(
        comb(counts[-1], a) * layer_count(counts[:-1], energy - a * weight)
        for a in range(min(counts[-1], energy // weight) + 1)
    )


def layer_labels(groups: tuple[tuple[int, ...], ...], energy: int):
    if energy < 0:
        return []
    if not groups:
        return [0] if energy == 0 else []
    weight = 1 << (len(groups) - 1)
    result = []
    for a in range(min(len(groups[-1]), energy // weight) + 1):
        lower = layer_labels(groups[:-1], energy - a * weight)
        if not lower:
            continue
        for chosen in combinations(groups[-1], a):
            mask = sum(1 << i for i in chosen)
            result.extend(mask | low for low in lower)
    return result


def rank_two(rows):
    pivots = {}
    for value in rows:
        while value:
            p = value.bit_length() - 1
            if p in pivots:
                value ^= pivots[p]
            else:
                pivots[p] = value
                break
    return len(pivots)


@lru_cache(None)
def q_power(exponent):
    result, bit = 1, 0
    while exponent:
        weight = 1 << bit
        if exponent & 1:
            result ^= (result << weight) ^ (result << (2 * weight))
        exponent >>= 1
        bit += 1
    return result


@lru_cache(maxsize=500000)
def primitive_polynomial(value, R):
    h = 1 << R
    W = sum((s + 1) << s for s in range(R + 1))
    divisor = q_power(W - 2 * h)
    quotient = 0
    while value and value.bit_length() >= divisor.bit_length():
        shift = value.bit_length() - divisor.bit_length()
        quotient ^= 1 << shift
        value ^= divisor << shift
    assert value == 0
    result = 0
    while quotient:
        bit = quotient & -quotient
        degree = bit.bit_length() - 1
        assert degree % 2 == 0
        result |= 1 << (degree // 2)
        quotient ^= bit
    return result


@lru_cache(maxsize=500000)
def primitive_entry(base, u_power, q_exponent):
    """z^base (1+z)^u_power (1+z+z^2)^q_exponent."""
    assert min(base, u_power, q_exponent) >= 0
    result, bit = 1 << base, 0
    while u_power or q_exponent:
        weight = 1 << bit
        if u_power & 1:
            result ^= result << weight
        if q_exponent & 1:
            result ^= (result << weight) ^ (result << (2 * weight))
        u_power >>= 1
        q_exponent >>= 1
        bit += 1
    return result


def projected_entry(a, b, shift, dg, do):
    on_difference = sum(c << s for s, c in enumerate(dg))
    off_difference = sum(c << s for s, c in enumerate(do))
    assert on_difference % 2 == off_difference % 2 == 0
    i, j = a - on_difference // 2, b - off_difference // 2
    return primitive_entry(shift // 2 + a, i + b - j, j)


@lru_cache(None)
def tilting_weights(n):
    value = n + 1
    top = value.bit_length() - 1
    return tuple([1 << i for i in range(top)] +
                 [1 << i for i in range(top) if value >> i & 1])


@lru_cache(None)
def tilting_character(n):
    result = {0: 1}
    for weight in tilting_weights(n):
        following = defaultdict(int)
        for k, count in result.items():
            following[k - weight] += count
            following[k + weight] += count
        result = dict(following)
    return result


@lru_cache(None)
def tensor_tiltings(m):
    character = {m - 2 * j: comb(m, j) for j in range(m + 1)}
    answer = []
    for n in range(m, -1, -1):
        mu = character.get(n, 0)
        if mu:
            assert mu > 0
            answer.append((n, mu))
            for k, count in tilting_character(n).items():
                character[k] = character.get(k, 0) - mu * count
                assert character[k] >= 0
    assert not any(character.values())
    assert sum(mu * 2 ** len(tilting_weights(n)) for n, mu in answer) == 2 ** m
    return tuple(answer)


@lru_cache(None)
def local_bicolored_tiltings(m):
    result = defaultdict(int)
    for k in range(m + 1):
        for n, mu in tensor_tiltings(k):
            for p, nu in tensor_tiltings(m - k):
                result[n, p] += comb(m, k) * mu * nu
    assert sum(mu * 2 ** (len(tilting_weights(n)) + len(tilting_weights(p)))
               for (n, p), mu in result.items()) == 4 ** m
    return tuple((n, p, mu) for (n, p), mu in result.items())


def tilting_blocks(R):
    """Actual direct-sum blocks after both energy constraints.

    The local multiplicity sums all XOR counts k.  Its determinant twist
    depends on m-n-p, so that k no longer changes the resulting matrix.
    """
    h = 1 << R
    W = sum((s + 1) << s for s in range(R + 1))
    length = R + (R + 1).bit_length() + 1
    aggregate = defaultdict(int)

    def walk(s, budget, on, off, mu):
        if s < 0:
            Won = sum(c << i for i, c in enumerate(on))
            Woff = sum(c << i for i, c in enumerate(off))
            if Won & 1 or (W - Woff) & 1:
                return
            a, b = Won // 2, h - (W - Woff) // 2
            if b < 0:
                return
            dimension = layer_count(tuple(on), a) * layer_count(tuple(off), b)
            if dimension:
                shift = 2 * (W - Won - Woff)
                aggregate[tuple(on), tuple(off), a, b, shift] += mu
            return
        m, weight = s + 1, 1 << s
        for n, p, multiplicity in local_bicolored_tiltings(m):
            cost = (m - p) * weight
            if cost > budget:
                continue
            on1, off1 = list(on), list(off)
            for w in tilting_weights(n):
                on1[s + w.bit_length() - 1] += 1
            for w in tilting_weights(p):
                off1[s + w.bit_length() - 1] += 1
            walk(s - 1, budget - cost, on1, off1, mu * multiplicity)

    walk(R, 2 * h, [0] * length, [0] * length, 1)
    result = []
    for (on, off, a, b, shift), mu in aggregate.items():
        dimension = layer_count(on, a) * layer_count(off, b)
        result.append((on, off, a, b, shift, mu, dimension))
    D = layer_count(tuple(range(1, R + 2)), h)
    assert sum(mu * dimension for *_, mu, dimension in result) == D * D
    return result


@lru_cache(None)
def reduced_entry_polynomial(on, off, distances_on, distances_off, shift):
    result = 1 << (shift + sum(c << s for s, c in enumerate(on)))
    for s, (k, c, dg, do) in enumerate(zip(on, off, distances_on, distances_off)):
        a, b = k - dg + do, c - do
        bit = 0
        while a or b:
            weight = 1 << (s + bit)
            if a & 1:
                result ^= result << weight
            if b & 1:
                result ^= (result << weight) ^ (result << (2 * weight))
            a >>= 1
            b >>= 1
            bit += 1
    return result


def reduced_block_ranks(on, off, a, b, shift, primitive_R=None):
    offset, on_groups, off_groups = 0, [], []
    for k, c in zip(on, off):
        on_groups.append(tuple(range(offset, offset + k)))
        off_groups.append(tuple(range(offset + k, offset + k + c)))
        offset += k + c
    left, right = layer_labels(tuple(on_groups), a), layer_labels(tuple(off_groups), b)
    labels = [x | y for x in left for y in right]
    on_masks = tuple(sum(1 << i for i in group) for group in on_groups)
    off_masks = tuple(sum(1 << i for i in group) for group in off_groups)
    dimension = len(labels)
    matrices = defaultdict(lambda: [0] * dimension)
    for row, x in enumerate(labels):
        for col in range(row, dimension):
            difference = x ^ labels[col]
            dg = tuple((difference & mask).bit_count() for mask in on_masks)
            do = tuple((difference & mask).bit_count() for mask in off_masks)
            polynomial = (
                projected_entry(a, b, shift, dg, do) if primitive_R is not None else
                reduced_entry_polynomial(on, off, dg, do, shift))
            while polynomial:
                bit = polynomial & -polynomial
                degree = bit.bit_length() - 1
                matrices[degree][row] |= 1 << col
                matrices[degree][col] |= 1 << row
                polynomial ^= bit
    return dimension, {degree: rank_two(rows) for degree, rows in matrices.items()}


def block_contains_central_monomial(on, off, a, b, shift, R):
    offset, on_groups, off_groups = 0, [], []
    for k, c in zip(on, off):
        on_groups.append(tuple(range(offset, offset + k)))
        off_groups.append(tuple(range(offset + k, offset + k + c)))
        offset += k + c
    left, right = layer_labels(tuple(on_groups), a), layer_labels(tuple(off_groups), b)
    on_masks = tuple(sum(1 << i for i in group) for group in on_groups)
    off_masks = tuple(sum(1 << i for i in group) for group in off_groups)
    on_distances = {
        tuple(((x ^ y) & mask).bit_count() for mask in on_masks)
        for i, x in enumerate(left) for y in left[i:]
    }
    off_distances = {
        tuple(((x ^ y) & mask).bit_count() for mask in off_masks)
        for i, x in enumerate(right) for y in right[i:]
    }
    basis = {}
    target = 1 << (1 << R)
    for dg, do in product(on_distances, off_distances):
        value = projected_entry(a, b, shift, dg, do)
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
        remainder = target
        while remainder:
            pivot = remainder.bit_length() - 1
            if pivot not in basis:
                break
            remainder ^= basis[pivot]
        if not remainder:
            return True, len(basis)
    return False, len(basis)


def central_detector_bound(R, progress=False):
    started = time.monotonic()
    D = layer_count(tuple(range(1, R + 2)), 1 << R)
    blocks = tilting_blocks(R)
    floor = detected_dimension = detected_blocks = 0
    for index, (on, off, a, b, shift, mu, dimension) in enumerate(blocks):
        detected, _ = block_contains_central_monomial(on, off, a, b, shift, R)
        if detected:
            floor += mu
            detected_dimension += mu * dimension
            detected_blocks += 1
        if progress and index % 100 == 0:
            print(json.dumps({"R": R, "block": index, "count": len(blocks),
                              "rank_floor": floor,
                              "elapsed_seconds": round(time.monotonic()-started, 2)}),
                  flush=True)
    return {"R": R, "label_count": D, "block_count": len(blocks),
            "detected_blocks": detected_blocks,
            "rank_floor": floor,
            "detected_total_dimension": detected_dimension,
            "rho_floor": log2(floor) / (2 * log2(D)) if floor else 0,
            "strict_seed_threshold": D ** sqrt(2),
            "strict_seed_obstruction_certified": floor ** 408 > D ** 577,
            "rational_upper_bound_on_sqrt2": [577, 408],
            "elapsed_seconds": round(time.monotonic()-started, 3)}


def asymptotic_type_witness(R, S=16):
    """Construct the single high-multiplicity detecting block in the proof.

    This verifies structural conditions at huge energy without enumerating
    the shell.  The mathematical proof supplies the polynomial basis.
    """
    K = R - (R + 1).bit_length() - 4
    assert K > S
    a, h = 1 << S, 1 << R
    types, multiplicity_lower = [], 1
    for s in range(R + 1):
        m = s + 1
        if s < S or s > K:
            n, p, lower = 0, m, 1
        elif s == S:
            n, p, lower = 2, 7, 1
            assert any(nn == n and pp == p and mu > 0
                       for nn, pp, mu in local_bicolored_tiltings(m))
        else:
            k = 2 * (m // 4)
            length = m - k - 7
            assert length >= 0
            n, mu = max(tensor_tiltings(k), key=lambda pair: pair[1])
            r, nu = max(tensor_tiltings(length), key=lambda pair: pair[1])
            p = r + 7
            lower = comb(m, k) * mu * nu
            assert mu * (k + 1) ** 3 >= 2 ** k
            assert nu * (length + 1) ** 3 >= 2 ** length
            assert lower * 128 * (m + 1) ** 9 >= 4 ** m
        assert n % 2 == 0 and n + p <= m and (m - n - p) % 2 == 0
        types.append((n, p))
        multiplicity_lower *= lower
    Delta = sum(((s + 1 - n - p) // 2) << s
                for s, (n, p) in enumerate(types))
    A = sum((n // 2) << s for s, (n, p) in enumerate(types))
    B = h - Delta - A
    assert all(value % a == 0 for value in (Delta, A, B))
    assert Delta + A <= h // 16 and Delta <= B <= h - a
    off_counts = defaultdict(int)
    for s, (n, p) in enumerate(types):
        for weight in tilting_weights(p):
            off_counts[s + weight.bit_length() - 1] += 1
    assert all(off_counts[s] >= 3 for s in range(S, R))
    ambient_bits = (R + 1) * (R + 2) // 2
    return {"R": R, "S": S, "K": K, "energy": h,
            "Delta": Delta, "A": A, "B": B,
            "multiplicity_lower_hex": hex(multiplicity_lower),
            "multiplicity_binary_log": log2(multiplicity_lower),
            "rho_lower_using_ambient_label_count":
                log2(multiplicity_lower) / (2 * ambient_bits),
            "types": types,
            "three_off_binary_chains_verified": True}


def calculate_tilting(R, progress=False, primitive=False):
    started = time.monotonic()
    h = 1 << R
    counts = tuple(range(1, R + 2))
    D = layer_count(counts, h)
    degree = 2 * h if primitive else 2 * sum(k << s for s, k in enumerate(counts))
    ranks = [0] * (degree + 1)
    blocks = tilting_blocks(R)
    for index, (on, off, a, b, shift, mu, dimension) in enumerate(blocks):
        actual_dimension, block_ranks = reduced_block_ranks(
            on, off, a, b, shift, R if primitive else None)
        assert dimension == actual_dimension
        for k, rank in block_ranks.items():
            ranks[k] += mu * rank
        if progress and index % 10 == 0:
            print(json.dumps({"R": R, "block": index, "count": len(blocks),
                              "dimension": dimension,
                              "elapsed_seconds": round(time.monotonic()-started, 2)}),
                  flush=True)
    products = [a * b for a, b in zip(ranks, reversed(ranks))]
    max_product = max(products)
    S = sum(sqrt(q) for q in products)
    return {
        "R": R, "energy": h, "label_count": D, "bit_count": sum(counts),
        "degree": degree, "method": "actual modular tilting direct sums",
        "primitive_normalization": primitive,
        "tilting_block_count": len(blocks),
        "max_block_dimension": max(block[-1] for block in blocks),
        "coefficient_ranks": ranks,
        "max_rank_product": max_product,
        "rho_raw_max": log2(max_product) / (4 * log2(D)),
        "raw_seed_S": S, "strict_seed_threshold": D ** sqrt(2),
        "rho_raw_seed": log2(S) / (2 * log2(D)),
        "elapsed_seconds": round(time.monotonic()-started, 3),
    }


@lru_cache(None)
def entry_polynomial(profile, distances_on, distances_off):
    """Rows (x,y), columns (u,v), with g=x xor y=u xor v.

    On g, the 2x2 factor is z^w(J+z^w I).  Outside g, it is
    (1+z^w)J+z^(2w)I.  The initial shift accounts for z^sum(w*g).
    """
    shift = sum(k << s for s, k in enumerate(profile))
    result = 1 << shift
    for s, (k, dg, do) in enumerate(zip(profile, distances_on, distances_off)):
        outside = s + 1 - k
        a, b = k - dg + do, outside - do
        bit = 0
        while a or b:
            weight = 1 << (s + bit)
            if a & 1:
                result ^= result << weight
            if b & 1:
                result ^= (result << weight) ^ (result << (2 * weight))
            a >>= 1
            b >>= 1
            bit += 1
    return result


def block_data(R, profile):
    energy = 1 << R
    offset = 0
    on, off = [], []
    for s, k in enumerate(profile):
        on.append(tuple(range(offset, offset + k)))
        off.append(tuple(range(offset + k, offset + s + 1)))
        offset += s + 1
    total_on = sum(k << s for s, k in enumerate(profile))
    left = layer_labels(tuple(on), total_on // 2)
    right = layer_labels(tuple(off), energy - total_on // 2)
    labels = [x | y for x in left for y in right]
    on_masks = tuple(sum(1 << i for i in group) for group in on)
    off_masks = tuple(sum(1 << i for i in group) for group in off)
    return labels, on_masks, off_masks


def profile_ranks(R, profile, primitive=False):
    labels, on, off = block_data(R, profile)
    dimension = len(labels)
    matrices = defaultdict(lambda: [0] * dimension)
    for row, x in enumerate(labels):
        for col in range(row, dimension):
            difference = x ^ labels[col]
            dg = tuple((difference & mask).bit_count() for mask in on)
            do = tuple((difference & mask).bit_count() for mask in off)
            polynomial = entry_polynomial(profile, dg, do)
            if primitive:
                polynomial = primitive_polynomial(polynomial, R)
            while polynomial:
                bit = polynomial & -polynomial
                degree = bit.bit_length() - 1
                matrices[degree][row] |= 1 << col
                matrices[degree][col] |= 1 << row
                polynomial ^= bit
    return dimension, {degree: rank_two(rows) for degree, rows in matrices.items()}


def calculate(R, progress=False, primitive=False):
    started = time.monotonic()
    counts = tuple(range(1, R + 2))
    energy = 1 << R
    D = layer_count(counts, energy)
    degree = 2 * energy if primitive else 2 * sum(k << s for s, k in enumerate(counts))
    ranks = [0] * (degree + 1)
    profiles = []
    weighted_dimension = 0
    for profile in product(*(range(c + 1) for c in counts)):
        total_on = sum(k << s for s, k in enumerate(profile))
        if total_on & 1 or total_on > 2 * energy:
            continue
        dimension = layer_count(profile, total_on // 2) * layer_count(
            tuple(c - k for c, k in zip(counts, profile)), energy - total_on // 2
        )
        if not dimension:
            continue
        multiplicity = prod(comb(c, k) for c, k in zip(counts, profile))
        actual_dimension, block_ranks = profile_ranks(R, profile, primitive)
        assert actual_dimension == dimension
        weighted_dimension += dimension * multiplicity
        for k, rank in block_ranks.items():
            ranks[k] += multiplicity * rank
        profiles.append({
            "xor_counts": profile,
            "dimension": dimension,
            "multiplicity": multiplicity,
            "coefficient_ranks": block_ranks,
        })
        if progress:
            print(json.dumps({"R": R, "profile": profile, "dimension": dimension,
                              "elapsed_seconds": round(time.monotonic()-started, 2)}),
                  flush=True)
    assert weighted_dimension == D * D
    products = [a * b for a, b in zip(ranks, reversed(ranks))]
    max_product = max(products)
    S = sum(sqrt(q) for q in products)
    return {
        "R": R,
        "energy": energy,
        "label_count": D,
        "bit_count": sum(counts),
        "degree": degree,
        "primitive_normalization": primitive,
        "xor_profile_count": len(profiles),
        "coefficient_ranks": ranks,
        "max_rank_product": max_product,
        "rho_raw_max": log2(max_product) / (4 * log2(D)),
        "raw_seed_S": S,
        "strict_seed_threshold": D ** sqrt(2),
        "rho_raw_seed": log2(S) / (2 * log2(D)),
        "elapsed_seconds": round(time.monotonic()-started, 3),
        "profiles": profiles,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--R", type=int, default=4)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--output")
    parser.add_argument("--tilting", action="store_true")
    parser.add_argument("--block-stats", action="store_true")
    parser.add_argument("--primitive", action="store_true")
    parser.add_argument("--central-bound", action="store_true")
    parser.add_argument("--asymptotic-witness", action="store_true")
    args = parser.parse_args()
    if args.asymptotic_witness:
        result = asymptotic_type_witness(args.R)
    elif args.central_bound:
        result = central_detector_bound(args.R, args.progress)
    elif args.block_stats:
        blocks = tilting_blocks(args.R)
        result = {"R": args.R, "count": len(blocks),
                  "max_dimension": max(block[-1] for block in blocks),
                  "sum_dimension": sum(block[-1] for block in blocks),
                  "sum_dimension_squared": sum(block[-1] ** 2 for block in blocks)}
    else:
        result = (calculate_tilting if args.tilting else calculate)(
            args.R, args.progress, args.primitive)
    if args.output:
        with open(args.output, "w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ("profiles", "coefficient_ranks", "types", "multiplicity_lower_hex")},
                     indent=2))


if __name__ == "__main__":
    main()
