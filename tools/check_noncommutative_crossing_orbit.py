#!/usr/bin/env python3
"""Exact orbit-kernel crossing pieces from noncommuting Heisenberg words.

This is a construction and bounded rank experiment, not a stronger factoring
theorem.  All arithmetic in the crossing tensors is over F_2.  Group coordinates
are (a,b,c) in F_(2^k)^3 with (a,b,c)(A,B,C)=(a+A,b+B,c+C+aB).
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import itertools
import json
import math
from collections import Counter
from pathlib import Path
import time


IRREDUCIBLE = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011}


class Heisenberg:
    def __init__(self, degree: int):
        self.k = degree
        self.q = 1 << degree
        self.mask = self.q - 1
        self.d = self.q ** 3
        self.width = 3 * degree
        self.modulus = IRREDUCIBLE[degree]
        self.generators = tuple(
            self.encode(1 << i, 0, 0) for i in range(degree)
        ) + tuple(self.encode(0, 1 << i, 0) for i in range(degree)) + tuple(
            self.encode(0, 0, 1 << i) for i in range(degree)
        )
        self.words = tuple(self.word_product(bits, range(self.width))
                           for bits in range(self.d))
        assert len(set(self.words)) == self.d
        self.coordinates = [0] * self.d
        for bits, element in enumerate(self.words):
            self.coordinates[element] = bits
        self.inverse = tuple(self.inv(x) for x in range(self.d))

    def field_mul(self, a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return result

    def encode(self, a: int, b: int, c: int) -> int:
        return a | (b << self.k) | (c << (2 * self.k))

    def decode(self, element: int) -> tuple[int, int, int]:
        return (element & self.mask, (element >> self.k) & self.mask,
                element >> (2 * self.k))

    def mul(self, left: int, right: int) -> int:
        a, b, c = self.decode(left)
        A, B, C = self.decode(right)
        return self.encode(a ^ A, b ^ B, c ^ C ^ self.field_mul(a, B))

    def inv(self, element: int) -> int:
        a, b, c = self.decode(element)
        return self.encode(a, b, c ^ self.field_mul(a, b))

    def word_product(self, bits: int, order: range | tuple[int, ...]) -> int:
        result = 0
        for j in order:
            if (bits >> j) & 1:
                result = self.mul(result, self.generators[j])
        return result


def polynomial_mul(a: int, b: int) -> int:
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def group_word_polynomial(group: Heisenberg, element: int,
                          order: tuple[int, ...],
                          weights: tuple[int, ...]) -> dict[int, int]:
    """Coefficients of product_j (g_j + z^w_j (1+g_j))^s_j."""
    bits = group.coordinates[element]
    state = {0: 1}
    for j in order:
        if not ((bits >> j) & 1):
            continue
        nxt: dict[int, int] = {}
        shift = weights[j]
        for value, polynomial in state.items():
            # z^w at the identity, (1+z^w) at g_j.
            nxt[value] = nxt.get(value, 0) ^ (polynomial << shift)
            other = group.mul(value, group.generators[j])
            nxt[other] = nxt.get(other, 0) ^ polynomial ^ (polynomial << shift)
        state = {value: polynomial for value, polynomial in nxt.items() if polynomial}
    return state


def group_word_polynomial_independent(group: Heisenberg, element: int,
                                     order: tuple[int, ...],
                                     weights: tuple[int, ...]) -> dict[int, int]:
    """Independent subset formula, including reordered noncommuting factors."""
    bits = group.coordinates[element]
    active = tuple(j for j in order if (bits >> j) & 1)
    answer: dict[int, int] = {}
    for subset in range(1 << len(active)):
        product = 0
        polynomial = 1
        for place, j in enumerate(active):
            if (subset >> place) & 1:
                product = group.mul(product, group.generators[j])
                polynomial = polynomial_mul(polynomial, 1 | (1 << weights[j]))
            else:
                polynomial <<= weights[j]
        answer[product] = answer.get(product, 0) ^ polynomial
    return {value: polynomial for value, polynomial in answer.items() if polynomial}


def binary_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        while row:
            lead = row.bit_length() - 1
            pivot = pivots.get(lead)
            if pivot is None:
                pivots[lead] = row
                break
            row ^= pivot
    return len(pivots)


def orbit_polynomial(group: Heisenberg, x: int, y: int, u: int, v: int,
                     kernels: list[dict[int, int]]) -> int:
    shift = group.mul(group.inverse[y], u)
    if group.mul(x, shift) != v:
        return 0
    quotient = group.mul(x, group.inverse[y])
    return kernels[quotient].get(shift, 0)


def profile(group: Heisenberg, order: tuple[int, ...], weights: tuple[int, ...],
            direct_row_check: bool = False,
            control_correction: tuple[int, ...] | None = None) -> dict:
    d = group.d
    degree = sum(weights)
    kernels = []
    for S in range(d):
        if control_correction is None:
            control = S
        else:
            a, b, c = group.decode(S)
            control = group.encode(a, b, c ^ control_correction[a * group.q + b])
        kernels.append(group_word_polynomial(group, control, order, weights))
    assert all(sum(poly.bit_count() & 1 for a, poly in kernel.items() if a == 0) == 1
               and all((poly.bit_count() & 1) == (a == 0) for a, poly in kernel.items())
               for kernel in kernels)
    # At fixed quotient S, the row block is a group convolution matrix indexed
    # by y and u: K(S,y^{-1}u).  Construct it independently of the dense tensor.
    row_ranks = [0] * (degree + 1)
    for S, kernel in enumerate(kernels):
        rows = [[0] * d for _ in range(degree + 1)]
        for y in range(d):
            for a, polynomial in kernel.items():
                u = group.mul(y, a)
                while polynomial:
                    bit = polynomial & -polynomial
                    rows[bit.bit_length() - 1][y] ^= 1 << u
                    polynomial ^= bit
        for j, matrix in enumerate(rows):
            row_ranks[j] += binary_rank(matrix)
    column_rows = [[0] * (d * d) for _ in range(degree + 1)]
    dense_row = [[0] * (d * d) for _ in range(degree + 1)] if direct_row_check else None
    for S, kernel in enumerate(kernels):
        for y in range(d):
            x = group.mul(S, y)
            for a, polynomial in kernel.items():
                u, v = group.mul(y, a), group.mul(x, a)
                assert orbit_polynomial(group, x, y, u, v, kernels) == polynomial
                while polynomial:
                    bit = polynomial & -polynomial
                    j = bit.bit_length() - 1
                    column_rows[j][x * d + u] ^= 1 << (y * d + v)
                    if dense_row is not None:
                        dense_row[j][x * d + y] ^= 1 << (u * d + v)
                    polynomial ^= bit
    column_ranks = [binary_rank(rows) for rows in column_rows]
    if dense_row is not None:
        assert row_ranks == [binary_rank(rows) for rows in dense_row]
    # An explicit sum of all column-cut matrices is the crossing permutation.
    for x in range(d):
        for u in range(d):
            value = 0
            for rows in column_rows:
                value ^= rows[x * d + u]
            assert value == 1 << (u * d + x)
    products = [r * c for r, c in zip(row_ranks, column_ranks)]
    widths = [((r - 1).bit_length() + (c - 1).bit_length()) if r and c else 0
              for r, c in zip(row_ranks, column_ranks)]
    if control_correction is not None:
        assert order == tuple(range(group.width))
        assert column_ranks[0] == central_correction_rank(group, control_correction)
    else:
        assert column_ranks[0] == compressed_constant_cut(group, order)[0]
    return {"field_degree": group.k, "labels": d, "binary_width": group.width,
            "order": list(order), "weights": list(weights), "piece_degree": degree,
            "control_correction": list(control_correction) if control_correction is not None else None,
            "row_ranks": row_ranks, "column_ranks": column_ranks,
            "rank_products": products, "max_rank_product": max(products),
            "total_widths": widths, "maximum_total_width": max(widths),
            "nonzero_piece_count": sum(r != 0 for r in row_ranks)}


def full_coefficient_span(group: Heisenberg,
                          control_correction: tuple[int, ...] | None = None) -> dict:
    """Every recombination of the unweighted coefficient family, exactly."""
    d, width = group.d, group.width
    order, weights = tuple(range(width)), (1,) * width
    kernels = []
    for S in range(d):
        if control_correction is None:
            control = S
        else:
            a, b, c = group.decode(S)
            control = group.encode(a, b, c ^ control_correction[a * group.q + b])
        kernels.append(group_word_polynomial(group, control, order, weights))
    assert binary_rank([poly for kernel in kernels for poly in kernel.values()]) == width + 1
    row_blocks = [[[0] * d for _ in range(d)] for _ in range(width + 1)]
    columns = [[0] * (d * d) for _ in range(width + 1)]
    for S, kernel in enumerate(kernels):
        for y in range(d):
            x = group.mul(S, y)
            for a, polynomial in kernel.items():
                u, v = group.mul(y, a), group.mul(x, a)
                while polynomial:
                    bit = polynomial & -polynomial
                    j = bit.bit_length() - 1
                    row_blocks[j][S][y] ^= 1 << u
                    columns[j][x * d + u] ^= 1 << (y * d + v)
                    polynomial ^= bit
    records = []
    for functional in range(1, 1 << (width + 1)):
        selected = [j for j in range(width + 1) if (functional >> j) & 1]
        r = 0
        for S in range(d):
            rows = [0] * d
            for j in selected:
                rows = [left ^ right for left, right in zip(rows, row_blocks[j][S])]
            r += binary_rank(rows)
        cols = [0] * (d * d)
        for j in selected:
            cols = [left ^ right for left, right in zip(cols, columns[j])]
        c = binary_rank(cols)
        if functional & 1:
            assert r == d * d  # Each group-algebra row block has augmentation one.
        records.append({"functional_mask": functional, "row_rank": r,
                        "column_rank": c, "rank_product": r * c})
    # An exact bottleneck span certificate.  Track the actual selected masks,
    # not just the dimension of the span, until their XOR is the crossing.
    target = (1 << (width + 1)) - 1
    pivots: dict[int, tuple[int, int]] = {}
    threshold = None
    witness = 0
    by_cost = sorted(records, key=lambda item: (item["rank_product"], item["functional_mask"]))
    for index, record in enumerate(by_cost):
        value, provenance = record["functional_mask"], 1 << index
        while value:
            lead = value.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = (value, provenance)
                break
            other, source = pivots[lead]
            value ^= other
            provenance ^= source
        value, provenance = target, 0
        while value and value.bit_length() - 1 in pivots:
            other, source = pivots[value.bit_length() - 1]
            value ^= other
            provenance ^= source
        if value == 0:
            threshold, witness = record["rank_product"], provenance
            break
    selected = [by_cost[j] for j in range(len(by_cost)) if (witness >> j) & 1]
    assert threshold is not None
    total_mask = 0
    for record in selected:
        total_mask ^= record["functional_mask"]
    assert total_mask == target
    lower_masks = [record["functional_mask"] for record in records
                   if record["rank_product"] < threshold]
    assert binary_rank(lower_masks + [target]) == binary_rank(lower_masks) + 1
    unit_slice = [record for record in records if record["functional_mask"] & 1]
    scale = 10 ** 12
    def cost_bounds(upper: bool) -> tuple[int, list[int]]:
        distances = [None] * (1 << (width + 1))
        parents = [None] * len(distances)
        distances[0] = 0
        queue = [(0, 0)]
        costs = []
        for record in records:
            squared = record["rank_product"] * scale * scale
            lower = math.isqrt(squared)
            costs.append(lower + int(upper and lower * lower != squared))
        while queue:
            distance, state = heapq.heappop(queue)
            if distance != distances[state]:
                continue
            if state == target:
                path = []
                while state:
                    prior, mask = parents[state]
                    path.append(mask)
                    state = prior
                return distance, path
            for record, cost in zip(records, costs):
                mask = record["functional_mask"]
                nxt, proposed = state ^ mask, distance + cost
                if distances[nxt] is None or proposed < distances[nxt]:
                    distances[nxt] = proposed
                    parents[nxt] = (state, mask)
                    heapq.heappush(queue, (proposed, nxt))
        raise AssertionError("crossing was not spanned")
    lower_sum, _ = cost_bounds(False)
    upper_sum, sum_witness = cost_bounds(True)
    assert lower_sum <= upper_sum
    return {"labels": d, "binary_width": width,
            "control_correction": list(control_correction) if control_correction is not None else None,
            "coefficient_dimension": width + 1, "all_recombinations": records,
            "minimum_bottleneck_rank_product": threshold, "witness": selected,
            "unit_slice_minimum_column_rank": min(record["column_rank"] for record in unit_slice),
            "minimum_rank_sum_interval": {"lower_numerator": lower_sum,
                                           "upper_numerator": upper_sum,
                                           "denominator": scale},
            "rank_sum_upper_witness_masks": sum_witness}


def compressed_constant_cut(group: Heisenberg, order: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    """A q^3 square matrix with exactly the constant coefficient's column rank.

    The uncompressed cut has q^6 rows/columns.  All of its nonzero rows have
    equal horizontal coordinates x and u, and identical central translates.
    Removing the zero rows and the duplicated central translates leaves this
    exact compressed matrix, not an associated-graded estimate.
    """
    q = group.q
    correction = tuple(
        group.decode(group.word_product(group.coordinates[group.encode(a, b, 0)], order))[2]
        for a in range(q) for b in range(q)
    )
    return central_correction_rank(group, correction), correction


def central_correction_matrix(group: Heisenberg, correction: tuple[int, ...]) -> list[int]:
    q = group.q
    rows = [0] * (q ** 3)
    for a, b, A, B in itertools.product(range(q), repeat=4):
        f = correction[(a ^ A) * q + (b ^ B)]
        delta = group.field_mul(a, B) ^ group.field_mul(A, b) ^ f
        epsilon = group.field_mul(a, b) ^ group.field_mul(A, B) ^ f
        rows[(a * q + b) * q + delta] ^= 1 << ((A * q + B) * q + epsilon)
    return rows


def central_correction_rank(group: Heisenberg, correction: tuple[int, ...]) -> int:
    return binary_rank(central_correction_matrix(group, correction))


def correction_symmetry_checks() -> dict:
    """Check the actual permutation/transpose identities, not just equal ranks."""
    group = Heisenberg(2)
    q, d = group.q, group.d
    state = 1729
    checked = 0
    for trial in range(8):
        correction = []
        for _ in range(q * q):
            state = (1103515245 * state + 12345) & 0x7fffffff
            correction.append((state >> 11) & (q - 1))
        correction = tuple(correction)
        rows = central_correction_matrix(group, correction)
        complementary = tuple(correction[a * q + b] ^ group.field_mul(a, b)
                              for a in range(q) for b in range(q))
        transposed = central_correction_matrix(group, complementary)
        for i in range(d):
            for j in range(d):
                assert ((rows[i] >> j) & 1) == ((transposed[j] >> i) & 1)
                checked += 1
        t_a, t_b, c = trial % q, (trial // 2) % q, (trial // 3) % q
        translated = tuple(correction[a * q + b] ^ group.field_mul(t_a, b)
                           ^ group.field_mul(t_b, a) ^ c
                           for a in range(q) for b in range(q))
        changed = central_correction_matrix(group, translated)
        permutation = []
        for a, b, delta in itertools.product(range(q), repeat=3):
            permutation.append(((a ^ t_a) * q + (b ^ t_b)) * q + (delta ^ c))
        for i in range(d):
            for j in range(d):
                assert ((changed[i] >> j) & 1) == ((rows[permutation[i]] >> permutation[j]) & 1)
                checked += 1
    return {"field_degree": 2, "arbitrary_corrections": 8,
            "matrix_entry_equalities": checked,
            "identities": ["f+ab gives transpose", "F_q-affine additions give row/column permutations"]}


def constant_order_screen(max_degree: int) -> list[dict]:
    output = []
    for k in range(1, max_degree + 1):
        group = Heisenberg(k)
        seen: dict[tuple[int, ...], tuple[int, tuple[int, ...]]] = {}
        histogram: Counter[int] = Counter()
        for noncentral_order in itertools.permutations(range(2 * k)):
            order = noncentral_order + tuple(range(2 * k, 3 * k))
            rank, correction = compressed_constant_cut(group, order)
            if correction not in seen:
                seen[correction] = (rank, order)
                histogram[rank] += 1
        canonical, _ = compressed_constant_cut(group, tuple(range(3 * k)))
        assert canonical == (group.q - 1) * 3 ** k + group.q
        minimum = min(rank for rank, order in seen.values())
        output.append({"field_degree": k, "labels": group.d,
                       "orders_examined": math.factorial(2 * k),
                       "distinct_central_corrections": len(seen),
                       "rank_histogram": {str(rank): count for rank, count in sorted(histogram.items())},
                       "minimum_column_rank": minimum,
                       "minimizing_orders": [list(order) for rank, order in seen.values() if rank == minimum],
                       "canonical_column_rank": canonical})
    return output


def checks(degree_two: bool = False, span_degree_two: bool = False,
           order_degree: int = 2) -> dict:
    start = time.monotonic()
    tested_groups = []
    for k in (1, 2):
        group = Heisenberg(k)
        # All group axioms by the actual implementation at k=1.  At k=2,
        # exact bilinearity of multiplication proves the associative formula;
        # enumerate field distributivity and all inverse checks.
        for a, b, c in itertools.product(range(group.q), repeat=3):
            assert group.field_mul(a, b ^ c) == group.field_mul(a, b) ^ group.field_mul(a, c)
            assert group.field_mul(a, group.field_mul(b, c)) == group.field_mul(group.field_mul(a, b), c)
        for x in range(group.d):
            assert group.mul(x, group.inverse[x]) == 0 == group.mul(group.inverse[x], x)
        if k == 1:
            for a, b, c in itertools.product(range(group.d), repeat=3):
                assert group.mul(a, group.mul(b, c)) == group.mul(group.mul(a, b), c)
        x, y = group.generators[0], group.generators[k]
        assert group.mul(x, y) != group.mul(y, x)
        order = tuple(range(group.width))
        for S in range(group.d):
            assert group_word_polynomial(group, S, order, (1,) * group.width) == \
                group_word_polynomial_independent(group, S, order, (1,) * group.width)
        tested_groups.append({"field_degree": k, "order": group.d,
                              "generator_count": group.width})
    small = Heisenberg(1)
    profiles = [profile(small, order, (1, 1, 1), True)
                for order in itertools.permutations(range(3))]
    profiles += [profile(small, (0, 1, 2), weights, True)
                 for weights in ((1, 1, 2), (1, 2, 1), (2, 1, 1), (1, 2, 3))]
    # Every function on F_2^2, not just a bilinear shear: compare the exact
    # compressed constant matrix with the independently materialized tensor.
    arbitrary_shears = [profile(small, (0, 1, 2), (1, 1, 1), True, correction)
                        for correction in itertools.product(range(2), repeat=4)]
    if degree_two:
        larger = Heisenberg(2)
        profiles.append(profile(larger, tuple(range(6)), (1,) * 6))
        profiles.append(profile(larger, (2, 3, 0, 1, 4, 5), (1,) * 6))
        # Cubic-in-bits input control, outside the bilinear correction model.
        correction = tuple(((a & 1) * (b & 1) * ((a >> 1) & 1))
                           for a in range(4) for b in range(4))
        arbitrary_shears.append(profile(larger, tuple(range(6)), (1,) * 6,
                                       control_correction=correction))
        # The first rank-28 normalized quadratic discovered by the independent
        # full native/Python screen.  All pieces, not just the constant one,
        # are constructed here to prevent a local improvement becoming a
        # stronger crossing claim by omission.
        winner = tuple((a >> 1) ^ ((a & 1) * (a >> 1))
                       ^ ((b & 1) * (b >> 1))
                       ^ (2 * (((a & 1) * (b & 1)) ^ ((a & 1) * (b >> 1))
                               ^ ((a >> 1) * (b & 1))))
                       for a in range(4) for b in range(4))
        improved = profile(larger, tuple(range(6)), (1,) * 6,
                           control_correction=winner)
        assert improved["column_ranks"][0] == 28
        assert improved["max_rank_product"] > profiles[-1]["max_rank_product"]
        arbitrary_shears.append(improved)
    spans = [full_coefficient_span(small)]
    if span_degree_two:
        spans.append(full_coefficient_span(Heisenberg(2)))
        winner_table = (0, 0, 0, 1, 0, 2, 2, 1, 1, 3, 1, 2, 0, 0, 2, 3)
        improved_span = full_coefficient_span(Heisenberg(2), winner_table)
        assert improved_span["minimum_bottleneck_rank_product"] == 1128444
        assert improved_span["unit_slice_minimum_column_rank"] == 28
        spans.append(improved_span)
    return {"status": "exact finite construction checks; no stronger exponent",
            "coefficient_field": "F_2", "groups": tested_groups,
            "profiles": profiles, "full_coefficient_spans": spans,
            "arbitrary_shear_profiles": arbitrary_shears,
            "central_correction_symmetries": correction_symmetry_checks(),
            "constant_coefficient_order_screen": constant_order_screen(order_degree),
            "elapsed_seconds": round(time.monotonic() - start, 4)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree-two", action="store_true",
                        help="also rank the 64-label, six-generator matrices")
    parser.add_argument("--span-degree-two", action="store_true",
                        help="also rank all 127 recombinations at 64 labels")
    parser.add_argument("--order-degree", type=int, choices=(1, 2, 3), default=2,
                        help="screen all noncentral orders through this field degree")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks(args.degree_two, args.span_degree_two, args.order_degree)
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
