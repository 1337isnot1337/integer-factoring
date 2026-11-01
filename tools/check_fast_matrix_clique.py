#!/usr/bin/env python3
"""Exact checks for the three-way clique and grouped-table contractions."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path
import random


@dataclass(frozen=True)
class Ring:
    name: str
    order: int
    dual: bool = False

    def add(self, a: int, b: int) -> int:
        return a ^ b if self.dual else (a + b) % self.order

    def sub(self, a: int, b: int) -> int:
        return a ^ b if self.dual else (a - b) % self.order

    def mul(self, a: int, b: int) -> int:
        if self.dual:
            return ((a & b) & 1) | ((((a & 1) * (b >> 1)) ^
                                      ((b & 1) * (a >> 1))) << 1)
        return (a * b) % self.order


RINGS = (Ring("F2", 2), Ring("F3", 3), Ring("Z/6Z", 6),
         Ring("F2[eps]/eps^2", 4, True))


def add_matrix(a, b, ring, subtract=False):
    op = ring.sub if subtract else ring.add
    return [[op(x, y) for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def strassen(a, b, ring, stats=None):
    """Literal seven-product recursion, with its scalar counts."""
    n = len(a)
    if n == 1:
        if stats is not None:
            stats["multiplications"] += 1
        return [[ring.mul(a[0][0], b[0][0])]]
    assert n > 0 and n & (n - 1) == 0
    h = n // 2

    def blocks(m):
        return ([r[:h] for r in m[:h]], [r[h:] for r in m[:h]],
                [r[:h] for r in m[h:]], [r[h:] for r in m[h:]])

    aa, ab, ac, ad = blocks(a)
    ba, bb, bc, bd = blocks(b)

    def plus(x, y):
        if stats is not None:
            stats["additions"] += h * h
        return add_matrix(x, y, ring)

    def minus(x, y):
        if stats is not None:
            stats["additions"] += h * h
        return add_matrix(x, y, ring, True)

    def mm(x, y):
        return strassen(x, y, ring, stats)

    p1 = mm(plus(aa, ad), plus(ba, bd))
    p2 = mm(plus(ac, ad), ba)
    p3 = mm(aa, minus(bb, bd))
    p4 = mm(ad, minus(bc, ba))
    p5 = mm(plus(aa, ab), bd)
    p6 = mm(minus(ac, aa), plus(ba, bb))
    p7 = mm(minus(ab, ad), plus(bc, bd))
    ca = plus(minus(plus(p1, p4), p5), p7)
    cb = plus(p3, p5)
    cc = plus(p2, p4)
    cd = plus(plus(minus(p1, p2), p3), p6)
    return [ca[i] + cb[i] for i in range(h)] + [cc[i] + cd[i] for i in range(h)]


def naive_product(a, b, ring):
    out = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            for j, y in enumerate(b[k]):
                out[i][j] = ring.add(out[i][j], ring.mul(x, y))
    return out


def trace_product(a, b, c, ring):
    p = strassen(a, b, ring)
    out = 0
    for i in range(len(p)):
        for j in range(len(p)):
            out = ring.add(out, ring.mul(p[i][j], c[j][i]))
    return out


def blocked_rectangular_trace(a, b, c, short, long, ring):
    """Stream the long dimension in square blocks, retaining the short A."""
    assert long % short == 0
    aa = [row[:short] for row in a[:short]]
    result = 0
    for first in range(0, long, short):
        bb = [row[first:first + short] for row in b[:short]]
        cc = [row[:short] for row in c[first:first + short]]
        result = ring.add(result, trace_product(aa, bb, cc, ring))
    return result


def recursive_conditioning(k, width, tables, ring):
    """One packed edge string; interval recursion with sequential children."""
    assert k > 0 and k & (k - 1) == 0
    labels = [0]
    queries = 0
    mask = (1 << width) - 1

    def edge_slot(i, j):
        if i > j:
            i, j = j, i
        return i * (2 * k - i - 1) // 2 + j - i - 1

    def evaluate(first, size):
        nonlocal queries
        if size == 1:
            index = position = 0
            for neighbor in range(k):
                if neighbor == first:
                    continue
                value = (labels[0] >> (width * edge_slot(first, neighbor))) & mask
                index |= value << (width * position)
                position += 1
            queries += 1
            return tables[first][index]
        half = size // 2
        result = 0
        for word in range(1 << (width * half * half)):
            for i in range(half):
                for j in range(half):
                    value = (word >> (width * (i * half + j))) & mask
                    offset = width * edge_slot(first + i, first + half + j)
                    labels[0] = (labels[0] & ~(mask << offset)) | (value << offset)
            left = evaluate(first, half)
            right = evaluate(first + half, half)
            result = ring.add(result, ring.mul(left, right))
        return result

    value = evaluate(0, k)
    assert queries == k * (1 << (((k * k - 1) // 3) * width))
    return value, queries


def split_three(k):
    q, r = divmod(k, 3)
    sizes = [q + (i < r) for i in range(3)]
    out, start = [], 0
    for size in sizes:
        out.append(tuple(range(start, start + size)))
        start += size
    return out


def clique_checks(k, width, ring, rng, group_sizes=None):
    alphabet = 1 << width
    edges = list(combinations(range(k), 2))
    table = [[rng.randrange(ring.order) for _ in range(alphabet ** (k - 1))]
             for _ in range(k)]
    adjacent = [[e for e in edges if v in e] for v in range(k)]

    def vertex(v, lab):
        ix, power = 0, 1
        for e in adjacent[v]:
            ix += power * lab[e]
            power *= alphabet
        return table[v][ix]

    direct = 0
    for word in product(range(alphabet), repeat=len(edges)):
        lab = dict(zip(edges, word))
        term = 1
        for v in range(k):
            term = ring.mul(term, vertex(v, lab))
        direct = ring.add(direct, term)

    if group_sizes is None:
        groups = split_three(k)
    else:
        assert sum(group_sizes) == k and len(group_sizes) == 3
        groups, first = [], 0
        for size in group_sizes:
            groups.append(tuple(range(first, first + size)))
            first += size
    blocks = {}
    for i, j in combinations(range(3), 2):
        blocks[i, j] = [e for e in edges if
                        (e[0] in groups[i] and e[1] in groups[j]) or
                        (e[1] in groups[i] and e[0] in groups[j])]
    dims = {key: alphabet ** len(es) for key, es in blocks.items()}
    d = max(dims.values())
    matrices = []
    query_count = 0
    for gi, row_key, col_key in ((0, (0, 2), (0, 1)),
                                 (1, (0, 1), (1, 2)),
                                 (2, (1, 2), (0, 2))):
        internal = [e for e in edges if e[0] in groups[gi] and e[1] in groups[gi]]
        matrix = [[0] * d for _ in range(d)]
        rows = list(product(range(alphabet), repeat=len(blocks[row_key])))
        cols = list(product(range(alphabet), repeat=len(blocks[col_key])))
        for ri, rw in enumerate(rows):
            for ci, cw in enumerate(cols):
                lab = dict(zip(blocks[row_key], rw))
                lab.update(zip(blocks[col_key], cw))
                for inside in product(range(alphabet), repeat=len(internal)):
                    lab.update(zip(internal, inside))
                    term = 1
                    for v in groups[gi]:
                        query_count += 1
                        term = ring.mul(term, vertex(v, lab))
                    matrix[ri][ci] = ring.add(matrix[ri][ci], term)
        matrices.append(matrix)
    fast = trace_product(*matrices, ring)
    assert fast == direct
    rectangular = len(groups[0]) <= len(groups[1]) == len(groups[2])
    if rectangular:
        assert blocked_rectangular_trace(*matrices, dims[0, 1], dims[1, 2], ring) == direct
    conditioned_queries = None
    if k & (k - 1) == 0:
        conditioned, conditioned_queries = recursive_conditioning(k, width, table, ring)
        assert conditioned == direct
    expected_queries = sum(len(g) * alphabet **
                           (len(g) * (k - len(g)) + len(g) * (len(g) - 1) // 2)
                           for g in groups)
    assert query_count == expected_queries
    return {"vertices": k, "edge_bits": width, "ring": ring.name,
            "group_sizes": list(map(len, groups)),
            "matrix_dimension": d, "table_queries": query_count,
            "streamed_rectangular_blocks": rectangular,
            "recursive_conditioning_queries": conditioned_queries,
            "value": direct}


def setting_checks(q, domain, ring, rng):
    """Independent generic two-variable factors, as in a grouped entry."""
    weights = [[rng.randrange(ring.order) for _ in range(domain)] for _ in range(q)]
    pairs = {(i, j): [[rng.randrange(ring.order) for _ in range(domain)]
                     for _ in range(domain)] for i, j in combinations(range(q), 2)}
    direct = 0
    for word in product(range(domain), repeat=q):
        term = 1
        for i, x in enumerate(word):
            term = ring.mul(term, weights[i][x])
        for (i, j), matrix in pairs.items():
            term = ring.mul(term, matrix[word[i]][word[j]])
        direct = ring.add(direct, term)

    groups = split_three(q)
    settings = [list(product(range(domain), repeat=len(g))) for g in groups]
    dims = list(map(len, settings))
    d = 1
    while d < max(dims):
        d *= 2
    cluster_weight = []
    for group, words in zip(groups, settings):
        values = []
        for word in words:
            lab = dict(zip(group, word))
            term = 1
            for i in group:
                term = ring.mul(term, weights[i][lab[i]])
            for i, j in combinations(group, 2):
                term = ring.mul(term, pairs[i, j][lab[i]][lab[j]])
            values.append(term)
        cluster_weight.append(values)

    matrices = []
    for ga, gb in ((0, 1), (1, 2), (2, 0)):
        matrix = [[0] * d for _ in range(d)]
        for ai, aw in enumerate(settings[ga]):
            for bi, bw in enumerate(settings[gb]):
                term = cluster_weight[ga][ai]
                for va, xa in zip(groups[ga], aw):
                    for vb, xb in zip(groups[gb], bw):
                        val = pairs[va, vb][xa][xb] if va < vb else pairs[vb, va][xb][xa]
                        term = ring.mul(term, val)
                matrix[ai][bi] = term
        matrices.append(matrix)
    assert trace_product(*matrices, ring) == direct
    return {"nodes": q, "domain": domain, "ring": ring.name,
            "matrix_dimension": d, "value": direct}


def run():
    rng = random.Random(20261030)
    multiplication_cases = []
    for ring in RINGS:
        for exponent in range(5):
            d = 1 << exponent
            a = [[rng.randrange(ring.order) for _ in range(d)] for _ in range(d)]
            b = [[rng.randrange(ring.order) for _ in range(d)] for _ in range(d)]
            stats = {"multiplications": 0, "additions": 0}
            assert strassen(a, b, ring, stats) == naive_product(a, b, ring)
            assert stats["multiplications"] == 7 ** exponent
            assert stats["additions"] == 6 * (7 ** exponent - 4 ** exponent)
            multiplication_cases.append({"ring": ring.name, "dimension": d, **stats})

    clique_cases = [clique_checks(k, width, ring, rng) for ring in RINGS
                    for k, width in ((2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (6, 1))]
    clique_cases.extend(clique_checks(5, 1, ring, rng, (1, 2, 2)) for ring in RINGS)
    setting_cases = [setting_checks(q, domain, ring, rng) for ring in RINGS
                     for q in range(1, 8) for domain in (2, 3)]

    schedule_cases = []
    for x0 in (6, 8, 18):
        k0 = 1 << x0
        c = k0 - 5
        k = k0
        for j in range(7):
            g = (1 << (3 * j + x0)) * k
            # omega0<3, rho<1: this is a rational upper bound for (12).
            eta_upper = 1 + Fraction(3, 1 << (j + 5))
            b_upper = (eta_upper + 1) / k + Fraction(2, g) + Fraction(c, g ** 4)
            assert b_upper < Fraction(4, k)
            assert k >= (1 << (j + 5))
            assert Fraction(4, k) <= Fraction(1, 1 << (j + 3))
            schedule_cases.append({"x0": x0, "level": j, "log2_K": k.bit_length() - 1,
                                   "construction_term_dominated": True})
            k = g * g
    return {"status": "pass", "scope": "exact finite identities and arithmetic counts; no improved logarithmic exponent",
            "matrix_multiplication_cases": multiplication_cases,
            "clique_cases": clique_cases,
            "grouped_entry_cases": setting_cases,
            "schedule_cases": schedule_cases,
            "space": "recursive conditioning: linear assignment storage; matrix variants: exponential ring storage"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"],
                      "matrix_cases": len(result["matrix_multiplication_cases"]),
                      "clique_cases": len(result["clique_cases"]),
                      "grouped_entry_cases": len(result["grouped_entry_cases"]),
                      "schedule_cases": len(result["schedule_cases"])}))


if __name__ == "__main__":
    main()
