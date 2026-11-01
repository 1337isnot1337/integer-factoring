"""Exact F2 pairing-type ranks using the proved constant local blocks.

Coefficients are indexed by monomials (r,h,c) in lexicographic order.
This is a small-seed research evaluator, not an asymptotic factoring claim.
"""

import math, time
from functools import lru_cache
from check_modular_crossing import rank_two
from check_modular_recombination import insert, represent


@lru_cache(None)
def monomials(t):
    return [(r, h, t - r - h) for r in range(t + 1) for h in range(t - r + 1)]


def multiply(p, q):
    out = set()
    for a in p:
        for b in q:
            x = tuple(u + v for u, v in zip(a, b))
            if x in out:
                out.remove(x)
            else:
                out.add(x)
    return out


def powpoly(p, n):
    out = {(0, 0, 0)}
    while n:
        if n & 1:
            out = multiply(out, p)
        p = multiply(p, p)
        n >>= 1
    return out


R = {(1, 0, 0)}
H = {(0, 1, 0)}
C = {(0, 0, 1)}
D = R | H | C


@lru_cache(None)
def block(t, a, b, c, odd):
    index = {m: j for j, m in enumerate(monomials(t))}
    forms = []
    es = []
    m = c if odd else a + c
    n = 1 << m
    for k in range(1 if odd else a + 1):
        for l in range(c + 1):
            p = multiply(powpoly(D, a) if odd else powpoly(R | H, k), powpoly(H | C, l))
            p = multiply(p, powpoly(C, b + c - l if odd else t - k - l))
            forms.append(sum(1 << index[x] for x in p))
            es.append([0] * n)
    for row in range(n):
        free = (n - 1) ^ row
        u = free
        while True:
            k = 0 if odd else (u & ((1 << a) - 1)).bit_count()
            l = u.bit_count() if odd else (u >> a).bit_count()
            es[k * (c + 1) + l][row] |= 1 << (row | u)
            if not u:
                break
            u = (u - 1) & free
    table = [0] * (1 << len(forms))
    rows = [0] * n
    for i in range(1, len(table)):
        bit = (i & -i).bit_length() - 1
        rows = [x ^ y for x, y in zip(rows, es[bit])]
        gray = i ^ (i >> 1)
        table[gray] = n if gray & 1 else rank_two(rows)
    return forms, table


def data(t, d):
    out = []
    for a in range(t + 1):
        for c in range(t - a + 1):
            b = t - a - c
            w = (
                math.factorial(t)
                // (math.factorial(a) * math.factorial(b) * math.factorial(c))
                * (d - 1 if d % 2 else d - 2) ** b
                * math.comb(d, 2) ** c
            )
            if w:
                forms, table = block(t, a, b, c, bool(d % 2))
                out.append((forms, table, w))
    return out


def table_ranks(t, d):
    types = data(t, d)
    n = len(monomials(t))
    ranks = [0] * (1 << n)
    codes = [0] * len(types)
    changes = [
        [
            sum(((p >> j) & 1) << k for k, p in enumerate(forms))
            for forms, table, w in types
        ]
        for j in range(n)
    ]
    for i in range(1, 1 << n):
        bit = (i & -i).bit_length() - 1
        for j, change in enumerate(changes[bit]):
            codes[j] ^= change
        ranks[i ^ (i >> 1)] = sum(
            table[code] * w for code, (forms, table, w) in zip(codes, types)
        )
    assert ranks[-1] == d ** (2 * t)
    return ranks


def optimize(t, d):
    start = time.monotonic()
    rs = table_ranks(t, d)
    n = len(monomials(t))
    ix = {m: j for j, m in enumerate(monomials(t))}
    rev = [ix[(c, h, r)] for r, h, c in monomials(t)]
    reverse = [0] * (1 << n)
    for mask in range(1, len(reverse)):
        bit = mask & -mask
        reverse[mask] = reverse[mask ^ bit] ^ (1 << rev[bit.bit_length() - 1])
    costs = sorted((rs[mask] * rs[reverse[mask]], mask) for mask in range(1, len(rs)))
    basis = {}
    selected = []
    target = len(rs) - 1
    for cost, mask in costs:
        if insert(basis, mask, 1 << len(selected)):
            selected.append(mask)
        s = represent(basis, target)
        if s is not None:
            witness = [x for j, x in enumerate(selected) if (s >> j) & 1]
            return {
                "seed": d,
                "word_length": t,
                "dimension": n,
                "cutoff": cost,
                "rho": math.log2(cost) / (4 * t * math.log2(d)),
                "witness": [
                    {"mask": hex(x), "row": rs[x], "column": rs[reverse[x]]}
                    for x in witness
                ],
                "seconds": time.monotonic() - start,
            }
