#!/usr/bin/env python3
"""Exact F_2 checks for self-dual nilpotent-cycle crossing labels."""

from __future__ import annotations

import argparse
import itertools
import json
from functools import lru_cache, reduce
from pathlib import Path

from check_finite_seed_amplification import (
    PIECES, add, column, multiply, rank,
)
from check_modular_crossing import inverse_two


def transpose(a):
    return tuple(sum(((r >> j) & 1) << i for i, r in enumerate(a))
                 for j in range(len(a)))


def kron(a, b):
    return tuple(sum(((x >> j) & 1) * y << (len(b) * j)
                     for j in range(len(a))) for x in a for y in b)


def identity(n):
    return tuple(1 << i for i in range(n))


def transformed_seed():
    s = (5, 1, 2)
    inv = inverse_two(list(s), 3)
    left = kron(transpose(s), transpose(s))
    right = kron(transpose(inv), transpose(inv))
    return tuple(multiply(left, multiply(a, right)) for a in PIECES)


def cycle_basis(t):
    out = []
    for word in itertools.product(range(3), repeat=t):
        active = [j for j, a in enumerate(word) if a != 1]
        if not len(active) % 2:
            continue
        j = active[0]
        if word[j] != 0:
            continue
        terms = [word]
        for k in active[1:]:
            if word[k] == 2:
                w = list(word)
                w[j], w[k] = 2, 0
                terms.append(tuple(w))
        out.append(tuple(terms))
    assert len(out) == (3**t - (-1)**t) // 4
    return out


def gram(basis):
    supports = [set(a) for a in basis]
    return tuple(sum((len(a & b) % 2) << j for j, b in enumerate(supports))
                 for a in supports)


def complete_profile(seed, t):
    basis = cycle_basis(t)
    d = len(basis)
    g = gram(basis)
    monomials = list(itertools.combinations_with_replacement(range(20), t))
    indices = {a: i for i, a in enumerate(monomials)}
    entries = {
        q: tuple(i for i, a in enumerate(seed)
                 if (a[3*q[0] + q[1]] >> (3*q[2] + q[3])) & 1)
        for q in itertools.product(range(3), repeat=4)
    }

    @lru_cache(None)
    def polynomial(key):
        if not all(key):
            return 0
        out = 0
        for term in itertools.product(*key):
            out ^= 1 << indices[tuple(sorted(term))]
        return out

    rows = {}
    for x, y, u, v in itertools.product(range(d), repeat=4):
        p = 0
        for xx, yy, uu, vv in itertools.product(
                basis[x], basis[y], basis[u], basis[v]):
            key = tuple(sorted(entries[q] for q in zip(xx, yy, uu, vv)))
            p ^= polynomial(key)
        while p:
            bit = p & -p
            p ^= bit
            k = bit.bit_length() - 1
            rows.setdefault(k, [0]*(d*d))[x*d+y] ^= 1 << (u*d+v)
    expected = tuple(
        sum((((g[x] >> v) & 1)*((g[y] >> u) & 1)) << (u*d+v)
            for u in range(d) for v in range(d))
        for x in range(d) for y in range(d)
    )
    assert reduce(add, rows.values(), (0,)*(d*d)) == expected
    profile = [(monomials[k], rank(a), rank(column(a, d)))
               for k, a in sorted(rows.items())]
    maximum = max(r*c for _, r, c in profile)
    return {
        "t": t, "labels": d, "pieces": len(profile),
        "maximum_rank_product": maximum,
        "maximizers": [a for a in profile if a[1]*a[2] == maximum],
        "profile": profile,
        "raw_gram_crossing_verified": True,
    }


def core_check(seed):
    words = [w for w in itertools.product(range(3), repeat=3)
             if w.count(1) == 1]
    idx = {w: i for i, w in enumerate(words)}
    d = len(words)
    diff = [0]*d
    for j, w in enumerate(words):
        for k, z in enumerate(w):
            if z == 2:
                v = list(w)
                v[k] = 0
                diff[idx[tuple(v)]] ^= 1 << j
    diff = tuple(diff)
    assert not any(multiply(diff, diff))

    orders = list(itertools.permutations((4, 5, 17)))
    b = []
    for x in words:
        for y in words:
            row = 0
            for ui, u in enumerate(words):
                for vi, v in enumerate(words):
                    entry = 0
                    for typ in orders:
                        term = 1
                        for j, i in enumerate(typ):
                            term &= (seed[i][3*x[j]+y[j]] >>
                                     (3*u[j]+v[j])) & 1
                        entry ^= term
                    row |= entry << (d*ui+vi)
            b.append(row)
    b = tuple(b)

    active = tuple(itertools.product((0, 2), repeat=4))
    for i in (4, 5, 17):
        assert not any((seed[i][3*x+y] >> (3*u+v)) & 1
                       for x, y, u, v in active)
    actual = {(x, y, u, v) for x, y, u, v in active
              if (seed[16][3*x+y] >> (3*u+v)) & 1}
    assert actual == {(x, y, x, y) for x in (0, 2) for y in (0, 2)}

    dx, dy = kron(diff, identity(d)), kron(identity(d), diff)

    def transform(a):
        a = add(a, multiply(transpose(dx), multiply(a, dx)))
        return add(a, multiply(transpose(dy), multiply(a, dy)))

    bp = transform(b)
    assert transform(bp) == b
    j = add(identity(d*d), kron(diff, diff))
    assert multiply(j, j) == identity(d*d)
    assert column(bp, d) == multiply(transpose(j), multiply(column(b, d), j))
    core_ranks = (rank(b), rank(column(b, d)))
    cycle_ranks = (rank(bp), rank(column(bp, d)))
    assert core_ranks == cycle_ranks == (48, 48)

    # Direct sparse substitution on all four legs, independently of the
    # matrix-operation implementation of equation (7).
    basis = []
    for w in words:
        terms = [(0, w)]
        for k, value in enumerate(w):
            if value == 2:
                a = list(w)
                a[k] = 0
                terms.append((2, tuple(a)))
        basis.append(terms)
    checked = 0
    for x, y, u, v in itertools.product(range(d), repeat=4):
        value = 0
        for xx, yy, uu, vv in itertools.product(
                basis[x], basis[y], basis[u], basis[v]):
            first = (seed[16][3*xx[0]+yy[0]] >> (3*uu[0]+vv[0])) & 1
            if first:
                value ^= (b[d*idx[xx[1]]+idx[yy[1]]] >>
                          (d*idx[uu[1]]+idx[vv[1]])) & 1
        assert value == ((bp[d*x+y] >> (d*u+v)) & 1)
        checked += 1
    return {
        "m": 1, "labels_in_block": d,
        "core_ranks": core_ranks, "cycle_ranks": cycle_ranks,
        "direct_sparse_entries_checked": checked,
        "first_position_forcing_verified": True,
        "four_term_involution_verified": True,
        "column_invertible_equivalence_verified": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks = []
    for t in range(1, 7):
        basis = cycle_basis(t)
        g = gram(basis)
        assert multiply(g, g) == identity(len(basis))
        checks.append({"t": t, "labels": len(basis), "gram_involution": True})
    seed = transformed_seed()
    result = {
        "field": 2, "source_commit": "e96ee50560284301428d72421330d8773155d25b",
        "gram_checks": checks,
        "complete_profile_t3": complete_profile(seed, 3),
        "surviving_core": core_check(seed),
    }
    assert result["complete_profile_t3"]["pieces"] == 205
    assert result["complete_profile_t3"]["maximum_rank_product"] == 144
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "gram_checks": len(checks),
        "t3_labels": 7, "t3_pieces": 205, "t3_maximum_product": 144,
        "surviving_core": result["surviving_core"],
        "all_checks_passed": True,
    }, indent=2))


if __name__ == "__main__":
    main()
