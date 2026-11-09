#!/usr/bin/env python3
"""Primitive central idempotents of the twelve-site coefficient algebra.

The exact ten-dimensional table from the coefficient-table certificate is
decomposed over F5.  The decomposition uses generic elements, split linear
minimal polynomials, and CRT idempotents; no dense cut operator is built.
It writes the primitive central idempotents, their regular ranks, and the
trace-form block structure to evidence.
"""
from __future__ import annotations

import argparse
import json
import random
import time
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIME = 5
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'
TABLE = ROOT / 'evidence/bounded-spin-twelve-coefficient-table.json'


def ptrim(a):
    a = [x % PRIME for x in a]
    while a and not a[-1]:
        a.pop()
    return a


def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % PRIME
    return ptrim(c)


def psub(a, b):
    n = max(len(a), len(b))
    c = [0] * n
    for i in range(n):
        c[i] = ((a[i] if i < len(a) else 0)
                - (b[i] if i < len(b) else 0)) % PRIME
    return ptrim(c)


def pdivmod(a, b):
    a = ptrim(a)
    b = ptrim(b)
    assert b, 'division by zero polynomial'
    q = [0] * max(0, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, PRIME)
    while len(a) >= len(b) and a:
        f = a[-1] * inv % PRIME
        sh = len(a) - len(b)
        q[sh] = f
        for i, x in enumerate(b):
            a[sh + i] = (a[sh + i] - f * x) % PRIME
        a = ptrim(a)
    return ptrim(q), a


def pmod(a, b):
    return pdivmod(a, b)[1]


def peval(p, x):
    v = 0
    for c in p[::-1]:
        v = (v * x + c) % PRIME
    return v


def pxgcd(a, b):
    old_r, r = ptrim(a), ptrim(b)
    old_s, s = [1], [0]
    old_t, t = [0], [1]
    while r:
        q, rem = pdivmod(old_r, r)
        old_r, r = r, rem
        old_s, s = s, psub(old_s, pmul(q, s))
        old_t, t = t, psub(old_t, pmul(q, t))
    return ptrim(old_r), ptrim(old_s), ptrim(old_t)


def pinv(a, m):
    g, s, _ = pxgcd(a, m)
    assert len(g) == 1 and g[0] % PRIME, 'nonunit polynomial'
    inv = pow(g[0], -1, PRIME)
    return ptrim([x * inv % PRIME for x in s])


def factor_linear(p):
    """Factor p into powers of linear factors; the algebra is split."""
    rem = ptrim(p)
    factors = []
    for root in range(PRIME):
        f = [(-root) % PRIME, 1]
        mult = 0
        while len(rem) >= 2 and peval(rem, root) == 0:
            q, r = pdivmod(rem, f)
            if r:
                break
            rem = q
            mult += 1
        if mult:
            factors.append((f, mult))
    assert rem == [1], f'minimal polynomial is not split: {rem}'
    return factors


def crt_idempotents(p, factors):
    gs = []
    for f, mult in factors:
        g = [1]
        for _ in range(mult):
            g = pmul(g, f)
        gs.append(g)
    out = []
    for i, g in enumerate(gs):
        m_i = [1]
        for j, h in enumerate(gs):
            if i != j:
                m_i = pmul(m_i, h)
        h_i = pmod(pmul(m_i, pinv(pmod(m_i, g), g)), p)
        out.append(h_i)
    return out


def rank(a):
    a = [row[:] for row in a]
    if not a:
        return 0
    rows, cols = len(a), len(a[0])
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c] % PRIME), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c] % PRIME, -1, PRIME)
        a[r] = [x * inv % PRIME for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c] % PRIME:
                f = a[i][c]
                a[i] = [(x - f * y) % PRIME for x, y in zip(a[i], a[r])]
        r += 1
    return r


def nullspace(a):
    a = [row[:] for row in a]
    if not a:
        return []
    rows, cols = len(a), len(a[0])
    pivots = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c] % PRIME), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c] % PRIME, -1, PRIME)
        a[r] = [x * inv % PRIME for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c] % PRIME:
                f = a[i][c]
                a[i] = [(x - f * y) % PRIME for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        v = [0] * cols
        v[f] = 1
        for row, c in enumerate(pivots):
            v[c] = (-a[row][f]) % PRIME
        basis.append(v)
    return basis


class Algebra:
    def __init__(self, table, unit):
        self.table = table
        self.unit = unit
        self.size = len(unit)

    def mul(self, a, b):
        n = self.size
        out = [0] * n
        for i in range(n):
            if a[i]:
                for j in range(n):
                    if b[j]:
                        f = a[i] * b[j] % PRIME
                        row = self.table[i][j]
                        for k in range(n):
                            out[k] = (out[k] + f * row[k]) % PRIME
        return out

    def minpoly(self, a):
        n = self.size
        vecs = [self.unit[:]]
        for k in range(1, n + 1):
            ak = self.unit[:]
            for _ in range(k):
                ak = self.mul(ak, a)
            cols = k
            aug = [[vecs[j][i] for j in range(cols)] + [ak[i]]
                   for i in range(n)]
            r = 0
            pivots = []
            for c in range(cols):
                pivot = next((i for i in range(r, n) if aug[i][c] % PRIME), None)
                if pivot is None:
                    continue
                aug[r], aug[pivot] = aug[pivot], aug[r]
                inv = pow(aug[r][c] % PRIME, -1, PRIME)
                aug[r] = [x * inv % PRIME for x in aug[r]]
                for i in range(n):
                    if i != r and aug[i][c] % PRIME:
                        f = aug[i][c]
                        aug[i] = [(x - f * y) % PRIME
                                  for x, y in zip(aug[i], aug[r])]
                pivots.append(c)
                r += 1
            consistent = True
            for i in range(r, n):
                if not any(aug[i][:cols]) and aug[i][-1] % PRIME:
                    consistent = False
                    break
            if consistent:
                coef = [0] * cols
                for i, c in enumerate(pivots):
                    coef[c] = aug[i][-1] % PRIME
                return [(-coef[j]) % PRIME for j in range(cols)] + [1]
            vecs.append(ak)
        return None

    def eval_poly(self, a, poly):
        out = [0] * self.size
        for j, c in enumerate(poly):
            if c:
                t = self.unit[:]
                for _ in range(j):
                    t = self.mul(t, a)
                for i in range(self.size):
                    out[i] = (out[i] + c * t[i]) % PRIME
        return out

    def left_rank(self, e):
        n = self.size
        matrix = [[sum(e[a] * self.table[a][b][c] for a in range(n)) % PRIME
                   for c in range(n)] for b in range(n)]
        return rank(matrix)

    def trace_left(self, e):
        n = self.size
        return sum(e[a] * sum(self.table[a][b][b] for b in range(n))
                   for a in range(n)) % PRIME

    def block_basis(self, e):
        n = self.size
        basis = []
        for j in range(n):
            v = self.mul(e, [1 if k == j else 0 for k in range(n)])
            if any(v) and rank(basis + [v]) > len(basis):
                basis.append(v)
        return basis


def decompose(alg, seed=0, tries=200):
    rng = random.Random(seed)
    primitives = []
    stack = [alg.unit[:]]
    while stack:
        e = stack.pop()
        dim = alg.left_rank(e)
        if dim == 1:
            primitives.append(e)
            continue
        split = False
        for _ in range(tries):
            a = alg.mul(e, [rng.randrange(PRIME) for _ in range(alg.size)])
            if not any(a):
                continue
            p = alg.minpoly(a)
            if p is None:
                continue
            factors = factor_linear(p)
            non_x = [(f, m) for f, m in factors
                     if not (len(f) == 2 and f[0] == 0 and f[1] == 1)]
            if len(non_x) < 2:
                continue
            hs = crt_idempotents(p, factors)
            subs = []
            for (f, _), h in zip(factors, hs):
                if len(f) == 2 and f[0] == 0 and f[1] == 1:
                    continue
                ei = alg.mul(e, alg.eval_poly(a, h))
                if not any(ei):
                    continue
                assert alg.mul(ei, ei) == ei
                if ei != e:
                    subs.append(ei)
            uniq = []
            for ei in subs:
                if ei not in uniq:
                    uniq.append(ei)
            if len(uniq) >= 2:
                total = [0] * alg.size
                for ei in uniq:
                    for k in range(alg.size):
                        total[k] = (total[k] + ei[k]) % PRIME
                if total == e:
                    stack.extend(uniq)
                    split = True
                    break
        if not split:
            primitives.append(e)
    primitives.sort()
    return primitives


def compute():
    start = time.time()
    raw = TABLE.read_text()
    data = json.loads(raw)
    table = data['table']
    unit = data['unit']
    alg = Algebra(table, unit)
    n = alg.size

    for i in range(n):
        e_i = [1 if k == i else 0 for k in range(n)]
        assert alg.mul(unit, e_i) == e_i
        assert alg.mul(e_i, unit) == e_i

    idempotents = decompose(alg)
    total = [0] * n
    for e in idempotents:
        for k in range(n):
            total[k] = (total[k] + e[k]) % PRIME
    assert total == unit

    for i, e in enumerate(idempotents):
        assert any(e)
        assert alg.mul(e, e) == e
        for j, f in enumerate(idempotents):
            want = e if i == j else [0] * n
            assert alg.mul(e, f) == want

    blocks = []
    for idx, e in enumerate(idempotents):
        basis = alg.block_basis(e)
        dim = len(basis)
        assert alg.left_rank(e) == dim
        form = [[alg.trace_left(alg.mul(basis[a], basis[b]))
                 for b in range(dim)] for a in range(dim)]
        form_rank = rank(form)
        if dim == 1:
            assert form == [[1]], form
            kind = 'scalar'
            nilpotent = None
        elif dim == 2:
            assert form_rank in (1, 2), form
            kind = 'dual' if form_rank == 1 else 'field'
            nilpotent = None
            if kind == 'dual':
                found = None
                for a in range(PRIME):
                    for b in range(PRIME):
                        if not a and not b:
                            continue
                        x = [(a * basis[0][k] + b * basis[1][k]) % PRIME
                             for k in range(n)]
                        if alg.mul(x, x) == [0] * n:
                            found = x
                            break
                    if found is not None:
                        break
                assert found is not None, 'dual block has no nilpotent'
                assert rank([basis[0], found]) == 2
                assert alg.mul(basis[0], found) == found
                assert alg.mul(found, found) == [0] * n
                nilpotent = found
        else:
            raise AssertionError(f'unexpected block dimension {dim}')
        blocks.append({
            'index': idx,
            'idempotent': e,
            'regular_rank': dim,
            'block_basis': basis,
            'trace_form': form,
            'trace_form_rank': form_rank,
            'kind': kind,
            'nilpotent': nilpotent,
        })

    assert sorted(b['regular_rank'] for b in blocks) == [1] * 6 + [2] * 2
    assert sum(1 for b in blocks if b['kind'] == 'scalar') == 6
    assert sum(1 for b in blocks if b['kind'] == 'dual') == 2
    assert not any(b['kind'] == 'field' for b in blocks)

    radicals = [b['nilpotent'] for b in blocks if b['kind'] == 'dual']
    assert len(radicals) == 2
    assert rank(radicals) == 2
    for x in radicals:
        assert alg.mul(x, x) == [0] * n
    assert alg.mul(radicals[0], radicals[1]) == [0] * n

    full_form = [[alg.trace_left(alg.mul(
        [1 if k == i else 0 for k in range(n)],
        [1 if k == j else 0 for k in range(n)])) for j in range(n)]
        for i in range(n)]
    full_rank = rank(full_form)
    kernel = nullspace(full_form)
    assert full_rank == 8
    assert len(kernel) == 2
    assert rank(kernel + radicals) == 2
    for r in radicals:
        assert all(alg.trace_left(alg.mul(r, [1 if k == i else 0
                                              for k in range(n)])) == 0
                    for i in range(n))

    result = {
        'status': 'pass',
        'prime': PRIME,
        'length': 12,
        'dimension': 10,
        'basis_indices': data['basis_indices'],
        'table_sha256': sha256(raw.encode()).hexdigest(),
        'structure': 'F5^6 x F5[e]/(e^2) x F5[e]/(e^2)',
        'regular_rank_partition': [b['regular_rank'] for b in blocks],
        'primitive_central_idempotents': blocks,
        'trace_form': {
            'rank': full_rank,
            'radical_dimension': len(kernel),
            'radical_basis': radicals,
            'square_zero': True,
        },
        'checks': {
            'idempotent_count': len(blocks),
            'sum_is_unit': True,
            'pairwise_orthogonal': True,
            'six_scalar_blocks': True,
            'two_dual_blocks': True,
            'regular_rank_sum': sum(b['regular_rank'] for b in blocks),
            'trace_form_rank_plus_radical': full_rank + len(kernel),
        },
        'seconds': round(time.time() - start, 3),
        'scope': ('Exact primitive central idempotents of the twelve-site '
                  'coefficient algebra over F5. The algebra is the product '
                  'of six copies of F5 and two dual-number blocks. Cut-module '
                  'ranks, the two-cut pairing, and any factoring score change '
                  'remain open.'),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = compute()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({
        'status': result['status'],
        'structure': result['structure'],
        'idempotents': len(result['primitive_central_idempotents']),
        'regular_rank_partition': result['regular_rank_partition'],
        'trace_form': result['trace_form']['rank'],
        'radical_dimension': result['trace_form']['radical_dimension'],
        'seconds': result['seconds'],
    }, indent=2))


if __name__ == '__main__':
    main()
