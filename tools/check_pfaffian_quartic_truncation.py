#!/usr/bin/env python3
"""Exact quartic Gaussian truncation identity and one stopped q8 completion."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb, lcm
from pathlib import Path
from random import Random

import pfaffian_quartic_truncation as model

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT/'evidence/pfaffian-quartic-truncation.json'
base = model.base


def universal_identity():
    random, count, digest = Random(20261103), 0, sha256()
    for size in (4, 6, 8):
        edges = list(combinations(range(size), 2))
        codes = range(1 << len(edges)) if size == 4 else [random.getrandbits(len(edges)) for _ in range(64)]
        for code in codes:
            pair = {(1 << a) | (1 << b): (code >> i) & 1 for i, (a, b) in enumerate(edges)}
            @cache
            def pf(mask):
                if not mask:
                    return 1
                if mask.bit_count() % 2:
                    return 0
                first, out = mask & -mask, 0
                rest = mask ^ first
                for j in range(size):
                    bit = 1 << j
                    if rest & bit:
                        out ^= pair.get(first | bit, 0)*pf(rest ^ bit)
                return out
            gaussian = [pf(mask) for mask in range(1 << size)]
            quartic = [0]*(1 << size)
            quartic[0] = 1
            for support in range(1 << size):
                if support.bit_count() == 4 and gaussian[support]:
                    for mask in range(1 << size):
                        if not mask & support:
                            quartic[mask | support] ^= quartic[mask]
            assert quartic == [g if mask.bit_count() % 4 == 0 else 0 for mask, g in enumerate(gaussian)]
            result = [0]*(1 << size)
            for a, x in enumerate(gaussian):
                if x:
                    for b, y in enumerate(quartic):
                        if y and not a & b:
                            result[a | b] ^= 1
            assert result == [1 if mask == 0 else pair.get(mask, 0) for mask in range(1 << size)]
            digest.update(json.dumps([size, code, gaussian, quartic]).encode())
            count += 1
    return dict(controls=count, sha256=digest.hexdigest())


def literal_controls():
    @cache
    def pf(active, axis):
        if not active:
            return 1
        if len(active) % 2:
            return 0
        (v, c), rest = active[0], active[1:]
        out = 0
        for i, (w, d) in enumerate(rest):
            if c == d and (v >> axis & 1) != (w >> axis & 1):
                out ^= model.multiply(model.kernel(axis)[v][w], pf(rest[:i]+rest[i+1:], axis))
        return out
    @cache
    def reference_coefficient(word):
        x, y = (pf(word, axis) & 65535 for axis in (0, 1))
        out = x
        for j in range(16):
            if y >> j & 1:
                out ^= 1 << (32-j)
        if (x.bit_count()+y.bit_count()) % 2:
            out ^= 1 << 16
        assert out == model.quartic(word)
        return out
    @cache
    def at_axis(word, axis):
        polynomial = reference_coefficient(word)
        return polynomial if axis == 0 else sum(((polynomial >> j) & 1) << (32-j) for j in range(33))
    records, digest, witness = [], sha256(), None
    for axis in (0, 1):
        left, right = base.cube.axis_parts(axis)
        pair = [[model.kernel(axis)[v][w] & base.MASK for w in right] for v in left]
        count, doubles = 0, 0
        for word in product(range(4), repeat=8):
            active = tuple((v, c) for v, c in enumerate(word) if c)
            aa = tuple((i, word[v]) for i, v in enumerate(left) if word[v])
            bb = tuple((i, word[v]) for i, v in enumerate(right) if word[v])
            value = pf(active, axis)
            if len(active) % 2 == 0:
                for support in combinations(range(len(active)), 4):
                    quartet = tuple(active[i] for i in support)
                    if len({v >> axis & 1 for v, c in quartet}) < 2:
                        continue
                    rest = tuple(a for i, a in enumerate(active) if i not in support)
                    value ^= model.multiply(at_axis(quartet, axis), pf(rest, axis))
                if len(active) == 8:
                    double = 0
                    for support in combinations(range(1, 8), 3):
                        indices = (0,)+support
                        quartet = tuple(active[i] for i in indices)
                        rest = tuple(a for i, a in enumerate(active) if i not in indices)
                        if (len({v >> axis & 1 for v, c in quartet}) == 2
                                and len({v >> axis & 1 for v, c in rest}) == 2):
                            double ^= model.multiply(at_axis(quartet, axis), at_axis(rest, axis))
                    value ^= double
                    if double:
                        doubles += 1
                        if witness is None:
                            witness = dict(axis=axis, word=word, polynomial=double)
            assert value.bit_count() % 2 == int(all(word[v] == word[v ^ 7] for v in range(8)))
            assert model.jet_entry(aa, bb, left, right, pair, axis) == value & base.MASK
            assert not value or value.bit_length()-1 <= 8*len(active)
            digest.update(json.dumps([axis, word, value]).encode())
            count += 1
        assert count == 65536 and doubles == 609
        records.append(dict(axis=axis, literal_full_support_entries=count, nonzero_two_quartic_entries=doubles))
    assert witness is not None
    return dict(cases=records, two_quartic_witness=witness, sha256=digest.hexdigest())


def pressure_certificate():
    tail, known_sum = 3759, 14*6+322*16
    mean = Fraction(56*4096-known_sum, tail)
    assert mean == Fraction(32020, 537) and mean < 60
    rho = Fraction(4, 5)
    assert 8/rho == 10 and 4*rho == Fraction(16, 5) and 4*rho*10 == 32
    polynomial = {0: 1}
    for _ in range(5):
        following = {}
        for i, c in polynomial.items():
            for j, d in ((0, 1), (3, 14), (8, 322), (30, 3759)):
                following[i+j] = following.get(i+j, 0)+c*d
        polynomial = following
    for j in range(17):
        polynomial[5*j] = polynomial.get(5*j, 0)-comb(16, j)*7**j
    polynomial = {j: c for j, c in polynomial.items() if c}
    assert min(polynomial) == 3
    polynomial = {j-3: c for j, c in polynomial.items()}
    assert max(polynomial) == 147 and polynomial[0] == 70
    degree = 512
    denominator = lcm(*(comb(degree, j) for j in polynomial))
    bernstein = [0]*(degree+1)
    for j, c in polynomial.items():
        weight, value = c*(denominator//comb(degree, j)), 1
        for k in range(j, degree+1):
            bernstein[k] += weight*value
            if k < degree:
                value = value*(k+1)//(k+1-j)
    assert min(bernstein) > 0
    row = bernstein[:]
    for j in range(degree+1):
        value = comb(degree, j)*row[0]
        assert value % denominator == 0 and value//denominator == polynomial.get(j, 0)
        row = [b-a for a, b in zip(row, row[1:])]
    return dict(dimension=4096, determinant_degree=56*4096, known_orders=337,
                known_sum=known_sum, tail_count=tail, tail_mean=str(mean), tail_cap=60,
                pressure_floor=str(rho), polynomial=sorted(polynomial.items()),
                actual_polynomial_degree=147, bernstein_degree=degree,
                positive_coefficients=len(bernstein), monomial_reconstructions=degree+1,
                minimum_index=bernstein.index(min(bernstein)), minimum=str(Fraction(min(bernstein), denominator)),
                integer_coefficient_sha256=sha256(json.dumps(bernstein, separators=(',', ':')).encode()).hexdigest())


def compute():
    identity = universal_identity()
    patterns = [(c,)*4 for c in range(1, 8)]
    for a, b in combinations(range(1, 8), 2):
        patterns += sorted(set(permutations((a, a, b, b))))
    invariants = 0
    for support in combinations(range(8), 4):
        for colors in patterns:
            model.quartic(tuple(zip(support, colors)))
            invariants += 1
    assert invariants == 9310
    matrices, records, combined, six = {}, [], {}, 0
    for axis in (0, 1):
        for target in (0, 1):
            blocks, controls = model.charge_block(axis, target)
            low_rank = base.column_rank(blocks, 16)
            assert low_rank == (16 if target == 0 else 20)
            matrices[axis, target] = blocks
            six += controls
    for (axis, target), blocks in matrices.items():
        result = base.ranks(blocks)
        assert result['smith_prefix'] == ({0: 1, 16: 42} if target == 0 else {6: 2, 16: 40})
        independent = base.column_rank(blocks, 19)
        assert independent == result['toeplitz_ranks'][19]
        for order, count in result['smith_prefix'].items():
            combined.setdefault(axis, Counter())[order] += (1 if target == 0 else 7)*count
        records.append(dict(axis=axis, charge=target, **result, independent_column_rank_19=independent))
    assert six == 115616
    assert all(prefix == Counter({0: 1, 6: 14, 16: 322}) for prefix in combined.values())
    literal = literal_controls()
    names = ['tools/pfaffian_quartic_truncation.py', 'tools/check_pfaffian_quartic_truncation.py',
             'tools/check_fixed_curve_quartic_floor.py', 'tools/cube_homogeneous_gaussian.py',
             'tools/check_quartic_charge_fusion.py']
    return dict(scope='Ordinary quartic Pfaffian identity, exact generic low jets and STOP certificate for one fixed q8 covariant completion. No full profile, independent mathematical review or stronger complete factoring bound.',
                registered_base='ccd12e7ac337dba2b1f2c7f30073db725f912481', precision=19,
                identity_controls=identity, invariant_quartets=invariants,
                leading_six_occupied_controls=six, cases=records,
                whole_prefixes={axis: dict(prefix) for axis, prefix in combined.items()},
                literal_controls=literal, pressure_certificate=pressure_certificate(),
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(compute()))
    if args.write:
        EXPECTED.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(EXPECTED.read_text())
    print('Pfaffian quartic truncation: exact identity, repeated-colour jets and fixed-q8 4/5 STOP PASS')
