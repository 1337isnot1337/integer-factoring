#!/usr/bin/env python3
"""Exact eight-site rank algebra, including every characteristic-five extension.

NumPy is optional for the repository, but required for this focused replay.
All arrays are reduced modulo five before bounded integer matrix products.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
from random import Random

import numpy as np

from bounded_spin_crossing import paths, path_entry, path_norm, polynomial_entry

ROOT = Path(__file__).resolve().parents[1]
P, D, N = 5, 6, 169


def rank(a):
    a = np.array(a, dtype=np.int64, copy=True) % P
    r = 0
    for j in range(a.shape[1]):
        hits = np.flatnonzero(a[r:, j])
        if not len(hits):
            continue
        h = r + int(hits[0])
        a[[r, h]] = a[[h, r]]
        a[r, j:] = a[r, j:] * pow(int(a[r, j]), -1, P) % P
        a[r+1:, j:] = (a[r+1:, j:] - a[r+1:, j, None]*a[r, None, j:]) % P
        r += 1
        if r == len(a):
            break
    return r


def rref(rows):
    a = [list(map(int, row)) for row in rows if any(row)]
    i = 0
    for j in range(D):
        k = next((k for k in range(i, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[i], a[k] = a[k], a[i]
        a[i] = [x*pow(a[i][j], -1, P) % P for x in a[i]]
        for k in range(len(a)):
            if k != i:
                c = a[k][j]
                a[k] = [(x-c*y) % P for x, y in zip(a[k], a[i])]
        i += 1
        if i == len(a):
            break
    return tuple(tuple(row) for row in a[:i])


def kernel(rows):
    a = rref(rows)
    pivots = [next(j for j, x in enumerate(row) if x) for row in a]
    result = []
    for j in range(D):
        if j in pivots:
            continue
        v = [int(i == j) for i in range(D)]
        for p, row in zip(pivots, a):
            v[p] = -row[j] % P
        result.append(v)
    return result


def dot(a, b):
    return sum(int(x)*int(y) for x, y in zip(a, b)) % P


def coefficient_tensor(length):
    """Compress exact coordinate signature blocks; no ambient four-leg tensor.

    A signature fixes x&y and x^y. In that block the local entry is
    z^b(1-z)^(s-b), b=|x xor u|, s=|x xor y|. Contract the two
    path frames on each side, then restore canonical x,y,u,v axes.
    """
    ps = paths(length)
    words = [w for w in range(1 << length) if w.bit_count() == length//2]
    frame = np.array([[path_entry(path, w) % P for w in words] for path in ps],
                     dtype=np.int64)
    indices = {w: i for i, w in enumerate(words)}
    groups = defaultdict(list)
    for x, y in product(words, repeat=2):
        groups[x & y, x ^ y].append((x, y))
    size = len(ps)**2
    result = np.zeros((length+1, size, size), dtype=np.int64)
    for (_, difference), pairs in groups.items():
        a = np.stack([np.outer(frame[:, indices[x]], frame[:, indices[y]]).ravel()
                      for x, y in pairs], axis=1) % P
        local = np.zeros((length+1, len(pairs), len(pairs)), dtype=np.int64)
        s = difference.bit_count()
        for i, (x, _) in enumerate(pairs):
            for h, (_, u) in enumerate(pairs):
                b = (x ^ u).bit_count()
                for j in range(s-b+1):
                    local[b+j, i, h] = comb(s-b, j)*(-1)**j % P
        for j in range(length+1):
            result[j] = (result[j] + (a @ local[j] % P) @ a.T) % P
    q = len(ps)
    return result.reshape(length+1, q, q, q, q).transpose(0, 1, 2, 4, 3)


def matrix_poly(t, coefficients):
    a = np.zeros_like(t)
    identity = np.eye(len(t), dtype=np.int64)
    for c in reversed(coefficients):
        a = (a @ t + int(c)*identity) % P
    return a


def normalized_operators(tensor):
    length, q = len(tensor)-1, tensor.shape[1]
    norms = [path_norm(path) % P for path in paths(length)]
    inverse = np.array([pow(x*y, -1, P) for x, y in product(norms, repeat=2)])
    # Swapping columns at the crossing endpoint makes both endpoints I.
    # Diagonal norm factors and column permutations preserve each cut rank.
    return [tensor.transpose(0, 1, 2, 4, 3).reshape(length+1, q*q, q*q)
            * inverse[None, None, :] % P,
            tensor.transpose(0, 1, 3, 4, 2).reshape(length+1, q*q, q*q)
            * inverse[None, None, :] % P]


def check():
    certificate = json.loads((ROOT/'evidence/bounded-spin-rank-model.json').read_text())
    chosen, models = certificate['basis_indices'], certificate['models']
    controls = 0
    random = Random(1337)
    for length in (4, 8):
        t = coefficient_tensor(length)
        ps = paths(length)
        selected = (list(product(range(len(ps)), repeat=4)) if length == 4 else
                    [tuple(random.randrange(len(ps)) for _ in range(4)) for _ in range(24)])
        for ids in selected:
            raw = t[(slice(None),)+ids]
            denominator = path_norm(ps[ids[2]])*path_norm(ps[ids[3]]) % P
            assert (raw*pow(denominator, -1, P) % P).tolist() == polynomial_entry(
                tuple(ps[i] for i in ids))
            controls += 1
    operators = normalized_operators(t)
    assert t.shape == (9, 13, 13, 13, 13)
    assert rank(t.reshape(9, -1)) == D
    assert not np.any((2*t[1]+8*t[0]) % P)
    assert not np.any((sum(2*j*t[j] for j in range(9))-8*t.sum(axis=0)) % P)
    assert not np.any((2*t[-2]+8*t[-1]) % P)
    target = np.array(certificate['target'], dtype=np.int64)
    forms, cuts, powers_models = [], [], []
    raw_pairs = []
    for a, model in zip(operators, models):
        b = a[chosen]
        assert rank(b.reshape(D, -1)) == D
        identity = np.eye(N, dtype=np.int64)
        assert np.array_equal(a.sum(axis=0) % P, identity)
        assert np.array_equal(np.tensordot(target, b, axes=(0, 0)) % P, identity)
        generator = np.tensordot(model['generator_coefficients'], b, axes=(0, 0)) % P
        minimal = model['minimal_polynomial']
        assert not np.any(matrix_poly(generator, minimal))
        powers = [identity]
        for _ in range(D):
            powers.append(powers[-1] @ generator % P)
        assert rank(np.array(powers[:D]).reshape(D, -1)) == D
        polynomials = model['polynomial_basis']
        for matrix, polynomial in zip(b, polynomials):
            assert np.array_equal(matrix_poly(generator, polynomial), matrix)
        # Degree-six minimal polynomial and six independent powers prove
        # this span is the complete commutative algebra F5[T]/(minimal).
        product_polynomial = [1]
        blocks, total_primary = [], 0
        for block in model['blocks']:
            factor, exponent, nullities = (block[key] for key in
                                           ('factor', 'exponent', 'nullities'))
            assert len(factor) == 2 and factor[1] == 1 and exponent in (1, 2)
            x = matrix_poly(generator, factor)
            power = identity
            assert nullities[0] == 0
            for k in range(1, exponent+1):
                power = power @ x % P
                assert N-rank(power) == nullities[k]
                out = [0]*(len(product_polynomial)+1)
                for j, c in enumerate(product_polynomial):
                    out[j] = (out[j]+c*factor[0]) % P
                    out[j+1] = (out[j+1]+c) % P
                product_polynomial = out
            total_primary += nullities[-1]
            root = -factor[0] % P
            forms.append(tuple(sum(row[j]*pow(root, j, P) for j in range(D)) % P
                               for row in polynomials))
            indices = [len(forms)-1]
            if exponent == 2:
                forms.append(tuple(sum(row[j]*j*pow(root, j-1, P) for j in range(1, D)) % P
                                   for row in polynomials))
                indices.append(len(forms)-1)
            blocks.append((indices, nullities))
        assert product_polynomial == minimal and total_primary == N
        cuts.append(blocks)
        powers_models.append((b, blocks))
        raw_pairs.append([rank(matrix) for matrix in a])

    # Every vector over every extension has one of these zero patterns:
    # its vanishing forms define a common kernel, and every further form
    # vanishing on that kernel also vanishes at the vector.
    flats = {()}
    for form in forms:
        flats |= {rref(list(flat)+[form]) for flat in list(flats)}
    separator = certificate['separator']
    threshold = certificate['threshold']
    assert dot(separator, target) == 1
    cheaper, ranks_by_pattern = [], {}
    for flat in sorted(flats):
        u = kernel(flat)
        vanishing = tuple(not any(dot(form, v) for v in u) for form in forms)
        ranks = []
        for blocks in cuts:
            nullity = 0
            for indices, nullities in blocks:
                order = 0
                for j in indices:
                    if not vanishing[j]:
                        break
                    order += 1
                nullity += nullities[order]
            ranks.append(N-nullity)
        ranks_by_pattern[vanishing] = ranks
        if ranks[0]*ranks[1] < threshold:
            assert all(dot(separator, v) == 0 for v in u)
            cheaper.extend(u)
    assert len(flats) == 759 and len(rref(cheaper)) == 5

    # An upper certificate over the base field attains the same minmax.
    pieces = certificate['pieces']
    assert [sum(piece[j] for piece in pieces) % P for j in range(D)] == target.tolist()
    piece_ranks = []
    for w in pieces:
        pattern = tuple(dot(form, w) == 0 for form in forms)
        expected = ranks_by_pattern[pattern]
        actual = [rank(np.tensordot(w, b, axes=(0, 0)) % P) for b, _ in powers_models]
        assert actual == expected and actual[0]*actual[1] <= threshold
        piece_ranks.append(actual)

    # Exhaust all base-field directions and compare the algebraic oracle
    # with direct elimination on deterministic additional directions.
    all_digits = list(product(range(P), repeat=D))
    patterns = {tuple(dot(form, w) == 0 for form in forms) for w in all_digits}
    assert patterns <= ranks_by_pattern.keys()
    for w in random.sample(all_digits, 16):
        predicted = ranks_by_pattern[tuple(dot(form, w) == 0 for form in forms)]
        assert predicted == [rank(np.tensordot(w, b, axes=(0, 0)) % P)
                             for b, _ in powers_models]
    return {'status': 'pass', 'length': 8, 'physical_labels': 13,
            'coefficient_span_dimension': D, 'ambient_cut_dimension': N,
            'independent_magnetic_entry_controls': controls,
            'complete_eight_site_jet_identities': 3,
            'coefficient_rank_pairs': [list(pair) for pair in zip(*raw_pairs)],
            'primary_component_dimensions': [20, 1, 43, 70, 35],
            'double_component_nilpotent_rank': 13, 'kernel_flats': len(flats),
            'base_field_vectors': len(all_digits), 'base_field_zero_patterns': len(patterns),
            'exact_minmax_rank_product': threshold, 'separator': separator,
            'attaining_piece_rank_pairs': piece_ranks,
            'scope': 'Eight-site full coefficient span over every characteristic-five field extension; no uniform-length theorem or improved factoring bound.',
            'source_sha256': {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ('tools/bounded_spin_crossing.py',
                              'tools/check_bounded_spin_rank_algebra.py',
                              'evidence/bounded-spin-rank-model.json')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
