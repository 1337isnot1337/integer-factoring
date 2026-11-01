#!/usr/bin/env python3
"""Exact exterior-operator, cut-rank, and occupation-graph checks."""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import random
import subprocess
import tempfile
import time

from check_matrix_unit_monodromy import rank, zip_bits


def states(D, r):
    return [sum(1 << a for a in subset) for subset in combinations(range(D), r)]


def exterior_operators(D, r):
    basis = states(D, r)
    indices = {state: i for i, state in enumerate(basis)}
    operators = {}
    for p, q in product(range(D), repeat=2):
        matrix = [[0]*len(basis) for _ in basis]
        # Column convention: annihilate p, then create q. Exterior signs
        # reduce to one in characteristic two.
        for column, state in enumerate(basis):
            if not (state >> p & 1):
                continue
            removed = state ^ (1 << p)
            if removed >> q & 1:
                continue
            matrix[indices[removed | (1 << q)]][column] = 1
        operators[p, q] = matrix
    return basis, indices, operators


def direct(t, p, q, initial, algebra):
    basis, indices, operators = algebra
    dp = [0]*len(basis)
    dp[indices[initial]] = 1
    for position in range(t):
        a = (p >> (2*position)) & 3
        b = (q >> (2*position)) & 3
        matrix = operators[a, b]
        following = dp.copy() if a == b else [0]*len(basis)
        for i in range(len(basis)):
            for j in range(len(basis)):
                if matrix[i][j]:
                    following[j] ^= dp[i] ^ (dp[i] << 1)
        dp = following
    return dp[indices[initial]]


def walk(t, p, q, initial):
    state, a, b = initial, 0, 0
    for position in range(t):
        u = (p >> (2*position)) & 3
        v = (q >> (2*position)) & 3
        if u == v:
            a += state >> u & 1
        elif not (state >> v & 1) or state >> u & 1:
            return 0
        else:
            state ^= (1 << u) | (1 << v)
            b += 1
    if state != initial:
        return 0
    result = 1 << a
    for _ in range(b):
        result ^= result << 1
    return result


def incoming_edges(d, S):
    D = d*d
    incoming = [0]*D
    for a, b in product(range(D), repeat=2):
        x, u, v, y = a % d, a // d, b % d, b // d
        p, q = x+d*y, v+d*u
        if a != b and S >> q & 1 and not (S >> p & 1):
            incoming[a] |= 1 << b
    return incoming


def topological_order(incoming, T):
    answer = []
    remaining = T
    while remaining:
        choices = [a for a in range(len(incoming))
                   if remaining >> a & 1 and not (incoming[a] & remaining)]
        if not choices:
            return None
        chosen = choices[0]
        answer.append(chosen)
        remaining ^= 1 << chosen
    return answer


def verify_topological(incoming, T, order):
    assert sorted(order) == [a for a in range(len(incoming)) if T >> a & 1]
    position = {a: i for i, a in enumerate(order)}
    for a, b in product(order, repeat=2):
        if incoming[a] >> b & 1:
            assert position[b] < position[a]


def ternary_witnesses():
    witnesses = []
    for subset in combinations(range(9), 4):
        S = sum(1 << a for a in subset)
        incoming = incoming_edges(3, S)
        found = None
        for keep_occupied in combinations(subset, 3):
            outside = [a for a in range(9) if not (S >> a & 1)]
            for keep_empty in combinations(outside, 4):
                T = sum(1 << a for a in keep_occupied+keep_empty)
                order = topological_order(incoming, T)
                if order is not None:
                    verify_topological(incoming, T, order)
                    found = {'occupation_mask': S, 'acyclic_mask': T,
                             'topological_order': order}
                    break
            if found is not None:
                break
        assert found is not None
        witnesses.append(found)
    assert len(witnesses) == comb(9, 4) == 126
    return witnesses


def full_profiles(t, S, algebra):
    n = 1 << t
    row = [[direct(t, p, q, S, algebra) for q in range(n*n)] for p in range(n*n)]
    column = [[direct(t, zip_bits(t, x, y), zip_bits(t, v, u), S, algebra)
               for v in range(n) for y in range(n)]
              for x in range(n) for u in range(n)]
    result = []
    for polynomials in (row, column):
        result.append([rank([sum(((value >> degree) & 1) << j
                                for j, value in enumerate(values))
                            for values in polynomials])
                       for degree in range(t+1)])
    return result


def verify_triangular(t, S, T, order):
    checks = 0
    row_letters = [a for a in range(4) if not (S >> a & 1)]
    row_letters += [a for a in range(4) if S >> a & 1]
    words = list(product(row_letters, repeat=t))
    encoded = [sum(a << (2*i) for i, a in enumerate(word)) for word in words]
    for i, p in enumerate(words):
        for j, q in enumerate(words):
            value = walk(t, encoded[i], encoded[j], S)
            if i == j:
                assert value == 1 << sum(S >> a & 1 for a in p)
            elif i > j:
                assert value == 0
            checks += 1
    words = list(product(order, repeat=t))
    for i, a in enumerate(words):
        x = sum((v & 1) << j for j, v in enumerate(a))
        u = sum((v >> 1) << j for j, v in enumerate(a))
        for j, b in enumerate(words):
            v = sum((w & 1) << k for k, w in enumerate(b))
            y = sum((w >> 1) << k for k, w in enumerate(b))
            value = walk(t, zip_bits(t, x, y), zip_bits(t, v, u), S)
            if i == j:
                assert value == 1 << sum(S >> letter & 1 for letter in a)
            elif i < j:
                assert value == 0
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--native', type=Path, default=None,
                        help='Use an existing native executable; otherwise compile the sibling C++ source in a temporary directory.')
    args = parser.parse_args()
    started = time.monotonic()
    native_workspace = tempfile.TemporaryDirectory(prefix='exterior-monodromy-')
    native_directory = Path(native_workspace.name)
    if args.native is None:
        native_path = native_directory / 'matrix_unit_monodromy'
        subprocess.run(['c++', '-O2', '-std=c++17',
                        str(Path(__file__).with_name('matrix_unit_monodromy.cpp')),
                        '-o', str(native_path)], check=True)
    else:
        native_path = args.native.resolve()
    rng = random.Random(96271)
    algebras = {r: exterior_operators(4, r) for r in range(1, 4)}
    entry_checks = endpoint_checks = profile_checks = triangular_checks = 0
    for S in range(1, 15):
        algebra = algebras[S.bit_count()]
        for t in range(1, 4):
            pairs = list(product(range(4**t), repeat=2)) if t <= 2 else [
                (rng.randrange(4**t), rng.randrange(4**t)) for _ in range(128)]
            for p, q in pairs:
                value = direct(t, p, q, S, algebra)
                assert value == walk(t, p, q, S)
                assert value.bit_count() % 2 == (p == q)
                entry_checks += 1
                endpoint_checks += 1
        incoming = incoming_edges(2, S)
        if S in (6, 9):
            T = next(T for T in range(16) if T.bit_count() == 3
                     and topological_order(incoming, T) is not None)
        else:
            T = 15
        order = topological_order(incoming, T)
        assert order is not None
        verify_topological(incoming, T, order)
        triangular_checks += verify_triangular(2, S, T, order)
        if S in (3, 9):
            triangular_checks += verify_triangular(3, S, T, order)
            for t in range(1, 4):
                rr, cc = full_profiles(t, S, algebra)
                temporary = native_directory / f'profile-t{t}-s{S}.json'
                native = json.loads(subprocess.check_output(
                    [str(native_path), str(t), str(temporary), str(S), 'particles'], text=True))
                assert rr == native['row_ranks'] and cc == native['column_ranks']
                profile_checks += 2
    witnesses = ternary_witnesses()
    # The finite binary seven-site rejection uses an integer square.
    lower_square = (comb(7, 3)*2**7)*(comb(7, 3)*2**4)
    assert lower_square == 2508800 > 2**20
    output = {'status': 'passed', 'exterior_operator_walk_comparisons': entry_checks,
              'endpoint_checks': endpoint_checks, 'full_cut_rank_profile_comparisons': profile_checks,
              'triangular_entry_checks': triangular_checks,
              'binary_seven_site_lower_square': lower_square,
              'ternary_four_particle_occupations_checked': len(witnesses),
              'ternary_witnesses': witnesses,
              'scope': 'The exact witnesses and arithmetic support the stated finite classes; all-length rank bounds follow from the companion triangular proof.',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seconds': time.monotonic()-started}
    print(json.dumps({key: value for key, value in output.items() if key != 'ternary_witnesses'}, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n')
    native_workspace.cleanup()


if __name__ == '__main__':
    main()
