#!/usr/bin/env python3
"""Exact boundary-stable matrix seeds, endpoint automata, and CAR halves."""
from __future__ import annotations

import argparse
from itertools import product
import hashlib
import json
from pathlib import Path
import random
import time


def identity(n):
    return tuple(1 << i for i in range(n))


def zero(n):
    return (0,)*n


def add(A, B):
    return tuple(a ^ b for a, b in zip(A, B))


def multiply(A, B):
    rows = []
    for row in A:
        value = 0
        while row:
            bit = row & -row
            value ^= B[bit.bit_length()-1]
            row ^= bit
        rows.append(value)
    return tuple(rows)


def apply(A, v):
    return sum(((row & v).bit_count() & 1) << i for i, row in enumerate(A))


def flat(A):
    n = len(A)
    return sum(row << (n*i) for i, row in enumerate(A))


def rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            top = value.bit_length()-1
            if top not in pivots:
                pivots[top] = value
                break
            value ^= pivots[top]
    return len(pivots)


def unit(n, i, j):
    return tuple((1 << j) if k == i else 0 for k in range(n))


def tensor(A, B):
    m, n = len(A), len(B)
    rows = []
    for i, k in product(range(m), range(n)):
        rows.append(sum(((A[i] >> j) & 1)*((B[k] >> l) & 1) << (n*j+l)
                        for j, l in product(range(m), range(n))))
    return tuple(rows)


def trace(A):
    return sum(row >> i & 1 for i, row in enumerate(A)) & 1


def rank_factor(rows):
    basis = []
    pivots = {}
    for original in rows:
        value = original
        while value:
            top = value.bit_length()-1
            if top not in pivots:
                pivots[top] = len(basis)
                basis.append(value)
                break
            value ^= basis[pivots[top]]
    coordinates = []
    for original in rows:
        value, coordinate = original, 0
        while value:
            j = pivots[value.bit_length()-1]
            coordinate ^= 1 << j
            value ^= basis[j]
        coordinates.append(coordinate)
    return coordinates, basis


def cut_factors(table, column=False):
    rows = []
    for left0, left1, a, b in product(range(2), repeat=4):
        value = 0
        for right0, right1, c, d in product(range(2), repeat=4):
            quad = ((left0, right0, left1, right1) if column
                    else (left0, left1, right0, right1))
            bit = table[quad][2*a+c] >> (2*b+d) & 1
            value |= bit << (8*right0+4*right1+2*c+d)
        rows.append(value)
    U, V = rank_factor(rows)
    r = len(V)
    left, right = [], []
    for external in range(4):
        left.append([tensor(tuple(sum(((U[4*external+2*a+b] >> s) & 1) << b
                                      for b in range(2)) for a in range(2)), identity(2))
                     for s in range(r)])
        right.append([tensor(identity(2), tuple(sum(((V[s] >> (4*external+2*c+d)) & 1) << d
                                                   for d in range(2)) for c in range(2)))
                      for s in range(r)])
    return r, left, right


def family():
    I = identity(4)
    noise = []
    for a in range(4):
        noise.extend([add(unit(4, a, 0), unit(4, a, 3)),
                      unit(4, a, 1), unit(4, a, 2)])
    C = {}
    index = 0
    for q in product(range(2), repeat=4):
        x, y, u, w = q
        if x == w and y == u:
            C[q] = I
        else:
            C[q] = noise[index]
            index += 1
    assert index == 12
    Q = unit(4, 1, 0)
    return C, [{q: Q for q in C}, {q: add(C[q], Q) for q in C}], noise


def endpoint_reachable(C, eta, phi, values):
    b = len(next(iter(C.values())))
    generators = {}
    for q, M in C.items():
        generators[q] = tuple(M)+(values[q] << b,)
    basis = []
    pivots = {}
    def insert(value):
        while value:
            p = value.bit_length()-1
            if p not in pivots:
                pivots[p] = value
                basis.append(value)
                return
            value ^= pivots[p]
    insert(eta | (1 << b))
    current = 0
    while current < len(basis):
        v = basis[current]
        for G in generators.values():
            insert(apply(G, v))
        current += 1
    annihilator = phi | (1 << b)
    return all((annihilator & v).bit_count() % 2 == 0 for v in basis), len(basis)


def car_operators(m):
    n = 1 << m
    annihilators, creators = [], []
    for i in range(m):
        a, c = [0]*n, [0]*n
        for state in range(n):
            if state >> i & 1:
                a[state ^ (1 << i)] |= 1 << state
            else:
                c[state | (1 << i)] |= 1 << state
        annihilators.append(tuple(a))
        creators.append(tuple(c))
    return annihilators, creators


def left_right(A):
    n = len(A)
    L, R = [0]*(n*n), [0]*(n*n)
    for i, j in product(range(n), repeat=2):
        column = n*i+j
        for k in range(n):
            if A[k] >> i & 1:
                L[n*k+j] |= 1 << column
            if A[j] >> k & 1:
                R[n*i+k] |= 1 << column
    return tuple(L), tuple(R)


def normal_basis(a, c):
    m, b = len(a), len(a[0])
    answer = []
    for P, Q in product(range(1 << m), repeat=2):
        value = identity(b)
        for i in range(m):
            if P >> i & 1:
                value = multiply(value, c[i])
        for i in range(m):
            if Q >> i & 1:
                value = multiply(value, a[i])
        answer.append(value)
    return answer


def check_car(m):
    a, c = car_operators(m)
    n, b = 1 << m, 1 << (2*m)
    I, O = identity(b), zero(b)
    eta = sum(1 << (n*i+i) for i in range(n))
    ap, cp, aq, cq = [], [], [], []
    for ai, ci in zip(a, c):
        La, Ra = left_right(ai)
        Lc, Rc = left_right(ci)
        ap.append(add(La, Ra)); cp.append(Lc)
        aq.append(add(Lc, Rc)); cq.append(Ra)
    relation_checks = 0
    for A, C in ((ap, cp), (aq, cq)):
        for i in range(m):
            assert multiply(A[i], A[i]) == O and multiply(C[i], C[i]) == O
            assert apply(A[i], eta) == 0
            relation_checks += 3
            for j in range(m):
                assert add(multiply(A[i], C[j]), multiply(C[j], A[i])) == (I if i == j else O)
                assert multiply(A[i], A[j]) == multiply(A[j], A[i])
                assert multiply(C[i], C[j]) == multiply(C[j], C[i])
                relation_checks += 3
    for U, V in product(ap+cp, aq+cq):
        assert multiply(U, V) == multiply(V, U)
        relation_checks += 1
    U, V = normal_basis(ap, cp), normal_basis(aq, cq)
    assert rank(map(flat, U)) == rank(map(flat, V)) == n*n
    assert rank(flat(multiply(u, v)) for u, v in product(U, V)) == b*b
    new_terminal = [sum((apply(multiply(u, v), eta) & 1) << j for j, v in enumerate(V)) for u in U]
    assert rank(new_terminal) == 1
    units = [unit(n, i, j) for i, j in product(range(n), repeat=2)]
    old_U = [left_right(T)[0] for T in units]
    old_V = [left_right(T)[1] for T in units]
    old_terminal = [sum((apply(multiply(u, v), eta) & 1) << j for j, v in enumerate(old_V)) for u in old_U]
    assert rank(old_terminal) == n
    wick_checks = 0
    for x, y, u, v in product(range(m), repeat=4):
        M = multiply(multiply(multiply(ap[x], ap[y]), cp[u]), cp[v])
        expected = ((x == u and y == v) ^ (x == v and y == u))
        assert apply(M, eta) == (eta if expected else 0)
        wick_checks += 1
    interleaved = multiply(multiply(multiply(ap[0], cp[0]), ap[0]), cp[0])
    gathered = multiply(multiply(ap[0], ap[0]), multiply(cp[0], cp[0]))
    assert apply(interleaved, eta) & 1 == 1
    assert apply(gathered, eta) == 0
    return {'modes': m, 'fock_dimension': n, 'auxiliary_dimension': b,
            'car_and_cross_commutation_checks': relation_checks,
            'wick_vector_identities': wick_checks,
            'original_terminal_bilinear_rank': n, 'new_terminal_bilinear_rank': 1,
            'invalid_car_gather_rejected': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    C, pieces, noise = family()
    eta, phi = 9, 1
    inputs = list(C)
    values = {q: int(q[0] == q[3] and q[1] == q[2]) for q in C}
    assert all(apply(C[q], eta) == values[q]*eta for q in C)
    assert apply(pieces[0][inputs[0]], eta) not in (0, eta)
    assert any(multiply(pieces[0][q], pieces[1][r]) != multiply(pieces[1][r], pieces[0][q])
               for q, r in product(inputs, repeat=2))
    factors = [[cut_factors(piece, column) for piece in pieces] for column in (False, True)]
    assert [[f[0] for f in cut] for cut in factors] == [[1, 12], [1, 12]]
    assert [cut_factors(C, column)[0] for column in (False, True)] == [11, 11]
    factor_checks = commutation_checks = 0
    for column, cut in enumerate(factors):
        left_values, right_values = set(), set()
        for piece, (r, L, R) in zip(pieces, cut):
            for q in inputs:
                x, y, u, w = q
                left = 2*x+(u if column else y)
                right = 2*(y if column else u)+w
                actual = zero(4)
                for s in range(r):
                    actual = add(actual, multiply(L[left][s], R[right][s]))
                assert actual == piece[q]
                factor_checks += 1
            left_values.update(value for row in L for value in row)
            right_values.update(value for row in R for value in row)
        for L, R in product(left_values, right_values):
            assert multiply(L, R) == multiply(R, L)
            commutation_checks += 1
    rng = random.Random(784392)
    gathering_checks = 0
    for column, cut in enumerate(factors):
        for k in (2, 3):
            for _ in range(64):
                all_left, all_right, interleaved = identity(4), identity(4), identity(4)
                for j in range(k):
                    piece = rng.randrange(2)
                    r, L, R = cut[piece]
                    s = rng.randrange(r)
                    left, right = L[rng.randrange(4)][s], R[rng.randrange(4)][s]
                    all_left = multiply(all_left, left)
                    all_right = multiply(all_right, right)
                    interleaved = multiply(interleaved, multiply(left, right))
                assert interleaved == multiply(all_left, all_right)
                scalar = apply(interleaved, eta) & 1
                assert scalar == ((all_left[0] & apply(all_right, eta)).bit_count() & 1)
                gathering_checks += 1
    amplified_entries = type_comparisons = 0
    for k in (1, 2, 3):
        for indices in product(range(16), repeat=k):
            dp = [identity(4)]
            target = 1
            for index in indices:
                q = inputs[index]
                following = [zero(4)]*(len(dp)+1)
                for degree, M in enumerate(dp):
                    following[degree] = add(following[degree], multiply(M, pieces[0][q]))
                    following[degree+1] = add(following[degree+1], multiply(M, pieces[1][q]))
                dp = following
                target *= values[q]
            total = 0
            for M in dp:
                total ^= apply(M, eta) & 1
            assert total == target
            amplified_entries += 1
            if k <= 2 or indices[0] == indices[1] == indices[2]:
                direct = [zero(4)]*(k+1)
                for word in product(range(2), repeat=k):
                    M = identity(4)
                    for index, piece in zip(indices, word):
                        M = multiply(M, pieces[piece][inputs[index]])
                    direct[sum(word)] = add(direct[sum(word)], M)
                assert direct == dp
                type_comparisons += k+1
    densities = 0
    two_by_two = [tuple((mask >> (2*i)) & 3 for i in range(2)) for mask in range(16)]
    for A, B in product(two_by_two, repeat=2):
        if trace(A) != 1 or trace(B) != 1:
            continue
        D = tensor(A, B)
        assert trace(D) == 1
        assert any(trace(multiply(D, N)) for N in noise)
        densities += 1
    assert densities == 64
    accepted, reach_dimension = endpoint_reachable(C, eta, phi, values)
    assert accepted
    # A left observable quotient can pass without preserving the supplied
    # right boundary.  These two generators have endpoint values 1 and 0.
    quotient = {0: (0, 2), 1: (3, 0)}
    accepted_quotient, quotient_dimension = endpoint_reachable(quotient, 2, 2, {0: 1, 1: 0})
    assert accepted_quotient and apply(quotient[1], 2) != 0
    # Both one-letter scalar values vanish, but E01 E10 returns to eta.
    bad = {0: unit(2, 0, 1), 1: unit(2, 1, 0)}
    assert all(apply(M, 1) & 1 == 0 for M in bad.values())
    accepted_bad, bad_dimension = endpoint_reachable(bad, 1, 1, {0: 0, 1: 0})
    assert not accepted_bad
    car = [check_car(m) for m in (1, 2)]
    output = {'status': 'passed', 'exact_cut_factor_identities': factor_checks,
              'one_piece_extended_cut_profile': [11, 11],
              'two_piece_extended_cut_profiles': [[1, 1], [12, 12]],
              'global_cross_half_commutators': commutation_checks,
              'gathered_channel_word_checks': gathering_checks,
              'amplified_endpoint_entries': amplified_entries,
              'ordered_type_coefficient_comparisons': type_comparisons,
              'normalized_product_densities_rejected': densities,
              'endpoint_reachable_dimension': reach_dimension,
              'observable_quotient_reachable_dimension': quotient_dimension,
              'length_one_only_negative_control_reachable_dimension': bad_dimension,
              'car_structures': car,
              'scope': 'Exact constructive interface and noncollapse control; the toy rank sum is thirteen and does not improve the factoring score.',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seconds': time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
