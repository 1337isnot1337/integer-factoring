#!/usr/bin/env python3
"""Exact controls for the all-field Gaussian partner-space injectivity proof."""

import argparse
import json
from pathlib import Path

import numpy as np

from coherent_physical_tail_core import finite_field_images
from coherent_physical_tail_low_band import MUL, mul, power
from check_coherent_physical_tail_core import inverse, matrix_product, compound, action
from check_cancellation_module_automata import pmul


def polynomial_matrix_product(A, B):
    result = [[0] * len(B[0]) for _ in A]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for c in range(len(B)):
                result[i][j] ^= pmul(A[i][c], B[c][j])
    return result


def exact_digit_identities():
    checked = 0
    for k in (1, 2, 4, 8):
        F = [[3 << (2 * (a & ~b)) for b in range(k)] for a in range(k)]
        H = [[3 << (2 * (a & b)) for b in range(k)] for a in range(k)]
        J = [[int(a + b == k - 1) for b in range(k)] for a in range(k)]
        assert polynomial_matrix_product(F, J) == H
        expected = [[(1 ^ (1 << (2 * k))) * J[a][b] for b in range(k)] for a in range(k)]
        assert polynomial_matrix_product(F, H) == expected
        assert polynomial_matrix_product(H, list(map(list, zip(*F)))) == expected
        checked += 3 * k * k
    return checked


def data(k, t, side):
    n = 6 * k
    x, sigma = power(t, 2 * k), power(1 ^ power(t, 2 * k), -1)
    lam = mul(x, sigma)
    F = np.array([[mul(1 ^ t, power(t, 2 * (a & ~b))) for b in range(k)]
                  for a in range(k)], dtype=np.uint8)
    H = np.array([[mul(1 ^ t, power(t, 2 * (a & b))) for b in range(k)]
                  for a in range(k)], dtype=np.uint8)
    J = np.array([[int(a + b == k - 1) for b in range(k)] for a in range(k)], dtype=np.uint8)
    # Old and new creation coordinates use block order h,p,q,r,v,s
    # and u,P,Q,R,v,S. On the right R means the combined P+r.
    C = np.eye(n, dtype=np.uint8)
    C[4 * k:5 * k, :k] = F
    C[5 * k:6 * k, :k] = MUL[x if side == 'row' else lam, F]
    C[4 * k:5 * k, k:2 * k] = MUL[sigma if side == 'row' else 1, F]
    C[4 * k:5 * k, 3 * k:4 * k] = F
    if side == 'column':
        C[k:2 * k, 3 * k:4 * k] = np.eye(k, dtype=np.uint8)
    U = inverse(C)
    groups = finite_field_images(k, t, 0, mul, power)
    groups = groups[:2] if side == 'row' else groups[2:]
    operators = []
    for group in groups:
        for creator, annihilator in group:
            if side == 'column':
                creator, annihilator = ((annihilator[1], annihilator[0]),
                                        (creator[1], creator[0]))
            assert all(not any(vector[n:]) for op in (creator, annihilator) for vector in op)
            operators.append((creator, annihilator))
    return C, U, operators, F, H, J, x, sigma, lam


def coordinate_controls():
    checked = 0
    for k in (1, 2, 4, 8):
        for t in (2, 3):
            for side in ('row', 'column'):
                C, U, operators, F, H, J, x, sigma, lam = data(k, t, side)
                n = 6 * k
                for physical, (creator, annihilator) in enumerate(operators):
                    group, a = divmod(physical, k)
                    target = a + 4 * k * group
                    for which, op in enumerate((creator, annihilator)):
                        actual_c = matrix_product(U, np.array(op[0][:n], dtype=np.uint8)[:, None])[:, 0]
                        actual_a = matrix_product(np.array(op[1][:n], dtype=np.uint8)[None, :], C)[0]
                        expected_c, expected_a = np.zeros(n, dtype=np.uint8), np.zeros(n, dtype=np.uint8)
                        if which == 0:
                            expected_c[target] = 1
                        else:
                            expected_a[target] = 1
                            if group == 0:
                                expected_c[4 * k:5 * k] = H[:, a]
                                expected_c[3 * k:4 * k] = J[:, a]
                            else:
                                expected_c[:k] = H[:, a]
                                expected_c[2 * k:3 * k] = MUL[1 if side == 'row' else sigma, H[:, a]]
                                expected_c[3 * k:4 * k] = H[:, a]
                                expected_c[5 * k:6 * k] = MUL[mul(x, 1 ^ x) if side == 'row' else x, J[:, a]]
                        assert np.array_equal(actual_c, expected_c)
                        assert np.array_equal(actual_a, expected_a)
                        checked += 2 * n
                # The original leading logical creators are h and p.
                for a in range(2 * k):
                    actual = U[:, a]
                    expected = np.zeros(n, dtype=np.uint8)
                    if a < k:
                        expected[a] = 1
                        expected[4 * k:5 * k] = F[:, a]
                        expected[5 * k:6 * k] = MUL[x if side == 'row' else lam, F[:, a]]
                    else:
                        expected[a] = 1
                        expected[4 * k:5 * k] = MUL[sigma if side == 'row' else 1, F[:, a - k]]
                    assert np.array_equal(actual, expected)
                    checked += n
    return checked


def encoded_matrix(k, t, side):
    """After the partner basis change: hot z, partner b, free f."""
    _, _, _, F, H, J, x, sigma, lam = data(k, t, side)
    n, logical = 6 * k, 2 * k
    G = np.zeros((1 << n, 1), dtype=np.uint8)
    G[0, 0] = 1
    zero = [0] * n

    def add_pair(i, j, coefficient):
        nonlocal G
        if not coefficient:
            return
        a, b = zero.copy(), zero.copy()
        a[i] = b[j] = 1
        G ^= MUL[int(coefficient), action(action(G, (a, zero)), (b, zero))]

    for i in range(logical):
        add_pair(i, logical + i, 1)
    for a in range(k):
        for b in range(k):
            add_pair(a, k + b, H[a, b])
            add_pair(k + a, 5 * k + b,
                     mul(mul(x, 1 ^ x) if side == 'row' else x, int(J[a, b])))
    forms = np.zeros((n, logical), dtype=np.uint8)
    forms[:k, :k] = np.eye(k, dtype=np.uint8)
    forms[k:2 * k, :k] = F
    forms[5 * k:6 * k, :k] = MUL[x if side == 'row' else lam, F]
    forms[4 * k:5 * k, k:2 * k] = np.eye(k, dtype=np.uint8)
    forms[k:2 * k, k:2 * k] = MUL[sigma if side == 'row' else 1, F]
    output = np.zeros((1 << n, 1 << logical), dtype=np.uint8)
    output[:, 0] = G[:, 0]
    for state in range(1, 1 << logical):
        bit = state & -state
        output[:, state:state + 1] = action(output[:, state ^ bit: (state ^ bit) + 1],
                                            (forms[:, bit.bit_length() - 1], zero))
    free = forms[4 * k:6 * k]
    assert np.array_equal(matrix_product(free, inverse(free)), np.eye(logical, dtype=np.uint8))
    return output, compound(free)


def highest_degree_controls():
    checked, cases = 0, 0
    for k in (1, 2):
        d = 2 * k
        for side in ('row', 'column'):
            E, free = encoded_matrix(k, 2, side)
            E = E.reshape(1 << (4 * k), 1 << d, 1 << d)
            for degree in range(d + 1):
                functional = np.array([((37 * i + 19) % 255) + 1 if i.bit_count() <= degree else 0
                                       for i in range(1 << d)], dtype=np.uint8)
                chosen = next(i for i in range(1 << d) if i.bit_count() == degree)
                conditional = np.bitwise_xor.reduce(MUL[functional[None, :, None], E], axis=1)
                selected = conditional[[chosen | (state << d) for state in range(1 << d)]]
                assert np.array_equal(selected, MUL[functional[chosen], free])
                checked += selected.size
                cases += 1
    return dict(joint_functional_cases=cases, exact_partner_coefficient_entries=checked)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='passed', factoring_exponent_improved=False,
                  exact_polynomial_digit_entries=exact_digit_identities(),
                  exact_coordinate_entries=coordinate_controls(),
                  highest_degree=highest_degree_controls(),
                  scope='Controls for the explicit all-field injectivity proof; no complete high-band profiles or all-scale pressure floor.')
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
