#!/usr/bin/env python3
"""Reconstruct and certify the PR57 two-place pencil's integer scan input.

Predecessor: PR57 at 8a3ea3e6aeacba7d0a53764de5579e3707bb2ca7.
The generator and exact factorization are retained with explicit provenance.
"""
from __future__ import annotations

from hashlib import sha256

from joint_energy_conjugations import N, binary_product, model
from search_f4_recentered_crossings import five_letter_source


def build_input():
    data = model()
    generator = sum(1 << (N*r+c) for pair in ((4, 20), (7, 11), (8, 16))
                    for r in pair for c in pair)
    allowed = {N*r+c for r, c in data['positions']}
    assert all(not (generator >> i & 1) or i in allowed for i in range(N*N))
    assert binary_product(generator, generator) == 0
    coefficients = [sum(((data['operator'][r][c][0] >> degree) & 1) << (N*r+c)
                        for r in range(N) for c in range(N)) for degree in range(19)]
    changes = [binary_product(generator, value) ^ binary_product(value, generator)
               for value in coefficients]
    assert all(not value for i, value in enumerate(changes) if i not in (8, 12, 16))
    assert changes[8] and changes[16] and changes[12] == changes[8] ^ changes[16]
    assert all(binary_product(binary_product(generator, value), generator) == 0
               for value in coefficients)
    # The weighted zero-cut perturbation is exactly t^8 times a constant.
    for change, degree_shift in ((changes[8], 0), (changes[16], 4)):
        assert all(data['energies'][position // N] + degree_shift == 8
                   for position in range(N*N) if change >> position & 1)
    assert not ((changes[8] | changes[16]) & 1)
    assert data['operator'][0][0][0] & (1 << 8)
    assert all(not (low << data['energies'][r]) & 255
               for r, row in enumerate(data['operator']) for low, high in row)
    source, heights, ex, ey = five_letter_source()
    parameters = [value for group in heights for value in group[1:]] + [1.]
    lines = ['PENCIL1']
    for index, (lo, hi) in enumerate(source):
        x, y, u, v = ((index // 5**i) % 5 for i in range(4))
        position = 25 * (x + 5*y) + v + 5*u
        lines.append(f'{lo} {hi} {(changes[8] >> position) & 1} {(changes[16] >> position) & 1}')
    lines.append(' '.join(map(str, parameters)))
    text = ('\n'.join(lines) + '\n').encode()
    certificate = {'predecessor_pr': 57,
                   'predecessor_commit': '8a3ea3e6aeacba7d0a53764de5579e3707bb2ca7',
                   'generator_square_zero': True, 'commutes_with_physical_energy': True,
                   'sandwich_coefficients_checked': len(coefficients),
                   'nonzero_commutator_degrees': [8, 12, 16],
                   'commutator_factor': 't^8*(t+1)^4*(K+t^4*L)',
                   'pencil': 'A_q=A+q(t)(K+t^4*L), deg(q)<=11, q(1)=0',
                   'coefficient_family': 'q=(1+t)r, deg(r)<=10 over F4',
                   'exact_number_of_coefficient_vectors': 4**11,
                   'physical_pair_energies': [list(ex), list(ey)],
                   'uniform_zero_minimum_order': 8,
                   'opposite_cut_total_height': 5 * sum(sum(group) for group in heights[:4]),
                   'sufficient_infinity_sum_for_asymmetric_improvement': 'strictly greater than 60',
                   'equivalent_determinant_degree_target': 'strictly less than 70',
                   'input_sha256': sha256(text).hexdigest()}
    return text, certificate
