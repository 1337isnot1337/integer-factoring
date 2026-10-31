#!/usr/bin/env python3
"""Small binary operator sources from scalar Pascal contraction at n=2^k."""
from __future__ import annotations

from functools import lru_cache

from coupled_identity_operator_filters import _geometry, states
from weyl_lattice_filter_synthesis import polynomial_multiply


@lru_cache(maxsize=32)
def operator_basis(n, ramification=1, modes=3, geometry='split:1:2'):
    if n <= 0 or n & (n-1):
        raise ValueError('The scalar Pascal shortcut requires a power-of-two interval.')
    if modes not in (2, 3):
        raise ValueError('This synthesis prototype supports two or three binary modes.')
    data = _geometry(2, n, modes, geometry)
    alphabet = states(2, modes)
    selected = data['selected']
    powers = [1]
    scalar = 1 | (1 << (n*ramification))
    for exponent in range(modes):
        powers.append(polynomial_multiply(powers[-1], scalar))
    weights = tuple(ramification*sum(j*digit for j, digit in zip(selected, state)) for state in alphabet)
    gaussian = tuple(weight+(modes-sum(state))*n*ramification for state, weight in zip(alphabet, weights))
    lookup = {state: i for i, state in enumerate(alphabet)}
    rows, columns, blocks, orders, occupations = {}, {}, [], [], []
    offset = 0
    for occupation in range(modes+1):
        original = tuple(i for i, state in enumerate(alphabet) if sum(state) == occupation)
        ordered = sorted(original, key=lambda i: (gaussian[i], i))
        block_orders = tuple(gaussian[i] for i in ordered)
        for position, i in enumerate(ordered):
            rows[i], columns[i] = offset+position, position
        blocks.append({'occupation': occupation, 'indices': original,
                       'row_offset': offset, 'orders': block_orders})
        orders.extend(block_orders)
        occupations.extend([occupation]*len(original))
        offset += len(original)
    coefficients = {}
    for a in alphabet:
        for b in alphabet:
            if sum(a) != sum(b):
                continue
            entries = {}
            for source in alphabet:
                # Multiplication by Y^(1-b), followed by derivative^(1-a).
                if any(s and not digit for s, digit in zip(source, b)):
                    continue
                if any(not aa and bb and not s for aa, bb, s in zip(a, b, source)):
                    continue
                target = tuple(s+aa-bb for s, aa, bb in zip(source, a, b))
                assert all(digit in (0, 1) for digit in target)
                i, j = lookup[target], lookup[source]
                entries[rows[i], columns[j]] = powers[sum(b)] << weights[i]
            coefficients[a, b] = entries
    infinity = tuple(ramification*sum((2*j+1-2*n)*digit
                                      for j, digit in zip(selected, state)) for state in alphabet)
    return {'states': alphabet, 'selected': selected, 'retained': data['retained'],
            'coefficients': coefficients, 'gaussian_orders': tuple(orders),
            'row_occupations': tuple(occupations), 'infinity_orders': infinity,
            'blocks': tuple(blocks)}
