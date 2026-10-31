#!/usr/bin/env python3
"""Native exact-field zero profiles for the binary coupled identity core.

This only changes the Smith backend. The literal active tensor comes
from coupled_identity_mode_filters. Every Laurent coefficient is packed
completely, embedding F2 in the first F4 bit plane. A missing promised
pivot increases the series window; truncated coefficient arrays are
never supplied to the native oracle.
"""
from __future__ import annotations

from coupled_identity_mode_filters import active_tensor, zero_charge_blocks
from native_matrix_smith_f4 import MatrixOracle


def zero_profiles(n, ramification, modes, alpha, tensor_builder=active_tensor,
                  **builder_options):
    tensor = tensor_builder(2, n, ramification, modes, alpha, **builder_options)
    return profiles_from_tensor(modes, tensor)


def profiles_from_tensor(modes, tensor):
    """Profile any complete binary active tensor in physical order V,Y,U,X."""
    result = []
    for charge, matrix in zero_charge_blocks(2, modes, tensor):
        packed, shifts, maximum_bits = [], [], 1
        for row in matrix:
            values, powers = [], []
            for entry in row:
                if entry:
                    first = min(entry)
                    bits = sum(1 << (power - first) for power, value in entry.items() if value % 2)
                    maximum_bits = max(maximum_bits, bits.bit_length())
                    values.append((bits, 0))
                    powers.append(first)
                else:
                    values.append((0, 0))
                    powers.append(0)
            packed.append(values)
            shifts.append(powers)
        precision = max(512, 1 << (maximum_bits - 1).bit_length())
        while precision <= 8192:
            try:
                oracle = MatrixOracle(packed, len(matrix), shifts, precision)
                orders = oracle.profile()['orders']
                break
            except ArithmeticError:
                precision *= 2
        else:
            raise ArithmeticError('The complete binary core needs a larger supported series window.')
        if any(value != round(value) for value in orders):
            raise ArithmeticError('An unweighted integer Laurent matrix returned a noninteger order.')
        result.append({'charge': charge, 'orders': tuple(int(value) for value in orders),
                       'exact_jet_precision': precision,
                       'backend': 'native F4 oracle on complete F2 coefficients',
                       'maximum_coefficient_bits': maximum_bits})
    return tuple(result)
