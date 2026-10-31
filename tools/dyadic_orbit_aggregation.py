#!/usr/bin/env python3
"""Polynomial-cost interval aggregation over complete dyadic unipotent orbits."""

from dyadic_hyperbola import affine_interval_count, residue_interval


def cyclic_pieces(start, length, modulus):
    if not length:
        return ()
    start %= modulus
    stop = start + length
    if stop <= modulus:
        return ((start, stop),)
    return ((start, modulus), (0, stop - modulus))


def affine_residue_pairs(modulus, slope, offset, x_interval, y_interval):
    """Count pairs in integer intervals with y=slope*x+offset mod modulus."""
    assert modulus >= 1 and modulus & (modulus - 1) == 0 and slope % 2
    x0, x1 = x_interval
    y0, y1 = y_interval
    assert 0 <= x0 <= x1 and 0 <= y0 <= y1
    qx, rx = divmod(x1 - x0, modulus)
    qy, ry = divmod(y1 - y0, modulus)
    total = qx * qy * modulus + qx * ry + qy * rx
    for start, stop in cyclic_pieces(x0, rx, modulus):
        for lower, upper in cyclic_pieces(y0, ry, modulus):
            total += affine_interval_count(modulus, slope, offset,
                                            start, stop, lower, upper)
    return total


def unipotent_interval_sum(bits, depth, x_interval, y_interval,
                            output_scale=1, input_scale=1, input_shift=0,
                            output_shift=0):
    """Sum rectangle counts for alpha*L_k(beta*t+gamma)+delta, over all k.

    L_k(z)=z/(1+2^depth*k*z) mod 2^bits, 0<=k<2^(bits-depth).
    Both scales are odd. No group parameter or interval element is visited.
    """
    assert 1 <= depth <= bits and output_scale % 2 and input_scale % 2
    modulus = 1 << bits
    x0, x1 = x_interval
    y0, y1 = y_interval
    assert 0 <= x0 <= x1 <= modulus and 0 <= y0 <= y1 <= modulus
    size = 1 << (bits - depth)
    zero_input = -input_shift * pow(input_scale, -1, modulus) % modulus
    zero_output = output_shift % modulus
    total = size * (x0 <= zero_input < x1 and y0 <= zero_output < y1)
    slope = output_scale * input_scale
    offset = output_scale * input_shift + output_shift
    for valuation in range(bits):
        stride = 1 << (valuation + 1)
        root_x = ((1 << valuation) - input_shift) * pow(input_scale, -1, stride) % stride
        root_y = (slope * root_x + offset) % stride
        x_reduced = residue_interval(x0, x1, root_x, stride, modulus // stride)
        y_reduced = residue_interval(y0, y1, root_y, stride, modulus // stride)
        precision = min(bits, depth + 2 * valuation)
        pair_modulus = 1 << (precision - valuation - 1)
        pair_offset = (slope * root_x + offset - root_y) // stride
        multiplicity = 1 << min(bits - depth, 2 * valuation)
        total += multiplicity * affine_residue_pairs(pair_modulus, slope, pair_offset,
                                                       x_reduced, y_reduced)
    return total
