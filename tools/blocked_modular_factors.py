#!/usr/bin/env python3
"""Blocked half-rate crossing queries with contiguous binary channels.

The optional block_size argument is for finite construction checks. The
proved nearly linear-time family uses the default integer block schedule.
"""

from math import isqrt
import linear_modular_factors as local


def chosen_block_size(length):
    if length < 1:
        raise ValueError('A positive label length is required')
    logarithm = (length + 1).bit_length()
    root = isqrt(logarithm)
    root += root * root < logarithm
    logarithm2 = (logarithm + 1).bit_length()
    return min(length, root * logarithm2)


def block_layout(length, block_size=None):
    size = chosen_block_size(length) if block_size is None else min(length, block_size)
    if size < 1:
        raise ValueError('A positive block size is required')
    previous_size = previous_weight = None
    remaining = length
    while remaining:
        current = min(remaining, size)
        if current != previous_size:
            previous_weight = local.label_weight(current)
            previous_size = current
        weight = previous_weight
        yield current, weight, (4 * weight).bit_length()
        remaining -= current


def piece_code(length, shifts, block_size=None):
    """Construct one packed piece index, without a per-block pointer table."""
    code = bytearray()
    iterator = iter(shifts)
    for _, weight, width in block_layout(length, block_size):
        try:
            shift = next(iterator)
        except StopIteration as error:
            raise ValueError('Too few shifts') from error
        if not -2 * weight <= shift <= 2 * weight:
            raise ValueError('Shift outside the local piece range')
        value = shift + 2 * weight
        for bit in range(width - 1, -1, -1):
            code.append(48 + ((value >> bit) & 1))
    if next(iterator, None) is not None:
        raise ValueError('Too many shifts')
    return code.decode('ascii')


def decode_piece(length, code, block_size=None):
    position = 0
    for size, weight, width in block_layout(length, block_size):
        part = code[position:position + width]
        if len(part) != width or any(bit not in '01' for bit in part):
            raise ValueError('Malformed piece code')
        value = int(part, 2)
        if value > 4 * weight:
            raise ValueError('Padded piece code')
        yield size, weight, value - 2 * weight
        position += width
    if position != len(code):
        raise ValueError('Extra piece bits')


def next_piece(length, code, block_size=None):
    """Packed odometer; scans locate digits instead of storing offsets."""
    end = 0
    selected_start = selected_end = selected_value = None
    for _, weight, width in block_layout(length, block_size):
        value = int(code[end:end + width], 2)
        if value < 4 * weight:
            selected_start, selected_end, selected_value = end, end + width, value + 1
        elif value > 4 * weight:
            raise ValueError('Padded piece code')
        end += width
    if end != len(code):
        raise ValueError('Malformed piece length')
    if selected_start is None:
        return None
    width = selected_end - selected_start
    return (code[:selected_start] + format(selected_value, '0%db' % width)
            + '0' * (len(code) - selected_end))


def local_width(size, weight, shift, row_cut):
    count = local.channel_count(weight, shift, row_cut)
    width = 0 if count <= 1 else (count - 1).bit_length()
    return min(2 * size, width), count


def widths(length, code, block_size=None):
    row = column = 0
    for size, weight, shift in decode_piece(length, code, block_size):
        row += local_width(size, weight, shift, True)[0]
        column += local_width(size, weight, shift, False)[0]
    return row, column


def local_factor(size, weight, shift, row_cut, first_side, a, b, channel, count=None):
    """Row pairs are (x,y)|(v,u); column pairs are (x,u)|(v,y)."""
    if count is None:
        _, count = local_width(size, weight, shift, row_cut)
    if count > 1 << (2 * size):
        if first_side:
            return int(channel == (a << size) | b)
        first, second = channel >> size, channel & ((1 << size) - 1)
        x = local.partition_unrank(weight, first)
        v = local.partition_unrank(weight, a)
        if row_cut:
            y = local.partition_unrank(weight, second)
            u = local.partition_unrank(weight, b)
        else:
            u = local.partition_unrank(weight, second)
            y = local.partition_unrank(weight, b)
        return local.piece_entry(x, y, u, v, shift)
    if channel >= count:
        return 0
    a = local.partition_unrank(weight, a)
    b = local.partition_unrank(weight, b)
    angle, state = local.channel_unrank(weight, shift, row_cut, channel)
    if first_side:
        if row_cut:
            return local.top_factor(a, angle, state, b)
        return local.top_factor(a, angle, b, state)
    if row_cut:
        return local.bottom_factor(a, b, state)
    return local.bottom_factor(a, state, b)


def factor(length, code, row_cut, first_side, a_bits, b_bits, channel_bits, block_size=None):
    if len(a_bits) != length or len(b_bits) != length:
        raise ValueError('The external words have the wrong length')
    position = channel_position = 0
    value = 1
    for size, weight, shift in decode_piece(length, code, block_size):
        width, count = local_width(size, weight, shift, row_cut)
        a = int(a_bits[position:position + size], 2)
        b = int(b_bits[position:position + size], 2)
        chunk = channel_bits[channel_position:channel_position + width]
        if len(chunk) != width:
            raise ValueError('Channel word is too short')
        channel = int(chunk, 2) if width else 0
        if value:
            value &= local_factor(size, weight, shift, row_cut, first_side, a, b, channel, count)
        position += size
        channel_position += width
    if channel_position != len(channel_bits):
        raise ValueError('Channel word is too long')
    return value


def piece_entry(length, code, x_bits, y_bits, u_bits, v_bits, block_size=None):
    if any(len(bits) != length for bits in (x_bits, y_bits, u_bits, v_bits)):
        raise ValueError('The external words have the wrong length')
    position = 0
    value = 1
    for size, weight, shift in decode_piece(length, code, block_size):
        labels = tuple(local.partition_unrank(weight, int(bits[position:position + size], 2))
                       for bits in (x_bits, y_bits, u_bits, v_bits))
        if value:
            value &= local.piece_entry(*labels, shift)
        position += size
    return value
