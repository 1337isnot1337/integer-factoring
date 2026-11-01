#!/usr/bin/env python3
"""Complete natural profiles of bit-deletion downsets via charge blocks.

The selected one-particle variable is t=w^2. Every returned order is
converted back to the original w valuation. Endpoint permutations prove
full generic rank separately in every charge block. A short series
window raises an error instead of returning an incomplete profile.
This module does not support arbitrary physical gauges.
"""
from __future__ import annotations

import ctypes
from functools import lru_cache
import hashlib
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import time


@lru_cache(None)
def _library(words):
    source = Path(__file__).with_suffix('.cpp')
    dependencies = [source, source.with_name('selected_gaussian_coefficients.cpp'),
                    source.with_name('native_rectangular_smith_f4.cpp')]
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The charge-block constructor requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={words}']
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = hashlib.sha256(b''.join(path.read_bytes() for path in dependencies)
                            + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-hot-compounds-' + str(os.getuid()))
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = cache / (digest + '.so')
    if not target.exists():
        fd, temporary = tempfile.mkstemp(prefix=digest + '-', suffix='.so', dir=cache)
        os.close(fd)
        try:
            run = subprocess.run(compiler + flags + [str(source), '-o', temporary],
                                 capture_output=True, text=True)
            if run.returncode:
                raise RuntimeError(run.stderr[-12000:])
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    library = ctypes.CDLL(str(target))
    u64, u32 = ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_uint32)
    i64 = ctypes.POINTER(ctypes.c_int64)
    library.hot_code_compound_block.argtypes = [ctypes.c_int, ctypes.c_int, u64, ctypes.c_int,
                                                u32, u32, i64, u64]
    library.hot_code_compound_block.restype = ctypes.c_int
    return library


def one_particle(modes):
    """B(t)=[[D,J],[J,D]], D_ab=(1+t)t^(2(a|b)+1), J_aa=t^(2a+1)."""
    if not isinstance(modes, int) or not 1 <= modes <= 16:
        raise ValueError('Use one through sixteen hot logical modes.')
    result = []
    for i in range(2 * modes):
        side, a = divmod(i, modes)
        row = []
        for j in range(2 * modes):
            other, b = divmod(j, modes)
            row.append(3 << (2 * (a | b) + 1) if side == other
                       else (1 << (2 * a + 1) if a == b else 0))
        result.append(row)
    return result


def validate_code(letters, modes):
    letters = tuple(letters)
    if (not letters or len(set(letters)) != len(letters)
            or any(not isinstance(x, int) or not 0 <= x < 1 << modes for x in letters)):
        raise ValueError('Use a nonempty set of distinct hot physical masks.')
    chosen = set(letters)
    for value in letters:
        remaining = value
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            if value ^ bit not in chosen:
                raise ValueError('The physical code must be closed under deleting occupied bits.')
    return letters


def hot_energy(mask):
    return sum(8 * a + 4 for a in range(mask.bit_length()) if mask >> a & 1)


def compound_block(matrix, row_masks, column_masks, *, coefficients=False, precision=None):
    """Return a complete selected-minor block and its ordinary t-orders.

    The endpoint must exchange the two mode sets. If coefficients=False,
    only the complete orders cross the native interface.
    """
    n = len(matrix)
    if n % 2 or not 2 <= n <= 32 or any(len(row) != n for row in matrix):
        raise ValueError('Use a square binary polynomial matrix of even order 2 through 32.')
    if any(not isinstance(value, int) or value < 0 for row in matrix for value in row):
        raise ValueError('Binary polynomial bit masks must be nonnegative integers.')
    modes = n // 2
    rows, columns = tuple(row_masks), tuple(column_masks)
    if (not rows or len(rows) != len(columns) or len(rows) > 4096
            or len(set(rows)) != len(rows) or len(set(columns)) != len(columns)
            or any(not isinstance(mask, int) or not 0 <= mask < 1 << n for mask in rows + columns)
            or len({mask.bit_count() for mask in rows + columns}) != 1):
        raise ValueError('Supply distinct equally charged row and column sets of equal size, at most 4096.')
    charge, count = rows[0].bit_count(), len(rows)
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    coefficient_bound = charge * max(0, degree)
    words = 2
    while 64 * words <= coefficient_bound:
        words *= 2
    if precision is not None:
        if (not isinstance(precision, int) or precision < 128 or precision > 8192
                or precision & (precision - 1) or precision <= coefficient_bound):
            raise ValueError('Precision must be a power of two from 128 through 8192 above the coefficient bound.')
        words = precision // 64
    input_words = max(1, (degree + 64) // 64)
    mask = (1 << 64) - 1
    packed = (ctypes.c_uint64 * (n * n * input_words))(
        *(value >> (64 * k) & mask for row in matrix for value in row for k in range(input_words)))
    rr, cc = (ctypes.c_uint32 * count)(*rows), (ctypes.c_uint32 * count)(*columns)
    attempts = []
    while words <= 128:
        output = (ctypes.c_int64 * count)()
        values = (ctypes.c_uint64 * (count * count * words))() if coefficients else None
        status = _library(words).hot_code_compound_block(modes, count, packed, input_words,
                                                         rr, cc, output, values)
        attempts.append({'precision': 64 * words, 'status': status})
        if status == 0:
            break
        if status != 6 or precision is not None:
            raise ArithmeticError(f'The complete charge block failed (status {status}, precision {64 * words}).')
        words *= 2
    else:
        raise ArithmeticError('The complete charge-block profile exceeded the series budget.')
    result = {'charge': charge, 'size': count, 'orders_t': sorted(map(int, output)),
              'coefficient_degree_bound_t': coefficient_bound, 'precision_t': 64 * words,
              'precision_attempts': attempts}
    if coefficients:
        result['matrix'] = [[sum(int(values[(i * count + j) * words + k]) << (64 * k)
                                  for k in range(words)) for j in range(count)] for i in range(count)]
        if max(value.bit_length() - 1 for row in result['matrix'] for value in row) > coefficient_bound:
            raise ArithmeticError('A complete determinant exceeded its proved degree bound.')
    return result


def natural_profiles(xletters, yletters=None, *, modes=None, keep_blocks=False):
    """Return scalar-cancelled natural w-orders for paired downset codes."""
    xletters = tuple(xletters)
    yletters = xletters if yletters is None else tuple(yletters)
    if modes is None:
        modes = max(1, max((value.bit_length() for value in xletters + yletters), default=1))
    xletters, yletters = validate_code(xletters, modes), validate_code(yletters, modes)
    matrix = one_particle(modes)
    pairs = [(x, y) for x in xletters for y in yletters]
    charges = sorted({x.bit_count() + y.bit_count() for x, y in pairs})
    blocks, orders = [], []
    began = time.monotonic()
    for charge in charges:
        chosen = [(x, y) for x, y in pairs if x.bit_count() + y.bit_count() == charge]
        rows = [x | (y << modes) for x, y in chosen]
        columns = [y | (x << modes) for x, y in chosen]
        current = compound_block(matrix, rows, columns, coefficients=keep_blocks)
        current['seconds_from_start'] = time.monotonic() - began
        blocks.append(current)
        orders.extend(2 * value for value in current['orders_t'])
    if len(orders) != len(xletters) * len(yletters):
        raise ArithmeticError('The charge blocks did not recover every channel.')
    return {'x_letters': list(xletters), 'y_letters': list(yletters), 'modes': modes,
            'energies_x': list(map(hot_energy, xletters)), 'energies_y': list(map(hot_energy, yletters)),
            'orders_w': sorted(orders), 'blocks': blocks, 'seconds': time.monotonic() - began,
            'scope': 'Complete natural profiles after cancelling opposite scalar shifts; arbitrary gauges are not supported.'}
