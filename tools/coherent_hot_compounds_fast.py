#!/usr/bin/env python3
"""Exact hot-kernel profiles with a tighter complete-coefficient bound.

The normalized one-particle kernel removes its row/column monomials. A
selected minor has degree at most its removed monomial exponent. This
reduces the starting series window; unresolved full profiles still retry.
"""
from __future__ import annotations

import ctypes
from functools import lru_cache
from hashlib import sha256
from itertools import product
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import time

from coherent_hot_compounds import height, support, validate_maps


@lru_cache(None)
def _library(words):
    source = Path(__file__).with_suffix('.cpp')
    dependencies = (source, source.with_name('selected_gaussian_coefficients.cpp'),
                    source.with_name('native_rectangular_smith_f4.cpp'))
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The hot compound constructor requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={words}']
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = sha256(b''.join(path.read_bytes() for path in dependencies)
                    + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-coherent-fast-' + str(os.getuid()))
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
    library.coherent_hot_compound_fast.argtypes = [ctypes.c_int, ctypes.c_int,
                                                   u32, u32, u32, u32,
                                                   ctypes.POINTER(ctypes.c_int64), u64]
    library.coherent_hot_compound_fast.restype = ctypes.c_int
    return library


def compound_block(modes, rows, columns, *, coefficients=False, precision=None):
    count = len(rows)
    if not 1 <= modes <= 16 or not 1 <= count <= 4096 or len(columns) != count:
        raise ValueError('Use at most sixteen modes and 4096 paired block vectors.')
    rows, columns = tuple(map(tuple, rows)), tuple(map(tuple, columns))
    if any(not row or len(set(row)) != len(row)
           or any(not isinstance(v, int) or not 0 <= v < 1 << (2 * modes) for v in row)
           for row in rows + columns):
        raise ValueError('Each sparse block vector must contain distinct occupation masks.')
    charge = {v.bit_count() for row in rows + columns for v in row}
    if len(charge) != 1:
        raise ValueError('Every vector must have the same occupation number.')
    charge = charge.pop()
    low = (1 << modes) - 1

    def total_height(word):
        return height(word & low) + height(word >> modes)

    heights = [{total_height(v) for v in row} for row in rows + columns]
    if any(len(hs) != 1 for hs in heights):
        raise ValueError('Every block vector must be height homogeneous.')
    coefficient_bound = max(next(iter(hs)) for hs in heights[:count]) + max(
        next(iter(hs)) for hs in heights[count:])
    words = 2
    while 64 * words <= coefficient_bound:
        words *= 2
    if precision is not None:
        if (not isinstance(precision, int) or precision < 128 or precision > 8192
                or precision & (precision - 1) or precision <= coefficient_bound):
            raise ValueError('Precision must be a supported power of two above the degree bound.')
        words = precision // 64

    def pack(vectors):
        offsets = [0]
        for vector in vectors:
            offsets.append(offsets[-1] + len(vector))
        return ((ctypes.c_uint32 * len(offsets))(*offsets),
                (ctypes.c_uint32 * offsets[-1])(*(v for row in vectors for v in row)))

    ro, rr = pack(rows)
    co, cc = pack(columns)
    attempts = []
    while words <= 128:
        output = (ctypes.c_int64 * count)()
        values = (ctypes.c_uint64 * (count * count * words))() if coefficients else None
        status = _library(words).coherent_hot_compound_fast(modes, count, ro, rr, co, cc,
                                                          output, values)
        attempts.append({'precision': 64 * words, 'status': status})
        if status == 0:
            break
        if status != 6 or precision is not None:
            raise ArithmeticError(f'Fast hot block failed (status {status}, precision {64 * words}).')
        words *= 2
    else:
        raise ArithmeticError('The complete hot profile exceeded the series budget.')
    result = {'charge': charge, 'size': count, 'orders_t': sorted(map(int, output)),
              'coefficient_degree_bound_t': coefficient_bound, 'precision_t': 64 * words,
              'precision_attempts': attempts, 'endpoint_entries': count * count}
    if coefficients:
        result['matrix'] = [[sum(int(values[(i * count + j) * words + k]) << (64 * k)
                                  for k in range(words)) for j in range(count)] for i in range(count)]
        assert max(v.bit_length() - 1 for row in result['matrix'] for v in row) <= coefficient_bound
    return result


def profiles(embedding, retraction, modes, *, keep_blocks=False):
    embedding, retraction, heights, charges, checks = validate_maps(embedding, retraction, modes)
    q = len(embedding)
    ps, ls = tuple(map(support, embedding)), tuple(map(support, retraction))
    cuts, started = [], time.monotonic()
    for cut in ('zero', 'infinity'):
        blocks, orders = [], []
        for charge in sorted({a + b for a in charges for b in charges}):
            pairs = [(x, y) for x in range(q) for y in range(q) if charges[x] + charges[y] == charge]
            row_vectors, column_vectors = [], []
            for x, y in pairs:
                row_vectors.append(tuple(a | (b << modes)
                                         for a, b in product(ls[x], ls[y] if cut == 'zero' else ps[y])))
                column_vectors.append(tuple(a | (b << modes)
                                            for a, b in product(ps[y] if cut == 'zero' else ls[y], ps[x])))
            block = compound_block(modes, row_vectors, column_vectors, coefficients=keep_blocks)
            block['pairs'] = pairs
            blocks.append(block)
            orders.extend(block['orders_t'])
        if len(orders) != q * q:
            raise ArithmeticError('Occupation blocks omitted a physical channel.')
        cuts.append({'cut': cut, 'orders_t': sorted(orders), 'blocks': blocks})
    return {'modes': modes, 'q': q, 'embedding': list(embedding), 'retraction': list(retraction),
            'heights': list(heights), 'charges': list(charges), 'deletion_checks': checks,
            'zero_orders_t': cuts[0]['orders_t'], 'infinity_orders_t': cuts[1]['orders_t'],
            'cuts': cuts, 'seconds': time.monotonic() - started,
            'scope': 'Complete natural two-cut profiles with exact hot-kernel degree bound.'}
