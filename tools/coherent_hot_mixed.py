#!/usr/bin/env python3
"""Complete natural cut profiles of height-homogeneous coherent retractions.

Physical vectors may mix occupation numbers. Both map spaces must remain
invariant under every bit deletion. Height fixes parity, so the two complete
cut matrices split into even and odd blocks, with no occupation omission.
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

from coherent_hot_compounds import _basis, _in_span, height, support


def validate_maps(embedding, retraction, modes):
    embedding, retraction = tuple(embedding), tuple(retraction)
    q = len(embedding)
    if not 1 <= modes <= 16 or not q or len(retraction) != q:
        raise ValueError('Use equal nonempty map sizes and one through sixteen modes.')
    if any(not isinstance(v, int) or v <= 0 or v.bit_length() > 1 << modes
           for v in embedding + retraction):
        raise ValueError('Physical vectors must be nonzero binary occupation masks.')
    heights = []
    for left, right in zip(retraction, embedding):
        hs = {height(word) for word in support(left | right)}
        if len(hs) != 1:
            raise ValueError('Each paired physical label must have one height.')
        heights.append(hs.pop())
    if any((left & right).bit_count() % 2 != int(i == j)
           for i, left in enumerate(retraction) for j, right in enumerate(embedding)):
        raise ValueError('The physical retraction identity LP=I failed.')
    checks = 0
    for vectors in (embedding, retraction):
        pivots = _basis(vectors)
        if len(pivots) != q:
            raise ValueError('Physical map vectors are dependent.')
        for vector in vectors:
            for bit in range(modes):
                deleted = 0
                for word in support(vector):
                    if word >> bit & 1:
                        deleted ^= 1 << (word ^ (1 << bit))
                if not _in_span(deleted, pivots):
                    raise ValueError('A physical map space is not closed under bit deletion.')
                checks += 1
    return embedding, retraction, tuple(heights), checks


@lru_cache(None)
def _library(words):
    source = Path(__file__).with_suffix('.cpp')
    dependencies = (source, source.with_name('coherent_hot_compounds_fast.cpp'),
                    source.with_name('selected_gaussian_coefficients.cpp'),
                    source.with_name('native_rectangular_smith_f4.cpp'))
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The mixed hot constructor requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={words}']
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = sha256(b''.join(path.read_bytes() for path in dependencies)
                    + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-coherent-mixed-' + str(os.getuid()))
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
    library.coherent_hot_mixed_block.argtypes = [ctypes.c_int, ctypes.c_int,
                                                u32, u32, u32, u32,
                                                ctypes.POINTER(ctypes.c_int64), u64]
    library.coherent_hot_mixed_block.restype = ctypes.c_int
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
    parities = {value.bit_count() % 2 for row in rows + columns for value in row}
    if len(parities) != 1:
        raise ValueError('A mixed block must have one occupation parity.')
    parity = parities.pop()
    low = (1 << modes) - 1
    heights = [{height(value & low) + height(value >> modes) for value in row}
               for row in rows + columns]
    if any(len(hs) != 1 for hs in heights):
        raise ValueError('Each sparse vector must have one height.')
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
        status = _library(words).coherent_hot_mixed_block(modes, count, ro, rr, co, cc,
                                                         output, values)
        attempts.append({'precision': 64 * words, 'status': status})
        if status == 0:
            break
        if status != 6 or precision is not None:
            raise ArithmeticError(f'Mixed block failed (status {status}, precision {64 * words}).')
        words *= 2
    else:
        raise ArithmeticError('The complete mixed profile exceeded the series budget.')
    result = {'parity': parity, 'size': count, 'orders_t': sorted(map(int, output)),
              'coefficient_degree_bound_t': coefficient_bound, 'precision_t': 64 * words,
              'precision_attempts': attempts, 'endpoint_entries': count * count}
    if coefficients:
        result['matrix'] = [[sum(int(values[(i * count + j) * words + k]) << (64 * k)
                                  for k in range(words)) for j in range(count)] for i in range(count)]
        assert max(v.bit_length() - 1 for row in result['matrix'] for v in row) <= coefficient_bound
    return result


def profiles(embedding, retraction, modes, *, keep_blocks=False):
    embedding, retraction, heights, checks = validate_maps(embedding, retraction, modes)
    q = len(embedding)
    ps, ls = tuple(map(support, embedding)), tuple(map(support, retraction))
    cuts, started = [], time.monotonic()
    for cut in ('zero', 'infinity'):
        blocks, orders = [], []
        for parity in (0, 1):
            pairs = [(x, y) for x in range(q) for y in range(q)
                     if (heights[x] + heights[y]) % 2 == parity]
            if not pairs:
                continue
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
            raise ArithmeticError('Parity blocks omitted a physical channel.')
        cuts.append({'cut': cut, 'orders_t': sorted(orders), 'blocks': blocks})
    return {'modes': modes, 'q': q, 'embedding': list(embedding), 'retraction': list(retraction),
            'heights': list(heights), 'deletion_checks': checks,
            'zero_orders_t': cuts[0]['orders_t'], 'infinity_orders_t': cuts[1]['orders_t'],
            'cuts': cuts, 'seconds': time.monotonic() - started,
            'scope': 'Complete natural two-cut profiles; height homogeneous, possibly mixed occupation.'}
