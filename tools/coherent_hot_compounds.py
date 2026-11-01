#!/usr/bin/env python3
"""Complete two-cut profiles of homogeneous deletion-closed coherent codes.

Binary masks encode physical vectors in the occupation basis. The two
physical spaces must be mutual retracts and invariant under every bit
deletion. No local profile is substituted for a missing opposite cut.
The variable is the deflated reflection variable t=w^2.
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


def support(vector):
    result = []
    while vector:
        bit = vector & -vector
        result.append(bit.bit_length() - 1)
        vector ^= bit
    return tuple(result)


def height(word):
    return sum(2 * bit + 1 for bit in range(word.bit_length()) if word >> bit & 1)


def _basis(vectors):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return pivots


def _in_span(value, pivots):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            return False
        value ^= pivots[pivot]
    return True


def validate_maps(embedding, retraction, modes):
    embedding, retraction = tuple(embedding), tuple(retraction)
    q = len(embedding)
    if not 1 <= modes <= 16 or not q or len(retraction) != q:
        raise ValueError('Use equal nonempty map sizes and one through sixteen modes.')
    if any(not isinstance(v, int) or v <= 0 or v.bit_length() > 1 << modes
           for v in embedding + retraction):
        raise ValueError('Physical vectors must be nonzero binary occupation masks.')
    heights, charges = [], []
    for left, right in zip(retraction, embedding):
        words = support(left | right)
        hs, ns = {height(w) for w in words}, {w.bit_count() for w in words}
        if len(hs) != 1 or len(ns) != 1:
            raise ValueError('Each paired physical label must have one height and one occupation.')
        heights.append(hs.pop())
        charges.append(ns.pop())
    if any((left & right).bit_count() % 2 != int(i == j)
           for i, left in enumerate(retraction) for j, right in enumerate(embedding)):
        raise ValueError('The physical retraction identity LP=I failed.')
    deletion_checks = 0
    for vectors in (embedding, retraction):
        pivots = _basis(vectors)
        if len(pivots) != q:
            raise ValueError('Physical map vectors are dependent.')
        for vector in vectors:
            words = support(vector)
            for bit in range(modes):
                deleted = 0
                for word in words:
                    if word >> bit & 1:
                        deleted ^= 1 << (word ^ (1 << bit))
                if not _in_span(deleted, pivots):
                    raise ValueError('A physical map space is not closed under bit deletion.')
                deletion_checks += 1
    return embedding, retraction, tuple(heights), tuple(charges), deletion_checks


def one_particle(modes):
    return [[3 << (2 * ((i % modes) | (j % modes)) + 1)
             if i // modes == j // modes else
             (1 << (2 * (i % modes) + 1) if i % modes == j % modes else 0)
             for j in range(2 * modes)] for i in range(2 * modes)]


@lru_cache(None)
def _library(words):
    source = Path(__file__).with_suffix('.cpp')
    dependencies = (source, source.with_name('selected_gaussian_coefficients.cpp'),
                    source.with_name('native_rectangular_smith_f4.cpp'))
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The coherent compound constructor requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={words}']
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = sha256(b''.join(path.read_bytes() for path in dependencies)
                    + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-coherent-compounds-' + str(os.getuid()))
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
    library.coherent_hot_compound_block.argtypes = [ctypes.c_int, ctypes.c_int, u64, ctypes.c_int,
                                                     u32, u32, u32, u32,
                                                     ctypes.POINTER(ctypes.c_int64), u64]
    library.coherent_hot_compound_block.restype = ctypes.c_int
    return library


def compound_block(matrix, row_vectors, column_vectors, *, coefficients=False, precision=None):
    n, count = len(matrix), len(row_vectors)
    if n % 2 or not 2 <= n <= 32 or any(len(row) != n for row in matrix):
        raise ValueError('Use an even square one-particle matrix of order at most 32.')
    if not 1 <= count <= 4096 or len(column_vectors) != count:
        raise ValueError('Use equal nonempty block sizes at most 4096.')
    if any(not isinstance(v, int) or v < 0 for row in matrix for v in row):
        raise ValueError('One-particle coefficients must be binary polynomial masks.')
    rows, cols = tuple(map(tuple, row_vectors)), tuple(map(tuple, column_vectors))
    if any(not row or len(set(row)) != len(row)
           or any(not isinstance(v, int) or not 0 <= v < 1 << n for v in row)
           for row in rows + cols):
        raise ValueError('Each sparse block vector must contain distinct valid occupation masks.')
    charges = {v.bit_count() for row in rows + cols for v in row}
    if len(charges) != 1:
        raise ValueError('A compound block must have one occupation number.')
    charge = charges.pop()
    degree = max(v.bit_length() - 1 for row in matrix for v in row)
    coefficient_bound = charge * max(0, degree)
    words = 2
    while 64 * words <= coefficient_bound:
        words *= 2
    if precision is not None:
        if (not isinstance(precision, int) or precision < 128 or precision > 8192
                or precision & (precision - 1) or precision <= coefficient_bound):
            raise ValueError('Precision must be a supported power of two above the complete degree bound.')
        words = precision // 64
    input_words, mask = max(1, (degree + 64) // 64), (1 << 64) - 1
    packed = (ctypes.c_uint64 * (n * n * input_words))(
        *(v >> (64 * k) & mask for row in matrix for v in row for k in range(input_words)))

    def pack_vectors(vectors):
        offsets = [0]
        for vector in vectors:
            offsets.append(offsets[-1] + len(vector))
        return ((ctypes.c_uint32 * len(offsets))(*offsets),
                (ctypes.c_uint32 * offsets[-1])(*(v for vector in vectors for v in vector)))

    ro, rr = pack_vectors(rows)
    co, cc = pack_vectors(cols)
    attempts = []
    while words <= 128:
        output = (ctypes.c_int64 * count)()
        values = (ctypes.c_uint64 * (count * count * words))() if coefficients else None
        status = _library(words).coherent_hot_compound_block(n // 2, count, packed, input_words,
                                                            ro, rr, co, cc, output, values)
        attempts.append({'precision': 64 * words, 'status': status})
        if status == 0:
            break
        if status != 6 or precision is not None:
            raise ArithmeticError(f'Coherent block failed (status {status}, precision {64 * words}).')
        words *= 2
    else:
        raise ArithmeticError('The complete coherent profile exceeded the series budget.')
    result = {'charge': charge, 'size': count, 'orders_t': sorted(map(int, output)),
              'coefficient_degree_bound_t': coefficient_bound, 'precision_t': 64 * words,
              'precision_attempts': attempts, 'endpoint_entries': count * count}
    if coefficients:
        result['matrix'] = [[sum(int(values[(i * count + j) * words + k]) << (64 * k)
                                  for k in range(words)) for j in range(count)] for i in range(count)]
        if max(v.bit_length() - 1 for row in result['matrix'] for v in row) > coefficient_bound:
            raise ArithmeticError('A complete compound coefficient exceeded its degree bound.')
    return result


def profiles(embedding, retraction, modes, *, keep_blocks=False):
    embedding, retraction, heights, charges, checks = validate_maps(embedding, retraction, modes)
    q, matrix = len(embedding), one_particle(modes)
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
            block = compound_block(matrix, row_vectors, column_vectors, coefficients=keep_blocks)
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
            'scope': 'Complete deflated natural two-cut profiles; opposite scalar shifts cancel.'}
