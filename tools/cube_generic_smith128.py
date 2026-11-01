#!/usr/bin/env python3
"""Exact generic cut orders using homogeneous finite-field specialization.

Input coefficients are in F_2. Row and column degrees are twice their
nonnegative homogeneous weights. Every minor has total degree below 128.
Under these checked conditions, setting one residual variable to 1 and
the other to a degree-128 field generator preserves all axis valuations.
"""
from __future__ import annotations

import ctypes
from functools import lru_cache
from hashlib import sha256
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile


FIELD_DEGREE = 128
FIELD_MODULUS = (1 << FIELD_DEGREE) | 0x87
WORD_MASK = (1 << 64) - 1


def field_reduce(value):
    while value.bit_length() > FIELD_DEGREE:
        value ^= FIELD_MODULUS << (value.bit_length() - FIELD_DEGREE - 1)
    return value


def field_product(first, second):
    result = 0
    while second:
        bit = second & -second
        result ^= first << (bit.bit_length() - 1)
        second ^= bit
    return field_reduce(result)


@lru_cache(None)
def native_library():
    source = Path(__file__).with_suffix('.cpp')
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('A C++17 compiler is required for the generic Smith engine.')
    flags = ['-std=c++17', '-O3', '-shared', '-fPIC']
    cpuinfo = Path('/proc/cpuinfo')
    if (platform.machine().lower() in ('x86_64', 'amd64') and cpuinfo.exists()
            and 'pclmulqdq' in cpuinfo.read_text()):
        flags += ['-mpclmul', '-msse4.1']
    digest = sha256(source.read_bytes() + repr(compiler + flags).encode()).hexdigest()[:24]
    target = Path(tempfile.gettempdir()) / ('integer-factoring-cube-generic-' + digest + '.so')
    if not target.exists():
        temporary = target.with_name(target.name + '.' + str(os.getpid()) + '.tmp')
        completed = subprocess.run(compiler + flags + [str(source), '-o', str(temporary)],
                                   text=True, capture_output=True)
        if completed.returncode:
            raise RuntimeError('The generic Smith engine did not compile: ' + completed.stderr)
        temporary.replace(target)
    library = ctypes.CDLL(str(target))
    pointer = ctypes.POINTER(ctypes.c_uint64)
    library.cube_gf128_smith.argtypes = [ctypes.c_int, ctypes.c_int, pointer,
                                       ctypes.POINTER(ctypes.c_int64)]
    library.cube_gf128_smith.restype = ctypes.c_int
    library.cube_gf128_product.argtypes = [pointer, pointer, pointer]
    library.cube_gf128_product.restype = None
    return library


def generic_smith_profile(matrix, axis, row_degrees, column_degrees):
    """Return all generic Smith orders, or reject an unsupported input.

    Degrees are doubled weights: a term (a,b,c) in entry (i,j) must
    satisfy 2(a+b+c)=row_degrees[i]+column_degrees[j]. All coefficients
    are binary, represented by sets/frozensets of monomial exponent triples.
    """
    order = len(matrix)
    if not 1 <= order <= 64 or any(len(row) != order for row in matrix):
        raise ValueError('Use a nonempty square matrix of order at most 64.')
    if axis not in (0, 1, 2):
        raise ValueError('The axis must be 0, 1, or 2.')
    if len(row_degrees) != order or len(column_degrees) != order:
        raise ValueError('Supply one doubled weight per row and column.')
    if any(not isinstance(value, int) or value < 0
           for value in tuple(row_degrees) + tuple(column_degrees)):
        raise ValueError('Doubled weights must be nonnegative integers.')
    bound = (sum(row_degrees) + sum(column_degrees)) // 2
    if bound >= FIELD_DEGREE:
        raise ValueError('The degree-128 extension must exceed every minor total degree.')
    precision = bound + 1
    data = (ctypes.c_uint64 * (order * order * precision * 2))()
    remaining = [index for index in range(3) if index != axis]
    for i, row in enumerate(matrix):
        for j, polynomial in enumerate(row):
            if not isinstance(polynomial, (set, frozenset)):
                raise TypeError('Use binary sets of monomial exponent triples.')
            for term in polynomial:
                if (not isinstance(term, tuple) or len(term) != 3
                        or any(not isinstance(power, int) or power < 0 for power in term)):
                    raise ValueError('Every monomial needs three nonnegative integer exponents.')
                if 2 * sum(term) != row_degrees[i] + column_degrees[j]:
                    raise ValueError('The additive homogeneous weights do not match an entry.')
                exponent = term[axis]
                value = field_reduce(1 << term[remaining[1]])
                offset = ((i * order + j) * precision + exponent) * 2
                data[offset] ^= value & WORD_MASK
                data[offset + 1] ^= value >> 64
    output = (ctypes.c_int64 * order)()
    status = native_library().cube_gf128_smith(order, precision, data, output)
    if status == 1:
        raise ArithmeticError('The generic matrix is singular or has an unsupported valuation.')
    if status:
        raise ArithmeticError('The generic Smith engine returned status ' + str(status))
    orders = list(output)
    if orders != sorted(orders) or min(orders) < 0 or sum(orders) > bound:
        raise ArithmeticError('The returned profile violates its complete determinant budget.')
    return {'orders': orders, 'axis': axis, 'field_degree': FIELD_DEGREE,
            'complete_minor_degree_bound': bound, 'absolute_series_precision': precision,
            'scope': 'Exact generic orders over F_2 in the remaining two variables.'}
