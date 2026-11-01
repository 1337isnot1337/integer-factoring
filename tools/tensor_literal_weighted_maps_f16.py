#!/usr/bin/env python3
"""Literal F16 profiles with termwise Laurent weights and explicit clearing.

Sparse terms are indexed by (occupation mask, integer monomial exponent).
The common clearing power is never discarded: returned Laurent orders
subtract it, and the receipt records both forms and the actual tail cap.
"""
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

import tensor_hot_compounds_f16 as field


@lru_cache(None)
def library(words):
    source = Path(__file__).with_suffix('.cpp')
    dependencies = [source] + [source.with_name(name) for name in
        ('tensor_literal_hot_pfaffians_f16.cpp', 'tensor_hot_compounds_f16_fast.cpp',
         'tensor_hot_compounds_f16.cpp', 'coherent_hot_compounds_f4.cpp',
         'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')]
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or not shutil.which(compiler[0]):
        raise RuntimeError('C++17 compiler required.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={words}']
    if platform.machine().lower() in ('x86_64', 'amd64') and 'pclmulqdq' in Path('/proc/cpuinfo').read_text():
        flags += ['-mpclmul', '-msse2']
    digest = sha256(b''.join(path.read_bytes() for path in dependencies)
                    + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-weighted-literal-f16-' + str(os.getuid()))
    cache.mkdir(mode=0o700, exist_ok=True)
    target = cache / (digest + '.so')
    if not target.exists():
        fd, temporary = tempfile.mkstemp(suffix='.so', dir=cache)
        os.close(fd)
        try:
            subprocess.run(compiler + flags + [str(source), '-o', temporary],
                           check=True, capture_output=True, text=True)
            os.replace(temporary, target)
        except subprocess.CalledProcessError as error:
            raise RuntimeError(error.stderr) from error
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    result = ctypes.CDLL(str(target))
    u32, u64 = ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint64)
    i32 = ctypes.POINTER(ctypes.c_int32)
    result.tensor_literal_weighted_block_f16.argtypes = [ctypes.c_int, ctypes.c_int,
        u32, u32, u32, i32, u32, u32, u32, i32, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int64), u64, ctypes.c_int, ctypes.c_uint64]
    result.tensor_literal_weighted_block_f16.restype = ctypes.c_int
    result.tensor_literal_pfaffian_cache_info.argtypes = [u64]
    result.tensor_literal_pfaffian_cache_info.restype = None
    return result


def cache_info(words):
    values = (ctypes.c_uint64 * 3)()
    library(words).tensor_literal_pfaffian_cache_info(values)
    return dict(zip(('entries', 'requests', 'hits'), map(int, values)))


def tensor_terms(a, b, modes, a_weight=None, b_weight=None):
    """Tensor two constant maps, optionally adding assigned-minus-word gauges."""
    out = {}
    for u, x in a.items():
        for v, y in b.items():
            shift = (0 if a_weight is None else a_weight-field.height(u))
            shift += 0 if b_weight is None else b_weight-field.height(v)
            key = (u | (v << modes), shift)
            out[key] = out.get(key, 0) ^ field.mul(x, y)
            if not out[key]:
                del out[key]
    return out


def polynomial_tensor_terms(a, b, modes, a_weight=None, b_weight=None, parameter_sign=1):
    """Tensor Laurent maps; reverse their parameter exponents at infinity."""
    if parameter_sign not in (-1, 1):
        raise ValueError('The parameter sign must identify one of the two places.')
    out = {}
    for (u, e), x in a.items():
        for (v, f), y in b.items():
            shift = parameter_sign*(e+f)
            shift += 0 if a_weight is None else a_weight-field.height(u)
            shift += 0 if b_weight is None else b_weight-field.height(v)
            key = (u | (v << modes), shift)
            out[key] = out.get(key, 0) ^ field.mul(x, y)
            if not out[key]:
                del out[key]
    return out


def clearing_power(rows, columns):
    return max(0, -min(shift for vector in rows for _, shift in vector)
                  -min(shift for vector in columns for _, shift in vector))


def block(modes, rows, columns, coefficients=False, order_cap=None,
          clearing=None, cache_budget=8_000_000):
    count = len(rows)
    if not 1 <= modes <= 16 or not 1 <= count <= 4608 or len(columns) != count:
        raise ValueError('Invalid weighted literal block dimensions.')
    if order_cap is not None and (type(order_cap) is not int or order_cap < 1):
        raise ValueError('Expected a positive integer Laurent order cap.')
    if any(not vector or any(type(w) is not int or not 0 <= w < 1 << (2*modes)
                             or type(shift) is not int or abs(shift) > 8192 or not 0 < c < 16
                             for (w, shift), c in vector.items()) for vector in rows + columns):
        raise ValueError('Invalid sparse weighted F16 vector.')
    if len({w.bit_count() % 2 for vector in rows + columns for w, _ in vector}) != 1:
        raise ValueError('Mixed occupation parity in one block.')
    required = clearing_power(rows, columns)
    clearing = required if clearing is None else clearing
    if type(clearing) is not int or not required <= clearing <= 8192:
        raise ValueError('The common clearing power does not clear every term.')
    low = (1 << modes)-1
    term_degree = lambda w, shift: field.height(w & low)+field.height(w >> modes)+shift
    bound = max(term_degree(w, shift) for vector in rows for w, shift in vector)
    bound += max(term_degree(w, shift) for vector in columns for w, shift in vector)+clearing
    normalized_cap = None if order_cap is None else order_cap+clearing
    needed = bound+1 if coefficients or order_cap is None else normalized_cap
    words = 1
    while 64*words < needed:
        words *= 2

    def pack(vectors):
        offsets, masks, scalars, shifts = [0], [], [], []
        for vector in vectors:
            for (mask, shift), scalar in sorted(vector.items()):
                masks.append(mask)
                scalars.append(scalar)
                shifts.append(shift)
            offsets.append(len(masks))
        return ((ctypes.c_uint32 * len(offsets))(*offsets),
                (ctypes.c_uint32 * len(masks))(*masks),
                (ctypes.c_uint32 * len(scalars))(*scalars),
                (ctypes.c_int32 * len(shifts))(*shifts))

    row_data, column_data = pack(rows), pack(columns)
    attempts = []
    while words <= 128:
        orders = (ctypes.c_int64 * count)()
        values = (ctypes.c_uint64 * (4 * count * count * words))() if coefficients else None
        status = library(words).tensor_literal_weighted_block_f16(modes, count, *row_data, *column_data,
                    clearing, orders, values, normalized_cap or 0, cache_budget)
        attempts.append([64*words, status])
        if not status:
            break
        if status != 6:
            raise ArithmeticError(f'Weighted literal status {status}; cache info {cache_info(words)}.')
        words *= 2
    else:
        raise ArithmeticError('Weighted literal profile exceeded the series budget.')
    raw = sorted(map(int, orders))
    result = {'orders_t': [x-clearing for x in raw], 'cleared_orders_t': raw,
              'clearing_power': clearing, 'size': count, 'cleared_polynomial_degree_bound': bound,
              'precision_attempts': attempts, 'endpoint_entries': count*count,
              'endpoint_source': 'exact t=1 kernel identity and actual sparse-map contraction',
              'binary_pfaffian_cache': cache_info(words)}
    if order_cap is not None:
        result.update(order_cap=order_cap, cleared_order_cap=normalized_cap,
                      orders_are_exact_below_cap=True, tail_orders_lower_bound=order_cap,
                      exact_orders_below_cap=sum(x < normalized_cap for x in raw),
                      capped_tail_count=sum(x == normalized_cap for x in raw))
    if coefficients:
        result['matrix'] = [[tuple(sum(int(values[(4*(i*count+j)+plane)*words+k]) << (64*k)
                                      for k in range(words)) for plane in range(4))
                            for j in range(count)] for i in range(count)]
    return result
