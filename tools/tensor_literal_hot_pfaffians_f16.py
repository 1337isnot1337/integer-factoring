#!/usr/bin/env python3
"""Literal hot-reflection Pfaffian profiles for arbitrary homogeneous maps.

Deletion closure is not required. Capped requests compute exact finite
jets, while generic full rank follows from the independently checked
crossing endpoint. A bounded cache shares binary Pfaffians across entries.
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
        ('tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.cpp',
         'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
         'native_rectangular_smith_f4.cpp')]
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or not shutil.which(compiler[0]):
        raise RuntimeError('C++17 compiler required.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', f'-DF4_SERIES_WORDS={words}']
    if platform.machine().lower() in ('x86_64', 'amd64') and 'pclmulqdq' in Path('/proc/cpuinfo').read_text():
        flags += ['-mpclmul', '-msse2']
    digest = sha256(b''.join(path.read_bytes() for path in dependencies)
                    + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ('integer-factoring-literal-hot-f16-' + str(os.getuid()))
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
    result.tensor_literal_hot_block_f16.argtypes = [ctypes.c_int, ctypes.c_int,
        u32, u32, u32, u32, u32, u32, ctypes.POINTER(ctypes.c_int64), u64,
        ctypes.c_int, ctypes.c_uint64]
    result.tensor_literal_hot_block_f16.restype = ctypes.c_int
    result.tensor_literal_pfaffian_control.argtypes = [ctypes.c_int, ctypes.c_uint64, u64, ctypes.c_uint64]
    result.tensor_literal_pfaffian_control.restype = ctypes.c_int
    result.tensor_literal_pfaffian_cache_info.argtypes = [u64]
    result.tensor_literal_pfaffian_cache_info.restype = None
    return result


def cache_info(words):
    values = (ctypes.c_uint64 * 3)()
    library(words).tensor_literal_pfaffian_cache_info(values)
    return dict(zip(('entries', 'requests', 'hits'), map(int, values)))


def literal_pfaffian(modes, mask, precision_words=None, cache_budget=8_000_000):
    if not 1 <= modes <= 16 or type(mask) is not int or not 0 <= mask < 1 << (4 * modes):
        raise ValueError('Invalid literal Pfaffian mask.')
    bound = sum(2 * (i % modes) + 1 for i in range(4 * modes) if mask >> i & 1)
    words = 1
    while 64 * words <= bound:
        words *= 2
    if precision_words is not None:
        words = precision_words
    if type(words) is not int or words < 1 or words > 128 or words & (words - 1):
        raise ValueError('Invalid binary coefficient precision.')
    output = (ctypes.c_uint64 * words)()
    status = library(words).tensor_literal_pfaffian_control(modes, mask, output, cache_budget)
    if status:
        raise ArithmeticError(f'Literal Pfaffian status {status}.')
    return sum(int(value) << (64 * i) for i, value in enumerate(output))


def compound_block(modes, rows, columns, coefficients=False, order_cap=None, cache_budget=8_000_000):
    count = len(rows)
    if not 1 <= modes <= 16 or not 1 <= count <= 4608 or len(columns) != count:
        raise ValueError('Invalid literal block dimensions.')
    if order_cap is not None and (type(order_cap) is not int or order_cap < 1):
        raise ValueError('Expected a positive integer order cap.')
    if any(not vector or any(type(w) is not int or not 0 <= w < 1 << (2*modes) or not 0 < c < 16
                             for w, c in vector.items()) for vector in rows + columns):
        raise ValueError('Invalid sparse F16 vector.')
    if len({w.bit_count() % 2 for vector in rows + columns for w in vector}) != 1:
        raise ValueError('Mixed occupation parity.')
    low = (1 << modes) - 1
    heights = [{field.height(w & low) + field.height(w >> modes) for w in vector}
               for vector in rows + columns]
    if any(len(hh) != 1 for hh in heights):
        raise ValueError('A physical vector mixes inherited heights.')
    bound = max(next(iter(hh)) for hh in heights[:count]) + max(next(iter(hh)) for hh in heights[count:])
    words = 1
    needed = bound + 1 if coefficients or order_cap is None else order_cap
    while 64 * words < needed:
        words *= 2

    def pack(vectors):
        offsets, masks, scalars = [0], [], []
        for vector in vectors:
            for mask, scalar in sorted(vector.items()):
                masks.append(mask)
                scalars.append(scalar)
            offsets.append(len(masks))
        return ((ctypes.c_uint32 * len(offsets))(*offsets),
                (ctypes.c_uint32 * len(masks))(*masks),
                (ctypes.c_uint32 * len(scalars))(*scalars))

    row_data, column_data = pack(rows), pack(columns)
    attempts = []
    while words <= 128:
        orders = (ctypes.c_int64 * count)()
        values = (ctypes.c_uint64 * (4 * count * count * words))() if coefficients else None
        status = library(words).tensor_literal_hot_block_f16(modes, count, *row_data, *column_data,
                    orders, values, order_cap or 0, cache_budget)
        attempts.append([64 * words, status])
        if not status:
            break
        if status != 6:
            raise ArithmeticError(f'Literal block status {status}; cache info {cache_info(words)}.')
        words *= 2
    else:
        raise ArithmeticError('Literal Smith profile exceeded the series budget.')
    result = {'orders_t': sorted(map(int, orders)), 'size': count,
              'literal_polynomial_degree_bound': bound, 'precision_attempts': attempts,
              'endpoint_entries': count * count,
              'endpoint_source': 'exact t=1 kernel identity and actual sparse-map contraction',
              'binary_pfaffian_cache': cache_info(words)}
    if order_cap is not None:
        result.update(order_cap=order_cap, orders_are_exact_below_cap=True,
                      tail_orders_lower_bound=order_cap,
                      exact_orders_below_cap=sum(int(x) < order_cap for x in orders),
                      capped_tail_count=sum(int(x) == order_cap for x in orders))
    if coefficients:
        result['matrix'] = [[tuple(sum(int(values[(4 * (i * count + j) + plane) * words + k]) << (64 * k)
                                      for k in range(words)) for plane in range(4))
                            for j in range(count)] for i in range(count)]
    return result
