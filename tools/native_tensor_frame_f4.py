#!/usr/bin/env python3
"""Exact truncated polynomial matrix sandwiches for tensor-product frames."""
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


@lru_cache(None)
def library(precision):
    if precision not in (128, 256, 512, 1024):
        raise ValueError('Use a supported power-of-two precision.')
    source = Path(__file__).with_suffix('.cpp')
    included = source.with_name('native_rectangular_smith_f4.cpp')
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The exact frame helper requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', '-DF4_SERIES_WORDS='+str(precision//64)]
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = sha256(source.read_bytes()+included.read_bytes()+repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir())/('integer-factoring-tensor-frame-f4-'+str(os.getuid()))
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = cache/(digest+'.so')
    if not target.exists():
        fd, temporary = tempfile.mkstemp(prefix=digest+'-', suffix='.so', dir=cache)
        os.close(fd)
        try:
            compiled = subprocess.run(compiler+flags+[str(source), '-o', temporary], capture_output=True, text=True)
            if compiled.returncode:
                raise RuntimeError(compiled.stderr[-12000:])
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    result = ctypes.CDLL(str(target))
    result.tensor_frame_f4.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_uint64), ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_uint64)]
    result.tensor_frame_f4.restype = ctypes.c_int
    return result


def sandwich(left, middle, right, precision=128):
    n = len(left)
    if not 1 <= n <= 32 or any(len(matrix) != n or any(len(row) != n for row in matrix)
                             for matrix in (left, middle, right)):
        raise ValueError('Use three equally sized square matrices, of order at most32.')
    values = [entry for matrix in (left, middle, right) for row in matrix for entry in row]
    if any(len(pair) != 2 or any(not isinstance(a, int) or a < 0 for a in pair) for pair in values):
        raise ValueError('Use polynomial binary coefficient planes.')
    bits = max((a.bit_length() for pair in values for a in pair), default=0)
    if bits > precision:
        raise ValueError('Input exceeds the declared absolute precision.')
    input_words = max(1, (bits+63)//64); output_words = precision//64; mask = (1 << 64)-1
    packed = (ctypes.c_uint64*(6*n*n*input_words))(*(a >> (64*k) & mask
               for pair in values for a in pair for k in range(input_words)))
    result = (ctypes.c_uint64*(2*n*n*output_words))()
    status = library(precision).tensor_frame_f4(n, packed, input_words, result)
    if status:
        raise ArithmeticError('Polynomial frame transform failed: '+str(status))
    entries = [(sum(int(result[(2*i)*output_words+k]) << (64*k) for k in range(output_words)),
                sum(int(result[(2*i+1)*output_words+k]) << (64*k) for k in range(output_words)))
               for i in range(n*n)]
    return [entries[i*n:(i+1)*n] for i in range(n)]
