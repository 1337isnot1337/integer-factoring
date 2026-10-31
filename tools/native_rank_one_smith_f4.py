#!/usr/bin/env python3
"""Integer-order D+uv^T profile accelerator, with full replay required for results.

Input vectors contain polynomial F4 bit-plane pairs. The output cutoff
must be smaller than the work precision. This structured implementation
does not by itself certify that intermediate truncated tails are irrelevant;
retained search results are compared with complete-matrix profiles.
"""
from __future__ import annotations

import ctypes
from hashlib import sha256
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile


LIBRARIES = {}


def library(precision):
    if precision in LIBRARIES:
        return LIBRARIES[precision]
    if precision not in (128, 256, 512, 1024, 2048):
        raise ValueError('Use a supported power-of-two work precision.')
    source = Path(__file__).with_suffix('.cpp')
    included = source.with_name('native_rectangular_smith_f4.cpp')
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional structured search engine requires C++17.')
    flags = ['-O3', '-std=c++17', '-shared', '-fPIC', '-DF4_SERIES_WORDS='+str(precision//64)]
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64', 'amd64', 'i386', 'i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul', '-msse2']
    digest = sha256(source.read_bytes()+included.read_bytes()+repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir())/('integer-factoring-rank-one-f4-'+str(os.getuid()))
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = cache/(digest+'.so')
    if not target.exists():
        fd, temporary = tempfile.mkstemp(prefix=digest+'-', suffix='.so', dir=cache)
        os.close(fd)
        try:
            result = subprocess.run(compiler+flags+[str(source), '-o', temporary], capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr[-12000:])
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    result = ctypes.CDLL(str(target))
    result.rank_one_f4_profile.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_uint64), ctypes.c_int,
                                          ctypes.c_int, ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)]
    result.rank_one_f4_profile.restype = ctypes.c_int
    LIBRARIES[precision] = result
    return result


def rank_one_profile(diagonal, left, right, cutoff=128, work_precision=512):
    n = len(diagonal)
    if not 1 <= n <= 4096 or len(left) != n or len(right) != n:
        raise ValueError('Supply three nonempty vectors of the same supported length.')
    if cutoff < 1 or cutoff >= work_precision:
        raise ValueError('The retained output cutoff must be below the work precision.')
    values = list(diagonal)+list(left)+list(right)
    if any(len(pair) != 2 or any(not isinstance(a, int) or a < 0 for a in pair) for pair in values):
        raise ValueError('Use nonnegative binary coefficient bit planes.')
    bits = max((a.bit_length() for pair in values for a in pair), default=0)
    if bits > work_precision:
        raise ValueError('Input polynomial exceeds the work precision.')
    input_words = max(1, (bits+63)//64); mask = (1 << 64)-1
    packed = (ctypes.c_uint64*(6*n*input_words))(*(a >> (64*k) & mask
               for pair in values for a in pair for k in range(input_words)))
    orders = (ctypes.c_int64*n)(); counts = (ctypes.c_int64*6)()
    status = library(work_precision).rank_one_f4_profile(n, packed, input_words, cutoff, orders, counts)
    if status:
        raise ArithmeticError('Structured profile unresolved (status '+str(status)+').')
    return {'orders': list(orders), 'pivot_counts': dict(zip(
        ('diagonal_update', 'offdiagonal_to_arrow', 'arrow_center', 'arrow_leaf_diagonal',
         'arrow_leaf_column', 'arrow_leaf_row'), counts))}
