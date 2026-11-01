#!/usr/bin/env python3
"""Exact F5 matrix multiplication and rank via the optional FLINT library."""
from __future__ import annotations

import ctypes
from functools import lru_cache
from hashlib import sha256
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np


@lru_cache(None)
def library():
    source = Path(__file__).with_suffix('.c')
    compiler = shlex.split(os.environ.get('CC','cc'))
    version = subprocess.check_output(['pkg-config','--modversion','flint'],text=True).strip()
    cflags = shlex.split(subprocess.check_output(['pkg-config','--cflags','flint'],text=True))
    libs = shlex.split(subprocess.check_output(['pkg-config','--libs','flint'],text=True))
    flags = ['-O3','-shared','-fPIC']+cflags
    digest = sha256(source.read_bytes()+repr((compiler,flags,libs,version)).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir())/('integer-factoring-gf5-flint-'+str(os.getuid()))
    cache.mkdir(mode=0o700,parents=True,exist_ok=True)
    target = cache/(digest+'.so')
    if not target.exists():
        fd,temporary = tempfile.mkstemp(prefix=digest+'-',suffix='.so',dir=cache)
        os.close(fd)
        try:
            subprocess.run(compiler+flags+[str(source),'-o',temporary]+libs,
                           check=True,capture_output=True,text=True)
            os.replace(temporary,target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    result = ctypes.CDLL(str(target))
    pointer = ctypes.POINTER(ctypes.c_uint8)
    result.gf5_rank.argtypes = [pointer,ctypes.c_int64,ctypes.c_int64]
    result.gf5_rank.restype = ctypes.c_int64
    result.gf5_rref.argtypes = [pointer,pointer,ctypes.c_int64,ctypes.c_int64]
    result.gf5_rref.restype = ctypes.c_int64
    result.gf5_multiply.argtypes = [pointer,pointer,pointer,ctypes.c_int64,ctypes.c_int64,ctypes.c_int64]
    result.gf5_multiply.restype = None
    return result,version


def matrix(values):
    a = np.asarray(values)
    if a.ndim != 2 or not np.issubdtype(a.dtype,np.integer):
        raise ValueError('Use a two-dimensional integer coefficient array.')
    return np.ascontiguousarray(a % 5,dtype=np.uint8)


def multiply(left,right):
    left,right = matrix(left),matrix(right)
    if left.shape[1] != right.shape[0]:
        raise ValueError('Matrix inner dimensions differ.')
    out = np.zeros((left.shape[0],right.shape[1]),dtype=np.uint8)
    if not out.size or not left.shape[1]:
        return out
    pointer = ctypes.POINTER(ctypes.c_uint8)
    library()[0].gf5_multiply(out.ctypes.data_as(pointer),left.ctypes.data_as(pointer),
                             right.ctypes.data_as(pointer),left.shape[0],left.shape[1],right.shape[1])
    return out


def rank(values):
    a = matrix(values)
    if not a.size:
        return 0
    pointer = ctypes.POINTER(ctypes.c_uint8)
    return int(library()[0].gf5_rank(a.ctypes.data_as(pointer),a.shape[0],a.shape[1]))


def rref(values):
    a = matrix(values)
    out = np.zeros_like(a)
    if not a.size:
        return out,0
    pointer = ctypes.POINTER(ctypes.c_uint8)
    r = library()[0].gf5_rref(out.ctypes.data_as(pointer),a.ctypes.data_as(pointer),a.shape[0],a.shape[1])
    return out,int(r)


def pivot_columns(values):
    reduced,r = rref(values)
    return [int(np.flatnonzero(row)[0]) for row in reduced[:r]]


def inverse(values):
    a = matrix(values)
    if a.shape[0] != a.shape[1]:
        raise ValueError('Use a square matrix.')
    identity = np.eye(len(a),dtype=np.uint8)
    reduced,_ = rref(np.concatenate((a,identity),axis=1))
    if not np.array_equal(reduced[:,:len(a)],identity):
        raise ArithmeticError('Matrix is singular over F5.')
    return reduced[:,len(a):]
