#!/usr/bin/env python3
"""Finite Smith orders of an F4 Laurent matrix with a supplied generic rank.

The rank is an external mathematical promise. The 512-coefficient series
window cannot certify that later pivots vanish. Missing a promised pivot
raises an error. Field arithmetic is exact, weights use floating point,
and a proposed asymptotic pressure improvement needs independent exact
replay. Polynomials are pairs of bit planes as in the rectangular engine.
"""
from __future__ import annotations

import ctypes
import hashlib
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile

from native_weighted_smith_f4 import _pack

_LIBRARY=None
_U64=ctypes.POINTER(ctypes.c_uint64)
_I64=ctypes.POINTER(ctypes.c_int64)
_U32=ctypes.POINTER(ctypes.c_uint32)
_DOUBLE=ctypes.POINTER(ctypes.c_double)


def _library():
    global _LIBRARY
    if _LIBRARY is not None:return _LIBRARY
    source=Path(__file__).with_suffix('.cpp')
    included=source.with_name('native_rectangular_smith_f4.cpp')
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional matrix profile engine requires C++17.')
    flags=['-O3','-std=c++17','-shared','-fPIC']
    try:cpu=Path('/proc/cpuinfo').read_text()
    except OSError:cpu=''
    if platform.machine().lower() in ('x86_64','amd64','i386','i686') and 'pclmulqdq' in cpu:
        flags+=['-mpclmul','-msse2']
    digest=hashlib.sha256(source.read_bytes()+included.read_bytes()+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-matrix-f4-'+str(os.getuid()))
    cache.mkdir(mode=0o700,parents=True,exist_ok=True);target=cache/(digest+'.so')
    if not target.exists():
        fd,temporary=tempfile.mkstemp(prefix=digest+'-',suffix='.so',dir=cache);os.close(fd)
        try:
            run=subprocess.run(compiler+flags+[str(source),'-o',temporary],capture_output=True,text=True)
            if run.returncode:raise RuntimeError(run.stderr[-12000:])
            os.replace(temporary,target)
        finally:
            if os.path.exists(temporary):os.unlink(temporary)
    lib=ctypes.CDLL(str(target))
    lib.matrix_f4_profile.argtypes=[ctypes.c_int,ctypes.c_int,_U64,ctypes.c_int,_I64,_DOUBLE,_DOUBLE,
                                   ctypes.c_int,ctypes.c_double,_DOUBLE,_I64,_U32,_U32]
    lib.matrix_f4_profile.restype=ctypes.c_int
    _LIBRARY=lib;return lib


class MatrixOracle:
    def __init__(self,matrix,rank,shifts=None):
        self.n=len(matrix)
        if not 1<=self.n<=1024 or any(len(row)!=self.n for row in matrix):
            raise ValueError('Supply a square matrix of order 1 through 1024.')
        if not isinstance(rank,int) or not 1<=rank<=self.n:
            raise ValueError('Supply the proved generic rank between 1 and the matrix order.')
        self.rank=rank
        self._packed,self._words=_pack([v for row in matrix for v in row],self.n*self.n)
        if shifts is None:shifts=[[0]*self.n for _ in range(self.n)]
        if len(shifts)!=self.n or any(len(row)!=self.n for row in shifts):
            raise ValueError('One integer Laurent shift is required per entry.')
        flat=[v for row in shifts for v in row]
        if any(not isinstance(v,int) or abs(v)>1000000 for v in flat):
            raise ValueError('Laurent shifts must be integers of absolute value at most one million.')
        self._shifts=(ctypes.c_int64*len(flat))(*flat)

    def profile(self,row_weights=None,column_weights=None,scale=1.,infinity=False):
        rw=[0.]*self.n if row_weights is None else list(map(float,row_weights))
        cw=[0.]*self.n if column_weights is None else list(map(float,column_weights))
        scale=float(scale)
        if len(rw)!=self.n or len(cw)!=self.n or not math.isfinite(scale) or scale<=0:
            raise ValueError('Supply finite weights of the matrix order and positive finite scale.')
        if any(not math.isfinite(v) or abs(v)/scale>1000000 for v in rw+cw):
            raise ValueError('Weights must be finite and within the normalized conversion bound.')
        rr=(ctypes.c_double*self.n)(*rw);cc=(ctypes.c_double*self.n)(*cw)
        orders=(ctypes.c_double*self.rank)();values=(ctypes.c_int64*self.rank)()
        rows=(ctypes.c_uint32*self.rank)();cols=(ctypes.c_uint32*self.rank)()
        status=_library().matrix_f4_profile(self.n,self.rank,self._packed,self._words,self._shifts,
                                            rr,cc,int(bool(infinity)),scale,orders,values,rows,cols)
        if status:raise ArithmeticError(f'Matrix profile did not resolve every promised pivot (status {status}).')
        result={'orders':list(orders),'base_prefix':list(values),'rows':list(rows),'cols':list(cols)}
        if not all(math.isfinite(v) for v in result['orders']):raise ArithmeticError('Profile weights overflowed.')
        return result
