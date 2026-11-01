#!/usr/bin/env python3
"""Adaptive regular LDU cones for complete F4 matrix profiles.

Cached diagonals can reorder. Each hit checks all row/column regularity
conditions; no unverified affine extrapolation is used. Unresolved tracked
coefficient knowledge falls back to the existing complete native oracle.
Returned row/column IDs label affine diagonal weights for gradients; after
reordering they are not advertised as source-minor certificates.
"""
from __future__ import annotations

import ctypes
from hashlib import sha256
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile

from native_matrix_smith_f4 import MatrixOracle

_LIBRARIES = {}
_U64 = ctypes.POINTER(ctypes.c_uint64)
_I64 = ctypes.POINTER(ctypes.c_int64)
_U32 = ctypes.POINTER(ctypes.c_uint32)
_DOUBLE = ctypes.POINTER(ctypes.c_double)


def library(precision):
    if precision in _LIBRARIES:
        return _LIBRARIES[precision]
    source = Path(__file__).with_suffix('.cpp')
    included = source.with_name('native_rectangular_smith_f4.cpp')
    compiler = shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('Regular Smith chambers require a C++17 compiler.')
    flags = ['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={precision//64}']
    try:
        cpu = Path('/proc/cpuinfo').read_text()
    except OSError:
        cpu = ''
    if platform.machine().lower() in ('x86_64','amd64','i386','i686') and 'pclmulqdq' in cpu:
        flags += ['-mpclmul','-msse2']
    digest = sha256(source.read_bytes()+included.read_bytes()+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir())/('integer-factoring-smith-chambers-f4-'+str(os.getuid()))
    cache.mkdir(mode=0o700,parents=True,exist_ok=True)
    target = cache/(digest+'.so')
    if not target.exists():
        descriptor,temporary = tempfile.mkstemp(prefix=digest+'-',suffix='.so',dir=cache)
        os.close(descriptor)
        try:
            result = subprocess.run(compiler+flags+[str(source),'-o',temporary],capture_output=True,text=True)
            if result.returncode:
                raise RuntimeError('Smith chamber compilation failed:\n'+result.stderr[-12000:])
            os.replace(temporary,target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    lib = ctypes.CDLL(str(target))
    lib.smith_chamber_create.argtypes = [ctypes.c_int,ctypes.c_int,_U64,ctypes.c_int,_I64]
    lib.smith_chamber_create.restype = ctypes.c_uint64
    lib.smith_chamber_destroy.argtypes = [ctypes.c_uint64]
    lib.smith_chamber_destroy.restype = None
    lib.smith_chamber_profile.argtypes = [ctypes.c_uint64,_DOUBLE,_DOUBLE,ctypes.c_int,
        ctypes.c_double,_DOUBLE,_I64,_U32,_U32,_U64]
    lib.smith_chamber_profile.restype = ctypes.c_int
    lib.smith_chamber_compare.argtypes = [ctypes.c_int64]+[ctypes.c_double]*5
    lib.smith_chamber_compare.restype = ctypes.c_int
    _LIBRARIES[precision] = lib
    return lib


class CachedMatrixOracle:
    def __init__(self,matrix,rank,shifts=None,precision=512):
        self._initialize(MatrixOracle(matrix,rank,shifts=shifts,precision=precision))

    def _initialize(self,original):
        self.original = original
        self.n,self.rank,self.precision = original.n,original.rank,original.precision
        self._lib = library(self.precision)
        self._handle = self._lib.smith_chamber_create(self.n,self.rank,original._packed,
                                                     original._words,original._shifts)
        if not self._handle:
            raise ArithmeticError('The chamber source matrix was rejected.')
        self.cache_stats = {'calls':0,'hits':0,'builds':0,'tracked_unresolved':0,'fallbacks':0,
                            'resumed_pivots':0}

    @classmethod
    def from_oracle(cls,original):
        if isinstance(original,cls):
            return original
        result = cls.__new__(cls)
        result._initialize(original)
        return result

    def close(self):
        if getattr(self,'_handle',0):
            self._lib.smith_chamber_destroy(self._handle)
            self._handle = 0

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def profile(self,row_weights=None,column_weights=None,scale=1.,infinity=False):
        rw = [0.]*self.n if row_weights is None else list(map(float,row_weights))
        cw = [0.]*self.n if column_weights is None else list(map(float,column_weights))
        scale = float(scale)
        if (not self._handle or len(rw)!=self.n or len(cw)!=self.n or not math.isfinite(scale)
                or scale<=0 or any(not math.isfinite(v) or abs(v)/scale>1e6 for v in rw+cw)):
            raise ValueError('Supply complete finite weights and a positive scale within the native range.')
        rows,columns = (ctypes.c_double*self.n)(*rw),(ctypes.c_double*self.n)(*cw)
        orders = (ctypes.c_double*self.rank)()
        prefixes = (ctypes.c_int64*self.rank)()
        chosen_rows,chosen_columns = (ctypes.c_uint32*self.rank)(),(ctypes.c_uint32*self.rank)()
        stats = (ctypes.c_uint64*8)()
        status = self._lib.smith_chamber_profile(self._handle,rows,columns,int(bool(infinity)),scale,
                                                orders,prefixes,chosen_rows,chosen_columns,stats)
        self.cache_stats['calls'] += 1
        self.cache_stats.update(hits=int(stats[1]),builds=int(stats[2]),tracked_unresolved=int(stats[3]),
                                resumed_pivots=int(stats[6]))
        if status:
            self.cache_stats['fallbacks'] += 1
            result = self.original.profile(rw,cw,scale,infinity)
            result['chamber'] = {'hit':False,'legacy_fallback':True,'tracked_status':status}
            return result
        result = {'orders':list(orders),'base_prefix':list(prefixes),'rows':list(chosen_rows),
                  'cols':list(chosen_columns),'chamber':{'hit':bool(stats[0]),
                      'legacy_fallback':False,'regularity_inequalities':int(stats[4]),
                      'resumed_prefix_pivots':int(stats[5]),
                      'row_column_ids_are_weight_labels':True}}
        if not all(math.isfinite(v) for v in result['orders']):
            raise ArithmeticError('The chamber profile overflowed.')
        return result


def cache_sector_pressure(model):
    """Replace existing block oracles without changing the tensor constructors."""
    for sector,_ in model.cuts:
        sector.blocks = [(charge,rows,columns,CachedMatrixOracle.from_oracle(oracle))
                         for charge,rows,columns,oracle in sector.blocks]
    return model


def sector_cache_statistics(model):
    totals = {'calls':0,'hits':0,'builds':0,'tracked_unresolved':0,'fallbacks':0,'resumed_pivots':0}
    for sector,_ in model.cuts:
        for _,_,_,oracle in sector.blocks:
            if isinstance(oracle,CachedMatrixOracle):
                for key in totals:
                    totals[key] += oracle.cache_stats[key]
    return totals
