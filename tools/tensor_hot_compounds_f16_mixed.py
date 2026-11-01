#!/usr/bin/env python3
"""Exact height-homogeneous F16 compounds with mixed occupation parity.

Each vector has one height. Terms of different particle number pair to
zero; terms of equal number use the existing binary K-minor constructor.
"""
import ctypes
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import os, platform, shlex, shutil, subprocess, tempfile
from coherent_hot_compounds import height
import tensor_hot_compounds_f16 as field

@lru_cache(None)
def library(words):
    source=Path(__file__).with_suffix('.cpp')
    dependencies=[source,source.with_name('tensor_hot_compounds_f16_fast.cpp'),source.with_name('tensor_hot_compounds_f16.cpp'),source.with_name('coherent_hot_compounds_f4.cpp'),source.with_name('selected_gaussian_coefficients.cpp'),source.with_name('native_rectangular_smith_f4.cpp')]
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or not shutil.which(compiler[0]):raise RuntimeError('C++17 compiler required.')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={words}']
    if platform.machine().lower() in ('x86_64','amd64') and 'pclmulqdq' in Path('/proc/cpuinfo').read_text():flags+=['-mpclmul','-msse2']
    digest=sha256(b''.join(p.read_bytes() for p in dependencies)+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-coherent-f16-mixed-'+str(os.getuid()));cache.mkdir(mode=0o700,exist_ok=True)
    target=cache/(digest+'.so')
    if not target.exists():
        fd,tmp=tempfile.mkstemp(suffix='.so',dir=cache);os.close(fd)
        try:
            subprocess.run(compiler+flags+[str(source),'-o',tmp],check=True,capture_output=True,text=True)
            os.replace(tmp,target)
        except subprocess.CalledProcessError as e:raise RuntimeError(e.stderr) from e
        finally:
            if os.path.exists(tmp):os.unlink(tmp)
    lib=ctypes.CDLL(str(target));u32=ctypes.POINTER(ctypes.c_uint32);u64=ctypes.POINTER(ctypes.c_uint64)
    lib.coherent_hot_compound_f16_mixed.argtypes=[ctypes.c_int,ctypes.c_int,u32,u32,u32,u32,u32,u32,ctypes.POINTER(ctypes.c_int64),u64]
    lib.coherent_hot_compound_f16_mixed.restype=ctypes.c_int
    lib.coherent_hot_compound_f16_mixed_capped.argtypes=lib.coherent_hot_compound_f16_mixed.argtypes+[ctypes.c_int]
    lib.coherent_hot_compound_f16_mixed_capped.restype=ctypes.c_int
    lib.tensor_f16_scalar_control.argtypes=[ctypes.c_uint,ctypes.c_uint]
    lib.tensor_f16_scalar_control.restype=ctypes.c_int
    lib.tensor_f16_fast_minor_control.argtypes=[ctypes.c_int,ctypes.c_uint32,ctypes.c_uint32,u64]
    lib.tensor_f16_fast_minor_control.restype=ctypes.c_int
    return lib


def compound_block(modes,rows,columns,coefficients=False,order_cap=None):
    count=len(rows)
    if order_cap is not None and (type(order_cap)is not int or order_cap<1):raise ValueError('Expected a positive integer Smith-order cap.')
    if not 1<=modes<=16 or not 1<=count<=4608 or len(columns)!=count:raise ValueError('Invalid block size.')
    if any(not v or any(not isinstance(w,int) or not 0<=w<1<<(2*modes) or not 1<=c<16 for w,c in v.items()) for v in rows+columns):raise ValueError('Invalid vector.')
    charges={w.bit_count()%2 for v in rows+columns for w in v}
    if len(charges)!=1:raise ValueError('Mixed occupation parity.')
    low=(1<<modes)-1
    hs=[{height(w&low)+height(w>>modes) for w in v} for v in rows+columns]
    if any(len(h)!=1 for h in hs):raise ValueError('Mixed height.')
    bound=max(next(iter(h)) for h in hs[:count])+max(next(iter(h)) for h in hs[count:])
    # K-minors have degree at most half the regular entry bound. A
    # capped jet below64 therefore needs only one binary word per plane;
    # the separate monomial base carries the remaining entry height.
    words=1 if order_cap is not None and order_cap<=64 and not coefficients and bound<128 else 2
    while 64*words<=(bound//2 if order_cap is not None and not coefficients else bound):words*=2
    def pack(vectors):
        offsets=[0];masks=[];scalars=[]
        for v in vectors:
            for w,c in sorted(v.items()):masks.append(w);scalars.append(c)
            offsets.append(len(masks))
        return ((ctypes.c_uint32*len(offsets))(*offsets),(ctypes.c_uint32*len(masks))(*masks),(ctypes.c_uint32*len(scalars))(*scalars))
    rp,cp=pack(rows),pack(columns);attempts=[]
    while words<=128:
        orders=(ctypes.c_int64*count)();values=(ctypes.c_uint64*(4*count*count*words))() if coefficients else None
        if order_cap is None:
            status=library(words).coherent_hot_compound_f16_mixed(modes,count,*rp,*cp,orders,values)
        else:
            status=library(words).coherent_hot_compound_f16_mixed_capped(modes,count,*rp,*cp,orders,values,order_cap)
        attempts.append([64*words,status])
        if status==0:break
        if status!=6:raise ArithmeticError(f'F16 compound error {status} at precision {64*words}.')
        words*=2
    else:raise ArithmeticError('Incomplete profile: precision exhausted.')
    result={'orders_t':sorted(map(int,orders)),'size':count,'degree_bound_t':bound,'precision_attempts':attempts,'endpoint_entries':count*count}
    if order_cap is not None:
        result.update(order_cap=order_cap,orders_are_exact_below_cap=True,
                      tail_orders_lower_bound=order_cap,
                      exact_orders_below_cap=sum(int(value)<order_cap for value in orders),
                      capped_tail_count=sum(int(value)==order_cap for value in orders))
    if coefficients:
        result['matrix']=[[tuple(sum(int(values[(4*(i*count+j)+plane)*words+k])<<(64*k) for k in range(words)) for plane in range(4)) for j in range(count)] for i in range(count)]
    return result

