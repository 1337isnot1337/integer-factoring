#!/usr/bin/env python3
"""Construct selected strict Gaussian coefficients without the full alphabet.

Each coefficient is a polynomial determinant. The native computation
uses formal series only after proving a strict determinant degree bound,
so the returned coefficients are complete exact polynomials. A missing
pivot then certifies a zero determinant, rather than hiding an unknown
tail. Every constructed tensor also checks its full strict endpoint.
"""
from __future__ import annotations
import ctypes
from functools import lru_cache
import hashlib
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
from check_schur_downset_crossing import sparse_kernel
from relabelled_gaussian_sums import binary_inverse


@lru_cache(None)
def _library(words):
    source=Path(__file__).with_suffix('.cpp');dependency=source.with_name('native_rectangular_smith_f4.cpp')
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:raise RuntimeError('A C++17 compiler is required for selected coefficient construction.')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={words}']
    try:cpu=Path('/proc/cpuinfo').read_text()
    except OSError:cpu=''
    if platform.machine().lower() in ('x86_64','amd64','i386','i686') and 'pclmulqdq' in cpu:flags+=['-mpclmul','-msse2']
    digest=hashlib.sha256(source.read_bytes()+dependency.read_bytes()+repr((compiler,flags)).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-selected-gaussian-'+str(os.getuid()));cache.mkdir(mode=0o700,parents=True,exist_ok=True);target=cache/(digest+'.so')
    if not target.exists():
        fd,temporary=tempfile.mkstemp(prefix=digest+'-',suffix='.so',dir=cache);os.close(fd)
        try:
            run=subprocess.run(compiler+flags+[str(source),'-o',temporary],capture_output=True,text=True)
            if run.returncode:raise RuntimeError(run.stderr[-8000:])
            os.replace(temporary,target)
        finally:
            if os.path.exists(temporary):os.unlink(temporary)
    lib=ctypes.CDLL(str(target));U64=ctypes.POINTER(ctypes.c_uint64);U32=ctypes.POINTER(ctypes.c_uint32)
    lib.selected_binary_exterior.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,U64,ctypes.c_int,U32,U32,U64,ctypes.c_int];lib.selected_binary_exterior.restype=ctypes.c_int
    return lib


def strict_one_particle(indices,ambient):
    indices=tuple(indices);n=len(indices);kernel=sparse_kernel(indices,ambient);one=[row[:2*n] for row in kernel[2*n:]]
    inverse=binary_inverse([[int(b&~a==0) for b in indices] for a in indices]);repaired=[]
    for row in inverse:
        value=[0]*(2*n)
        for i,coefficient in enumerate(row):
            if coefficient:value=[x^y for x,y in zip(value,one[i])]
        repaired.append(value)
    one[:n]=repaired
    return one


def selected_tensor(indices,ambient,xletters,yletters=None):
    indices=tuple(indices);n=len(indices);xletters=tuple(xletters);yletters=tuple(xletters if yletters is None else yletters)
    if not 1<=n<=16 or any(a<0 or a>=1<<n for a in xletters+yletters):raise ValueError('Invalid physical mode masks.')
    if not 1<=len(xletters)<=64 or not 1<=len(yletters)<=64 or len(set(xletters))!=len(xletters) or len(set(yletters))!=len(yletters):raise ValueError('Use distinct nonempty physical alphabets of size at most 64.')
    one=strict_one_particle(indices,ambient);degree=max(value.bit_length()-1 for row in one for value in row)
    rank_bound=max(x.bit_count() for x in xletters)+max(y.bit_count() for y in yletters)
    coefficient_degree_bound=rank_bound*degree
    words=1
    while 64*words<=coefficient_degree_bound:words*=2
    if words>64:raise ValueError('The requested determinant exceeds the selected coefficient budget.')
    input_words=max(1,(degree+64)//64);mask=(1<<64)-1
    packed=(ctypes.c_uint64*(4*n*n*input_words))(*(value>>(64*k)&mask for row in one for value in row for k in range(input_words)))
    xx=(ctypes.c_uint32*len(xletters))(*xletters);yy=(ctypes.c_uint32*len(yletters))(*yletters);count=len(xletters)**2*len(yletters)**2;out=(ctypes.c_uint64*(count*words))()
    status=_library(words).selected_binary_exterior(n,len(xletters),len(yletters),packed,input_words,xx,yy,out,1)
    if status:raise ArithmeticError(f'Selected exact exterior coefficient construction failed (status {status}).')
    values=[sum(int(out[i*words+k])<<(64*k) for k in range(words)) for i in range(count)]
    assert max(value.bit_length()-1 for value in values)<=coefficient_degree_bound
    return {'values':[(value,0) for value in values],'q_x':len(xletters),'q_u':len(yletters),'dimensions':[len(xletters),len(yletters),len(yletters),len(xletters)],'denominator_degree':0,'strict_endpoint':True,'indices':indices,'ambient':ambient,'x_letters':xletters,'y_letters':yletters,'coefficient_degree_bound':coefficient_degree_bound,'exact_precision':64*words,'degree':max(value.bit_length()-1 for value in values)}
