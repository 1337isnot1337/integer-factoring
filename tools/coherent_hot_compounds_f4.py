#!/usr/bin/env python3
"""Exact F4 coherent retractions and complete two-cut occupation profiles.

Vectors are sparse dictionaries {occupation_word: F4_coefficient}.
F4 uses omega^2=omega+1, encoded by 0,1,2,3. No conjugation is used.
"""
from __future__ import annotations
import ctypes
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import os, platform, shlex, shutil, subprocess, tempfile, time
from coherent_hot_compounds import height


def mul(a,b):
    z=0
    for _ in range(2):
        if b&1:z^=a
        b>>=1;a<<=1
        if a&4:a^=7
    return z


def axpy(x,a,y):
    z=dict(x)
    for k,v in y.items():
        z[k]=z.get(k,0)^mul(a,v)
        if not z[k]:del z[k]
    return z


def basis(vectors):
    pivots={}
    for v in vectors:
        v=dict(v)
        while v:
            k=max(v)
            if k in pivots:v=axpy(v,v[k],pivots[k])
            else:
                inv=mul(v[k],v[k]);pivots[k]={i:mul(inv,a) for i,a in v.items()};break
    return pivots


def in_span(v,pivots):
    v=dict(v)
    while v:
        k=max(v)
        if k not in pivots:return False
        v=axpy(v,v[k],pivots[k])
    return True


def dot(x,y):
    result=0
    for k,v in x.items():result^=mul(v,y.get(k,0))
    return result


def validate_maps(embedding,retraction,modes):
    if not 1<=modes<=16 or not embedding or len(embedding)!=len(retraction):raise ValueError('Invalid map dimensions.')
    q=len(embedding);hs=[];charges=[]
    for p,l in zip(embedding,retraction):
        if not p or not l or any(not isinstance(w,int) or not 0<=w<1<<modes or c not in (1,2,3) for v in (p,l) for w,c in v.items()):raise ValueError('Invalid sparse F4 vector.')
        words=set(p)|set(l);hh={height(w) for w in words};ch={w.bit_count() for w in words}
        if len(hh)!=1 or len(ch)!=1:raise ValueError('Maps must preserve height and occupation.')
        hs.append(hh.pop());charges.append(ch.pop())
    if any(dot(l,p)!=int(i==j) for i,l in enumerate(retraction) for j,p in enumerate(embedding)):raise ValueError('LP=I failed.')
    for vectors in (embedding,retraction):
        pivots=basis(vectors)
        if len(pivots)!=q:raise ValueError('Dependent map.')
        for v in vectors:
            for bit in range(modes):
                d={w^(1<<bit):c for w,c in v.items() if w>>bit&1}
                if not in_span(d,pivots):raise ValueError('Deletion closure failed.')
    return hs,charges


@lru_cache(None)
def library(words):
    source=Path(__file__).with_suffix('.cpp')
    dependencies=[source,source.with_name('selected_gaussian_coefficients.cpp'),source.with_name('native_rectangular_smith_f4.cpp')]
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or not shutil.which(compiler[0]):raise RuntimeError('C++17 compiler required.')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={words}']
    if platform.machine().lower() in ('x86_64','amd64') and 'pclmulqdq' in Path('/proc/cpuinfo').read_text():flags+=['-mpclmul','-msse2']
    digest=sha256(b''.join(p.read_bytes() for p in dependencies)+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-coherent-f4-'+str(os.getuid()));cache.mkdir(mode=0o700,exist_ok=True)
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
    lib.coherent_hot_compound_f4.argtypes=[ctypes.c_int,ctypes.c_int,u32,u32,u32,u32,u32,u32,ctypes.POINTER(ctypes.c_int64),u64]
    lib.coherent_hot_compound_f4.restype=ctypes.c_int
    return lib


def compound_block(modes,rows,columns,coefficients=False):
    count=len(rows)
    if not 1<=modes<=16 or not 1<=count<=4096 or len(columns)!=count:raise ValueError('Invalid block size.')
    if any(not v or any(not isinstance(w,int) or not 0<=w<1<<(2*modes) or c not in (1,2,3) for w,c in v.items()) for v in rows+columns):raise ValueError('Invalid vector.')
    charges={w.bit_count() for v in rows+columns for w in v}
    if len(charges)!=1:raise ValueError('Mixed occupation.')
    low=(1<<modes)-1
    hs=[{height(w&low)+height(w>>modes) for w in v} for v in rows+columns]
    if any(len(h)!=1 for h in hs):raise ValueError('Mixed height.')
    bound=max(next(iter(h)) for h in hs[:count])+max(next(iter(h)) for h in hs[count:])
    words=2
    while 64*words<=bound:words*=2
    def pack(vectors):
        offsets=[0];masks=[];scalars=[]
        for v in vectors:
            for w,c in sorted(v.items()):masks.append(w);scalars.append(c)
            offsets.append(len(masks))
        return ((ctypes.c_uint32*len(offsets))(*offsets),(ctypes.c_uint32*len(masks))(*masks),(ctypes.c_uint32*len(scalars))(*scalars))
    rp,cp=pack(rows),pack(columns);attempts=[]
    while words<=128:
        orders=(ctypes.c_int64*count)();values=(ctypes.c_uint64*(2*count*count*words))() if coefficients else None
        status=library(words).coherent_hot_compound_f4(modes,count,*rp,*cp,orders,values);attempts.append([64*words,status])
        if status==0:break
        if status!=6:raise ArithmeticError(f'F4 compound error {status} at precision {64*words}.')
        words*=2
    else:raise ArithmeticError('Incomplete profile: precision exhausted.')
    result={'orders_t':sorted(map(int,orders)),'size':count,'degree_bound_t':bound,'precision_attempts':attempts,'endpoint_entries':count*count}
    if coefficients:
        result['matrix']=[[(sum(int(values[(2*(i*count+j))*words+k])<<(64*k) for k in range(words)),sum(int(values[(2*(i*count+j)+1)*words+k])<<(64*k) for k in range(words))) for j in range(count)] for i in range(count)]
    return result


def tensor_vector(x,y,modes):return {a|(b<<modes):mul(c,d) for a,c in x.items() for b,d in y.items()}


def profiles(embedding,retraction,modes,keep_blocks=False):
    hs,charges=validate_maps(embedding,retraction,modes);q=len(hs);cuts=[];start=time.monotonic()
    for cut in ('zero','infinity'):
        blocks=[];orders=[]
        for charge in sorted({a+b for a in charges for b in charges}):
            pairs=[(x,y) for x in range(q) for y in range(q) if charges[x]+charges[y]==charge]
            rows=[tensor_vector(retraction[x],retraction[y] if cut=='zero' else embedding[y],modes) for x,y in pairs]
            cols=[tensor_vector(embedding[y] if cut=='zero' else retraction[y],embedding[x],modes) for x,y in pairs]
            block=compound_block(modes,rows,cols,keep_blocks);blocks.append(block);orders+=block['orders_t']
        if len(orders)!=q*q:raise ArithmeticError('Incomplete cut.')
        cuts.append({'cut':cut,'orders_t':sorted(orders),'blocks':blocks})
    return {'q':q,'modes':modes,'heights':hs,'charges':charges,'zero_orders_t':cuts[0]['orders_t'],'infinity_orders_t':cuts[1]['orders_t'],'cuts':cuts,'seconds':time.monotonic()-start}
