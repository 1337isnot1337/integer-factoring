#!/usr/bin/env python3
"""Exact F16 coherent retractions and complete two-cut occupation profiles.

Vectors are sparse dictionaries {occupation_word: F16_coefficient}.
F16 is F4[z]/(z^2+z+omega), with four binary coefficient planes.
The embedded F4 uses omega^2=omega+1. No conjugation is used.
Derived from coherent_hot_compounds_f4.py; binary minors reuse that engine.
"""
from __future__ import annotations
import ctypes
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import os, platform, shlex, shutil, subprocess, tempfile, time
from coherent_hot_compounds import height


from tensor_orthogonal_scaling import TowerField

FIELD = TowerField(2)
MUL = [[FIELD.mul(a,b) for b in range(16)] for a in range(16)]
INV = [0] + [FIELD.inv(a) for a in range(1,16)]


def mul(a,b):
    return MUL[a][b]


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
                inv=INV[v[k]];pivots[k]={i:mul(inv,a) for i,a in v.items()};break
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
        if not p or not l or any(not isinstance(w,int) or not 0<=w<1<<modes or not 1<=c<16 for v in (p,l) for w,c in v.items()):raise ValueError('Invalid sparse F16 vector.')
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
    dependencies=[source,source.with_name('coherent_hot_compounds_f4.cpp'),source.with_name('selected_gaussian_coefficients.cpp'),source.with_name('native_rectangular_smith_f4.cpp')]
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or not shutil.which(compiler[0]):raise RuntimeError('C++17 compiler required.')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={words}']
    if platform.machine().lower() in ('x86_64','amd64') and 'pclmulqdq' in Path('/proc/cpuinfo').read_text():flags+=['-mpclmul','-msse2']
    digest=sha256(b''.join(p.read_bytes() for p in dependencies)+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-coherent-f16-'+str(os.getuid()));cache.mkdir(mode=0o700,exist_ok=True)
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
    lib.coherent_hot_compound_f16.argtypes=[ctypes.c_int,ctypes.c_int,u32,u32,u32,u32,u32,u32,ctypes.POINTER(ctypes.c_int64),u64]
    lib.coherent_hot_compound_f16.restype=ctypes.c_int
    lib.coherent_hot_compound_f16_capped.argtypes=lib.coherent_hot_compound_f16.argtypes+[ctypes.c_int]
    lib.coherent_hot_compound_f16_capped.restype=ctypes.c_int
    lib.tensor_f16_scalar_control.argtypes=[ctypes.c_uint,ctypes.c_uint]
    lib.tensor_f16_scalar_control.restype=ctypes.c_int
    return lib


def compound_block(modes,rows,columns,coefficients=False,order_cap=None):
    count=len(rows)
    if order_cap is not None and (type(order_cap)is not int or order_cap<1):raise ValueError('Expected a positive integer Smith-order cap.')
    if not 1<=modes<=16 or not 1<=count<=4096 or len(columns)!=count:raise ValueError('Invalid block size.')
    if any(not v or any(not isinstance(w,int) or not 0<=w<1<<(2*modes) or not 1<=c<16 for w,c in v.items()) for v in rows+columns):raise ValueError('Invalid vector.')
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
        orders=(ctypes.c_int64*count)();values=(ctypes.c_uint64*(4*count*count*words))() if coefficients else None
        if order_cap is None:
            status=library(words).coherent_hot_compound_f16(modes,count,*rp,*cp,orders,values)
        else:
            status=library(words).coherent_hot_compound_f16_capped(modes,count,*rp,*cp,orders,values,order_cap)
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


def tensor_vector(x,y,modes):return {a|(b<<modes):mul(c,d) for a,c in x.items() for b,d in y.items()}


def profiles(embedding,retraction,modes,keep_blocks=False,order_cap=None):
    hs,charges=validate_maps(embedding,retraction,modes);q=len(hs);cuts=[];start=time.monotonic()
    for cut in ('zero','infinity'):
        blocks=[];orders=[]
        for charge in sorted({a+b for a in charges for b in charges}):
            pairs=[(x,y) for x in range(q) for y in range(q) if charges[x]+charges[y]==charge]
            rows=[tensor_vector(retraction[x],retraction[y] if cut=='zero' else embedding[y],modes) for x,y in pairs]
            cols=[tensor_vector(embedding[y] if cut=='zero' else retraction[y],embedding[x],modes) for x,y in pairs]
            block=compound_block(modes,rows,cols,keep_blocks,order_cap);blocks.append(block);orders+=block['orders_t']
        if len(orders)!=q*q:raise ArithmeticError('Incomplete cut.')
        cuts.append({'cut':cut,'orders_t':sorted(orders),'blocks':blocks})
    result={'q':q,'modes':modes,'heights':hs,'charges':charges,'zero_orders_t':cuts[0]['orders_t'],'infinity_orders_t':cuts[1]['orders_t'],'cuts':cuts,'seconds':time.monotonic()-start}
    if order_cap is not None:
        result.update(order_cap=order_cap,profile_kind='exact capped orders min(actual_order, order_cap)',
                      pressure_interpretation='upper bound from pessimistic tail orders',
                      computed_smith_orders=sum(sum(value<order_cap for value in cut['orders_t']) for cut in cuts),
                      capped_tail_count=sum(sum(value==order_cap for value in cut['orders_t']) for cut in cuts))
    return result



from collections import defaultdict


def apply_operator(vector,operator):
    out={}
    for word,c in vector.items():
        for removed,a in operator.items():
            if word&removed==removed:
                target=word^removed;out[target]=out.get(target,0)^mul(c,a)
                if not out[target]:del out[target]
    return out


def differential_kernel(operators,modes,height_cap):
    if not 1<=modes<=16 or not operators or type(height_cap)is not int or height_cap<0:raise ValueError('Invalid kernel dimensions.')
    for op in operators:
        if not op or any(type(w)is not int or not 0<=w<1<<modes or not 1<=c<16 for w,c in op.items()):raise ValueError('Invalid F16 differential operator.')
        if len({height(w) for w in op})!=1 or len({w.bit_count() for w in op})!=1:raise ValueError('Operator must be homogeneous in height and occupation.')
    cells=defaultdict(list)
    for word in range(1<<modes):
        if height(word)<=height_cap:cells[height(word),word.bit_count()].append(word)
    generators=[];ambient=0;image_ranks=[]
    for cell,words in sorted(cells.items()):
        pivots={};ambient+=len(words)
        for word in words:
            pre={word:1};image={}
            for oi,op in enumerate(operators):
                for target,c in apply_operator(pre,op).items():image[(oi,target)]=c
            while image:
                k=max(image);factor=image[k]
                if k in pivots:
                    old_image,old_pre=pivots[k];image=axpy(image,factor,old_image);pre=axpy(pre,factor,old_pre)
                else:
                    inv=INV[factor];image={k:mul(inv,v) for k,v in image.items()};pre={k:mul(inv,v) for k,v in pre.items()};pivots[k]=(image,pre);break
            if not image:generators.append(pre)
        image_ranks.append({'height':cell[0],'charge':cell[1],'ambient':len(words),'image_rank':len(pivots)})
    pivots=basis(generators)
    for pivot in sorted(pivots):
        for other in pivots:
            if other!=pivot and pivot in pivots[other]:pivots[other]=axpy(pivots[other],pivots[other][pivot],pivots[pivot])
    words=sorted(pivots,key=lambda w:(height(w),w));p=[pivots[w] for w in words];l=[{w:1} for w in words]
    hs,charges=validate_maps(p,l,modes)
    for v in p:
        for op in operators:assert not apply_operator(v,op)
    assert len(p)==ambient-sum(c['image_rank'] for c in image_ranks)
    return {'modes':modes,'height_cap':height_cap,'operators':operators,'ambient_dimension':ambient,'q':len(p),'embedding':p,'retraction':l,'pivot_words':words,'heights':hs,'charges':charges,'cell_image_ranks':image_ranks}



def inverse(matrix):
    n=len(matrix)
    if any(len(row)!=n for row in matrix):raise ValueError('Expected a square matrix.')
    A=[row[:]+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return None
        A[pivot],A[j]=A[j],A[pivot];c=INV[A[j][j]];A[j]=[mul(c,x) for x in A[j]]
        for i in range(n):
            if i!=j:
                c=A[i][j];A[i]=[x^mul(c,y) for x,y in zip(A[i],A[j])]
    return [row[n:] for row in A]


def rank(matrix):
    return len(basis([{j:x for j,x in enumerate(row) if x} for row in matrix]))


def subset_gram_certificate(coefficients):
    if not coefficients or not coefficients[0]:raise ValueError('Expected nonempty coefficients.')
    m=len(coefficients[0])
    if any(len(row)!=m or any(x not in range(16) for x in row) for row in coefficients):raise ValueError('Expected a rectangular F16 matrix.')
    records=[]
    for mask in range(1<<m):
        A=[[x for j,x in enumerate(row) if mask>>j&1] for row in coefficients]
        G=[[0 for _ in A] for _ in A]
        for i,row in enumerate(A):
            for j,other in enumerate(A):
                for a,b in zip(row,other):G[i][j]^=mul(a,b)
        r=rank(A);g=rank(G);records.append({'active_mask':mask,'rank':r,'gram_rank':g})
    return {'all_subsets_nondegenerate':all(r['rank']==r['gram_rank'] for r in records),'subsets':records}


def orthogonal_retraction(code):
    """Return L=(P^T P)^(-1)P^T, or None for a degenerate code."""
    p=code['embedding'];l=[{} for _ in p];cells={}
    for i,cell in enumerate(zip(code['heights'],code['charges'])):cells.setdefault(cell,[]).append(i)
    for indices in cells.values():
        G=[[dot(p[i],p[j]) for j in indices] for i in indices];V=inverse(G)
        if V is None:return None
        for i,row in zip(indices,V):
            for c,j in zip(row,indices):l[i]=axpy(l[i],c,p[j])
    validate_maps(p,l,code['modes'])
    return l

