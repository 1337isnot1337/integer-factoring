#!/usr/bin/env python3
"""Certified low-order coherent profiles, with explicit unknown-tail counts.

Every returned lower-order bound is <= its true Smith order. Replacing the
unknown tail by cutoff+1 gives an UPPER channel partition function, not a
complete Smith profile. The coefficient and elimination arithmetic is exact.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import ctypes,os,platform,shlex,shutil,subprocess,tempfile,time
from coherent_hot_compounds_f4 import height,validate_maps,tensor_vector

@lru_cache(None)
def library(words,cutoff):
    source=Path(__file__).with_suffix('.cpp')
    dependencies=[source,source.with_name('selected_gaussian_coefficients.cpp'),source.with_name('native_rectangular_smith_f4.cpp')]
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or not shutil.which(compiler[0]):raise RuntimeError('C++17 compiler required.')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={words}',f'-DHOT_PROFILE_CUTOFF={cutoff}']
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



def partial_block(modes,rows,columns,cutoff):
    count=len(rows)
    if type(cutoff)is not int or cutoff<0:raise ValueError('Use a nonnegative integer cutoff.')
    if not 1<=modes<=16 or not count or len(columns)!=count:raise ValueError('Invalid block size.')
    if any(not v or any(type(w)is not int or not 0<=w<1<<(2*modes) or c not in (1,2,3) for w,c in v.items()) for v in rows+columns):raise ValueError('Invalid vector.')
    if len({w.bit_count() for v in rows+columns for w in v})!=1:raise ValueError('Mixed occupation.')
    low=(1<<modes)-1
    hs=[{height(w&low)+height(w>>modes) for w in v} for v in rows+columns]
    if any(len(h)!=1 for h in hs):raise ValueError('Mixed height.')
    hs=[next(iter(h)) for h in hs]
    minimum=(min(hs[:count])+min(hs[count:]))//2
    if minimum>cutoff:return {'size':count,'cutoff':cutoff,'retained_orders_t':[],'tail_count':count,'skipped_by_minimum_entry_order':minimum}
    if count>4096:raise ValueError('Unskipped block exceeds supported order4096.')
    words=1
    while 64*words<=cutoff:words*=2
    def pack(vectors):
        offsets=[0];masks=[];scalars=[]
        for v in vectors:
            for w,c in sorted(v.items()):masks.append(w);scalars.append(c)
            offsets.append(len(masks))
        return ((ctypes.c_uint32*len(offsets))(*offsets),(ctypes.c_uint32*len(masks))(*masks),(ctypes.c_uint32*len(scalars))(*scalars))
    rp,cp=pack(rows),pack(columns);orders=(ctypes.c_int64*count)()
    status=library(words,cutoff).coherent_hot_compound_f4(modes,count,*rp,*cp,orders,None)
    if status:raise ArithmeticError(f'Partial F4 compound status {status}.')
    retained=sorted(int(v) for v in orders if v<=cutoff)
    assert all(0<=v<=cutoff for v in retained)
    return {'size':count,'cutoff':cutoff,'retained_orders_t':retained,'tail_count':count-len(retained),'precision_bits':64*words}


def partial_profiles(p,l,modes,cutoff,progress=None):
    hs,charges=validate_maps(p,l,modes);q=len(hs);cuts=[];start=time.monotonic()
    for cut in ('zero','infinity'):
        blocks=[];orders=[];tail=0
        for charge in sorted({a+b for a in charges for b in charges}):
            pairs=[(x,y) for x in range(q) for y in range(q) if charges[x]+charges[y]==charge]
            rows=[tensor_vector(l[x],l[y] if cut=='zero' else p[y],modes) for x,y in pairs]
            cols=[tensor_vector(p[y] if cut=='zero' else l[y],p[x],modes) for x,y in pairs]
            block=partial_block(modes,rows,cols,cutoff);blocks.append(block);orders+=block['retained_orders_t'];tail+=block['tail_count']
            if progress:progress({'cut':cut,'charge':charge,'size':len(pairs),'retained':len(block['retained_orders_t']),'seconds':time.monotonic()-start})
        assert len(orders)+tail==q*q
        cuts.append({'cut':cut,'retained_orders_t':sorted(orders),'tail_count':tail,'blocks':blocks})
    return {'q':q,'modes':modes,'heights':hs,'charges':charges,'cutoff':cutoff,'cuts':cuts,'seconds':time.monotonic()-start,'complete_profile_claimed':False}
