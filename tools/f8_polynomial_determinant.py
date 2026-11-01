"""Exact ordinary polynomial determinant over F8, with degree guards."""
import ctypes,os,subprocess,tempfile,platform,shlex,shutil
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import numpy as np
@lru_cache(None)
def library(precision):
    root=Path(__file__).parent
    names=['f8_polynomial_determinant.cpp','coherent_hot_compounds_f8.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp']
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or not shutil.which(compiler[0]):raise RuntimeError('A C++17 compiler is required')
    flags=['-O3','-std=c++17','-shared','-fPIC',f'-DF4_SERIES_WORDS={precision//64}']
    try:cpu=Path('/proc/cpuinfo').read_text()
    except OSError:cpu=''
    if platform.machine().lower() in ('x86_64','amd64') and 'pclmulqdq' in cpu:flags+=['-mpclmul','-msse2']
    h=sha256(b''.join((root/n).read_bytes()for n in names)+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-f8-poly-'+str(os.getuid()));cache.mkdir(mode=0o700,exist_ok=True)
    path=cache/(h+'.so')
    if not path.exists():
        fd,temp=tempfile.mkstemp(suffix='.so',dir=cache);os.close(fd)
        try:
            r=subprocess.run(compiler+flags+[str(root/names[0]),'-o',temp],capture_output=True,text=True)
            if r.returncode:raise RuntimeError(r.stderr)
            os.replace(temp,path)
        finally:
            if os.path.exists(temp):os.unlink(temp)
    lib=ctypes.CDLL(str(path));lib.determinant_f8.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_uint64),ctypes.c_int,ctypes.POINTER(ctypes.c_uint64)];lib.determinant_f8.restype=ctypes.c_int;return lib

def determinant(C):
    C=np.asarray(C)
    if C.ndim!=3 or C.shape[1]!=C.shape[2] or C.dtype.kind not in 'ui' or np.any(C<0)or np.any(C>7):raise ValueError('Use an integer F8 polynomial matrix')
    D=C.shape[0]-1;n=C.shape[1]
    if D<0 or not 1<=n<=128:raise ValueError('Invalid dimensions')
    precision=128
    while precision<=2*n*D:precision*=2
    if precision>32768:raise ValueError('Degree guard exceeds supported precision')
    nw=max(1,(D+64)//64);packed=np.zeros((n,n,3,nw),dtype=np.uint64)
    for k in range(D+1):
        for a in range(3):packed[:,:,a,k//64]|=((C[k]>>a)&1).astype(np.uint64)<<np.uint64(k%64)
    words=precision//64;out=(ctypes.c_uint64*(3*words))();status=library(precision).determinant_f8(n,D,packed.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),nw,out)
    if status==3:return []
    if status:raise ArithmeticError('Polynomial determinant status '+str(status))
    planes=[sum(int(out[a*words+w])<<(64*w)for w in range(words))for a in range(3)]
    return [sum(((p>>k)&1)<<a for a,p in enumerate(planes))for k in range(max(p.bit_length()for p in planes))]
