#!/usr/bin/env python3
"""Optional F4 rectangular weighted Smith search; Python uses only the stdlib.

Field entries are pairs of polynomial bit masks for 1 and omega, with
omega**2=omega+1. Dense entries have index x+A*(y+B*(u+B*v)), with alphabets A,B,B,A.
Parameters contain infinity heights, zero heights, and log(s)>0. Within
each end, list legs x,y,u,v and every nonzero letter; letter0 has height0.
Heights have already been multiplied by log(s). Cut1 uses C at infinity
and R at zero; cut0 reverses them. evaluate() takes independent X/V and
Y/U probability distributions. types() numerically optimizes both.

The adjacent C++17 source is compiled to a content-addressed temporary
cache on first use. Polynomial arithmetic is exact, but real pivot weights
and pressures are floating point: this is a search tool, not a certificate.
A missing pivot in the 512-coefficient window raises an error. A proposed
improvement needs independent exact replay and its endpoint proof.

A common denominator of degree d that is a unit at zero is restored once
with denominator_degree=d. This arithmetic does not establish rational
path admissibility. For typed rational paths set endpoint_check=False and
supply polynomiality and endpoint proofs separately.

One native seed is active at a time. Constructing or activating another
oracle makes previous instances inactive; call activate() before reuse.
"""
from __future__ import annotations
from native_weighted_smith_f4 import _pack
import argparse
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import threading

_PRECISION=512
_LIBRARY=None
_GENERATION=0
_LOCK=threading.RLock()
_DOUBLE=ctypes.POINTER(ctypes.c_double)
_U64=ctypes.POINTER(ctypes.c_uint64)
_U32=ctypes.POINTER(ctypes.c_uint32)
_I64=ctypes.POINTER(ctypes.c_int64)


def _library():
    global _LIBRARY
    if _LIBRARY is not None:return _LIBRARY
    source=Path(__file__).with_suffix('.cpp')
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional rectangular oracle needs a C++17 compiler; set CXX if needed.')
    flags=['-O3','-std=c++17','-shared','-fPIC']
    try:cpu=Path('/proc/cpuinfo').read_text()
    except OSError:cpu=''
    if platform.machine().lower() in ('x86_64','amd64','i386','i686') and 'pclmulqdq' in cpu:
        flags+=['-mpclmul','-msse2']
    digest=hashlib.sha256(source.read_bytes()+repr((compiler,flags,platform.machine())).encode()).hexdigest()[:24]
    cache=Path(tempfile.gettempdir())/('integer-factoring-rectangular-f4-'+str(os.getuid()))
    cache.mkdir(mode=0o700,parents=True,exist_ok=True);target=cache/(digest+'.so')
    if not target.exists():
        fd,tmp=tempfile.mkstemp(prefix=digest+'-',suffix='.so',dir=cache);os.close(fd)
        try:
            run=subprocess.run(compiler+flags+[str(source),'-o',tmp],capture_output=True,text=True)
            if run.returncode:raise RuntimeError('Rectangular oracle compilation failed:\n'+run.stderr[-12000:])
            os.replace(tmp,target)
        finally:
            if os.path.exists(tmp):os.unlink(tmp)
    lib=ctypes.CDLL(str(target))
    lib.rectangular_f4_set_seed.argtypes=[ctypes.c_int,ctypes.c_int,_U64,ctypes.c_int,ctypes.c_int]
    lib.rectangular_f4_set_seed.restype=ctypes.c_int
    lib.rectangular_f4_cost_gradient.argtypes=[_DOUBLE,_DOUBLE,_DOUBLE,ctypes.c_int,_DOUBLE,_DOUBLE,_DOUBLE]
    lib.rectangular_f4_cost_gradient.restype=ctypes.c_int
    lib.rectangular_f4_prefixes.argtypes=[_DOUBLE,ctypes.c_int,_DOUBLE,_I64,_U32,_U32]
    lib.rectangular_f4_prefixes.restype=ctypes.c_int
    _LIBRARY=lib;return lib


def _logsumexp(values):
    top=max(values)
    return top+math.log(sum(math.exp(v-top) for v in values))


class RectangularOracle:
    def __init__(self,values,A,B,denominator_degree=0,endpoint_check=True):
        if any(not isinstance(q,int) or not 2<=q<=64 for q in (A,B)):
            raise ValueError('A and B must be integers from 2 through 64.')
        if not isinstance(denominator_degree,int) or denominator_degree<0:
            raise ValueError('denominator_degree must be a nonnegative integer.')
        values=list(values)
        self.A=A;self.B=B;self.n=A*B;self.alphabets=[A,B,B,A]
        self.offsets=[];self.stride=0
        for q in self.alphabets:self.offsets.append(self.stride);self.stride+=q-1
        self.parameters=2*self.stride+1;self.denominator_degree=denominator_degree
        self._packed,self._words=_pack(values,A*A*B*B)
        self._endpoint_check=bool(endpoint_check);self.activate()

    def activate(self):
        global _GENERATION
        with _LOCK:
            _GENERATION+=1;self._generation=None
            status=_library().rectangular_f4_set_seed(self.A,self.B,self._packed,self._words,int(self._endpoint_check))
            if status:raise ArithmeticError(f'Rectangular seed rejected (status {status}).')
            self._generation=_GENERATION
        return self

    def _params(self,params,cut):
        x=list(map(float,params))
        if len(x)!=self.parameters or cut not in(0,1):raise ValueError('Wrong parameter length or cut.')
        if any(not math.isfinite(v) for v in x) or x[-1]<=0:raise ValueError('Finite parameters and log(s)>0 are required.')
        if max(map(abs,x[:-1]),default=0)/x[-1]>1_000_000:raise ValueError('Normalized heights exceed the safe conversion range.')
        return x,(ctypes.c_double*len(x))(*x)

    def evaluate(self,params,piX,piY,cut=1,want_spectra=False):
        x,native_x=self._params(params,cut);distributions=[]
        for p,q in ((piX,self.A),(piY,self.B)):
            p=list(map(float,p))
            if len(p)!=q or any(not math.isfinite(v) or v<0 for v in p) or abs(sum(p)-1)>1e-9:
                raise ValueError('Each paired type must be a finite probability distribution.')
            distributions.append((ctypes.c_double*q)(*p))
        value=ctypes.c_double();gradient=(ctypes.c_double*self.parameters)();orders=(ctypes.c_double*(2*self.n))()
        with _LOCK:
            if self._generation!=_GENERATION:raise RuntimeError('Another seed replaced this oracle; call activate() before reuse.')
            status=_library().rectangular_f4_cost_gradient(native_x,*distributions,cut,ctypes.byref(value),gradient,orders)
        if status:raise ArithmeticError(f'Full rectangular profile unresolved (status {status}).')
        cost=value.value-self.denominator_degree*x[-1];gradient=list(gradient);orders=list(orders)
        gradient[-1]-=self.denominator_degree
        orders[:self.n]=[v+self.denominator_degree*x[-1] for v in orders[:self.n]]
        if not all(math.isfinite(v) for v in [cost]+gradient+orders):raise ArithmeticError('Floating pressure overflowed.')
        return (cost,gradient,orders) if want_spectra else (cost,gradient)

    def prefixes(self,params,cut=1):
        """Attaining nested minors, encoded by newly selected row/column indices."""
        x,native_x=self._params(params,cut)
        orders=(ctypes.c_double*(2*self.n))()
        bases=(ctypes.c_int64*(2*self.n))()
        rows=(ctypes.c_uint32*(2*self.n))()
        cols=(ctypes.c_uint32*(2*self.n))()
        with _LOCK:
            if self._generation!=_GENERATION:raise RuntimeError('Another seed replaced this oracle; call activate() before reuse.')
            status=_library().rectangular_f4_prefixes(native_x,cut,orders,bases,rows,cols)
        if status:raise ArithmeticError(f'Prefix profile unresolved (status {status}).')
        orders,bases=list(orders),list(bases)
        for i in range(self.n):
            orders[i]+=self.denominator_degree*x[-1]
            bases[i]+=self.denominator_degree*(i+1)
        return {'orders':orders,'base_prefix':bases,'rows':list(rows),'cols':list(cols)}

    def types(self,params,cut=1):
        x,_=self._params(params,cut);GX=[0.]*self.A;GY=[0.]*self.B
        for side in range(2):
            for leg in range(4):
                g=GX if leg in(0,3) else GY
                for a in range(1,self.alphabets[leg]):g[a]+=x[side*self.stride+self.offsets[leg]+a-1]
        _,_,orders=self.evaluate(x,[1/self.A]*self.A,[1/self.B]*self.B,cut,True)
        minimumX=min(GX);minimumY=min(GY)
        constant=minimumX+minimumY+_logsumexp([-v for v in orders[:self.n]])+_logsumexp([-v for v in orders[self.n:]])
        def equation(R):
            return constant-R*(_logsumexp([-(v-minimumX)/R for v in GX])+_logsumexp([-(v-minimumY)/R for v in GY]))
        low,high=1e-12,8.
        if equation(low)<=0:raise ArithmeticError('Nonpositive paired Gibbs root; check exact profiles and endpoint hypotheses.')
        while equation(high)>0:
            high*=2
            if high>1e8:raise ArithmeticError('Paired Gibbs root was not bracketed.')
        for _ in range(90):
            middle=(low+high)/2
            if equation(middle)>0:low=middle
            else:high=middle
        R=(low+high)/2
        px=[math.exp(-v/R-_logsumexp([-u/R for u in GX])) for v in GX]
        py=[math.exp(-v/R-_logsumexp([-u/R for u in GY])) for v in GY]
        return R/2,px,py,orders,GX,GY


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed',type=Path);parser.add_argument('--params',type=Path)
    parser.add_argument('--log-s',type=float,default=.1);parser.add_argument('--cut',type=int,choices=(0,1),default=1)
    parser.add_argument('--output',type=Path);args=parser.parse_args()
    data=json.loads(args.seed.read_text());A=data.get('q_x',data.get('dimensions',[0])[0]);B=data.get('q_u',data.get('dimensions',[0,0])[1])
    oracle=RectangularOracle(data['values'],A,B,data.get('denominator_degree',0),data.get('strict_endpoint',True))
    params=json.loads(args.params.read_text()) if args.params else [0.]*(oracle.parameters-1)+[args.log_s]
    rho,px,py,orders,gx,gy=oracle.types(params,args.cut)
    result={'exploratory':True,'dimensions':[A,B,B,A],'rho':rho,'type_x':px,'type_y':py,
            'infinity_orders_times_log_s':orders[:oracle.n],'zero_orders_times_log_s':orders[oracle.n:],
            'G_X_times_log_s':gx,'G_Y_times_log_s':gy,'cut':args.cut}
    text=json.dumps(result,indent=2)
    if args.output:args.output.write_text(text+'\n')
    print(text)

if __name__=='__main__':main()
