"""Coordinate alphabet restrictions of the strict reflection Gaussian."""
import argparse,itertools,json,math,random,time
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.special import logsumexp
from reflection_dual_rail_gaussian import tensor,physical_bit_energies,logical_heights
from native_matrix_smith_f4 import MatrixOracle


class Search:
    def __init__(self,M):
        self.M=M;self.seed=tensor(M);self.q=self.seed['q_x'];self.m=M//2
        self.energy=[sum(e for i,e in enumerate(physical_bit_energies(M))if state>>i&1)for state in range(self.q)]
    def profile(self,C,cut):
        q=self.q;values=self.seed['values'];pairs=list(itertools.product(C,repeat=2));matrix=[]
        for a,b in pairs:
            row=[]
            for c,d in pairs:
                x,y,u,v=(a,b,c,d)if cut==0 else(a,c,b,d)
                row.append(values[x+q*y+q*q*u+q*q*q*v])
            matrix.append(row)
        oracle=MatrixOracle(matrix,len(pairs),precision=512)
        rowweights=[(self.energy[a]+self.energy[b])/2 if cut==0 else self.energy[b]/2 for a,b in pairs]
        colweights=[0 if cut==0 else self.energy[b]/2 for a,b in pairs]
        result=oracle.profile(rowweights,colweights,infinity=bool(cut))
        orders=[round(v)for v in result['orders']]
        assert max(abs(a-b)for a,b in zip(orders,result['orders']))<1e-7
        return sorted(orders)
    def evaluate(self,C,name='subset'):
        C=tuple(sorted(C));nu=self.profile(C,0);mu=self.profile(C,1);E=np.array([self.energy[c]for c in C],float)
        aa=np.array(nu,float);bb=np.array(mu,float)
        def rho(t):
            raw=logsumexp(-t*aa)+logsumexp(-t*bb);gaps=t*(E-min(E));constant=raw+2*t*min(E)
            if constant<1e-12:return 1.
            def diff(R):return constant-2*R*logsumexp(-gaps/R)
            high=2.
            while diff(high)>0:high*=2
            return brentq(diff,1e-6,high,xtol=1e-11)/2
        result=minimize_scalar(rho,bounds=(1e-5,2),method='bounded',options={'xatol':1e-9})
        return {'name':name,'alphabet':C,'rho':float(result.fun),'log_scale':float(result.x),'profiles':[dict(Counter(nu)),dict(Counter(mu))],'energies':list(E.astype(int).tolist())}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--interval',type=int,default=8);ap.add_argument('--random',type=int,default=0);ap.add_argument('--out',required=True);args=ap.parse_args()
    start=time.time();S=Search(args.interval);rng=random.Random(20261031);m=S.m;q=S.q;results=[];seen=set()
    candidates=[(tuple(range(q)),'full')]
    for parity in(0,1):candidates.append((tuple(x for x in range(q)if x.bit_count()%2==parity),'parity'))
    for weight in range(1,m):candidates.append((tuple(x for x in range(q)if x.bit_count()<=weight),'weight-ball'))
    for cap in range(min(S.energy),max(S.energy)+1):
        C=tuple(x for x in range(q)if S.energy[x]<=cap)
        if 2<=len(C)<q:candidates.append((C,'energy-ball'))
    for mask in range(1,q):
        candidates.append((tuple(x for x in range(q)if(x&mask).bit_count()%2==0),'hyperplane'))
    for size in range(5,q):
        for _ in range(args.random):candidates.append(((0,)+tuple(rng.sample(range(1,q),size-1)),'random'))
    for C,name in candidates:
        C=tuple(sorted(C))
        if C in seen or len(C)<2:continue
        seen.add(C);result=S.evaluate(C,name);results.append(result)
        best=min(x['rho']for x in results)
        print(json.dumps({'number':len(results),'elapsed':time.time()-start,'best':best,**{k:v for k,v in result.items()if k not in('profiles','energies')}}),flush=True)
        sources=[Path(__file__),Path(__file__).with_name('reflection_dual_rail_gaussian.py'),
                 Path(__file__).with_name('native_matrix_smith_f4.py'),Path(__file__).with_name('native_matrix_smith_f4.cpp'),
                 Path(__file__).with_name('native_rectangular_smith_f4.cpp')]
        Path(args.out).write_text(json.dumps({'args':vars(args),'elapsed':time.time()-start,'results':results,
                                            'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest()for p in sources},
                                            'scope':'finite exact binary polynomial profiles with exploratory numerical type-optimized pressure; no exclusion outside the listed alphabets and natural gauges'},indent=2)+'\n')


if __name__=='__main__':main()
