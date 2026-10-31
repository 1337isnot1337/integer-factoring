#!/usr/bin/env python3
"""Finite natural-gauge screen for growing physical frames; no global exclusion."""
import json,math,time,argparse
from pathlib import Path
from itertools import combinations,product
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from reflection_gaussian_physical_frames import frame_matrix,matrix_oracle,logical_heights

def complete_factors(M,gates):
    profiles=[]
    for side in(0,1):
        b,rows,cols=frame_matrix(M,side,gates);h=np.array(logical_heights(M,'zero' if side==0 else 'infinity'))
        data=matrix_oracle(b).profile(np.r_[h[rows],-h[rows]],np.r_[h[cols],-h[cols]],infinity=bool(side))
        values=data['orders'];assert all(x==round(x) for x in values)
        values=list(map(round,values));assert values==[-v for v in values[::-1]]
        profiles.append(values[M:])
    return profiles

def root_for(M,profiles,ell):
    factors=np.array(profiles,dtype=float).reshape(-1)
    physical=np.tile(np.array([4*(M-1-2*i) for i in range(M//2)]),2)*ell
    log_b=-2*M*M*ell+np.logaddexp(factors*ell/2,-factors*ell/2).sum()
    def equation(r):return log_b-r*np.logaddexp(0,-physical/r).sum()
    if equation(1e-9)<=0:return 0.,float(equation(math.sqrt(2)))
    upper=4.
    while equation(upper)>0:upper*=2
    rate=brentq(equation,1e-9,upper,xtol=1e-12)
    return rate/2,float(equation(math.sqrt(2)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--M',type=int,default=16);p.add_argument('--trials',type=int,default=1500);p.add_argument('--seed',type=int,default=75013);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();M=args.M;rng=np.random.default_rng(args.seed);start=time.monotonic()
    baseline=complete_factors(M,[]);opt=minimize_scalar(lambda x:root_for(M,baseline,math.exp(x))[0],bounds=(-9,1),method='bounded')
    ell=math.exp(opt.x);base=np.sort(np.array(baseline).reshape(-1))[::-1];prefix=np.cumsum(base)
    basephi=float(np.logaddexp(base*ell/2,-base*ell/2).sum())
    seen={tuple(base)};records=[];best={'rho':float(opt.fun),'scale':ell,'gates':[],'profiles':baseline};failures=0;dominating=0;neutral=0
    candidates=[]
    powers=sorted(set((-4*M,-2*M,-M,-4,-2,0,2,4,M,2*M,4*M)))
    for leg,kind,matching,pair,gamma in product(range(4),('create','annihilate'),('reverse','adjacent'),combinations(powers,2),(1,2)):
        candidates.append([{'leg':leg,'kind':kind,'matching':matching,'powers':pair,'gamma':gamma}])
    rng.shuffle(candidates)
    def random_gate():
        i,j=rng.choice(len(powers),2,replace=False)
        return {'leg':int(rng.integers(4)),'kind':('create','annihilate')[rng.integers(2)],
                'matching':('reverse','adjacent','all')[rng.integers(3)],'powers':(powers[i],powers[j]),'gamma':int(rng.integers(1,4))}
    def save(count):
        report={'status':'exploratory exact-profile screen at natural physical heights',
                'M':M,'baseline_rho':float(opt.fun),'best':best,'cases':count,'unique_combined_profiles':len(seen),
                'majorization_dominating':dominating,'neutral':neutral,'failures':failures,'nondominating_profiles':records,
                'scope':'Only the tested finite words at the natural physical heights; no improved bound or universal obstruction.',
                'seconds':time.monotonic()-start,'seed':args.seed}
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    for trial in range(args.trials):
        if trial<args.trials//2:gates=candidates[trial%len(candidates)]
        else:
            gates=[random_gate() for _ in range(int(rng.integers(2,5)))]
            if trial%3==0:
                # A natural-gauge integral prefix changes the flag before
                # the remaining nonregular frame acts.
                leg=int(rng.integers(4));sgn=-1 if leg<2 else 1
                gates=[{'leg':leg,'kind':'create','matching':'reverse','powers':(0,sgn*2*M),'gamma':1}]+gates
        try:profiles=complete_factors(M,gates)
        except ArithmeticError:
            failures+=1;continue
        current=np.sort(np.array(profiles).reshape(-1))[::-1];key=tuple(current)
        same=np.array_equal(current,base);neutral+=int(same)
        dominant=bool(np.all(np.cumsum(current)>=prefix));dominating+=int(dominant)
        if key not in seen:
            seen.add(key)
            if not dominant:
                fit=minimize_scalar(lambda x:root_for(M,profiles,math.exp(x))[0],bounds=(-9,1),method='bounded')
                scale=math.exp(fit.x);rate,gap=root_for(M,profiles,scale)
                record={'rho':rate,'scale':scale,'critical_pressure':gap,'gates':gates,'profiles':profiles}
                records.append(record)
                if rate<best['rho']:best=record;print({'trial':trial,'new_best':rate,'gates':gates},flush=True)
        if trial%100==0:
            save(trial+1);print({'trial':trial,'best':best['rho'],'unique':len(seen),'nondominating':len(records),'neutral':neutral,'seconds':time.monotonic()-start},flush=True)
    save(args.trials)

if __name__=='__main__':main()
