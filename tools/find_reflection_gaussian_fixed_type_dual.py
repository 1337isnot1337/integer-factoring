#!/usr/bin/env python3
"""Discover and exactly validate a global fixed-type Gaussian gauge dual.

Floating linear programs select minors; rational row reduction and the
independent determinant/entropy checker validate the resulting certificate.
The construction covers the returned type only, not all physical types.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math
import numpy as np
from scipy.optimize import linprog
from scipy.special import expit
from reflection_unrestricted_gaussian import ReflectionGaussianPressure

def rref(matrix):
    work=[[F(value) for value in row] for row in matrix]
    pivots=[];row=0
    for col in range(len(work[0])):
        selected=next((i for i in range(row,len(work)) if work[i][col]),None)
        if selected is None:continue
        work[row],work[selected]=work[selected],work[row]
        value=work[row][col];work[row]=[x/value for x in work[row]]
        for i in range(len(work)):
            if i!=row and work[i][col]:
                value=work[i][col];work[i]=[x-value*y for x,y in zip(work[i],work[row])]
        pivots.append(col);row+=1
        if row==len(work):break
    return work,pivots

def collect(M,trials=300):
    P=ReflectionGaussianPressure(M);nat=P.natural(1.);rng=np.random.default_rng(98135)
    active={}
    for trial in range(trials):
        x=nat.copy()
        if trial:x[:-1]+=rng.normal(size=len(x)-1)*1e-4
        for side,(oracle,rows,cols) in enumerate(P.cuts):
            h=x[side*2*M:(side+1)*2*M]
            rw=np.r_[h[rows],-h[rows]];cw=np.r_[h[cols],-h[cols]]
            data=oracle.profile(rw,cw,1.,infinity=bool(side))
            coefficients=np.zeros(4*M,dtype=int)
            rowmask=colmask=0
            for k in range(1,M+1):
                r,c=data['rows'][k-1],data['cols'][k-1]
                rowmask|=1<<r;colmask|=1<<c
                coefficients[side*2*M+rows[r%M]]+=1 if r<M else -1
                coefficients[side*2*M+cols[c%M]]+=1 if c<M else -1
                base=data['base_prefix'][k-1]
                actual=int(base+coefficients@nat[:-1])
                expected=-sum(4*M-2-4*j for j in range(k))
                assert actual==expected,(side,k,actual,expected)
                key=(side,k,rowmask,colmask)
                active[key]={'side':side,'size':k,'rows':rowmask,'columns':colmask,
                             'valuation':base,'coefficients':coefficients.tolist()}
    return P,list(active.values())

def system(M,u,active):
    N=len(active);rows=6*M+1
    A=np.zeros((rows,N+M),dtype=int);b=[F(0)]*rows
    c=[(1-u**(4*M-2-4*j))/(2*(1+u**(4*M-2-4*j))) for j in range(M)]+[F(0)]
    for side in range(2):
        for k in range(M):b[side*M+k]=c[k]-c[k+1]
    for j,minor in enumerate(active):
        A[minor['side']*M+minor['size']-1,j]=1
        A[2*M:6*M,j]=minor['coefficients']
        A[-1,j]=minor['valuation']
    for side in range(2):
        for leg in range(4):
            for i in range(M//2):
                h=side*2*M+leg*(M//2)+i
                typ=(0 if leg in(0,3) else 1)*(M//2)+i
                A[2*M+h,N+typ]=-1
                b[2*M+h]=F(-1,2)
    b[-1]=F(-M*M)
    return A,b

def add_direction(M,P,active,direction,dual=None):
    nat=P.natural(1.);x=nat+1e-5*np.asarray(direction)
    known={(a['side'],a['size'],a['rows'],a['columns']) for a in active}
    added=0;largest=0.
    for side,(oracle,rows,cols) in enumerate(P.cuts):
        h=x[side*2*M:(side+1)*2*M]
        data=oracle.profile(np.r_[h[rows],-h[rows]],np.r_[h[cols],-h[cols]],x[-1],infinity=bool(side))
        coefficients=np.zeros(4*M,dtype=int);rowmask=colmask=0
        for k in range(1,M+1):
            r,c=data['rows'][k-1],data['cols'][k-1]
            rowmask|=1<<r;colmask|=1<<c
            coefficients[side*2*M+rows[r%M]]+=1 if r<M else -1
            coefficients[side*2*M+cols[c%M]]+=1 if c<M else -1
            base=data['base_prefix'][k-1]
            assert base+coefficients@nat[:-1]==-sum(4*M-2-4*j for j in range(k))
            key=(side,k,rowmask,colmask)
            reduced=(dual[side*M+k-1]+coefficients@dual[2*M:6*M]+base*dual[6*M]) if dual is not None else 1.
            largest=max(largest,float(reduced))
            if key not in known and reduced>1e-9:
                active.append({'side':side,'size':k,'rows':rowmask,'columns':colmask,
                               'valuation':base,'coefficients':coefficients.tolist()})
                known.add(key);added+=1
    return added,largest

def solve(M,u,active,P):
    ell=-math.log(float(u));x=P.natural(ell);rho=P.evaluate(x)
    energy=np.tile([4*(M-1-2*i) for i in range(M//2)],2)
    target=expit(-energy*ell/(2*rho))
    for round_index in range(240):
        N=len(active);A,b=system(M,u,active)
        symmetric=np.zeros((M//2,N+M),dtype=int)
        for i in range(M//2):symmetric[i,N+i]=1;symmetric[i,N+M//2+i]=-1
        A=np.r_[A,symmetric];b=b+[F(0)]*(M//2)
        ae=np.c_[A,np.zeros((len(b),M))]
        au=np.zeros((2*M,N+2*M));bu=np.r_[target,-target]
        for i in range(M):
            au[i,N+i]=1;au[i,N+M+i]=-1
            au[M+i,N+i]=-1;au[M+i,N+M+i]=-1
        obj=np.r_[np.zeros(N+M),np.ones(M)]
        fit=linprog(obj,A_ub=au,b_ub=bu,A_eq=ae,b_eq=np.array(list(map(float,b))),
                    bounds=[(0,None)]*N+[(0,1)]*M+[(0,None)]*M,method='highs',
                    options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
        if fit.success:
            y=fit.eqlin.marginals
            direction=-np.r_[y[2*M:6*M],y[6*M]]
            added,violation=add_direction(M,P,active,direction,y)
            print({'target_round':round_index,'distance':fit.fun,'added':added,'violation':violation},flush=True)
            if not added:break
            continue
        slack=np.zeros((len(b),4*M+1))
        for i in range(4*M+1):slack[2*M+i,i]=1
        elastic=linprog(np.r_[np.zeros(N+M),np.ones(8*M+2)],
                        A_eq=np.c_[A,slack,-slack],b_eq=np.array(list(map(float,b))),
                        bounds=[(0,None)]*N+[(0,1)]*M+[(0,None)]*(8*M+2),method='highs',
                        options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
        assert elastic.success,elastic.message
        y=elastic.eqlin.marginals
        direction=-np.r_[y[2*M:6*M],y[6*M]]
        added,violation=add_direction(M,P,active,direction,y)
        print({'dual_round':round_index,'residual':elastic.fun,'added':added},flush=True)
        if not added:break
    print({'M':M,'u':str(u),'minors':N,'feasible':fit.success,'message':fit.message},flush=True)
    if not fit.success:return None
    support=[i for i in range(N) if fit.x[i]>1e-11]+list(range(N,N+M))
    reduced,pivots=rref([list(row)+[value] for row,value in zip(A[:,support].tolist(),b)])
    assert len(support) not in pivots,'Exact active system inconsistent'
    pivots=[p for p in pivots if p<len(support)]
    free=[i for i in range(len(support)) if i not in pivots]
    solution=[F(0)]*len(support)
    for i in free:solution[i]=F(float(fit.x[support[i]])).limit_denominator(10**12)
    for row,pivot in enumerate(pivots):
        value=reduced[row][-1]
        for i in free:value-=reduced[row][i]*solution[i]
        solution[pivot]=value
    for row in range(len(b)):
        assert sum(F(int(A[row,j]))*v for j,v in zip(support,solution))==b[row]
    assert all(v>=0 for j,v in zip(support,solution) if j<N),min(v for j,v in zip(support,solution) if j<N)
    probs=[solution[support.index(N+i)] for i in range(M)]
    assert all(0<p<1 for p in probs)
    h=sum(-float(p)*math.log(float(p))-(1-float(p))*math.log1p(-float(p)) for p in probs)
    cost=2*sum(math.log1p(float(u**(4*i+2))) for i in range(M))+ell*sum(float(p)*e for p,e in zip(probs,energy))
    print({'p':[float(p) for p in probs],'target':target.tolist(),'support':len(support)-M,
           'rho_fixed':cost/(2*h),'target_distance':fit.fun},flush=True)
    return {'interval':M,'u':str(u),'probabilities':list(map(str,probs)),
            'dual':[{**active[j],'multiplier':str(v)} for j,v in zip(support,solution) if j<N and v],
            'rho_fixed_floating':cost/(2*h),'exact_stationarity':True,'all_collected_minors':N}

if __name__=='__main__':
    from check_reflection_gaussian_fixed_type_duals import validate_certificate
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval',dest='M',type=int,default=8)
    parser.add_argument('--trials',type=int,default=100)
    parser.add_argument('--exp-negative-scale',dest='u',default='29/36')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    model,active=collect(args.M,args.trials)
    result=solve(args.M,F(args.u),active,model)
    if result is None:
        raise ArithmeticError('No exact certificate was found in the collected minor family.')
    certificate={'interval':args.M,'exp_negative_scale':result['u'],
                 'probabilities':result['probabilities'],
                 'dual_minors':[{'side':minor['side'],'size':minor['size'],
                                'row_mask':hex(minor['rows']),'column_mask':hex(minor['columns']),
                                'valuation':minor['valuation'],'multiplier':minor['multiplier']}
                               for minor in result['dual']]}
    verification=validate_certificate(certificate)
    report={'status':'Exact rational fixed-type dual; independent determinant and entropy validation passed.',
            'certificates':[certificate],'verification':[verification],
            'discovery':{'random_trials':args.trials,'random_seed':98135,
                         'collected_minors':result['all_collected_minors']}}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
