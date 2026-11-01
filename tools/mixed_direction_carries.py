"""Polynomial Teichmuller sums with forward and reverse Frobenius rows.

Uses the canonical unramified ring and scalar moments. This experimental
extension has full rows (one relative column), explicit polynomial factors,
and no denominator or arbitrary circuit interface.
"""
from itertools import product
from math import prod


def polynomial_sum_mixed(ring,width,variables,factors,directions,nonzero=False):
    """At row rho, variable i is evaluated at U_i^(2^(epsilon_i*rho mod w))."""
    if width<1 or ring.degree%width or variables<1:
        raise ValueError('Incompatible ring, width or variable count')
    if len(directions)!=variables or any(s not in (-1,1) for s in directions):
        raise ValueError('One sign per summation variable is required')
    mask=[nonzero]*variables if isinstance(nonzero,bool) else list(nonzero)
    if len(mask)!=variables:raise ValueError('Bad punctured-domain mask')
    zero=(0,)*variables
    bounds=[0]*variables;degree=0
    for row,terms in factors:
        if not 0<=row<width:raise ValueError('Row outside the Frobenius cycle')
        for c,counts in terms:
            if len(counts)!=variables or any(n<0 for n in counts):
                raise ValueError('Invalid ordinary row monomial')
        degree+=max([0]+[sum(n) for c,n in terms if c])
        for i in range(variables):bounds[i]+=max([0]+[n[i] for c,n in terms if c])
    local=[]
    for rho in range(width):
        values={zero:1}
        for row,terms in factors:
            if row!=rho:continue
            updated={}
            for old,a in values.items():
                for c,counts in terms:
                    out=tuple(x+y for x,y in zip(old,counts))
                    updated[out]=ring.add(updated.get(out,0),ring.mul(a,c))
            values={k:c for k,c in updated.items() if c}
        local.append(values)
    answer=0;transitions=0;peak=0;visited=0
    q=1 << width
    for initial in product(*(range(b+1) for b in bounds)):
        weight=1
        for i,a in enumerate(initial):
            moment=q if not mask[i] and a==0 else q-1
            weight=ring.mul(weight,ring.integer(moment))
        if not weight:continue
        visited+=1
        states={initial:weight}
        for rho in range(width):
            updated={}
            for incoming,a in states.items():
                for counts,c in local[rho].items():
                    out=[]
                    for i,(carry,n,sign) in enumerate(zip(incoming,counts,directions)):
                        if sign==1:
                            if (carry+n)&1:break
                            target=(carry+n)//2
                        else:target=2*carry-n
                        if not 0<=target<=bounds[i]:break
                        out.append(target)
                    else:
                        target=tuple(out)
                        updated[target]=ring.add(updated.get(target,0),ring.mul(a,c))
                        transitions+=1
            states={k:c for k,c in updated.items() if c}
            peak=max(peak,len(states))
            if not states:break
        answer=ring.add(answer,states.get(initial,0))
    return answer,dict(width=width,variables=variables,directions=list(directions),
                       precision=ring.precision,coefficient_degree=ring.degree,
                       degree_budget=degree,carry_bounds=bounds,
                       carry_box_size=prod(b+1 for b in bounds),
                       nonzero_initial_states=visited,peak_states=peak,
                       accepted_transitions=transitions,local_terms=list(map(len,local)))
