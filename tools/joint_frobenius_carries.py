"""Exact joint coefficients and traces of synchronized Frobenius products.

A factor is (row, terms). Each term is (coefficient, counts), with counts
indexed by variable then relative column. Field operations are explicit.
"""
from itertools import product
from math import gcd


def synchronization(width, variable_count, factor_count, incidences):
    """Maximum t|width with row[f]-shift[v] = exponent (mod t).

    Incidences are (factor, variable, binary_exponent); parallel edges matter.
    Returns t, variable shifts, factor rows, and integer cycle residuals.
    """
    if width < 1: raise ValueError('Positive field width required')
    size=variable_count+factor_count; edges=[[] for _ in range(size)]
    for f,v,e in incidences:
        if not (0<=f<factor_count and 0<=v<variable_count):raise ValueError('Bad incidence')
        a=v;b=variable_count+f;edges[a].append((b,e));edges[b].append((a,-e))
    potentials=[None]*size
    for root in range(size):
        if potentials[root] is not None:continue
        potentials[root]=0;stack=[root]
        while stack:
            a=stack.pop()
            for b,e in edges[a]:
                if potentials[b] is None:potentials[b]=potentials[a]+e;stack.append(b)
    residuals=[potentials[variable_count+f]-potentials[v]-e for f,v,e in incidences]
    t=width
    for e in residuals:t=gcd(t,e)
    shifts=[x%t for x in potentials[:variable_count]]
    rows=[x%t for x in potentials[variable_count:]]
    assert all((e+shifts[v]-rows[f])%t==0 for f,v,e in incidences)
    return t,shifts,rows,residuals


def coefficient(field, variables, t, factors, targets):
    if variables<1 or t<1 or len(targets)!=variables or min(targets)<0:raise ValueError('Bad dimensions')
    dimension=next((len(h) for _,terms in factors for _,h in terms),variables)
    if dimension%variables:raise ValueError('Bad count vector')
    D=dimension//variables
    S=sum(max([0]+[sum(h)for c,h in terms if c]) for _,terms in factors)
    degrees=[sum(max([0]+[sum(h[v*D+j]<<(row+t*j)for j in range(D))for c,h in terms if c])for row,terms in factors)for v in range(variables)]
    if any(k>d for k,d in zip(targets,degrees)):return 0
    J=max(1,(max(degrees+[1]).bit_length()+t-1)//t)
    zero=(0,)*(variables*J);local=[]
    for rho in range(t):
        values={zero:1}
        for row,terms in factors:
            if not 0<=row<t:raise ValueError('Unsynchronized row')
            if row!=rho:continue
            nxt={}
            for old,a in values.items():
                for b,h in terms:
                    if len(h)!=dimension or min(h)<0:raise ValueError('Bad term')
                    if not b:continue
                    out=tuple(old[v*J+j]+(h[v*D+j]if j<D else 0)for v in range(variables)for j in range(J))
                    nxt[out]=nxt.get(out,0)^field.mul(a,b)
            values={h:c for h,c in nxt.items()if c}
        local.append(values)
    # Every carry is bounded by S. At the first row of each column use
    # the sharper bound from all lower-position contributions before carries.
    bounds=[]
    for v in range(variables):
        for j in range(1,J):
            lower=sum(max([0]+[sum(h[v*D+k]<<(row+t*k)for k in range(min(D,j)))for c,h in terms if c])for row,terms in factors)
            bounds.append(min(S,lower>>(t*j)))
    answer=0;cache={}
    for tail in product(*(range(b+1)for b in bounds)):
        initial=tuple(a for v in range(variables)for a in ((0,)+tail[v*(J-1):(v+1)*(J-1)]))
        final=tuple(a for v in range(variables)for a in (tail[v*(J-1):(v+1)*(J-1)]+(0,)))
        states={initial:1}
        for rho in range(t):
            nxt={}
            for incoming,a in states.items():
                key=(rho,incoming)
                if key not in cache:
                    transitions=[]
                    for h,b in local[rho].items():
                        total=tuple(x+y for x,y in zip(incoming,h))
                        if all((total[v*J+j]&1)==((targets[v]>>(rho+t*j))&1)for v in range(variables)for j in range(J)):
                            outgoing=tuple(x>>1 for x in total);assert max(outgoing)<=S
                            transitions.append((outgoing,b))
                    cache[key]=transitions
                for outgoing,b in cache[key]:nxt[outgoing]=nxt.get(outgoing,0)^field.mul(a,b)
            states={h:c for h,c in nxt.items()if c}
        answer^=states.get(final,0)
    return answer


def frobenius_factor(field,width,t,variables,factor,power):
    row,terms=factor;D=width//t;new=[]
    for c,h in terms:
        out=[0]*(variables*D)
        for v in range(variables):
            for j in range(D):
                exponent=(row+t*j+power)%width
                out[v*D+exponent//t]+=h[v*D+j]
        new.append((field.power(c,1<<power),tuple(out)))
    return ((row+power)%t,new)


def joint_trace_coefficients(field,field_degree,width,t,variables,numerators,denominators,nonzero=False):
    if width<1 or field_degree%width or width%t:raise ValueError('Incompatible fields or row modulus')
    factors=list(numerators)
    for factor in denominators:
        for h in [0]+list(range(2,field_degree+1)):
            factors.append(frobenius_factor(field,width,t,variables,factor,h))
    D=width//t;b=1<<width
    degrees=[sum(max([0]+[sum(h[v*D+j]<<(row+t*j)for j in range(D))for c,h in terms if c])for row,terms in factors)for v in range(variables)]
    punctured=[nonzero]*variables if isinstance(nonzero,bool)else list(nonzero)
    if len(punctured)!=variables:raise ValueError('Bad domain mask')
    answer=0
    for target in product(*(range(0 if nz else b-1,d+1,b-1)for d,nz in zip(degrees,punctured))):
        answer^=coefficient(field,variables,t,factors,target)
    return answer


def joint_trace(field,field_degree,width,t,variables,numerators,denominators,nonzero=False):
    """Direct cyclic carries: one unique cycle for each retained monomial."""
    if width<1 or field_degree%width or width%t:raise ValueError('Incompatible fields or row modulus')
    D=width//t;dimension=variables*D;b=1<<width;factors=list(numerators)
    for f in denominators:
        for h in [0]+list(range(2,field_degree+1)):
            factors.append(frobenius_factor(field,width,t,variables,f,h))
    punctured=[nonzero]*variables if isinstance(nonzero,bool)else list(nonzero)
    if len(punctured)!=variables:raise ValueError('Bad domain mask')
    zero=(0,)*dimension;local=[]
    for rho in range(t):
        values={zero:1}
        for row,terms in factors:
            if not 0<=row<t:raise ValueError('Unsynchronized row')
            if row!=rho:continue
            nxt={}
            for h,a in values.items():
                for c,k in terms:
                    if len(k)!=dimension or min(k)<0:raise ValueError('Bad monomial')
                    out=tuple(x+y for x,y in zip(h,k));nxt[out]=nxt.get(out,0)^field.mul(a,c)
            values={k:c for k,c in nxt.items()if c}
        local.append(values)
    bounds=[]
    for v in range(variables):
        for j in range(D):
            rotated=sum(max([0]+[sum(k[v*D+a]<<( (row+t*a-t*j)%width)for a in range(D))for c,k in terms if c])for row,terms in factors)
            bounds.append(rotated//(b-1))
    ranges=[range(1 if j==0 and not punctured[v]else 0,bounds[v*D+j]+1)for v in range(variables)for j in range(D)]
    answer=0;cache={}
    for initial in product(*ranges):
        final=tuple(initial[v*D+(j+1)%D]for v in range(variables)for j in range(D));states={initial:1}
        for rho in range(t):
            nxt={}
            for incoming,a in states.items():
                key=(rho,incoming)
                if key not in cache:
                    transitions=[]
                    for k,c in local[rho].items():
                        total=tuple(x+y for x,y in zip(incoming,k))
                        if all(x%2==0 for x in total):transitions.append((tuple(x//2 for x in total),c))
                    cache[key]=transitions
                for outgoing,c in cache[key]:nxt[outgoing]=nxt.get(outgoing,0)^field.mul(a,c)
            states={k:c for k,c in nxt.items()if c}
        answer^=states.get(final,0)
    return answer
