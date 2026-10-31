#!/usr/bin/env python3
"""A complete degree-five crossing and its common-cofactor exclusion problem."""
from __future__ import annotations

from itertools import combinations,product

from check_joint_determinant_barrier import determinant,polynomial_gcd
from check_multipoint_hermite_crossings import rank
from matching_curve import matching_curve,cut_matrices

KEYS = tuple(product(range(2),repeat=4))


def seed_data():
    return matching_curve(5,(0,2,3))


def matrices(seed):
    return cut_matrices(seed)


def affine_product(a,b):
    out = {}
    for i,x in enumerate(a):
        if not x:
            continue
        for j,y in enumerate(b):
            if y:
                monomial = tuple(sorted(k-1 for k in (i,j) if k))
                out[monomial] = (out.get(monomial,0)+x*y) % 5
    return {monomial:c for monomial,c in out.items() if c}


def cofactor_problem(seed):
    mats = matrices(seed)
    cofactors = [[[determinant([[m[a][b] for b in range(4) if b != j]
                               for a in range(4) if a != i],5)
                  for j in range(4)] for i in range(4)] for m in mats]
    rows = [[cofactors[0][i][j][k] for i,j in product(range(4),repeat=2)]
            +[-cofactors[1][i][j][k] % 5 for i,j in product(range(4),repeat=2)]
            for k in range(16)]
    rows.append([sum(cofactors[0][i][j]) % 5 for i,j in product(range(4),repeat=2)]+[0]*16)
    work = [row+[int(i == 16)] for i,row in enumerate(rows)]
    pivots = []
    r = 0
    for c in range(32):
        hit = next((i for i in range(r,len(work)) if work[i][c]),None)
        if hit is None:
            continue
        work[r],work[hit] = work[hit],work[r]
        inverse = pow(work[r][c],-1,5)
        work[r] = [x*inverse % 5 for x in work[r]]
        for i in range(len(work)):
            if i != r and work[i][c]:
                scalar = work[i][c]
                work[i] = [(a-scalar*b) % 5 for a,b in zip(work[i],work[r])]
        pivots.append(c)
        r += 1
        if r == len(work):
            break
    assert all(any(row[:32]) or not row[-1] for row in work)
    free = [j for j in range(32) if j not in pivots]
    affine = [[0]*(len(free)+1) for _ in range(32)]
    for i,j in enumerate(free):
        affine[j][1+i] = 1
    for row,j in zip(work,pivots):
        affine[j] = [row[-1]]+[-row[c] % 5 for c in free]
    assert rank(rows,5) == r
    for i,row in enumerate(rows):
        assert [sum(row[j]*affine[j][k] for j in range(32)) % 5
                for k in range(len(free)+1)] == [int(i == 16)]+[0]*len(free)
    equations = []
    for offset in (0,16):
        for i,j in combinations(range(4),2):
            for k,l in combinations(range(4),2):
                a = affine_product(affine[offset+4*i+k],affine[offset+4*j+l])
                b = affine_product(affine[offset+4*i+l],affine[offset+4*j+k])
                eq = {m:(a.get(m,0)-b.get(m,0)) % 5 for m in set(a)|set(b)}
                eq = {m:c for m,c in eq.items() if c}
                if eq:
                    equations.append([[c,list(m)] for m,c in sorted(eq.items())])
    return {'linear_rank':r,'free_indices':free,'affine_coordinates':affine,
            'quadratic_equations':equations}


def candidate_problem():
    seed,values = seed_data()
    for point in range(5):
        assert all(sum(x*pow(point,k,5) for k,x in enumerate(seed[key])) % 5 == values[point][i]
                   for i,key in enumerate(KEYS))
    assert [seed[key][-1] for key in KEYS] == list(values[-1])
    ds = [determinant(m,5) for m in matrices(seed)]
    return {'prime':5,'degree':5,'keys':KEYS,'seed':[seed[k] for k in KEYS],
            'values':values,'generator_seed':1337,'entry_span_dimension':rank(list(seed.values()),5),
            'determinants':ds,'determinant_gcd_degree':len(polynomial_gcd(*ds,5))-1,
            **cofactor_problem(seed)}


def singular_script(problem):
    def polynomial(equation):
        return '+'.join(str(c)+''.join('*x'+str(i+1) for i in monomial)
                        for c,monomial in equation) or '0'
    n = len(problem['free_indices'])
    count = len(problem['quadratic_equations'])
    lines = ['ring R=5,('+','.join('x'+str(i+1) for i in range(n))+'),dp;',
             'option(redSB);',
             'ideal I='+','.join(polynomial(eq) for eq in problem['quadratic_equations'])+';',
             'matrix T;','ideal G=liftstd(I,T);',
             'print("BASIS_SIZE="+string(size(G)));','print("DIMENSION="+string(dim(G)));',
             'if(size(G)==1 && G[1]==1){','matrix MI[1]['+str(count)+']=I;',
             'matrix CK=MI*T;','print("CAS_IDENTITY="+string(CK[1,1]));',
             'print("CERTIFICATE_BEGIN");','int ci; poly cp;',
             'for(ci=1;ci<=nrows(T);ci++){cp=T[ci,1];while(cp!=0){print(string(ci-1)+":"+string(leadcoef(cp))+":"+string(leadexp(cp)));cp=cp-lead(cp);}}',
             'print("CERTIFICATE_END");','}','quit;']
    return '\n'.join(lines)+'\n'
