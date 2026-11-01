"""Joint geometric-box sums via norm-one parametrization and carry contraction."""
from itertools import product
from joint_frobenius_carries import joint_trace


def network_sum(field,width,left_count,right_count,edges,left_charges=None,right_charges=None):
    """edges=(left,right,tag); field has degree2*width. Includes all poles."""
    if field.degree!=2*width:raise ValueError('Require the explicit quadratic extension')
    b=1<<width
    left_charges=[0]*left_count if left_charges is None else left_charges
    right_charges=[0]*right_count if right_charges is None else right_charges
    if len(left_charges)!=left_count or len(right_charges)!=right_count:raise ValueError('Bad boundary charges')
    omega=next(1<<i for i in range(field.degree)if field.power(1<<i,b)!=(1<<i))
    conjugate=field.power(omega,b);answer=0
    for infinity in product((False,True),repeat=right_count):
        live=[j for j in range(right_count)if not infinity[j]]
        variables=left_count+len(live);indices={j:left_count+k for k,j in enumerate(live)}
        if variables==0:raise ValueError('Require at least one left variable')
        numerators=[];denominators=[]
        def polynomial(terms):
            merged={}
            for c,positions in terms:
                h=[0]*variables
                for v in positions:h[v]+=1
                h=tuple(h);merged[h]=merged.get(h,0)^c
            return (0,[(c,h)for h,c in merged.items()if c])
        for i,charge in enumerate(left_charges):
            charge%=b-1
            for h in range(width):
                if charge>>h&1:
                    _,terms=polynomial([(1,(i,))]);numerators.append((h,terms))
        for j,charge in enumerate(right_charges):
            if infinity[j]:continue
            charge%=b+1;z=indices[j]
            for h in range(width+1):
                if charge>>h&1:
                    _,terms=polynomial([(1,(z,)),(field.power(omega,1<<h),())]);numerators.append((h%width,terms))
                    _,terms=polynomial([(1,(z,)),(field.power(conjugate,1<<h),())]);denominators.append((h%width,terms))
        for i,j,gamma in edges:
            if not (0<=i<left_count and 0<=j<right_count):raise ValueError('Bad edge')
            gb=field.power(gamma,b)
            if infinity[j]:
                numerators.append(polynomial([(1,()),(gb,(i,))]))
                denominators.append(polynomial([(1,()),(gamma,(i,))]))
            else:
                z=indices[j]
                numerators.extend([polynomial([(1,(z,)),(omega,()),(gb,(i,z)),(field.mul(gb,conjugate),(i,))]),polynomial([(1,(z,)),(conjugate,())])])
                denominators.extend([polynomial([(1,(z,)),(omega,())]),polynomial([(1,(z,)),(conjugate,()),(gamma,(i,z)),(field.mul(gamma,omega),(i,))])])
        answer^=joint_trace(field,field.degree,width,width,variables,numerators,denominators,[True]*left_count+[False]*len(live))
    return answer
