#!/usr/bin/env python3
"""Literal certificates for a complete rational retraction row module."""
import json
from pathlib import Path
from hashlib import sha256
from joint_hot_retraction_module import saturated_basis,source_rows,divide_linear,polynomial_right_inverse,polynomial_multiply,dependency,degree,row_module,ambient,deflate
from coherent_hot_compounds_f4 import mul
from coherent_hot_compounds import height


def product(A,B):
    result=[]
    for row in A:
        output=[]
        for j in range(len(B[0])):
            v=0
            for x,r in zip(row,B):v^=polynomial_multiply(x,r[j])
            output.append(v)
        result.append(output)
    return result


def at_two(A):
    out=[]
    for row in A:
        values=[]
        for p in row:
            v=0
            for i in range(p.bit_length()-1,-1,-1):v=mul(v,2)^((p>>i)&1)
            values.append(v)
        out.append(values)
    return out


def rank_f4(A):
    A=[r[:] for r in A];r=0
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A)) if A[i][j]),None)
        if i is None:continue
        A[r],A[i]=A[i],A[r];inv=mul(A[r][j],A[r][j]);A[r]=[mul(inv,x) for x in A[r]]
        for i in range(r+1,len(A)):
            c=A[i][j]
            if c:A[i]=[x^mul(c,y) for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r


def main():
    words,A,_=ambient();selected=[x+9*y for y in range(8) for x in range(8)]
    hs=[height(words[x])+height(words[y]) for y in range(8) for x in range(8)]
    zero_values=[];infinity_values=[]
    for i,oldi in enumerate(selected):
        for j,oldj in enumerate(selected):
            v=deflate(A[oldi][oldj][0])
            if v:
                zero_values.append((v&-v).bit_length()-1+hs[i]);infinity_values.append(hs[j]-(v.bit_length()-1))
    assert min(zero_values)==16 and min(infinity_values)==-16
    strict,steps=saturated_basis();S=[[divide_linear(v,1) for v in row] for row in strict];R,bound=polynomial_right_inverse(S);B=source_rows();I=product(S,R)
    assert I==[[int(i==j) for j in range(17)] for i in range(17)]
    U=product(B,R);assert product(U,S)==B;assert rank_f4(at_two(U))==17
    leading=[sum(((v>>degree(row))&1)<<j for j,v in enumerate(row)) for row in S];assert not dependency(leading)
    dims=[]
    for protected in (4,8,9,10):
        counts=[]
        for h in range(17):
            m=row_module(h,protected);counts.append(m['row_image_dimension'])
            for row in m['rows_t']:
                assert all(v.bit_count()%2==0 for v in row)
        dims.append({'protected_t_jets':protected,'row_dimensions_by_height_0_through_16':counts})
    names=[Path(__file__).name,'joint_hot_retraction_module.py','restrict_implicit_reflection.py','reflection_dual_rail_gaussian.py','check_schur_gaussian_reduction.py']
    report={'natural_regularization_minima':[min(zero_values),min(infinity_values)],'literal_nonzero_entry_valuation_controls':len(zero_values)+len(infinity_values),'status':'passed','basis':S,'polynomial_right_inverse':R,'source_to_basis_matrix':U,'right_inverse_degree_bound':bound,'literal_right_inverse_entries':289,'literal_source_containment_entries':1088,'source_rank_at_F4_two':17,'basis_rank':17,'saturation_and_reduction':steps,'bounded_rows':dims,'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},'scope':'Exact polynomial normal form for fixed-coordinate rational left retractions; not a rank-minimization or factoring bound.'}
    Path('evidence/joint-hot-retraction-module.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('basis','polynomial_right_inverse','source_to_basis_matrix','saturation_and_reduction','source_sha256','bounded_rows')}))

if __name__=='__main__':main()
