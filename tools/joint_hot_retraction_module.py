#!/usr/bin/env python3
"""A complete polynomial model of rational joint left retractions.

Initial saturation uses 0 and 1. The checked polynomial right inverse then
certifies the complete rational row space, without restricting pole locations.
Polynomial output and the strict endpoint give a finite bounded-jet model.
"""
from functools import lru_cache
from restrict_implicit_reflection import restricted_tensor
from coherent_hot_compounds import height


@lru_cache(None)
def ambient():
    words=tuple(sorted((w for w in range(16) if height(w)<=8),key=lambda w:(height(w),w)))
    assert len(words)==9
    labels=[sum(1<<(7-i) for i in range(4) if w>>i&1) for w in words]
    seed=restricted_tensor(16,labels);Q=len(words)
    shift=min((v&-v).bit_length()-1 for v,_ in seed['values'] if v)
    values=[(v>>shift,0) for v,_ in seed['values']]
    matrix=[[values[x+Q*y+Q*Q*u+Q**3*v] for u in range(Q) for v in range(Q)] for y in range(Q) for x in range(Q)]
    for i,row in enumerate(matrix):
        for j,(v,_) in enumerate(row):assert v.bit_count()%2==int(i==j)
    return words,matrix,shift


def deflate(p):
    assert not any(p>>i&1 for i in range(1,p.bit_length(),2))
    return sum(((p>>(2*i))&1)<<i for i in range((p.bit_length()+1)//2))


def polynomial_multiply(a,b):
    result=0
    while b:
        bit=b&-b;b^=bit;result^=a<<(bit.bit_length()-1)
    return result


def source_rows():
    words,A,_=ambient();selected=[x+9*y for y in range(8) for x in range(8)];excluded=[i for i in range(81) if i not in selected]
    return [[deflate(A[e][j][0]) for j in selected] for e in excluded]


def dependency(vectors):
    pivots={}
    for i,v in enumerate(vectors):
        c=1<<i
        while v:
            k=v.bit_length()-1
            if k in pivots:v^=pivots[k][0];c^=pivots[k][1]
            else:pivots[k]=(v,c);break
        if not v:return c
    return 0


def degree(row):return max(v.bit_length()-1 for v in row)


def divide_linear(v,root):
    if not root:
        assert not v&1;return v>>1
    q=0
    while v.bit_length()>1:
        s=v.bit_length()-2;q^=1<<s;v^=3<<s
    assert not v;return q


@lru_cache(None)
def saturated_basis():
    words,A,_=ambient();selected=[x+9*y for y in range(8) for x in range(8)];excluded=[i for i in range(81) if i not in selected]
    rows=[[deflate(A[e][j][0]) for j in selected] for e in excluded];initial_potential=sum(degree(r) for r in rows);steps=[]
    for root in (0,1):
        while True:
            evaluations=[sum(((v&1) if root==0 else v.bit_count()%2)<<j for j,v in enumerate(row)) for row in rows]
            relation=dependency(evaluations)
            if not relation:break
            selected_rows=[i for i in range(len(rows)) if relation>>i&1];pivot=max(selected_rows,key=lambda i:degree(rows[i]));old=degree(rows[pivot]);new=[0]*64
            for i in selected_rows:new=[x^y for x,y in zip(new,rows[i])]
            if not any(new):rows.pop(pivot);steps.append({'root':root,'removed_dependency':selected_rows});continue
            new=[divide_linear(v,root) for v in new];assert degree(new)<old;rows[pivot]=new;steps.append({'root':root,'rows':selected_rows,'pivot':pivot})
    reductions=0
    while True:
        ds=[degree(r) for r in rows];leading=[sum(((v>>d)&1)<<j for j,v in enumerate(row)) for row,d in zip(rows,ds)];rel=dependency(leading)
        if not rel:break
        selected_rows=[i for i in range(len(rows)) if rel>>i&1];pivot=max(selected_rows,key=lambda i:ds[i]);d=ds[pivot];new=[0]*64
        for i in selected_rows:new=[x^(v<<(d-ds[i])) for x,v in zip(new,rows[i])]
        if not any(new):rows.pop(pivot)
        else:assert degree(new)<d;rows[pivot]=new
        reductions+=1
    for root in (0,1):assert not dependency([sum(((v&1) if root==0 else v.bit_count()%2)<<j for j,v in enumerate(row)) for row in rows])
    # Coefficients must vanish at 1 to make the correction vanish there.
    strict=[[v^(v<<1) for v in row] for row in rows]
    return strict,{'rank':len(rows),'initial_degree_potential':initial_potential,'saturation_steps':steps,'leading_reductions':reductions,'strict_basis_degrees':[degree(r) for r in strict]}


def row_module(row_height,protected=4,shift_radius=None):
    source,certificate=saturated_basis();words,_,_=ambient();h=[height(w) for w in words[:8]];colh=[h[v]+h[u] for u in range(8) for v in range(8)]
    lower=16+protected-row_height;uppers=[x+16-protected for x in colh];width=max(uppers)+1
    allowed=[((1<<(u+1))-1)^((1<<lower)-1) if u>=lower else 0 for u in uppers]
    pivots={};generators=[];variables=[]
    for i,row in enumerate(source):
        for shift in range(max(uppers)-degree(row)+1):
            parts=[v<<shift for v in row];image=sum((v&~mask)<<(j*width) for j,(v,mask) in enumerate(zip(parts,allowed)));combo=1<<len(variables);variables.append((i,shift))
            while image:
                k=image.bit_length()-1
                if k in pivots:image^=pivots[k][0];combo^=pivots[k][1]
                else:pivots[k]=(image,combo);break
            if not image:generators.append(combo)
    output={}
    for c in generators:
        parts=[0]*64
        while c:
            bit=c&-c;c^=bit;i,s=variables[bit.bit_length()-1];parts=[a^(v<<s) for a,v in zip(parts,source[i])]
        assert all(not(v&~m) for v,m in zip(parts,allowed))
        packed=sum(v<<(j*width) for j,v in enumerate(parts))
        while packed:
            k=packed.bit_length()-1
            if k in output:packed^=output[k]
            else:output[k]=packed;break
    rows=[[(v>>(j*width))&((1<<width)-1) for j in range(64)] for v in output.values()]
    assert all(p.bit_count()%2==0 for row in rows for p in row)
    return {'row_height':row_height,'protected_jets_t':protected,'coefficient_variables':len(variables),'row_image_dimension':len(rows),'rows_t':rows,'saturated_at':[0,1]}


def polynomial_right_inverse(S,maximum_degree=32):
    rows=len(S);columns=len(S[0]);degree_bound=max(v.bit_length()-1 for row in S for v in row)
    for bound in range(maximum_degree+1):
        width=degree_bound+bound+1;pivots={}
        for j in range(columns):
            for shift in range(bound+1):
                image=sum((S[i][j]<<shift)<<(i*width) for i in range(rows));combo=1<<(j*(bound+1)+shift)
                while image:
                    k=image.bit_length()-1
                    if k in pivots:image^=pivots[k][0];combo^=pivots[k][1]
                    else:pivots[k]=(image,combo);break
        solutions=[]
        for i in range(rows):
            image=1<<(i*width);combo=0
            while image:
                k=image.bit_length()-1
                if k not in pivots:break
                image^=pivots[k][0];combo^=pivots[k][1]
            if image:break
            solutions.append(combo)
        if len(solutions)==rows:
            R=[[sum(((solutions[c]>>(j*(bound+1)+k))&1)<<k for k in range(bound+1)) for c in range(rows)] for j in range(columns)]
            return R,bound
    raise ArithmeticError('No polynomial right inverse found within the declared bound.')
