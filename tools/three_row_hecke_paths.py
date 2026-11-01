#!/usr/bin/env python3
"""Three-row seminormal Hecke paths over F8 at a seventh root."""
from functools import lru_cache
from check_fourier_multiplication_gates import BinaryField
F=BinaryField(0b1011);M=tuple(tuple(F.multiply(a,b) for b in range(8)) for a in range(8));q=2

def paths(shape,level=None):
    out=[]
    def walk(counts,path):
        if tuple(counts)==tuple(shape):out.append(tuple(path));return
        for r in range(len(shape)):
            if counts[r]>=shape[r] or (r and counts[r]>=counts[r-1]):continue
            counts[r]+=1
            if level is None or counts[0]-counts[-1]<=level:walk(counts,path+[r])
            counts[r]-=1
    walk([0]*len(shape),[]);return out


def local(path,site,index,prefix=None):
    counts=list(prefix or (0,0,0));contents=[]
    for r in path:contents.append(counts[r]-r);counts[r]+=1
    if path[site]==path[site+1]:return ((index[path],q),)
    d=contents[site+1]-contents[site];den=1^F.phase(-d)
    if not den:raise ArithmeticError('Singular seminormal denominator.')
    a=M[q^1][F.power(den,6)];b=q^a
    other=path[:site]+(path[site+1],path[site])+path[site+2:]
    if other not in index:
        assert a==1,(path,site,d,a)
        return ((index[path],a),)
    return ((index[path],a),(index[other],b))


def operators(shape,level=None):
    ps=paths(shape,level);ix={p:i for i,p in enumerate(ps)}
    return ps,[[local(p,i,ix) for p in ps] for i in range(sum(shape)-1)]


def apply(v,op):
    out={}
    for i,c in v.items():
        for j,a in op[i]:out[j]=out.get(j,0)^M[c][a]
    return {i:c for i,c in out.items() if c}


def check_relations(shape,level=None):
    ps,ops=operators(shape,level);counts={'quadratic':0,'braid':0,'commuting':0}
    for a in range(len(ps)):
        v={a:1}
        for i,T in enumerate(ops):
            x=apply(v,T);z=apply(x,T)
            for j,c in x.items():z[j]=z.get(j,0)^M[q^1][c]
            z[a]=z.get(a,0)^q;assert not any(z.values()),('quadratic',a,i);counts['quadratic']+=1
            if i+1<len(ops):
                U=ops[i+1];assert apply(apply(apply(v,T),U),T)==apply(apply(apply(v,U),T),U),('braid',a,i);counts['braid']+=1
            for j in range(i+2,len(ops)):
                assert apply(apply(v,T),ops[j])==apply(apply(v,ops[j]),T),('commuting',a,i,j);counts['commuting']+=1
    return {'shape':shape,'level':level,'dimension':len(ps),'checks':counts}


def local_motif_certificate():
    count=0;records=[];witness=None
    for a in range(5):
        for b in range(5-a):
            prefix=(a+b,b,0);ps=[]
            def walk(c,p):
                if len(p)==3:ps.append(tuple(p));return
                for r in range(3):
                    if r and c[r]>=c[r-1]:continue
                    c[r]+=1
                    if c[0]-c[2]<=4:walk(c,p+[r])
                    c[r]-=1
            walk(list(prefix),[]);ix={p:i for i,p in enumerate(ps)};ops=[[local(p,i,ix,prefix) for p in ps] for i in range(2)]
            for i in range(len(ps)):
                v={i:1};T,U=ops
                assert apply(apply(apply(v,T),U),T)==apply(apply(apply(v,U),T),U)
                for op in ops:
                    x=apply(v,op);z=apply(x,op)
                    for j,c in x.items():z[j]=z.get(j,0)^M[q^1][c]
                    z[i]=z.get(i,0)^q;assert not any(z.values())
                count+=1
            def e(v,op):
                z=apply(v,op)
                for i,c in v.items():z[i]=z.get(i,0)^M[q][c]
                return {i:c for i,c in z.items() if c}
            for i in range(len(ps)):
                v={i:1};E=e(v,ops[0]);three=e(e(E,ops[1]),ops[0]);expect={j:M[q][c] for j,c in E.items()}
                if three!=expect and witness is None:witness={'prefix':prefix,'path':ps[i],'triple_product':three,'q_times_e':expect}
            records.append({'prefix':prefix,'paths':ps})
    assert count==147 and witness is not None
    return {'alcove_weights':15,'three_step_states':count,'quadratic_checks':2*count,'braid_checks':count,'local_records':records,'non_TL_witness':witness}
