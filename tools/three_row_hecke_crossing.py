#!/usr/bin/env python3
"""Finite three-row level-four Hecke block crossings, with checked endpoints.

No uniform block-swap theorem is assumed: every returned finite tensor has
its full endpoint checked. Packed construction is optional NumPy research
code, not a practical factoring algorithm.
"""
from three_row_hecke_paths import paths,local,M,q,F
import hashlib,time


def example_paths(alphabet):
    base=paths((2,2,2),4)
    if alphabet==5:return base
    if alphabet==9:
        one=(0,1,2);result=sorted(set([one+p for p in base]+[p+one for p in base]));assert len(result)==9;return result
    raise ValueError('The retained examples have alphabet5 or9.')


def actions(global_paths):
    import numpy as np
    index={p:i for i,p in enumerate(global_paths)};out=[]
    for site in range(len(global_paths[0])-1):
        diagonal=np.zeros(len(global_paths),dtype=np.uint8);source=[];target=[]
        for i,p in enumerate(global_paths):
            action=local(p,site,index);diagonal[i]=q
            for j,c in action:
                if j==i:diagonal[i]^=c
            for j,c in action:
                if j!=i and c:
                    assert c==diagonal[i];source.append(i);target.append(j)
        assert len(set(target))==len(target)
        out.append((diagonal,np.array(source,dtype=np.int64),np.array(target,dtype=np.int64)))
    return out


def coefficient_tensor(physical_paths,progress=None,max_primary_buffer_bytes=2**31):
    import numpy as np
    start=time.monotonic();ps=[tuple(p) for p in physical_paths]
    if not ps or len(set(ps))!=len(ps):raise ValueError('Use distinct nonempty physical paths.')
    length=len(ps[0])
    if not length or length%3 or any(len(p)!=length for p in ps):raise ValueError('Use one positive length divisible by3.')
    m=length//3;allowed=set(paths((m,m,m),4))
    if any(p not in allowed for p in ps):raise ValueError('A physical path leaves the alcove or is not closed.')
    qs=paths((2*m,2*m,2*m),4);ix={p:i for i,p in enumerate(qs)};N=len(ps);S=len(qs)
    boundary=np.array([ix[p+r] for p in ps for r in ps]);word=[i for j in range(length) for i in range(length+j-1,j-1,-1)];D=len(word);words=D//64+1
    primary_bytes=3*S*N*N*3*words*8
    if primary_bytes>max_primary_buffer_bytes:raise ValueError('The three primary packed buffers exceed the requested memory guard.')
    ops=actions(qs)
    states=np.zeros((S,N*N,3,words),dtype=np.uint64);states[boundary,np.arange(N*N),0,0]=1
    following=np.empty_like(states);E=np.empty_like(states)
    bits=np.array([[[(M[c][1<<b]>>a)&1 for b in range(3)] for a in range(3)] for c in range(8)],dtype=np.uint64)
    for step,site in enumerate(word):
        diag,src,dst=ops[site]
        following[:,:,0,:]=states[:,:,2,:]
        following[:,:,1,:]=states[:,:,0,:]^states[:,:,2,:]
        following[:,:,2,:]=states[:,:,1,:]
        for a in range(3):
            E[:,:,a,:]=0
            for b in range(3):E[:,:,a,:]^=states[:,:,b,:]&(-bits[diag,a,b])[:,None,None]
        for w in range(words-1,0,-1):E[:,:,:,w]=(E[:,:,:,w]<<np.uint64(1))|(E[:,:,:,w-1]>>np.uint64(63))
        E[:,:,:,0]<<=np.uint64(1)
        following^=E;following[dst]^=E[src];states,following=following,states
        if progress and (step+1)%length==0:progress({'step':step+1,'total_steps':D,'seconds':time.monotonic()-start})
    selected=states[boundary];T=np.zeros((D+1,N*N,N*N),dtype=np.uint8)
    for k in range(D+1):
        for a in range(3):T[k]^=(((selected[:,:,a,k//64]>>np.uint64(k%64))&np.uint64(1)).astype(np.uint8)<<a)
    T=T.transpose(0,2,1).reshape(D+1,N,N,N,N);endpoint=np.bitwise_xor.reduce(T,axis=0);scalar=int(endpoint[0,0,0,0]);assert scalar
    expected=np.zeros((N,N,N,N),dtype=np.uint8)
    for i in range(N):
        for j in range(N):expected[i,j,j,i]=scalar
    if not np.array_equal(endpoint,expected):raise ArithmeticError('Projected block braid is not a scalar swap.')
    assert scalar==F.phase(6*m*m)
    T=np.array(M,dtype=np.uint8)[F.power(scalar,6),T]
    report={'physical_paths':ps,'physical_q':N,'global_dimension':S,'field_polynomial':11,'hecke_parameter':q,'braid_scalar':scalar,'endpoint_entries':N**4,'degree_bound':D,'actual_degree':max(k for k in range(D+1) if np.any(T[k])),'primary_buffer_bytes':primary_bytes,'coefficient_shape':list(T.shape),'coefficient_sha256':hashlib.sha256(T.tobytes()).hexdigest(),'seconds':time.monotonic()-start,'scope':'Finite exact construction with checked swap endpoint; no uniform endpoint theorem or pressure improvement claimed.'}
    return T,report
