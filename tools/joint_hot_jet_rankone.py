#!/usr/bin/env python3
"""All-field first-residual rank-one constraints for rational joint retractions.

The polynomial row model is complete for the fixed embedding. Earlier RAW
jets at both natural places are fixed; this is not a profile-only condition.
"""
from functools import lru_cache
from joint_hot_retraction_module import ambient,saturated_basis,divide_linear,degree
from coherent_hot_compounds import height

def inflate(p):return sum(1<<(2*i) for i in range(p.bit_length()) if p>>i&1)

def build():
    words,A,_=ambient();selected=[x+9*y for y in range(8) for x in range(8)]
    return list(words[:8]),[[A[i][j] for j in selected] for i in selected]

def regular_matrices():
    letters,A=build();q=8;h=[2*height(w) for w in letters];hs=[h[x]+h[y] for y in range(q) for x in range(q)];A=[[v[0] for v in row] for row in A]
    M=[]
    for i,row in enumerate(A):
        assert all(not((v<<hs[i])&((1<<32)-1)) for v in row)
        M.append([(v<<hs[i])>>32 for v in row])
    def rev(v,d):
        assert v.bit_length()-1<=d
        return sum(1<<(d-j) for j in range(v.bit_length()) if v>>j&1)
    N=[[rev(A[x+q*y][v+q*u],h[v]+h[u]+32) for v in range(q) for y in range(q)] for u in range(q) for x in range(q)]
    return letters,hs,A,M,N


@lru_cache(None)
def source_basis():
    strict,_=saturated_basis();S=[[divide_linear(v,1) for v in row] for row in strict]
    return [[inflate(v)^(inflate(v)<<1) for v in row] for row in S]


def row_module(row_height,protected=8,shift_radius=None):
    source=source_basis();words,_,_=ambient();h=[2*height(w) for w in words[:8]];hs=[h[v]+h[u] for u in range(8) for v in range(8)]
    lower=32+protected-row_height;uppers=[x+32-protected for x in hs];width=max(uppers)+1
    allowed=[((1<<(u+1))-1)^((1<<lower)-1) if u>=lower else 0 for u in uppers]
    pivots={};variables=[];generators=[]
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
        row=[0]*64
        while c:
            bit=c&-c;c^=bit;i,s=variables[bit.bit_length()-1];row=[a^(v<<s) for a,v in zip(row,source[i])]
        assert all(not(v&~m) for v,m in zip(row,allowed));packed=sum(v<<(j*width) for j,v in enumerate(row))
        while packed:
            k=packed.bit_length()-1
            if k in output:packed^=output[k]
            else:output[k]=packed;break
    rows=[[(v>>(j*width))&((1<<width)-1) for j in range(64)] for v in output.values()]
    assert all(v.bit_count()%2==0 for row in rows for v in row)
    return {'row_height':row_height,'protected_jets_w':protected,'row_image_dimension':len(rows),'rows_t':rows}

def pmul(a,b,mask):
    r=0
    while b:
        bit=b&-b;b^=bit;r^=a<<(bit.bit_length()-1)
    return r&mask


def inverse(a,n):
    assert a&1;r=1
    for d in range(1,n):
        v=0
        for j in range(1,d+1):v^=((a>>j)&1)*((r>>(d-j))&1)
        r|=v<<d
    return r


def reverse(a,d):
    assert a.bit_length()-1<=d
    return sum(1<<(d-i) for i in range(a.bit_length()) if a>>i&1)


def leading_residual(matrix,k):
    n=len(matrix);mask=(1<<(k+1))-1;M=[[v&mask for v in row] for row in matrix];U=[1<<i for i in range(n)];V=U.copy();rank=0
    while rank<n:
        best=None
        for i in range(rank,n):
            for j in range(rank,n):
                v=M[i][j]
                if v:
                    order=(v&-v).bit_length()-1
                    if best is None or order<best[0]:best=(order,i,j)
        if best is None or best[0]>=k:break
        order,i,j=best;M[rank],M[i]=M[i],M[rank];U[rank],U[i]=U[i],U[rank]
        for row in M:row[rank],row[j]=row[j],row[rank]
        V[rank],V[j]=V[j],V[rank]
        inv=inverse(M[rank][rank]>>order,k+1)
        for i in range(rank+1,n):
            if not M[i][rank]:continue
            f=pmul(M[i][rank]>>order,inv,mask)
            if f&1:U[i]^=U[rank]
            for j in range(rank,n):M[i][j]^=pmul(f,M[rank][j],mask)
        for j in range(rank+1,n):
            if not M[rank][j]:continue
            f=pmul(M[rank][j]>>order,inv,mask)
            if f&1:V[j]^=V[rank]
            for i in range(n):M[i][j]^=pmul(f,M[i][rank],mask)
        rank+=1
    residual=[sum(((M[i][j]>>k)&1)<<(j-rank) for j in range(rank,n)) for i in range(rank,n)]
    return U[rank:],V[rank:],residual,rank


def project(rows,U,V):
    out=[]
    for u in U:
        row=0
        while u:
            bit=u&-u;u^=bit;row^=rows[bit.bit_length()-1]
        out.append(sum(((row&v).bit_count()%2)<<j for j,v in enumerate(V)))
    return out


def pack(rows,width):return sum(row<<(i*width) for i,row in enumerate(rows))


def lift(k,return_data=True):
    letters,hs,A,M,N=regular_matrices();U0,V0,C0,r0=leading_residual(M,k);Ui,Vi,Ci,ri=leading_residual(N,k);n0,ni=len(C0),len(Ci)
    variables=[];pivots={};modules={h:row_module(h,protected=k) for h in sorted(set(hs))}
    for target,h in enumerate(hs):
        for correction in modules[h]['rows_t']:
            rowbits=sum(((v>>(32+k-h))&1)<<j for j,v in enumerate(correction))
            zrows=[0]*64;zrows[target]=rowbits;delta0=pack(project(zrows,U0,V0),n0)
            x,y=target%8,target//8;irows=[0]*64
            for j,v in enumerate(correction):
                u,vv=j//8,j%8
                if (v>>(hs[j]+32-k))&1:irows[x+8*u]|=1<<(y+8*vv)
            deltai=pack(project(irows,Ui,Vi),ni)
            image=delta0|(deltai<<(n0*n0));combo=1<<len(variables);variables.append((target,correction))
            while image:
                bit=image.bit_length()-1
                if bit in pivots:image^=pivots[bit][0];combo^=pivots[bit][1]
                else:pivots[bit]=(image,combo);break
    return {'variables':variables,'pivots':pivots,'C0':C0,'Ci':Ci,'n0':n0,'ni':ni,'letters':letters,'A':A}

def reduce(v,P):
 c=0
 for k in sorted(P,reverse=True):
  if v>>k&1:v^=P[k][0];c^=P[k][1]
 return v,c

def analyze(jet):
 d=lift(jet,True);n=d['n0'];ni=d['ni'];bits=n*n;mask=(1<<bits)-1
 v=d['C0'][0];assert v&1
 assert all(not(H&((1<<n)-1)) and not((H>>bits)&1) for H,_ in d['pivots'].values())
 assert all(r in (0,v) for r in d['C0'])
 constraints={};W=[]
 for H,pre in d['pivots'].values():
  Z=H&mask;forced=sum((v if Z>>(i*n)&1 else 0)<<(i*n) for i in range(n));bad=Z^forced;h=H;c=pre
  while bad:
   k=bad.bit_length()-1
   if k not in constraints:constraints[k]=(bad,h,c);break
   b,hh,cc=constraints[k];bad^=b;h^=hh;c^=cc
  if not bad:W.append((h,c))
 # Project to infinity while retaining one original parameter preimage.
 V={}
 for h,c in W:
  hi=h>>bits
  while hi:
   k=hi.bit_length()-1
   if k not in V:V[k]=(hi,c);break
   hi^=V[k][0];c^=V[k][1]
 # Split the allowed infinity image by its anchor column and row.
 anchors={};kernel={}
 def anchor(H):return sum(((H>>(i*ni))&1)<<i for i in range(ni)) | ((H&((1<<ni)-1))<<ni)
 for H,c in V.values():
  a=anchor(H)
  while a:
   k=a.bit_length()-1
   if k not in anchors:anchors[k]=(a,H,c);break
   aa,hh,cc=anchors[k];a^=aa;H^=hh;c^=cc
  if not a:
   H,c0=reduce(H,kernel);c^=c0
   if H:kernel[H.bit_length()-1]=(H,c)
 A=list(anchors.values());rank=len(A);C=pack(d['Ci'],ni);u=sum(((C>>(i*ni))&1)<<i for i in range(ni));v=C&((1<<ni)-1)
 assert u&1 and v&1 and C==sum(v<<(i*ni) for i in range(ni) if u>>i&1)
 def outer(x,y):return sum(y<<(i*ni) for i in range(ni) if x>>i&1)
 du=[a&((1<<ni)-1) for a,H,c in A];dv=[a>>ni for a,H,c in A]
 assert all(not x&1 for x in du+dv)
 L=[];Q=[];squares=[]
 for i,(a,H,c) in enumerate(A):
  L.append(reduce(outer(du[i],v)^outer(u,dv[i])^H,kernel)[0]);squares.append(reduce(outer(du[i],dv[i]),kernel)[0])
 for i in range(rank):
  for j in range(i+1,rank):Q.append((i,j,reduce(outer(du[i],dv[j])^outer(du[j],dv[i]),kernel)[0]))
 # Compress the output space, preserving every quadratic zero exactly.
 span={}
 def encode(v):
  code=0
  while v:
   k=v.bit_length()-1
   if k not in span:span[k]=(v,1<<len(span))
   vv,cc=span[k];v^=vv;code^=cc
  return code
 LL=[encode(v) for v in L];QQ=[(i,j,encode(v)) for i,j,v in Q];SS=[encode(v) for v in squares]
 qspan={}
 for vv in SS+[v for _,_,v in QQ]:
  vv,_=reduce(vv,qspan)
  if vv:qspan[vv.bit_length()-1]=(vv,0)
 lspan={}
 for vv in LL:
  vv,_=reduce(vv,qspan);vv,_=reduce(vv,lspan)
  if vv:lspan[vv.bit_length()-1]=(vv,0)
 # New span entries never alter old coordinates: verify every decoding.
 def decode(c):
  v=0
  for raw,unit in span.values():
   if c&unit:v^=raw
  return v
 assert all(decode(c)==v for c,v in zip(LL,L))
 assert all(decode(c)==v for c,v in zip(SS,squares))
 assert all(decode(c)==v for (_,_,c),(_,_,v) in zip(QQ,Q))
 result={'jet':jet,'zero_rankone_image_dimension':len(W),'infinity_image_rank':len(V),'anchor_dimension':rank,'anchor_kernel_dimension':len(kernel),'constraint_output_dimension':len(span),'linear':LL,'quadratic':QQ,'squares':SS,'quadratic_output_rank':len(qspan),'linear_rank_mod_quadratics':len(lspan),'du':du,'dv':dv}
 return result
