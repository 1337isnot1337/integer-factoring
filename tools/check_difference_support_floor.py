"""Check active-difference rank floors using multivariate coefficients over F2."""
import json, itertools

def rank_mod(rows,p):
 rows=[row[:] for row in rows];r=0
 for c in range(len(rows[0]) if rows else 0):
  j=next((i for i in range(r,len(rows)) if rows[i][c]%p),None)
  if j is None:continue
  rows[r],rows[j]=rows[j],rows[r];inv=pow(rows[r][c],-1,p);rows[r]=[(v*inv)%p for v in rows[r]]
  for i in range(r+1,len(rows)):
   a=rows[i][c]
   if a:rows[i]=[(x-a*y)%p for x,y in zip(rows[i],rows[r])]
  r+=1
 return r

def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k not in p:p[k]=x;break
   x^=p[k]
 return len(p)

def matrices(C,m,B=None):
 if B is None:B=C
 a,b=len(C),len(B);idx={x:i for i,x in enumerate(B)};jdx={x:i for i,x in enumerate(C)}
 R=[[0]*(a*b) for _ in range(1<<m)];S=[[0]*(a*b) for _ in range(1<<m)]
 for xi,x in enumerate(C):
  for yi,y in enumerate(B):
   D=x^y;I=x&y;sub=D
   while True:
    u=I|sub;v=I|(D^sub)
    if u in idx and v in jdx:
     ui,vi=idx[u],jdx[v];swap=x^u;free=D^swap;t=free
     while True:
      mon=swap|t;R[mon][xi*b+yi]|=1<<(ui*a+vi);S[mon][xi*b+ui]|=1<<(yi*a+vi)
      if not t:break
      t=(t-1)&free
    if not sub:break
    sub=(sub-1)&D
 return R,S

def mul(a,b,k,mod):
 r=0
 while b:
  if b&1:r^=a
  b>>=1;a<<=1
  if a>>k:a^=mod
 return r

def main():
 cases=queries=0
 for m in (1,2,3):
  for cmask in range(1,1<<(1<<m)):
   C=[x for x in range(1<<m) if cmask>>x&1];q=len(C);diff={x^y for x in C for y in C}-{0}
   weights={}
   for D in diff|{0}:
    fiber=[x for x in C if x^D in C]
    h=len({x&~D for x in fiber});g=len({x&D for x in fiber})
    assert h*g>=len(fiber)
    weights[D]=(h,g,len(fiber))
   assert sum(v[2] for v in weights.values())==q*q
   R,S=matrices(C,m);rr=[0]*(q*q);ss=rr[:];old=0
   for t in range(1<<(1<<m)):
    mask=t^(t>>1)
    if t:
     j=(mask^old).bit_length()-1;rr=[a^b for a,b in zip(rr,R[j])];ss=[a^b for a,b in zip(ss,S[j])]
    active=sum((mask>>d)&1 for d in diff)
    assert rank(rr)>=active and rank(ss)>=2*active
    hr=sum(h for D,(h,g,n) in weights.items() if mask>>D&1)
    gc=sum(g for D,(h,g,n) in weights.items() if mask>>D&1)
    assert rank(rr)>=hr and rank(ss)>=gc
    old=mask;queries+=1
   cases+=1
 paired_cases=paired_queries=0
 for amask in range(1,16):
  for bmask in range(1,16):
   A=[x for x in range(4) if amask>>x&1];B=[x for x in range(4) if bmask>>x&1]
   R,S=matrices(A,2,B);a,b=len(A),len(B);diff={x^y for x in A for y in B};weights={}
   for D in diff:
    fiber=[x for x in A if x^D in B];h=len({x&~D for x in fiber});g=len({x&D for x in fiber})
    assert h*g>=len(fiber) and len(fiber)<=min(a,b)
    weights[D]=(h,g,len(fiber))
   assert sum(v[2] for v in weights.values())==a*b
   for mask in range(16):
    rr=[0]*(a*b);ss=rr[:]
    for j in range(4):
     if mask>>j&1:rr=[x^y for x,y in zip(rr,R[j])];ss=[x^y for x,y in zip(ss,S[j])]
    assert rank(rr)>=sum(h for D,(h,g,n) in weights.items() if mask>>D&1)
    assert rank(ss)>=sum(g for D,(h,g,n) in weights.items() if mask>>D&1)
    paired_queries+=1
   paired_cases+=1
 ternary_cases=ternary_queries=0
 for cmask in range(1,16):
  C=[x for x in range(4) if cmask>>x&1];q=len(C);Ds={x^y for x in C for y in C};weights={}
  for D in Ds:
   fiber=[x for x in C if x^D in C];weights[D]=(len({x&~D for x in fiber}),len({x&D for x in fiber}))
  entries=[]
  for xi,x in enumerate(C):
   for yi,y in enumerate(C):
    for ui,u in enumerate(C):
     for vi,v in enumerate(C):
      if (x^y)!=(u^v) or (x&y)!=(u&v):continue
      swap=x^u;free=(x^y)^swap;monomials=[];t=free
      while True:
       monomials.append((swap|t,(-1)**t.bit_count()))
       if not t:break
       t=(t-1)&free
      entries.append((xi*q+yi,ui*q+vi,xi*q+ui,yi*q+vi,monomials))
  for coeff in itertools.product(range(3),repeat=4):
   R=[[0]*(q*q) for _ in range(q*q)];S=[[0]*(q*q) for _ in range(q*q)]
   for i,j,k,l,mon in entries:R[i][j]=S[k][l]=sum(coeff[d]*sgn for d,sgn in mon)%3
   assert rank_mod(R,3)>=sum(h for D,(h,g) in weights.items() if coeff[D])
   assert rank_mod(S,3)>=sum(g for D,(h,g) in weights.items() if coeff[D])
   ternary_queries+=1
  ternary_cases+=1
 fib=[]
 for m in range(1,15):
  C=[x for x in range(1<<m) if not(x&(x>>1))];D={x^y for x in C for y in C}
  assert len(D)==1<<m
  fib.append([m,len(C),len(D)])
 graphs=[]
 for k,mod in [(1,3),(2,7),(3,11),(4,19)]:
  C=[a|(b<<k)|(mul(a,b,k,mod)<<(2*k)) for a in range(1<<k) for b in range(1<<k)]
  D={x^y for x in C for y in C}-{0};q=len(C)
  assert len(D)==((1<<(2*k))-1)*(1<<k)
  assert 2*len(D)**2>q**3
  graphs.append([k,q,len(D)])
 # Rational proof that phi^10=55phi+34 < 128 via phi < 17/10.
 assert 5*100<24**2 and 55*17+340<1280
 return dict(scope='Finite checks of the ordinary active-difference theorem; no factoring improvement',
             exhaustive_codes=cases,multivariate_functionals=queries,paired_codes=paired_cases,
             paired_multivariate_functionals=paired_queries,ternary_codes=ternary_cases,
             ternary_functionals=ternary_queries,hard_core_sumsets=fib,multiplication_graph_sumsets=graphs)

if __name__=='__main__':print(json.dumps(main(),indent=2))
