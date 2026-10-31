"""Check active-difference rank floors using multivariate coefficients over F2."""
import json

def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k not in p:p[k]=x;break
   x^=p[k]
 return len(p)

def matrices(C,m):
 q=len(C);idx={x:i for i,x in enumerate(C)};R=[[0]*(q*q) for _ in range(1<<m)];S=[[0]*(q*q) for _ in range(1<<m)]
 for xi,x in enumerate(C):
  for yi,y in enumerate(C):
   D=x^y;I=x&y;a=D
   while True:
    u=I|a;v=I|(D^a)
    if u in idx and v in idx:
     ui,vi=idx[u],idx[v];swap=x^u;free=D^swap;b=free
     while True:
      mon=swap|b;R[mon][xi*q+yi]|=1<<(ui*q+vi);S[mon][xi*q+ui]|=1<<(yi*q+vi)
      if not b:break
      b=(b-1)&free
    if not a:break
    a=(a-1)&D
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
   R,S=matrices(C,m);rr=[0]*(q*q);ss=rr[:];old=0
   for t in range(1<<(1<<m)):
    mask=t^(t>>1)
    if t:
     j=(mask^old).bit_length()-1;rr=[a^b for a,b in zip(rr,R[j])];ss=[a^b for a,b in zip(ss,S[j])]
    active=sum((mask>>d)&1 for d in diff)
    assert rank(rr)>=active and rank(ss)>=2*active
    old=mask;queries+=1
   cases+=1
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
             exhaustive_codes=cases,multivariate_functionals=queries,hard_core_sumsets=fib,multiplication_graph_sumsets=graphs)

if __name__=='__main__':print(json.dumps(main(),indent=2))
