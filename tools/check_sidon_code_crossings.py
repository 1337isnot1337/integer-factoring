"""Exact F2 crossing ranks for quadratic graph codes; exploratory, no rate claim."""
import json, math, heapq, random

def mul(a,b,k,mod):
    r=0
    while b:
        if b&1:r^=a
        b>>=1;a<<=1
        if a>>k:a^=mod
    return r

def rank(rows):
    piv={}
    for x in rows:
        while x:
            i=x.bit_length()-1
            if i not in piv:piv[i]=x;break
            x^=piv[i]
    return len(piv)

def probe(code,m,label):
    q=len(code);idx={x:i for i,x in enumerate(code)}
    diffs=[x^y for i,x in enumerate(code) for y in code[i+1:]]
    sidon=len(diffs)==len(set(diffs))
    R=[[0]*(q*q) for _ in range(m+1)]
    C=[[0]*(q*q) for _ in range(m+1)]
    entries=0
    for xi,x in enumerate(code):
      for yi,y in enumerate(code):
        D=x^y;I=x&y;s=D.bit_count();sub=D
        while True:
          u=I|sub;v=I|(D^sub)
          if u in idx and v in idx:
            ui,vi=idx[u],idx[v];b=(x^u).bit_count();a=s-b;t=a
            while True:
              j=b+t
              R[j][xi*q+yi]^=1<<(ui*q+vi)
              C[j][xi*q+ui]^=1<<(yi*q+vi)
              if t==0:break
              t=(t-1)&a
            entries+=1
          if sub==0:break
          sub=(sub-1)&D
    # Endpoint must be literal crossing in both flattenings.
    def combine(arr,mask):
      out=[0]*(q*q)
      for j in range(m+1):
        if mask>>j&1:out=[a^b for a,b in zip(out,arr[j])]
      return out
    target=(1<<(m+1))-1
    assert combine(R,target)==[1<<(yi*q+xi) for xi in range(q) for yi in range(q)]
    assert combine(C,target)==[1<<(ui*q+xi) for xi in range(q) for ui in range(q)]
    # Independent dense constructor, using binomial coefficients rather than
    # subset iteration, verifies every polynomial entry in these small cases.
    for xi,x in enumerate(code):
      for yi,y in enumerate(code):
       for ui,u in enumerate(code):
        for vi,v in enumerate(code):
          good=(x&y)==(u&v) and (x^y)==(u^v)
          b=(x^u).bit_count();s=(x^y).bit_count()
          for j in range(m+1):
            expected=int(good and 0<=j-b<=s-b and math.comb(s-b,j-b)%2==1)
            assert ((R[j][xi*q+yi]>>(ui*q+vi))&1)==expected
            assert ((C[j][xi*q+ui]>>(yi*q+vi))&1)==expected
    pairs=[]
    for mask in range(target+1):pairs.append([rank(combine(R,mask)),rank(combine(C,mask))])
    for mask,(r,c) in enumerate(pairs):
        assert rank(combine(R,mask)[::-1])==r
        assert rank(combine(C,mask)[::-1])==c
        if sidon:
            active=sum((mask>>((x^y).bit_count()))&1 for i,x in enumerate(code) for y in code[i+1:])
            assert r>=active and c>=2*active
    scale=1<<24
    weights=[math.isqrt(a*b*scale*scale) for a,b in pairs]
    dist=[10**100]*(target+1);dist[0]=0;prev={};heap=[(0,0)]
    while heap:
      d,u=heapq.heappop(heap)
      if d!=dist[u]:continue
      for v,w in enumerate(weights):
        z=u^v;nd=d+w
        if nd<dist[z]:dist[z]=nd;prev[z]=(u,v);heapq.heappush(heap,(nd,z))
    path=[];v=target
    while v:u,w=prev[v];path.append(w);v=u
    # A shortest-path dual checks every edge without trusting path optimality.
    assert all(dist[u^v]<=dist[u]+weights[v] for u in range(target+1) for v in range(target+1))
    return dict(label=label,bits=m,q=q,entries=entries,ranks=pairs,rank_sum_lower_numerator=dist[target],denominator=scale,
                lower_approx=dist[target]/scale,upper_approx=sum(math.sqrt(pairs[v][0]*pairs[v][1]) for v in path),
                q_power_sqrt2=q**math.sqrt(2),exceeds_q_power_1_5=dist[target]**2>q**3*scale**2,
                witness_masks=path,sidon=sidon,scope='Complete coefficient-span F2 finite control; not a factoring improvement')

if __name__=='__main__':
  results=[]
  for k,mod in [(1,3),(2,7)]:
    code=[a|(b<<k)|(mul(a,b,k,mod)<<(2*k)) for a in range(1<<k) for b in range(1<<k)]
    results.append(probe(code,3*k,'multiplication-graph-k'+str(k)))
  for k,mod in [(2,7),(3,11)]:
    code=[a|(mul(mul(a,a,k,mod),a,k,mod)<<k) for a in range(1<<k)]
    results.append(probe(code,2*k,'cubic-graph-k'+str(k)))
  rng=random.Random(20261031)
  for m in range(3,6):
    for t in range(4):
      candidates=list(range(1<<m));rng.shuffle(candidates);code=[];seen=set()
      for x in candidates:
        ds={x^y for y in code}
        if not ds&seen:code.append(x);seen.update(ds)
      results.append(probe(code,m,f'sidon-control-{m}-{t}'))
  print(json.dumps(results,indent=2))
