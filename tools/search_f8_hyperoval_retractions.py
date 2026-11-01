import sys,json,itertools,time
from pathlib import Path
sys.path.insert(0,'tools')
from coherent_f8_maps import mul,height,validate_maps,tensor_vector
from coherent_hot_compounds_f8 import partial_block
from collections import Counter
columns=[(1,t,mul(t,t))for t in range(8)]+[(0,1,0),(0,0,1)]
lower=[w for w in range(4096)if height(w)<25];boundary=sorted(w for w in range(4096)if height(w)==25 and w.bit_count()==3)
def inverse(A):
 n=len(A);a=[list(row)+[int(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n)if a[i][j]);a[i],a[j]=a[j],a[i];inv=(0,1,5,6,7,2,3,4)[a[j][j]];a[j]=[mul(inv,v)for v in a[j]]
  for i in range(n):
   if i!=j:
    c=a[i][j];a[i]=[v^mul(c,w)for v,w in zip(a[i],a[j])]
 return [row[n:]for row in a]
trials=list(itertools.combinations(range(10),3));trials.remove((0,8,9));trials=[(0,8,9)]+trials;records=[];start=time.monotonic();path=Path('evidence/f8-hyperoval-charge3-retraction-screen.json')
for S in trials:
 inv=inverse([[columns[j][i]for j in S]for i in range(3)]);rows=[]
 for row in inv:
  coeff=[]
  for c in columns:
   v=0
   for a,b in zip(row,c):v^=mul(a,b)
   coeff.append(v)
  rows.append(coeff)
 p=[{w:1}for w in lower+[31]]+[{boundary[j]:rows[i][j]for j in range(10)if rows[i][j]}for i in range(3)];l=[{w:1}for w in lower+[31]]+[{boundary[j]:1}for j in S];_,ch=validate_maps(p,l,12)
 pairs=[(x,y)for x in range(96)for y in range(96)if ch[x]+ch[y]==3];cuts=[]
 for cut in ('zero','infinity'):
  rr=[tensor_vector(l[x],l[y]if cut=='zero'else p[y],12)for x,y in pairs];cc=[tensor_vector(p[y]if cut=='zero'else l[y],p[x],12)for x,y in pairs]
  r=partial_block(12,rr,cc,63);cuts.append({'cut':cut,'counts':dict(Counter(r['retained_orders_t'])),'tail_count':r['tail_count']})
 rec={'pivot_columns':S,'cuts':cuts,'orders63_total':sum(c['counts'].get(63,0)for c in cuts),'seconds':time.monotonic()-start};records.append(rec)
 path.write_text(json.dumps({'completed':len(records)==120,'registered_trials':trials,'completed_trials':len(records),'records':records,'scope':'All120 coordinate retractions of the fixed F8 hyperoval rowspace, but ONLY charge3 orders through63. Not full crossing profiles or pressure bounds.'},indent=2)+'\n');print(json.dumps({'trial':len(records),'S':S,'orders63_total':rec['orders63_total'],'seconds':rec['seconds']}),flush=True)
