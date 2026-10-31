#!/usr/bin/env python3
"""Exact rectangular residual repair and independent-type endpoint checks.

The local numerator is rational with a recorded common denominator. These
checks establish admissibility and do not claim a pressure improvement.
"""
import argparse
import json
from pathlib import Path
from functools import lru_cache
from itertools import product,combinations_with_replacement,permutations
from check_coordinate_core_crossings import one_particle,mul,power,exact_divide,at_one,c_order,substitute_ramified,digest
from check_residual_quartic_codes import residual_operator

ROOT=Path(__file__).resolve().parent.parent
F=[1153,66,4,8,528,32,64,1152,256,512,1025,2048]
Finv=[129,66,4,8,528,32,64,1025,256,512,1153,2048]
hX=[0,0,1,1];hU=[0,0,1]

def apply(A,v):return sum(((row&v).bit_count()&1)<<i for i,row in enumerate(A))
assert all(apply(F,apply(Finv,1<<i))==1<<i for i in range(12))
assert all(apply(Finv,apply(F,1<<i))==1<<i for i in range(12))

def matrix_multiply(left,right):
 result=[]
 for row in left:
  value=0
  for j in range(12):
   if row>>j&1:value^=right[j]
  result.append(value)
 return result

def matrix_add(*matrices):
 result=[0]*12
 for matrix in matrices:
  result=[a^b for a,b in zip(result,matrix)]
 return result

def matrix_units(entries):
 result=[0]*12
 for i,j in entries:result[i]^=1<<j
 return result

IDENTITY=[1<<i for i in range(12)]
F_ALPHA=matrix_units([(0,9),(3,10),(11,8)])
F_BETA=matrix_units([(8,11),(9,10),(10,3)])
F_MIXED=matrix_units([(8,8),(9,9),(10,10)])
A=matrix_multiply(Finv,F_ALPHA)
B=matrix_multiply(Finv,F_BETA)
A2=matrix_multiply(A,A)
B2=matrix_multiply(B,B)
PARTIAL=matrix_add(IDENTITY,matrix_units([(10,0),(10,7)]))

def mixed_corner(alpha=0,beta=0):
 assert alpha in (0,1) and beta in (0,1)
 f=matrix_add(F,*([F_ALPHA] if alpha else []),*([F_BETA] if beta else []),
              *([F_MIXED] if alpha and beta else []))
 ia=matrix_add(IDENTITY,A,A2) if alpha else IDENTITY
 ib=matrix_add(IDENTITY,B,B2) if beta else IDENTITY
 inverse=matrix_multiply(matrix_multiply(ia,ib),Finv)
 assert matrix_multiply(f,inverse)==IDENTITY
 assert matrix_multiply(inverse,f)==IDENTITY
 return f,inverse

def check_parameter_factorization():
 assert not any(matrix_multiply(A2,A))
 assert not any(matrix_multiply(B2,B))
 assert matrix_multiply(Finv,F_MIXED)==matrix_multiply(B,A)
 assert matrix_multiply(F,A)==F_ALPHA
 assert matrix_multiply(F,B)==F_BETA
 assert matrix_multiply(PARTIAL,PARTIAL)==IDENTITY
 triangular=matrix_multiply(PARTIAL,F)
 heights=[hX[i%4]+hU[i//4] for i in range(12)]
 for row,mask in enumerate(triangular):
  assert mask>>row&1
  for column in range(12):
   if row!=column and mask>>column&1:assert heights[row]<heights[column]
 return {'factorization':'F(alpha,beta)=F0*(I+beta*B)*(I+alpha*A)',
         'A_rows':A,'B_rows':B,'A_cubed_zero':True,'B_cubed_zero':True,
         'determinant_polynomial':1,'partial_repair_rows':PARTIAL,
         'partial_repaired_residual_rows':triangular}

def check_literal_corner(alpha=0,beta=0):
 residual=residual_operator()
 xmap=(1,0,2,3)
 umap=((1,),(0,),(2,3) if alpha else (2,))
 ymap=((1,),(0,),(2,3) if beta else (2,))
 rows=[0]*12
 for x,y,u,v in product(range(4),range(3),range(3),range(4)):
  value=0
  for yy,uu in product(ymap[y],umap[u]):
   value^=residual[xmap[v]+4*yy]>>(xmap[x]+4*uu)&1
  rows[v+4*y]|=value<<(x+4*u)
 assert rows==mixed_corner(alpha,beta)[0]

def build(indices,M,check_types=True,alpha=0,beta=0,repair='full'):
 n=len(indices);b=1<<(n-2);qx=4*b;qu=3*b;top=1<<(n-1);T=indices[-1]
 current_F,current_inverse=mixed_corner(alpha,beta)
 if repair=='partial':
  assert alpha==beta==0
  current_repair=PARTIAL
 else:
  assert repair=='full'
  current_repair=current_inverse
 endpoint=matrix_multiply(current_repair,current_F)
 mat=one_particle(indices,M)
 @lru_cache(None)
 def minor(rr,cc):
  if rr.bit_count()!=cc.bit_count():return 0
  if not rr:return 1
  rb=rr&-rr;i=rb.bit_length()-1;value=0;remaining=cc
  while remaining:
   cb=remaining&-remaining;remaining^=cb;j=cb.bit_length()-1
   if mat[i][j]:value^=mul(mat[i][j],minor(rr^rb,cc^cb))
  return value
 mx=[((a%b)<<1)|(1,0,top,1|top)[a//b] for a in range(qx)]
 mu=[[((a%b)<<1)|s for s in ((1,) if a//b==0 else (0,) if a//b==1
       else (top,1|top) if alpha else (top,))] for a in range(qu)]
 my=[[((a%b)<<1)|s for s in ((1,) if a//b==0 else (0,) if a//b==1
       else (top,1|top) if beta else (top,))] for a in range(qu)]
 def pos(x,y,u,v):return x+qx*(y+qu*(u+qu*v))
 N=[0]*(qx*qx*qu*qu)
 for x,y,u,v in product(range(qx),range(qu),range(qu),range(qx)):
  p=0
  for yy,uu in product(my[y],mu[u]):
   xx,vv=mx[x],mx[v]
   rr=vv|(yy<<n);cc=uu|(xx<<n);p^=mul(5,minor(rr,cc))
   if xx&1 and uu&1 and yy&top and vv&top:p^=minor((vv^top)|((yy^top)<<n),(uu^1)|((xx^1)<<n))
  N[pos(x,y,u,v)]=p
 # c^2 H_out^-1 F^-1 H_out, with exact core heights in T units.
 gate=[]
 for row,mask in enumerate(current_repair):
  vt,yt=row%4,row//4;ht=hX[vt]+hU[yt]
  edges=[]
  for col in range(12):
   if mask>>col&1:
    vs,ys=col%4,col//4;hs=hX[vs]+hU[ys]
    assert 2+hs-ht>=0
    edges.append((vs,ys,power(3,2+hs-ht)))
  gate.append(edges)
 values=[0]*len(N)
 for x,y,u,v in product(range(qx),range(qu),range(qu),range(qx)):
  ans=0
  for vs,ys,p in gate[v//b+4*(y//b)]:ans^=mul(p,N[pos(x,ys*b+y%b,u,vs*b+v%b)])
  values[pos(x,y,u,v)]=ans
 cancel=min(c_order(p) for p in values if p)
 den=4-cancel
 if cancel:values=[exact_divide(p,power(3,cancel)) for p in values]
 hx=[sum(indices[j+1] for j in range(n-2) if a>>j&1)+T*hX[a//b] for a in range(qx)]
 hu=[sum(indices[j+1] for j in range(n-2) if a>>j&1)+T*hU[a//b] for a in range(qu)]
 for x,y,u,v in product(range(qx),range(qu),range(qu),range(qx)):
  p=substitute_ramified(values[pos(x,y,u,v)],T);shift=hu[y]+hx[v]-hx[x]-hu[u]-den*T
  if shift<0:
   assert not p&((1<<(-shift))-1);p>>=-shift
  else:p<<=shift
  expected=(endpoint[v//b+4*(y//b)]>>(x//b+4*(u//b)))&1
  expected*=int(x%b==v%b and y%b==u%b)
  assert p&1==expected
 checks=0
 if check_types:
  denominator=power(3,2*den)
  xtypes=[sorted(set(permutations(t))) for t in combinations_with_replacement(range(qx),2)]
  utypes=[sorted(set(permutations(t))) for t in combinations_with_replacement(range(qu),2)]
  for X,U in product(xtypes,utypes):
   for x,v,y,u in product(X,X,U,U):
    num=mul(values[pos(x[0],y[0],u[0],v[0])],values[pos(x[1],y[1],u[1],v[1])])
    p=exact_divide(num,denominator)
    assert at_one(p)==int(x==v and y==u)
    checks+=1
 extras=0
 if den==0:
  for x,y,u,v in product(range(qx),range(qu),range(qu),range(qx)):
   value=at_one(values[pos(x,y,u,v)])
   if value!=int(x==v and y==u):
    assert value and hu[y]+hx[v]-hx[x]-hu[u]>0
    extras+=1
 data={'dimensions':[qx,qu,qu,qx],'q_x':qx,'q_u':qu,'indices':list(indices),'ambient':M,'alpha':alpha,'beta':beta,'repair_kind':repair,'core_sizes':[4,3,3,4],'bulk_size':b,'degree':max(p.bit_length()-1 for p in values),'denominator_degree':den,'cancelled_c_order':cancel,'strict_endpoint':False,'typed_endpoint':'independent common X/V and U/Y full-letter types','logical_heights_x_v':hx,'logical_heights_u_y':hu,'residual_rows':current_F,'inverse_rows':current_inverse,'repair_rows':current_repair,'ramified_residual_rows':endpoint,'ramified_endpoint_checks':len(values),'two_copy_type_checks':checks,'polynomial_endpoint_positive_extras':extras if den==0 else None,'coefficient_sha256':digest(values),'values':values}
 return data

def check_no_bulk_wick(data):
 """A necessary Gaussian relation involving only retained physical states."""
 assert data['indices']==[0,1] and data['repair_kind']=='partial'
 assert data['alpha']==data['beta']==data['denominator_degree']==0
 xmap=(1,0,2,3);umap=(1,0,2)
 coefficients={}
 for x,y,u,v in product(range(4),range(3),range(3),range(4)):
  mask=xmap[x]|(umap[y]<<2)|(umap[u]<<4)|(xmap[v]<<6)
  coefficients[mask]=data['values'][x+4*(y+3*(u+3*v))]
 assert coefficients[0]==1
 failures=[]
 checks=0
 witness=0
 for mask,value in coefficients.items():
  if mask.bit_count()!=4:continue
  a,b,c,d=[1<<i for i in range(8) if mask>>i&1]
  pairs=[(a|b,c|d),(a|c,b|d),(a|d,b|c)]
  assert all(first in coefficients and second in coefficients for first,second in pairs)
  residual=mul(coefficients[0],value)
  for first,second in pairs:residual^=mul(coefficients[first],coefficients[second])
  checks+=1
  if residual:failures.append(mask)
  if mask==153:witness=residual
 assert checks==41 and len(failures)==5
 assert witness==17==power(3,4)
 retained=[0,153,9,144,17,136,129,24]
 return {'variable_order':['x0','xT','y0','yT','u0','uT','v0','vT'],
         'quartic_mask':153,'vacuum':1,'wick_residual_polynomial_bits':witness,
         'wick_residual_factor':'(1+w)^4',
         'relation_coefficient_bits':{str(mask):coefficients[mask] for mask in retained},
         'accessible_quartic_relations_checked':checks,
         'accessible_quartic_failure_masks':failures}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path)
 parser.add_argument('--seed-directory',type=Path)
 args=parser.parse_args()
 factorization=check_parameter_factorization()
 for alpha,beta in product((0,1),repeat=2):check_literal_corner(alpha,beta)
 records=[]
 wick=None
 for S,M,tag in [((0,1),2,'01'),((0,1,2),4,'012'),((0,1,2,3),4,'0123'),((0,2,5),8,'025')]:
  d=build(S,M)
  records.append({k:v for k,v in d.items() if k!='values'})
  if args.seed_directory:
   args.seed_directory.mkdir(parents=True,exist_ok=True)
   (args.seed_directory/f'rectangular-corrected-{tag}.json').write_text(json.dumps(d)+'\n')
 for S,M,tag in [((0,1),2,'01'),((0,1,2),4,'012'),((0,1,2,3),4,'0123')]:
  d=build(S,M,alpha=1,beta=1)
  assert d['denominator_degree']==0
  records.append({k:v for k,v in d.items() if k!='values'})
  if args.seed_directory:
   (args.seed_directory/f'rectangular-mixed-11-{tag}.json').write_text(json.dumps(d)+'\n')
 for S,M,tag in [((0,1),2,'01'),((0,1,2),4,'012'),((0,1,2,3),4,'0123')]:
  d=build(S,M,repair='partial')
  if tag=='01':wick=check_no_bulk_wick(d)
  records.append({k:v for k,v in d.items() if k!='values'})
  if args.seed_directory:
   (args.seed_directory/f'rectangular-partial-{tag}.json').write_text(json.dumps(d)+'\n')
 result={'status':'exact residual repair and paired-type polynomial crossings verified; no exponent improvement claimed','parameter_factorization':factorization,'non_gaussian_partial_control':wick,'cases':records}
 if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))

if __name__=='__main__':main()
