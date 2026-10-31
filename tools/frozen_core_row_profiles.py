#!/usr/bin/env python3
"""Compact exact active profiles for the polynomial frozen-core row family.

Each cut has a 16-by-16 active matrix at M=4 and a 64-by-64 active
matrix at M>=8. Twelve source matrices are compiled together; arbitrary
F4 rows are then combined without rerunning the Fock program. The fixed
spectator factors are supplied by the existing full-core basis theorem.

The optional simultaneous binary polynomial unit divisions preserve
all row combinations and all real weighted Smith profiles. They are
performed only on complete polynomial arrays, never on truncated data.
"""
from __future__ import annotations
import ctypes,time
from check_rectangular_defect_core_basis import core_basis_data
from rectangular_defect_profiles import Program,Expr,library,MASK,precision_sizes,lcm2,quotient2,gcd2,retained_known_below
from check_rectangular_f4_repair import scale,add,ZERO

def source_expressions(program,initial,beta):
 data=program.data;M=data['M'];T=M-1;infinity=data['place']=='infinity'
 u0=('in_' if infinity else 'right_')+'u0';uT=('in_' if infinity else 'right_')+f'u{T}'
 y0=('out_' if infinity else 'left_')+'y0';yT=('out_' if infinity else 'left_')+f'y{T}'
 v0=('out_' if infinity else 'right_')+'v0';vT=('out_' if infinity else 'right_')+f'v{T}'
 core=(1,0,2,3);hx=(0,0,1,1);hu=(0,0,1)
 def middle(a):
  a=program.add(a,program.number(program.number(a,uT),u0))
  word=[vT+'_create',yT+'_create',u0+'_annihilate','in_x0_annihilate'] if infinity else ['left_x0_create',yT+'_create',u0+'_annihilate',vT+'_annihilate']
  quartic=program.word(a,word);a=program.add(a,program.scale(quartic,1,z=2 if infinity else 0,u=-2))
  if beta:
   b=program.car(program.number(a,yT),y0+'_annihilate');a=program.add(a,program.scale(b,beta,z=M if infinity else -1))
  return a
 def physical(v,y):
  if infinity:return (1-M)*bool(core[v]&1)+(M-1)*bool(core[v]&2)-M*(core[y].bit_count())
  return bool(core[y]&1)+(2*M-1)*bool(core[y]&2)
 base=middle(initial) if infinity else None;out=[]
 for source in range(12):
  vp,yp=source%4,source//4
  if infinity:
   a=program.unit(base,[v0,vT],core[0],core[vp]);a=program.unit(a,[y0,yT],core[0],core[yp])
  else:
   a=program.unit(initial,[v0,vT],core[vp],core[0]);a=middle(a);a=program.unit(a,[y0,yT],core[0],core[yp])
  dh=hx[vp]+hu[yp];dz=physical(0,0)-physical(vp,yp)-(dh if infinity else 0)
  out.append(program.scale(a,1,z=dz,u=dh))
 return out

def compile_components(data,beta=1,profile_precision=None):
 if not data.get('full_core',False):raise ValueError('Frozen rows require full physical core coordinates.')
 if beta not in (0,1,2,3):raise ValueError('The output parameter must lie in F4.')
 if profile_precision is None:profile_precision=min(8192,max(512,64*data['M']))
 start=time.monotonic();bases=[core_basis_data(data,side) for side in ('left','right')];extended=dict(data);extended['active_images']=dict(data['active_images']);names=[]
 for side,basis in zip(('left','right'),bases):
  prefix='basis_'+side+'_';names.append([prefix+name for name in basis['core_names']]);extended['active_images'].update({prefix+name:v for name,v in basis['regular_images'].items()})
 program=Program(extended);k=program.k;n=1<<k;left,right=bases;ln,rn=names
 diagonal=[sum(s for j,s in enumerate(data['active_orders']) if mask>>j&1) for mask in range(n)]
 def vacuum(a,cn):
  for name in cn:a=program.add(a,program.number(a,name))
  return a
 value=vacuum(Expr(0),rn)
 for bit,name in enumerate(rn):
  if 5>>bit&1:value=program.car(value,name+'_create')
 value=program.diagonal(value,diagonal);values=source_expressions(program,value,beta);outputs=[]
 for value in values:
  for bit,name in enumerate(ln):
   if 5>>bit&1:value=program.car(value,name+'_annihilate')
  outputs.append(vacuum(value,ln))
 base=min(v.base for v in outputs);denominator=1
 for value in outputs:denominator=lcm2(denominator,value.denominator)
 output_ids=[]
 for value in outputs:
  factor=quotient2(denominator,value.denominator)<<(value.base-base);output_ids.append(program.polynomial(value,factor).register)
 bound=max(program.degree_bounds[v] for v in output_ids);window,builder=precision_sizes(program,base,data['M'],profile_precision);words=builder//64;output_words=window//64
 instructions=(ctypes.c_int*(5*len(program.instructions)))(*(v for row in program.instructions for v in row));packed=(ctypes.c_uint64*(2*words*len(program.coefficients)))(*(v>>(64*w)&MASK for pair in program.coefficients for v in pair for w in range(words)));terms=[];offsets=[0]
 for car in program.cars:terms.extend(v for term in car for v in term);offsets.append(len(terms)//3)
 cc=(ctypes.c_int*len(terms))(*terms);co=(ctypes.c_int*len(offsets))(*offsets);dd=(ctypes.c_int*n)(*(v-min(diagonal) for v in diagonal));ii=(ctypes.c_int*len(left['vacuum_rows']))(*left['vacuum_rows']);jj=(ctypes.c_int*len(right['vacuum_columns']))(*right['vacuum_columns'])
 ls=n//left['multiplicity'];rs=n//right['multiplicity'];assert ls>=12
 # Batch the twelve components through the existing coordinate exporter.
 # Left blocks label components; only the first duplicated right block is read.
 ids=[output_ids[i%12] for i in range(ls) for j in range(rs)];rr=(ctypes.c_int*len(ids))(*ids);output=(ctypes.c_uint64*(2*output_words*n*n))();shifts=(ctypes.c_int64*(n*n))();meta=(ctypes.c_int64*4)()
 lib=library(builder);I=ctypes.POINTER(ctypes.c_int);U=ctypes.POINTER(ctypes.c_uint64);L=ctypes.POINTER(ctypes.c_int64);lib.fock_f4_coordinates_v2.argtypes=[ctypes.c_int,ctypes.c_int,I,ctypes.c_int,U,ctypes.c_int,I,I,I,ctypes.c_int,I,ctypes.c_int,I,I,ctypes.c_int,ctypes.c_int,U,L,L];lib.fock_f4_coordinates_v2.restype=ctypes.c_int
 status=lib.fock_f4_coordinates_v2(k,len(program.instructions),instructions,len(program.coefficients),packed,len(program.cars),co,cc,dd,left['multiplicity'],ii,right['multiplicity'],jj,rr,base,output_words,output,shifts,meta)
 if status:raise ArithmeticError((status,list(meta)))
 nn=left['multiplicity'];assert nn==right['multiplicity'];components=[];common_shifts=[];maximum_bits=0
 for source in range(12):
  matrix=[]
  for i in range(nn):
   row=[]
   for j in range(nn):
    index=(source*nn+i)*n+j;planes=tuple(sum(int(output[(2*index+p)*output_words+w])<<(64*w) for w in range(output_words)) for p in range(2));row.append((planes,int(shifts[index])))
   matrix.append(row)
  components.append(matrix)
 aligned=[[[ZERO]*nn for _ in range(nn)] for _ in range(12)]
 for i in range(nn):
  row=[]
  for j in range(nn):
   shift=min((components[s][i][j][1] for s in range(12) if components[s][i][j][0]!=ZERO),default=0);row.append(shift)
   for s in range(12):
    p,v=components[s][i][j]
    if p!=ZERO:
     q=tuple(a<<(v-shift) for a in p);aligned[s][i][j]=q;maximum_bits=max(maximum_bits,max(q).bit_length())
  common_shifts.append(row)
 md={'M':data['M'],'place':data['place'],'matrix_order':nn,'exact_polynomial_matrix':bound<builder and int(meta[3])<window,'known_below':retained_known_below(meta,window),'builder_known_below':int(meta[1]),'minimum_valuation':int(meta[0]),'program_degree_bound':bound,'profile_precision':max(window,1<<max(0,(maximum_bits-1).bit_length())),'output_window':window,'builder_precision':builder,'row_charges':left['coordinate_row_charges'][5*nn:6*nn],'column_charges':right['coordinate_column_charges'][5*nn:6*nn],'matrix_seconds':time.monotonic()-start,'unit_denominator_hex':hex(denominator),'instructions':len(program.instructions),'maximum_normalized_degree':int(meta[3]),'maximum_aligned_bits':maximum_bits}
 return aligned,common_shifts,md

def combine_components(comp,g,require_crossing=True):
 if len(g)!=12 or any(c not in (0,1,2,3) for c in g):raise ValueError('Supply twelve F4 row coefficients.')
 if require_crossing and (g[0]^g[10])!=1:raise ValueError('The frozen row must satisfy g[0]+g[10]=1.')
 n=len(comp[0]);out=[[ZERO]*n for _ in range(n)]
 for s,c in enumerate(g):
  if c:
   for i in range(n):
    for j in range(n):out[i][j]=add(out[i][j],scale(comp[s][i][j],c))
 return out

def reduce_components(comp,sh,md):
 if not md['exact_polynomial_matrix']:raise ValueError('Unit reduction requires complete polynomial coefficient arrays.')
 comp=[[row[:] for row in matrix] for matrix in comp]
 n=len(comp[0])
 def divide_group(positions):
  positions=list(positions);unit=0
  for s,i,j in positions:
   for q in comp[s][i][j]:
    if q:unit=gcd2(unit,q>>((q&-q).bit_length()-1))
   if unit==1:break
  if unit>1:
   for s,i,j in positions:comp[s][i][j]=tuple(quotient2(q,unit) for q in comp[s][i][j])
  return unit or 1
 common=divide_group((s,i,j) for s in range(12) for i in range(n) for j in range(n));rows=[divide_group((s,i,j) for s in range(12) for j in range(n)) for i in range(n)];cols=[divide_group((s,i,j) for s in range(12) for i in range(n)) for j in range(n)]
 bits=max(q.bit_length() for mat in comp for row in mat for p in row for q in p);out=dict(md);out['common_unit_removed']=hex(common);out['row_units_removed']=list(map(hex,rows));out['column_units_removed']=list(map(hex,cols));out['unit_reduced_maximum_bits']=bits;out['profile_precision']=max(512,32*md['M'],1<<max(0,(bits-1).bit_length()))
 out['source_unit_denominator_hex']=out.pop('unit_denominator_hex')
 out['normalized_by_regular_units']=True
 return comp,sh,out


def compile_interval(M,beta=1,profile_precision=None,normalize_units=True):
 if M<4 or M&(M-1):raise ValueError('Use a power-of-two interval M>=4.')
 from check_rectangular_defect_core_basis import data_for
 data=[data_for(place,M,full_core=True) for place in ('infinity','zero')]
 compiled=[compile_components(d,beta,profile_precision) for d in data]
 if normalize_units:compiled=[reduce_components(*value) for value in compiled]
 return data,compiled
