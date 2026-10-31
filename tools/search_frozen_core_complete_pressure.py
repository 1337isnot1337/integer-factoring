#!/usr/bin/env python3
"""Exhaust a frozen-row family at fixed rational gauges and scale.

Filtered linear coefficient flags partition every endpoint-valid F4 row
by its exact least entry order at the two cuts. A single-entry lower
bound can reject a whole cell. Every remaining row receives complete
integer-weight Smith profiles and an outward-rounded, then exact if
needed, comparison with the Gaussian row at the same gauges.

This is a finite fixed-gauge result. It does not bound other gauges,
other scales, or a larger crossing family.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import heapq
import itertools
import json
import math
from pathlib import Path
import time

from check_rectangular_f4_repair import fmul,finv,ZERO
from check_three_mode_pressure_certificate import DyadicIntervals
from frozen_core_row_profiles import compile_interval,combine_components
from native_matrix_smith_f4 import MatrixOracle
from search_frozen_core_critical_pressure import AffineSpace


def dot(first,second):
    value=0
    for a,b in zip(first,second):value^=fmul(a,b)
    return value


def add_equations(space,rows):
    out=AffineSpace();out.rows={p:list(row) for p,row in space.rows.items()}
    try:
        for row in rows:out.add(row)
    except ArithmeticError:return None
    return out


def cardinality(space):
    return 0 if space is None else 4**(12-len(space.rows))


def coefficient_section(flags):
    """Choose one row per common coefficient fiber and record its size."""
    basis={}
    for flag in flags:
        for entry in flag:
            row=list(entry['form'])
            for pivot,old in sorted(basis.items()):
                if row[pivot]:
                    c=row[pivot];row=[x^fmul(c,y) for x,y in zip(row,old)]
            pivot=next((i for i,c in enumerate(row) if c),None)
            if pivot is not None:
                c=finv(row[pivot]);basis[pivot]=[fmul(c,x) for x in row]
    endpoint=[int(i in (0,10)) for i in range(12)]
    for pivot,row in sorted(basis.items()):
        if endpoint[pivot]:
            c=endpoint[pivot];endpoint=[x^fmul(c,y) for x,y in zip(endpoint,row)]
    assert not any(endpoint),'The coefficient kernel must preserve the endpoint equation.'
    free=[i for i in range(12) if i not in basis]
    section=add_equations(AffineSpace(),[[int(j==i) for j in range(12)] for i in free])
    assert section is not None
    return section,4**len(free),{'coefficient_rank':len(basis),'fixed_coordinates':free}


def entry_flag(comp,shifts,metadata,gauge,M):
    """Exact coefficient row-space filtration, stopping at full rank."""
    n=len(comp[0]);heap=[];rows=metadata['row_charges'];cols=metadata['column_charges']
    for i in range(n):
        for j in range(n):
            bits=0
            for source in comp:bits|=source[i][j][0]|source[i][j][1]
            if not bits:continue
            base=Fraction(shifts[i][j])+gauge*M*(rows[i]-cols[j])
            if base.denominator!=1:raise ValueError('These gauges require a further fixed ramification.')
            bit=bits&-bits;exponent=bit.bit_length()-1
            heapq.heappush(heap,(int(base)+exponent,i,j,int(base),bits))
    basis={};flag=[];terms=0
    while heap:
        order,i,j,base,bits=heapq.heappop(heap);bit=bits&-bits;exponent=bit.bit_length()-1
        vec=tuple(((source[i][j][0]>>exponent)&1)|(((source[i][j][1]>>exponent)&1)<<1) for source in comp)
        reduced=list(vec)
        for pivot,row in sorted(basis.items()):
            if reduced[pivot]:
                coefficient=reduced[pivot];reduced=[x^fmul(coefficient,y) for x,y in zip(reduced,row)]
        pivot=next((k for k,x in enumerate(reduced) if x),None)
        if pivot is not None:
            coefficient=finv(reduced[pivot]);basis[pivot]=[fmul(coefficient,x) for x in reduced]
            flag.append({'order':order,'form':vec})
            if len(basis)==12:break
        terms+=1;bits^=bit
        if bits:
            exponent=(bits&-bits).bit_length()-1
            heapq.heappush(heap,(base+exponent,i,j,base,bits))
    return flag,{'independent_forms':len(flag),'terms_examined':terms+1}


class ProductComparison:
    def __init__(self,scale,baseline,precision=192):
        self.scale=scale;self.baseline=tuple(tuple(x) for x in baseline)
        self.arithmetic=DyadicIntervals(precision);self.root=self.arithmetic.rational(scale.numerator,scale.denominator)
        self.baseline_interval=self.product(self.baseline)
        self.exact_baseline=None;self.exact_fallbacks=0

    @lru_cache(None)
    def power(self,exponent):return self.arithmetic.power(self.root,exponent)

    def product(self,profiles):
        sums=[]
        for profile in profiles:
            terms=[self.power(-int(order)) for order in profile]
            sums.append((sum(x[0] for x in terms),sum(x[1] for x in terms)))
        return self.arithmetic.product(sums)

    def ground_rejects(self,order,multiplicity=1):
        # This lower bound is a single actual term in each partition.
        return multiplicity*self.power(-order)[0]>=self.baseline_interval[1]

    def compare(self,profiles):
        profiles=tuple(tuple(x) for x in profiles)
        if profiles==self.baseline:return 0
        interval=self.product(profiles)
        if interval[0]>self.baseline_interval[1]:return 1
        if interval[1]<self.baseline_interval[0]:return -1
        self.exact_fallbacks+=1
        if self.exact_baseline is None:
            self.exact_baseline=math.prod(sum((self.scale**(-x) for x in profile),Fraction()) for profile in self.baseline)
        value=math.prod(sum((self.scale**(-x) for x in profile),Fraction()) for profile in profiles)
        return (value>self.exact_baseline)-(value<self.exact_baseline)


def leading_rank(matrix,shifts,metadata,gauge,M,ground):
    """Rank over F4 of the exact first weighted coefficient matrix."""
    n=len(matrix);rows=[]
    for i in range(n):
        planes=[0,0]
        for j in range(n):
            exponent=Fraction(ground-shifts[i][j])-gauge*M*(metadata['row_charges'][i]-metadata['column_charges'][j])
            assert exponent.denominator==1
            exponent=int(exponent)
            if exponent>=0:
                for side in range(2):planes[side]|=((matrix[i][j][side]>>exponent)&1)<<j
        rows.append(tuple(planes))
    basis={}
    def multiply(row,c):
        a,b=row
        return row if c==1 else ((b,a^b) if c==2 else (a^b,a))
    for row in rows:
        while row[0]|row[1]:
            bit=(row[0]|row[1])&-(row[0]|row[1]);pivot=bit.bit_length()-1
            c=int(bool(row[0]&bit))+2*int(bool(row[1]&bit))
            if pivot not in basis:
                basis[pivot]=multiply(row,finv(c));break
            other=multiply(basis[pivot],c);row=(row[0]^other[0],row[1]^other[1])
    return len(basis)


def profiles(compiled,g,gauges,M,matrices=None):
    result=[]
    for side,((comp,shifts,metadata),gauge) in enumerate(zip(compiled,gauges)):
        if not metadata['exact_polynomial_matrix']:raise ArithmeticError('Complete polynomial coefficients are required.')
        rw=[gauge*M*c for c in metadata['row_charges']];cw=[-gauge*M*c for c in metadata['column_charges']]
        if any(x.denominator!=1 for x in rw+cw):raise ValueError('Use integral weighted row and column orders.')
        matrix=combine_components(comp,g) if matrices is None else matrices[side]
        oracle=MatrixOracle(matrix,len(comp[0]),shifts,precision=metadata['profile_precision'])
        raw=oracle.profile(list(map(int,rw)),list(map(int,cw)))['orders']
        if any(x!=int(x) for x in raw):raise ArithmeticError('An integer-weight profile returned a noninteger order.')
        result.append(tuple(map(int,raw)))
    return tuple(result)


def screen(compiled,M,gauges,scale,deadline,cap):
    flags=[];flag_controls=[]
    for (comp,shifts,metadata),gauge in zip(compiled,gauges):
        flag,control=entry_flag(comp,shifts,metadata,gauge,M);flags.append(flag);flag_controls.append(control)
    gaussian=profiles(compiled,[1]+[0]*11,gauges,M);compare=ProductComparison(scale,gaussian)
    section,fiber_multiplicity,section_metadata=coefficient_section(flags)
    cells=[];covered=pruned=rank_pruned=queried=unresolved=equal=larger=0;improvements=[];digest=hashlib.sha256()
    for i,f in enumerate(flags[0]):
        for j,h in enumerate(flags[1]):
            prefix=[x['form'] for x in flags[0][:i]]+[x['form'] for x in flags[1][:j]]
            space=add_equations(section,prefix)
            if space is None:continue
            left=add_equations(space,[f['form']]);right=add_equations(space,[h['form']]);both=add_equations(space,[f['form'],h['form']])
            count=cardinality(space)-cardinality(left)-cardinality(right)+cardinality(both)
            if not count:continue
            covered+=count;cell={'first_nonzero_flags':[i,j],'ground_orders':[f['order'],h['order']],'rows':count,'kernel_dimension':12-len(space.rows)}
            if compare.ground_rejects(f['order']+h['order']):
                cell['status']='certified single-entry lower bound';pruned+=count
            elif time.monotonic()>=deadline:
                cell['status']='not exhausted within the requested budget';unresolved+=count
            else:
                made=0;cell_queries=0;cell_rank_pruned=0;cell_unresolved=0;free=12-len(space.rows)
                for assignment in itertools.product(range(4),repeat=free):
                    g=space.solution(assignment)
                    if not dot(f['form'],g) or not dot(h['form'],g):continue
                    if time.monotonic()>=deadline:
                        unresolved+=count-made;cell_unresolved+=count-made;break
                    matrices=[combine_components(comp,g) for comp,_,_ in compiled]
                    ranks=[leading_rank(matrix,shifts,metadata,gauge,M,entry['order'])
                           for matrix,(_,shifts,metadata),gauge,entry in zip(matrices,compiled,gauges,(f,h))]
                    assert min(ranks)>0
                    if compare.ground_rejects(f['order']+h['order'],math.prod(ranks)):
                        rank_pruned+=1;cell_rank_pruned+=1;made+=1;continue
                    if queried>=cap:
                        unresolved+=1;cell_unresolved+=1;made+=1;continue
                    try:
                        actual=profiles(compiled,g,gauges,M,matrices)
                    except ArithmeticError as error:
                        unresolved+=1;cell_unresolved+=1;made+=1
                        cell.setdefault('profile_errors',[])
                        if len(cell['profile_errors'])<3:cell['profile_errors'].append(str(error))
                        continue
                    # The filtered flag must equal the true first Smith order.
                    assert actual[0][0]==f['order'] and actual[1][0]==h['order']
                    assert all(profile.count(entry['order'])==rank for profile,entry,rank in zip(actual,(f,h),ranks))
                    relation=compare.compare(actual);made+=1;queried+=1;cell_queries+=1
                    digest.update(bytes(g));digest.update(json.dumps(actual,separators=(',',':')).encode())
                    if relation==0:equal+=1
                    elif relation>0:larger+=1
                    else:improvements.append({'g':g,'orders':actual})
                cell['rows_examined']=made;cell['profile_queries']=cell_queries;cell['leading_rank_pruned']=cell_rank_pruned;cell['unresolved_rows']=cell_unresolved;cell['status']='exhausted' if not cell_unresolved else 'unresolved within the requested budget or profile precision'
            cells.append(cell)
    assert covered*fiber_multiplicity==4**11,(covered,fiber_multiplicity,4**11)
    assert pruned+rank_pruned+queried+unresolved==covered
    for cell in cells:
        cell['distinct_tensor_rows']=cell['rows'];cell['rows']*=fiber_multiplicity
        if 'unresolved_rows' in cell:cell['unresolved_rows']*=fiber_multiplicity
    return {'gauges':list(map(str,gauges)),'scale':str(scale),'total_affine_rows':covered*fiber_multiplicity,'distinct_tensor_rows':covered,'parameter_fiber_multiplicity':fiber_multiplicity,'coefficient_section':section_metadata,'certifiably_pruned':(pruned+rank_pruned)*fiber_multiplicity,'entry_bound_pruned':pruned*fiber_multiplicity,'leading_rank_pruned':rank_pruned*fiber_multiplicity,'profile_queries':queried,'profile_covered_rows':queried*fiber_multiplicity,'unresolved_rows':unresolved*fiber_multiplicity,'equal_to_gaussian':equal*fiber_multiplicity,'strictly_larger_than_gaussian':larger*fiber_multiplicity,'improvements':improvements,'complete':not unresolved,'profile_digest_sha256':digest.hexdigest(),'exact_comparison_fallbacks':compare.exact_fallbacks,'dyadic_precision':compare.arithmetic.precision,'gaussian_profiles':gaussian,'flags':flags,'flag_controls':flag_controls,'cells':cells}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--M',type=int,default=32);parser.add_argument('--beta',type=int,default=1);parser.add_argument('--grid',default='0');parser.add_argument('--scale',type=Fraction);parser.add_argument('--seconds',type=float,default=3600);parser.add_argument('--cap',type=int,default=100000);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    started=time.monotonic();_,compiled=compile_interval(args.M,args.beta)
    scale=args.scale or Fraction(round(math.exp(math.sqrt(2)*math.log(args.M)/args.M)*10**8),10**8)
    if scale<=1:raise ValueError('The exact profile scale must exceed one.')
    grid=[Fraction(x) for x in args.grid.split(',')];result={'status':'exact finite comparison with the Gaussian background; no claim for other gauges or scales','M':args.M,'beta':args.beta,'scale':str(scale),'runs':[],'constructor_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('frozen_core_row_profiles.py',Path(__file__).name)}}
    for gauges in itertools.product(grid,repeat=2):
        run=screen(compiled,args.M,gauges,scale,started+args.seconds,args.cap);result['runs'].append(run);result['seconds']=time.monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:run[k] for k in ['gauges','total_affine_rows','certifiably_pruned','profile_queries','unresolved_rows','complete']}|{'improvements':len(run['improvements']),'seconds':result['seconds']}),flush=True)
        if time.monotonic()-started>=args.seconds:break
if __name__=='__main__':main()
