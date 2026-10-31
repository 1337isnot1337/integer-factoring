#!/usr/bin/env python3
"""Exhaust an explicit F4 grid of endpoint-valid hot-core perturbations.

The eight coordinates give six quadratic coefficients and two quartic
coefficients. Their prescribed Laurent polynomials vanish at w=1.
Every active Smith profile is exact; pressures guide further work and
are floating numerical values, not strict rate certificates.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
from reflection_hot_core_profiles import compile_interval,profile_data,optimize_scale,pressure
from check_rectangular_f4_repair import fmul
from native_matrix_smith_f4 import MatrixOracle

QUADRATICS=(3,5,6,9,10,12)


def polynomial_specs(P):
    return [(mask,{-4*P:1,-2*P:1} if mask==3 else {0:1,2*P:1}) for mask in QUADRATICS]+[(15,{-4*P:1,0:1}),(15,{0:1,4*P:1})]


def factor_from_row(row,P):
    factor={0:{0:1}}
    for coefficient,(mask,poly) in zip(row,polynomial_specs(P)):
        if not coefficient:continue
        target=factor.setdefault(mask,{})
        for power,value in poly.items():
            target[power]=target.get(power,0)^fmul(coefficient,value)
            if not target[power]:del target[power]
    return {mask:poly for mask,poly in factor.items() if poly}


def compile_sources(M,P,place):
    data=compile_interval(M,place,P);sign=1 if place=='zero' else -1
    sources=[data['components'][0]]
    for mask,poly in polynomial_specs(P):
        matrix=[]
        for row in data['components'][mask]:
            output=[]
            for entry in row:
                value={}
                for a,x in entry.items():
                    for b,y in poly.items():
                        degree=a+sign*b;value[degree]=value.get(degree,0)^fmul(x,y)
                        if not value[degree]:del value[degree]
                output.append(value)
            matrix.append(output)
        sources.append(matrix)
    shifts=[[min((power for source in sources for power in source[i][j]),default=0) for j in range(8)] for i in range(8)]
    packed=[[[tuple(sum(1<<(power-shifts[i][j]) for power,value in entry.items() if value>>bit&1) for bit in range(2))
               for j,entry in enumerate(row)] for i,row in enumerate(source)] for source in sources]
    width=max(1,max(value.bit_length() for source in packed for row in source for pair in row for value in pair))
    precision=max(512,1<<(2*width-1).bit_length())
    if precision>8192:raise ArithmeticError('The compiled source exceeds the supported complete profile budget.')
    return packed,shifts,precision


def combine(packed,coefficients):
    matrix=[]
    for i in range(8):
        row=[]
        for j in range(8):
            a,b=packed[0][i][j]
            for coefficient,source in zip(coefficients,packed[1:]):
                x,y=source[i][j]
                if coefficient==1:a^=x;b^=y
                elif coefficient==2:a^=y;b^=x^y
                elif coefficient==3:a^=x^y;b^=x
            row.append((a,b))
        matrix.append(row)
    return matrix


def exact_profiles(compiled,row):
    result=[]
    for packed,shifts,base_precision in compiled:
        matrix=combine(packed,row);precision=base_precision
        while precision<=8192:
            try:orders=MatrixOracle(matrix,8,shifts,precision).profile()['orders'];break
            except ArithmeticError:precision*=2
        else:raise ArithmeticError('The complete active profile remained unresolved.')
        if any(x!=int(x) for x in orders):raise ArithmeticError('Expected integer Smith orders.')
        result.append(tuple(map(int,orders)))
    return tuple(result)


def log_partition(orders,ell):
    ground=min(orders);shifted=[order-ground for order in orders];shifted.remove(0)
    return -ground*ell+math.log1p(math.fsum(math.exp(-ell*order) for order in shifted))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--interval',type=int,default=32);parser.add_argument('--ramification',type=int,default=1)
    parser.add_argument('--seconds',type=float,default=3600);parser.add_argument('--top',type=int,default=40);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    M,P=args.interval,args.ramification;start=time.monotonic();deadline=start+args.seconds
    compiled=[compile_sources(M,P,place) for place in ('zero','infinity')]
    control=profile_data(M,P,{0:{0:1}});baseline=optimize_scale(control);scales=[baseline['log_s']*factor for factor in (.5,.75,1,1.5,2)]
    # Independent literal coefficient expansion controls the compiled grid.
    controls=0
    for row in ((0,)*8,(1,0,2,0,3,0,1,2),(3,2,1,3,2,1,3,2)):
        direct=profile_data(M,P,factor_from_row(row,P))
        assert exact_profiles(compiled,row)==tuple(tuple(x) for x in direct['profiles']);controls+=2
    base_log=[sum(log_partition(orders,ell) for orders in control['profiles']) for ell in scales]
    result={'status':'finite F4 coefficient grid, exact profiles, exploratory floating pressures','interval':M,'ramification':P,
            'coefficient_grid':'F4^8','planned_candidates':4**8,'scales':scales,'gaussian':baseline,
            'independent_literal_profile_controls':controls,'profile_classes':[],'optimized':[],
            'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in (Path(__file__).name,'reflection_hot_core_profiles.py')}}
    classes={};completed=0;unresolved=[];best=0.;best_row=None
    def save():
        result['completed_candidates']=completed;result['complete']=completed==4**8 and not unresolved
        result['unresolved']=unresolved;result['minimum_sampled_relative_log_cost']=best;result['best_sampled_row']=best_row
        result['same_gauge_gain_below_1e_minus_9']=best<-1e-9;result['seconds']=time.monotonic()-start
        result['profile_classes']=sorted(classes.values(),key=lambda record:min(record['relative_log_costs']))
        tmp=args.output.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(args.output)
    for row in itertools.product(range(4),repeat=8):
        if time.monotonic()>=deadline:break
        try:profiles=exact_profiles(compiled,row)
        except ArithmeticError as error:unresolved.append({'row':row,'reason':str(error)});completed+=1;continue
        if profiles not in classes:
            costs=[sum(log_partition(orders,ell) for orders in profiles)-base for ell,base in zip(scales,base_log)]
            classes[profiles]={'representative':row,'profiles':profiles,'rows':0,'relative_log_costs':costs}
            if min(costs)<best:best=min(costs);best_row=row
        classes[profiles]['rows']+=1;completed+=1
        if completed%2048==0:
            save();print(json.dumps({'completed':completed,'classes':len(classes),'minimum_relative':best,'seconds':result['seconds']}),flush=True)
    save()
    for record in result['profile_classes'][:args.top]:
        if time.monotonic()>=deadline:break
        data=profile_data(M,P,factor_from_row(record['representative'],P));opt=optimize_scale(data)
        result['optimized'].append({'row':record['representative'],'optimization':opt,'factor':data['factor'],'degree_bound':data['polynomial_degree_bound']})
    save()

if __name__=='__main__':main()
