#!/usr/bin/env python3
"""Complete F4^8 coefficient grids on fixed compensated negative-height rays."""
from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_hot_core_profiles import profile_data,optimize_scale
from reflection_hot_state_coordinates import gaussian_coordinates
from reflection_hot_state_gauges import defect_block
from reflection_hot_compensated_gauges import (
    CompensatedHotStatePressure,retained_internal_factor,exterior_product,physical_bulks,occupations)
from search_reflection_hot_core_filters import polynomial_specs,factor_from_row
from search_reflection_hot_positive_rays import row_oracles,relative_cost


RAYS=(
    ((1,0),(0,1)),((0,1),(1,0)),((1,1),(0,0)),
    ((0,0),(1,1)),((1,1),(1,1)),((2,0),(0,2)),((2,2),(2,2)),
)


def compile_sources(interval):
    answer=[]
    for place in ('zero','infinity'):
        coordinates=gaussian_coordinates(interval,place,1,1,512,True)
        retained=retained_internal_factor(interval,place)
        factors=[retained]+[exterior_product({mask:poly},retained) for mask,poly in polynomial_specs(1)]
        blocks=[]
        for sector in range(2):
            components=[defect_block(coordinates,factor,sector) for factor in factors]
            assert all(md['exact_polynomial_matrix'] for _,_,md in components)
            shifts=[[min((offset[i][j] for matrix,offset,_ in components if any(matrix[i][j])),default=0)
                     for j in range(8)] for i in range(8)]
            sources=[[[(a<<(offset[i][j]-shifts[i][j]),b<<(offset[i][j]-shifts[i][j])) if a or b else (0,0)
                       for j,(a,b) in enumerate(row)] for i,row in enumerate(matrix)] for matrix,offset,_ in components]
            bits=max(value.bit_length() for source in sources for row in source for entry in row for value in entry)
            blocks.append((sources,shifts,max(512,1<<(bits-1).bit_length())))
        answer.append((coordinates,blocks))
    return answer


def profiles(compiled,oracles,ray):
    bulk=np.asarray(ray);physical=physical_bulks(bulk);answer=[]
    for place,((coordinates,_),side) in enumerate(zip(compiled,oracles)):
        left,right=coordinates['left'],coordinates['right']
        rows=bulk[place,0]*occupations(left,'left');columns=bulk[place,1]*occupations(right,'right')
        for letter,labels in left['labels'].items():
            rows-=2*physical[place,'xyuv'.index(letter)]*np.asarray(labels)
        for letter,labels in right['labels'].items():
            columns-=2*physical[place,'xyuv'.index(letter)]*np.asarray(labels)
        orders=[]
        for sector,oracle in enumerate(side):
            while True:
                try:
                    part=oracle.profile(rows[coordinates['row_groups'][sector]],columns[coordinates['column_groups'][sector]])['orders']
                    break
                except ArithmeticError:
                    if oracle.precision>=8192:
                        raise
                    oracle.precision*=2
            assert all(value==int(value) for value in part)
            orders+=list(map(int,part))
        answer.append(tuple(sorted(orders)))
    return tuple(answer)


def physical_partition(interval,ray,ell,rate):
    bulk=sum(map(sum,ray));natural=np.asarray(physical_bit_energies(interval),dtype=float)
    return 2*float(np.logaddexp(0,-ell*(natural[-1]-bulk)/rate)
                   +np.logaddexp(0,-ell*(natural[:-1]+bulk)/rate).sum())


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--interval',type=int,default=64)
    parser.add_argument('--seconds',type=float,default=2400.);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic();compiled=compile_sources(args.interval)
    baseline=optimize_scale(profile_data(args.interval,1,{0:{0:1}}));scales=[baseline['log_s']*x for x in (.5,1.,2.)]
    spectators=np.asarray(compiled[0][0]['data']['spectator_orders'])
    zero=row_oracles(compiled,(0,)*8);gaussian=[profiles(compiled,zero,ray) for ray in RAYS]
    controls=0
    for row in ((0,)*8,(0,1,1,1,1,0,0,1),(2,3,1,0,2,1,3,1)):
        direct=CompensatedHotStatePressure(args.interval,factor_from_row(row,1));oracles=row_oracles(compiled,row)
        for ray in RAYS:
            offsets=np.zeros((2,4,2));offsets[:,:,1]=-2*physical_bulks(np.asarray(ray))
            _,detail=direct.evaluate(1.,offsets,ray)
            assert profiles(compiled,oracles,ray)==tuple(tuple(sorted(map(int,part))) for part in detail['active_orders'])
            controls+=2
    classes=[{} for _ in RAYS];digest=sha256();completed=0
    result={'status':'complete finite coefficient grid at fixed compensated integer heights; numerical partition comparisons',
            'interval':args.interval,'planned_candidates':4**8,'boundary_bulk_rays':RAYS,
            'hot_relative_height_rule':'minus twice the bulk of its boundary','scales':scales,'gaussian_natural':baseline,
            'independent_compiled_profile_controls':controls,'minimum_same_height_relative_log_cost':0.,
            'best_sampled_baseline_excess':math.inf,'best_sampled_critical_excess':math.inf,'best_sampled_row':None,'complete':False}
    root=Path(__file__).resolve().parent.parent
    files=['tools/search_reflection_hot_compensated_rays.py','tools/reflection_hot_compensated_gauges.py',
           'tools/search_reflection_hot_positive_rays.py','tools/reflection_hot_state_coordinates.py']
    def save():
        result['completed_candidates']=completed;result['class_counts_by_ray']=list(map(len,classes))
        result['profile_digest_sha256']=digest.hexdigest();result['seconds']=time.monotonic()-start
        result['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest() for name in files}
        temporary=args.output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for row in itertools.product(range(4),repeat=8):
        if time.monotonic()-start>=args.seconds:
            break
        oracles=row_oracles(compiled,row)
        for index,ray in enumerate(RAYS):
            orders=profiles(compiled,oracles,ray)
            digest.update(json.dumps((row,index,orders),separators=(',',':')).encode())
            if orders in classes[index]:
                continue
            classes[index][orders]=row
            for ell in scales:
                active=relative_cost(orders,ell)
                relative=active-relative_cost(gaussian[index],ell)
                result['minimum_same_height_relative_log_cost']=min(result['minimum_same_height_relative_log_cost'],relative)
                full=active+sum(float(np.logaddexp(0,-ell*(spectators+sum(side))).sum()) for side in ray)
                critical=full-math.sqrt(2)*physical_partition(args.interval,ray,ell,math.sqrt(2))
                rate=2*baseline['rho'];excess=full-rate*physical_partition(args.interval,ray,ell,rate)
                result['best_sampled_critical_excess']=min(result['best_sampled_critical_excess'],critical)
                if excess<result['best_sampled_baseline_excess']:
                    result['best_sampled_baseline_excess']=excess;result['best_sampled_row']={'row':row,'ray_index':index,'log_s':ell,'profiles':orders}
        completed+=1
        if completed%1024==0:
            save();print(json.dumps({'completed':completed,'minimum_relative':result['minimum_same_height_relative_log_cost'],
                                    'baseline_excess':result['best_sampled_baseline_excess'],'seconds':result['seconds']}),flush=True)
    result['complete']=completed==4**8;save()


if __name__=='__main__':
    main()
