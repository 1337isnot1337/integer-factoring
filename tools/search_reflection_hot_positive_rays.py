#!/usr/bin/env python3
"""Exact F4^8 hot-defect profiles at prescribed positive physical heights."""
from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp

from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_hot_core_profiles import profile_data, optimize_scale
from reflection_hot_state_coordinates import gaussian_coordinates
from reflection_hot_state_gauges import HotStatePressure, defect_block
from search_reflection_hot_core_filters import polynomial_specs, factor_from_row, combine


RAYS = (
    ((2,2,2,2),(0,0,0,0)),
    ((0,0,0,0),(2,2,2,2)),
    ((2,2,2,2),(2,2,2,2)),
    ((2,0,0,2),(2,0,0,2)),
    ((0,2,2,0),(0,2,2,0)),
)


def compile_sources(interval, ramification=1, window=512):
    answer = []
    for place in ('zero','infinity'):
        coordinates = gaussian_coordinates(interval,place,ramification,1,window,True)
        compiled = []
        for sector in range(2):
            components = [defect_block(coordinates,{0:{0:1}},sector)]
            components += [defect_block(coordinates,{mask:poly},sector) for mask,poly in polynomial_specs(ramification)]
            assert all(md['exact_polynomial_matrix'] for _,_,md in components)
            shifts = [[min((offsets[i][j] for matrix,offsets,_ in components if any(matrix[i][j])),default=0)
                       for j in range(8)] for i in range(8)]
            sources = [[[(entry[0] << (offsets[i][j]-shifts[i][j]),entry[1] << (offsets[i][j]-shifts[i][j]))
                         if any(entry) else (0,0) for j,entry in enumerate(row)] for i,row in enumerate(matrix)]
                       for matrix,offsets,_ in components]
            bits = max(value.bit_length() for source in sources for row in source for entry in row for value in entry)
            precision = max(window,1 << (bits-1).bit_length())
            assert precision <= 8192
            compiled.append((sources,shifts,precision))
        answer.append((coordinates,compiled))
    return answer


def row_oracles(compiled,row):
    return [[MatrixOracle(combine(sources,row),8,shifts,precision) for sources,shifts,precision in blocks]
            for _,blocks in compiled]


def profiles(compiled,oracles,ray):
    result = []
    for place,((coordinates,_),side) in enumerate(zip(compiled,oracles)):
        orders = []
        left,right = coordinates['left'],coordinates['right']
        rw = sum((ray[place]['xyuv'.index(letter)]*np.asarray(labels) for letter,labels in left['labels'].items()),np.zeros(16))
        cw = sum((ray[place]['xyuv'.index(letter)]*np.asarray(labels) for letter,labels in right['labels'].items()),np.zeros(16))
        for sector,oracle in enumerate(side):
            while True:
                try:
                    part = oracle.profile(rw[coordinates['row_groups'][sector]],cw[coordinates['column_groups'][sector]])['orders']
                    break
                except ArithmeticError:
                    if oracle.precision >= 8192:
                        raise
                    oracle.precision *= 2
            assert all(value == int(value) for value in part)
            orders += list(map(int,part))
        result.append(tuple(sorted(orders)))
    return tuple(result)


def relative_cost(orders,ell):
    return sum(float(logsumexp(-ell*np.asarray(part))) for part in orders)


def physical_log_partition(interval,ramification,ray,ell,rate):
    natural = ramification*np.asarray(physical_bit_energies(interval))
    gx = natural[-1]+sum(ray[place][leg] for place in range(2) for leg in (0,3))
    gy = natural[-1]+sum(ray[place][leg] for place in range(2) for leg in (1,2))
    return float(np.logaddexp(0,-ell*gx/rate)+np.logaddexp(0,-ell*gy/rate)
                 +2*np.logaddexp(0,-ell*natural[:-1]/rate).sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval',type=int,default=16)
    parser.add_argument('--ramification',type=int,default=1)
    parser.add_argument('--seconds',type=float,default=1800.)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    start = time.monotonic();compiled = compile_sources(args.interval,args.ramification)
    baseline = optimize_scale(profile_data(args.interval,args.ramification,{0:{0:1}}))
    scales = [baseline['log_s']*factor for factor in (.5,1.,2.)]
    spectators = np.asarray(compiled[0][0]['data']['spectator_orders'])
    base_oracles = row_oracles(compiled,(0,)*8)
    gaussian_profiles = [profiles(compiled,base_oracles,ray) for ray in RAYS]
    # This complete coefficient-path replay is independent of the compiled
    # shared shifts and F4 combination routine.
    controls = 0
    for row in ((0,)*8,(0,1,1,1,1,0,0,1),(2,3,1,0,2,1,3,1)):
        model = HotStatePressure(args.interval,factor_from_row(row,args.ramification),args.ramification,internal_regular=True)
        oracles = row_oracles(compiled,row)
        for ray in RAYS:
            offsets = np.zeros((2,4,2));offsets[:,:,1] = ray
            _,detail = model.evaluate(1.,offsets)
            assert profiles(compiled,oracles,ray) == tuple(tuple(sorted(map(int,part))) for part in detail['active_orders'])
            controls += 2
    classes = [dict() for _ in RAYS]
    digest = sha256();completed = 0
    result = {'status':'exact finite F4^8 profiles at prescribed integer heights and sampled positive scales',
              'interval':args.interval,'ramification':args.ramification,'planned_candidates':4**8,
              'rays':RAYS,'scales':scales,'gaussian_natural':baseline,
              'independent_compiled_profile_controls':controls,
              'minimum_same_height_relative_log_cost':0.,'best_sampled_baseline_excess':math.inf,
              'best_sampled_critical_excess':math.inf,'best_sampled_row':None,'complete':False}
    root = Path(__file__).resolve().parent.parent
    files = ['tools/search_reflection_hot_positive_rays.py','tools/reflection_hot_state_gauges.py',
             'tools/reflection_hot_state_coordinates.py','tools/reflection_hot_state_support.py',
             'tools/search_reflection_hot_core_filters.py']
    def save():
        result['completed_candidates'] = completed
        result['class_counts_by_ray'] = list(map(len,classes))
        result['profile_digest_sha256'] = digest.hexdigest()
        result['seconds'] = time.monotonic()-start
        result['source_sha256'] = {name:sha256((root/name).read_bytes()).hexdigest() for name in files}
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for row in itertools.product(range(4),repeat=8):
        if time.monotonic()-start >= args.seconds:
            break
        oracles = row_oracles(compiled,row)
        for ray_index,ray in enumerate(RAYS):
            orders = profiles(compiled,oracles,ray)
            digest.update(json.dumps((row,ray_index,orders),separators=(',',':')).encode())
            if orders not in classes[ray_index]:
                classes[ray_index][orders] = row
                for ell in scales:
                    cost = relative_cost(orders,ell)
                    relative = cost-relative_cost(gaussian_profiles[ray_index],ell)
                    result['minimum_same_height_relative_log_cost'] = min(result['minimum_same_height_relative_log_cost'],relative)
                    full = cost+2*float(np.logaddexp(0,-ell*spectators).sum())
                    critical = full-math.sqrt(2)*physical_log_partition(args.interval,args.ramification,ray,ell,math.sqrt(2))
                    natural_rate = 2*baseline['rho']
                    excess = full-natural_rate*physical_log_partition(args.interval,args.ramification,ray,ell,natural_rate)
                    result['best_sampled_critical_excess'] = min(result['best_sampled_critical_excess'],critical)
                    if excess < result['best_sampled_baseline_excess']:
                        result['best_sampled_baseline_excess'] = excess
                        result['best_sampled_row'] = {'row':row,'ray_index':ray_index,'log_s':ell,'profiles':orders}
        completed += 1
        if completed%1024==0:
            save();print(json.dumps({'completed':completed,'classes':result['class_counts_by_ray'],
                                    'minimum_relative':result['minimum_same_height_relative_log_cost'],
                                    'baseline_excess':result['best_sampled_baseline_excess'],'seconds':result['seconds']}),flush=True)
    result['complete'] = completed == 4**8
    save()


if __name__=='__main__':
    main()
