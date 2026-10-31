#!/usr/bin/env python3
"""Explore exact hot-core defects with all monotone physical bit gauges."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from reflection_hot_core_profiles import profile_data, optimize_scale
from reflection_hot_state_gauges import HotStatePressure
from search_reflection_hot_core_filters import factor_from_row


EXCEPTIONAL = (
    (0,1,1,1,1,0,0,1), (0,1,1,1,1,1,0,1),
    (2,1,1,1,1,3,1,1), (1,1,1,1,1,1,1,1),
    (0,1,1,1,1,0,0,0), (0,2,2,2,2,0,0,0),
    (0,1,1,1,1,0,0,2))


def optimize_candidate(interval, ramification, row, baseline, rng, starts, maxfun, unrestricted=False):
    factor = factor_from_row(row, ramification)
    model = HotStatePressure(interval, factor, ramification, internal_regular=not unrestricted)
    initial = np.zeros(9);initial[0] = math.log(baseline['log_s'])
    points = [initial]
    if starts > 1:
        point = initial.copy();point[1:] = (-1 if unrestricted else 1)
        points.append(point)
    while len(points) < starts:
        point = initial.copy();point[0] += rng.uniform(-.4,.4)
        point[1:] = rng.uniform(-6 if unrestricted else 0,6,size=8)
        points.append(point)
    best = None
    def objective(parameters):
        ell = math.exp(parameters[0])
        heights = parameters[1:] if unrestricted else np.maximum(parameters[1:],0.)
        offsets = np.zeros((2,4,2));offsets[:,:,1] = heights.reshape(2,4)
        value, (gradient, ds), detail = model.evaluate(ell, offsets, True)
        nonlocal best
        if best is None or value < best['rho']:
            best = {'rho':value,'log_s':ell,'relative_bit_heights':heights.reshape(2,4).tolist(),
                    'active_orders':detail['active_orders'],
                    'core_type_x':detail['type_x_core'],'core_type_y':detail['type_y_core']}
        return value, np.r_[ell*ds,gradient[:,:,1].reshape(-1)]
    runs = []
    maximum_height = min(256.,8.*interval*ramification)
    for point in points:
        run = minimize(objective,point,jac=True,method='L-BFGS-B',
                       bounds=[(math.log(baseline['log_s']/8),math.log(baseline['log_s']*8))]
                              + [(-maximum_height if unrestricted else 0.,maximum_height)]*8,
                       options={'maxfun':maxfun,'maxiter':maxfun,'ftol':1e-13,'gtol':2e-9,'maxls':20})
        runs.append({'rho':float(run.fun),'evaluations':int(run.nfev),'success':bool(run.success),
                     'message':str(run.message)})
    return {'row':list(row),'factor':factor,'best':best,'starts':runs,'queries':model.queries,
            'active_modes':model.records[0][0]['data']['active_indices'],
            'spectator_orders':model.records[0][0]['data']['spectator_orders'],
            'coefficient_metadata':model.records[0][1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval',type=int,default=16)
    parser.add_argument('--ramification',type=int,default=1)
    parser.add_argument('--random',type=int,default=0)
    parser.add_argument('--starts',type=int,default=3)
    parser.add_argument('--maxfun',type=int,default=120)
    parser.add_argument('--seed',type=int,default=13371337)
    parser.add_argument('--representatives',type=Path)
    parser.add_argument('--unrestricted',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    start = time.monotonic();rng = np.random.default_rng(args.seed)
    rows = [(0,)*8]+list(EXCEPTIONAL)
    if args.representatives:
        source = json.loads(args.representatives.read_text())
        rows = list(dict.fromkeys([(0,)*8]+[tuple(record['representative']) for record in source['profile_classes']]))
    base_candidates = len(rows)
    while len(rows) < base_candidates+args.random:
        row = tuple(map(int,rng.integers(0,4,size=8)))
        if row not in rows:
            rows.append(row)
    baseline = optimize_scale(profile_data(args.interval,args.ramification,{0:{0:1}}))
    result = {'status':'exact coefficient matrices; exploratory floating physical-state optimization',
              'interval':args.interval,'ramification':args.ramification,'seed':args.seed,
              'unrestricted_state_heights':args.unrestricted,
              'representative_source':str(args.representatives) if args.representatives else None,
              'gaussian_natural':baseline,'planned_candidates':len(rows),'complete':False,'records':[]}
    if args.resume and args.output.exists():
        previous = json.loads(args.output.read_text())
        assert previous['interval'] == args.interval and previous['ramification'] == args.ramification
        assert previous.get('unrestricted_state_heights',False) == args.unrestricted
        assert [record['row'] for record in previous['records']] == [list(row) for row in rows[:len(previous['records'])]]
        result['records'] = previous['records']
        result['prior_seconds'] = previous.get('prior_seconds',0.) + previous.get('seconds',0.)
        if 'random_generator_state' in previous:
            rng.bit_generator.state = previous['random_generator_state']
        else:
            for _ in previous['records']:
                for _ in range(max(0,args.starts-2)):
                    rng.uniform(-.4,.4);rng.uniform(-6 if args.unrestricted else 0,6,size=8)
    root = Path(__file__).resolve().parent.parent
    files = ['tools/search_reflection_hot_state_gauges.py','tools/reflection_hot_state_support.py',
             'tools/reflection_hot_state_coordinates.py','tools/reflection_hot_state_gauges.py',
             'tools/native_fock_f4_parity.cpp']
    def save():
        result['seconds'] = time.monotonic()-start
        result['completed_candidates'] = len(result['records'])
        result['random_generator_state'] = rng.bit_generator.state
        result['source_sha256'] = {name:sha256((root/name).read_bytes()).hexdigest() for name in files}
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(args.output)
    for row in rows[len(result['records']):]:
        record = optimize_candidate(args.interval,args.ramification,row,baseline,rng,args.starts,args.maxfun,args.unrestricted)
        result['records'].append(record)
        save()
        print(json.dumps({'completed':len(result['records']),'planned':len(rows),'row':row,
                          'rho':record['best']['rho'],'queries':record['queries'],'seconds':result['seconds']}),flush=True)
    result['complete'] = True
    result['best'] = min(result['records'],key=lambda record:record['best']['rho'])
    save()


if __name__=='__main__':
    main()
