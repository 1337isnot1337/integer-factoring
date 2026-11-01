#!/usr/bin/env python3
"""Registered complete-profile and free-height pilot on cyclic Pauli words."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from random import Random
import subprocess
import time

import numpy as np
from scipy.optimize import minimize

from bounded_complete_support_fit import fit
from cyclic_pauli_crossings import construct, literal_seed, profiles
from cyclic_pauli_pressure import ScaledPauliPressure
from mixed_graph_retractions import FixedProfileSupport
from normalized_support_gap import normalized_gap
from search_f4_recentered_crossings import root_for_profiles


def make_loops(case):
    if 'loops' in case:
        return case['loops']
    rng, k, loops = Random(case['seed']), case['k'], []
    for _ in range(case['loop_count']):
        while True:
            shift = [rng.randrange(3) for _ in range(k)]
            a, b = shift+shift, [rng.randrange(3) for _ in range(2*k)]
            if any(a+b):
                break
        loops.append({'position': rng.randrange(k+1), 'a': a, 'b': b,
                      'coefficient': rng.randrange(1,4)})
    return loops


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refusing to overwrite an earlier run.')
    manifest = json.loads(args.manifest.read_text())
    root = Path(__file__).resolve().parents[1]
    for name, expected in manifest['source_sha256'].items():
        if sha256((root/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Registered source changed: {name}.')
    source = subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    if subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True):
        raise ValueError('Run from a clean frozen source.')
    begun = time.monotonic()
    report = {'status':'running_incomplete','source_commit':source,'manifest':manifest,'cases':[],
              'scope':'Finite literal qutrit sources, complete weighted block cuts and bounded numerical '
                       'whole-height fits. No global minimum, rigorous support or stronger score.'}

    def save():
        report['seconds'] = time.monotonic()-begun
        p = args.output.with_suffix('.tmp')
        p.write_text(json.dumps(report,indent=2)+'\n')
        p.replace(args.output)

    def best(candidate, metadata, planes):
        if candidate is None or candidate['rate'] >= report.get('best_support',{}).get('rate',math.inf):
            return
        report['best_support'] = dict(candidate, metadata=metadata)
        # Exact coefficient planes remain compact and portable in an NPZ artifact.
        temporary = args.output.with_suffix('.best-planes.tmp.npz')
        np.savez_compressed(temporary,low=planes[0],high=planes[1])
        temporary.replace(args.output.with_suffix('.best-planes.npz'))
        print(json.dumps({'event':'best','rate':candidate['rate'],'seconds':time.monotonic()-begun}),flush=True)

    save()
    for case in manifest['cases']:
        if time.monotonic()-begun >= manifest['total_seconds']:
            break
        planes, metadata = construct(case['k'], make_loops(case))
        if 'expected_planes_sha256' in case:
            assert metadata['planes_sha256'] == case['expected_planes_sha256']
        row = {'label':case['label'],'metadata':metadata,'status':'running_incomplete','fits':[]}
        report['cases'].append(row)
        try:
            unit, blocks = profiles(planes, metadata, precision=manifest['precision'])
        except (ArithmeticError,ValueError) as error:
            row.update(status='unit_profile_rejected',error=str(error))
            save()
            continue
        q,n,d = metadata['q'],metadata['q']**2,metadata['degree']
        row.update(unit_profiles=unit,unit_blocks=blocks,
            bare_mean_sum=sum(unit[0])+sum(unit[1]),
            bare_Jensen_excludes_all_angles_and_scales=sum(unit[0])+sum(unit[1])<=0)
        # raw infinity orders already include -degree; nonpositive mean sum
        # gives C >= (a+b) log(n) by the two Gibbs Jensen inequalities.
        fixed = FixedProfileSupport(unit,[0]*q,[0]*q)
        closed = []
        if case.get('bare_exclusion') or row['bare_Jensen_excludes_all_angles_and_scales']:
            row['bare_search_status'] = 'skipped_proved_bare_obstruction'
            row['bare_obstruction'] = case.get('bare_exclusion','nonpositive_raw_mean_Jensen')
            bare = {'scale':.05*math.sqrt(2),'rho':None}
            point = np.asarray((math.log(.05),0.))
            rate = fixed.evaluate(point,False)
            closed.append({'variables':point.tolist(),'objective':normalized_gap(rate,point),
                           'rate':rate,'status':'one_starting_point_only_no_bare_optimization'})
        else:
            bare = root_for_profiles(unit,[0]*q,[0]*q)
            row['bare_balanced_fit'] = bare
            for y in (-1.,0.,1.):
                def objective(point):
                    rate, gradient = fixed.evaluate(point)
                    return normalized_gap(rate,point,gradient)
                result = minimize(objective,(math.log(.08),y),jac=True,method='SLSQP',
                    bounds=((-9.,math.log(4.)),(-2.,2.)),options={'maxiter':120,'ftol':1e-11})
                closed.append({'variables':result.x.tolist(),'objective':float(result.fun),
                               'rate':fixed.evaluate(result.x,False),'success':bool(result.success)})
        row['bare_normalized_fits'] = closed
        best_closed = min(closed,key=lambda item:item['objective'])
        for item in closed:
            best({'rate':item['rate'],'parameters':[0.]*(8*(q-1))+item['variables'],
                  'from':'fixed_complete_unit_profiles'},metadata,planes)
        model = ScaledPauliPressure(planes,metadata,precision=manifest['precision'])
        local = np.random.default_rng(case['seed'])
        starts = []
        starts.append(('balanced-zero',np.r_[local.normal(0,.05,8*(q-1)),
            math.log(max(.05,bare['scale']/math.sqrt(2))),0.],True))
        heights = np.asarray([sum((2*i+1)*((word//3**(metadata['k']-1-i))%3)
                                   for i in range(metadata['k'])) for word in range(q)],float)[1:]
        zero = np.zeros(q-1)
        graded = np.r_[zero,zero,heights,heights,heights,heights,zero,zero]
        starts.append(('balanced-graded',np.r_[graded+local.normal(0,.05,len(graded)),math.log(.15),0.],True))
        starts.append(('free-angle-zero',np.r_[local.normal(0,.05,8*(q-1)),
            max(math.log(.05),best_closed['variables'][0]),best_closed['variables'][1]],False))
        for label, point, balanced in starts:
            remaining = manifest['total_seconds']-(time.monotonic()-begun)
            if remaining <= 0:
                break
            attempt = fit(None,point,model=model,seconds=min(manifest['trial_seconds'],remaining),
                          maxiter=manifest['maxiter'],balanced=balanced,height_bound=128.)
            attempt['start_label'] = label
            row['fits'].append(attempt)
            best(attempt['best'],metadata,planes)
            save()
        row['status'] = 'complete_attempts' if len(row['fits']) == 3 else 'budget_reached_incomplete'
        save()
        print(json.dumps({'label':case['label'],'degree':d,'bare_rho':bare['rho'],
            'zero_leading_vector_ranks':metadata['zero_R']['vector_matrix_ranks'],
            'fit_supports':[f['best']['rate'] if f['best'] else None for f in row['fits']]}),flush=True)
    report['status'] = ('complete_registered_attempts' if len(report['cases'])==len(manifest['cases']) and
        all(c['status'] in ('complete_attempts','unit_profile_rejected') for c in report['cases'])
        else 'budget_reached_incomplete')
    save()
    print(json.dumps({'status':report['status'],'seconds':report['seconds']}),flush=True)


if __name__=='__main__':
    main()
