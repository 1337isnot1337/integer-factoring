#!/usr/bin/env python3
"""Independent exact matrix and weighted-profile controls for full alphabets."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np

from check_full_alphabet_quartic_repair import build
from check_rectangular_defect_coordinates import point_control, original_cut, coordinate_weights
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from rectangular_defect_coordinates import coordinate_matrix
from check_rectangular_defect_core_basis import data_for
from search_rectangular_core_pressure import CoreDefectPressure


def profile_controls(data, compiled, seed, cases=3):
    matrix, shifts, metadata = compiled
    reduced = MatrixOracle(matrix, len(matrix), shifts, precision=metadata['profile_precision'])
    original, _, _ = original_cut(seed, data['place'], (0, 0, [0]*4, [0]*4, [0]*4))
    direct = MatrixOracle(original, len(original))
    rng = random.Random(7913+(data['place']=='zero'))
    records = []
    for case in range(cases):
        values = [rng.randrange(-8, 9)/4 for _ in range(18)] if case else [0]*18
        gauges = (values[0], values[1], values[2:6], values[6:10], values[10:14])
        if data.get('full_core', False):
            gauges += (values[14:18],)
        rw, cw = coordinate_weights(metadata, gauges)
        expected = reduced.profile(rw, cw)['orders']
        inactive = [i for i in range(2*data['M']) if i not in data['active_indices']]
        for index, order in zip(inactive, data['spectator_orders']):
            sign = 1 if data['place']=='infinity' or index<data['M'] else -1
            spectator = order+sign*(gauges[0]+gauges[1])
            expected = sorted(expected+[value+spectator for value in expected])
        _, rw, cw = original_cut(seed, data['place'], gauges)
        actual = direct.profile(rw, cw, infinity=data['place']=='infinity')['orders']
        if data['place']=='infinity':
            actual = [value+seed['denominator_degree'] for value in actual]
        assert len(actual)==len(expected)
        error = max(abs(a-b) for a, b in zip(sorted(actual), sorted(expected)))
        assert error<1e-9, (data['place'], case, error)
        records.append({'case':case, 'orders_checked':len(actual), 'maximum_error':error})
    return records


def pressure_control(full_core):
    M = 4
    seed = build(tuple(range(M)), M, gamma=2, check_types=False)
    direct = RectangularOracle(seed['values'], 16, 16, seed['denominator_degree'], False)
    reduced = CoreDefectPressure(M, full_core=full_core, full_alphabet_gamma=2)
    rng = np.random.default_rng(9871+full_core)
    theta = np.r_[rng.uniform(-.19, .19, reduced.parameters-1), .171]
    masks = [((letter % 4)<<1)|(1,0,8,9)[letter//4] for letter in range(16)]
    parameters = []
    for side in range(2):
        weights = ([M]*M, [-M]*M, [2*a-M+1 for a in range(M)], [2*a-M+1 for a in range(M)]) if side==0 else ([-2*M+1+2*a for a in range(M)], [1+2*a for a in range(M)], [0]*M, [0]*M)
        extras = (theta[1],theta[0],theta[1],theta[0]) if side==0 else (theta[2],-theta[2],-theta[3],theta[3])
        offset = 4+reduced.core_stride*side
        core = {}
        for letter in ('x','y','u','v'):
            if letter in reduced.leg_slices:
                start,count = reduced.leg_slices[letter]
                core[letter] = [0,*theta[offset+start:offset+start+count]]
            else:
                core[letter] = [0]*4
        for letter, ws, extra in zip(('x','y','u','v'),weights,extras):
            values = [sum(ws[a]*theta[-1]+extra for a in range(M) if mask>>a&1)
                      +core[letter][(mask&1)|((mask>>(M-2))&2)] for mask in masks]
            parameters.extend(value-values[0] for value in values[1:])
    parameters.append(theta[-1])
    rho, gradient, detail = reduced.evaluate(theta)
    reference = direct.types(parameters)[0]
    assert abs(rho-reference)<1e-10
    errors = []
    for j in range(reduced.parameters):
        a,b = theta.copy(),theta.copy()
        a[j]+=1e-6;b[j]-=1e-6
        errors.append(abs((reduced.evaluate(a)[0]-reduced.evaluate(b)[0])/2e-6-gradient[j]))
    assert max(errors)<1e-6
    return {'full_core':full_core,'parameters':theta.tolist(),'rho':rho,
            'direct_rho':reference,'pressure_error':abs(rho-reference),
            'derivative_coordinates':reduced.parameters,'maximum_derivative_error':max(errors),
            'entropy':detail['entropy']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    entries,profiles=[],[]
    for full_core,gamma in ((False,0),(False,2),(True,2)):
        seed=build(tuple(range(4)),4,gamma,check_types=False)
        for place in ('infinity','zero'):
            data=data_for(place,4,full_core=full_core)
            compiled=coordinate_matrix(data,full_alphabet_gamma=gamma)
            for z in (2,3):
                entry=point_control(data,compiled,z,seed)
                entries.append({'full_core':full_core,'gamma':gamma,**entry})
            profiles.append({'full_core':full_core,'gamma':gamma,'place':place,
                             'cases':profile_controls(data,compiled,seed)})
            print(json.dumps({'full_core':full_core,'gamma':gamma,'place':place,'status':'pass'}),flush=True)
    pressure=[pressure_control(full_core) for full_core in (False,True)]
    sources=['check_full_alphabet_quartic_repair.py','full_alphabet_quartic_profiles.py',
             'rectangular_defect_coordinates.py','search_rectangular_core_pressure.py',Path(__file__).name]
    result={'status':'pass; exact finite matrix/profile controls and floating Gibbs derivative controls; no pressure certificate',
            'entrywise_cases':entries,'weighted_profiles':profiles,'pressure_controls':pressure,
            'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources}}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
