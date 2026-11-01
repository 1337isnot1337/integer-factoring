#!/usr/bin/env python3
"""Free-scale chain rule, balanced normalization and exact pencil profiles."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from check_joint_energy_conjugations import cut_matrix
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure
from two_place_f4_census import build_input
from weighted_smith_oracle_f4 import WeightedSmithOracle


def inputs(root):
    manifest = json.loads((root/'evidence/two-place-free-gauge-asymmetric-manifest.json').read_text())
    text,_ = build_input()
    fields = text.decode().split()
    entries = [tuple(map(int,fields[1+4*i:1+4*(i+1)])) for i in range(625)]
    result = []
    for case in manifest['cases']:
        low,high = case['q_planes']
        values = [(a^(low if k else 0)^(low<<4 if l else 0),
                   b^(high if k else 0)^(high<<4 if l else 0)) for a,b,k,l in entries]
        assert sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest() == case['values_sha256']
        result.append((case,values))
    return manifest,result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest,cases = inputs(root)
    rng,records,orders_checked,derivatives = np.random.default_rng(5051091),[],0,0
    for position in (0,33,66):
        case,values = cases[position]
        model = ScaledAsymmetricPressure({'q_x':5,'q_u':5,'values':values})
        natural = np.array(manifest['natural_heights_and_scale'],float)
        raw = model.base.oracle.prefixes(natural,1)['orders']
        assert raw == case['natural_orders']
        exact_parameters = natural.copy()
        exact_parameters[:-1] += rng.integers(-2,3,len(natural)-1)
        exact = []
        for side in range(2):
            h = [[0]+list(map(int,exact_parameters[side*16+4*leg:side*16+4*(leg+1)])) for leg in range(4)]
            legs = (0,2,1,3) if side == 0 else (0,1,2,3)
            rw = [h[legs[0]][i%5]+h[legs[1]][i//5] for i in range(25)]
            column_legs = (1,3) if side == 0 else (3,2)
            cw = [h[column_legs[0]][i%5]+h[column_legs[1]][i//5] for i in range(25)]
            profile,_ = WeightedSmithOracle(cut_matrix(values,side==0),infinity=side==0).profile(rw,cw)
            exact += list(profile)
        assert exact == model.base.oracle.prefixes(exact_parameters,1)['orders']
        orders_checked += 50
        ordinary = model.base.oracle.types(natural*.21,1)[0]
        balanced = np.r_[natural[:-1],math.log(.21/math.sqrt(2)),0.]
        rate,_ = model.evaluate(balanced,gradient=False)
        assert abs(rate-math.sqrt(2)*ordinary)<1e-9
        controls = []
        for log_tan in (-2.,.7):
            point = np.r_[natural[:-1]+rng.normal(size=32)*.1,math.log(.21),log_tan]
            value,gradient,data = model.evaluate(point)
            errors = []
            for coordinate in range(len(point)):
                step = 1e-6
                delta = np.zeros_like(point)
                delta[coordinate] = step
                finite = (model.evaluate(point+delta,False)[0]-model.evaluate(point-delta,False)[0])/(2*step)
                errors.append(abs(float(finite-gradient[coordinate])))
                derivatives += 1
            assert max(errors)<1e-5
            controls.append({'point':point.tolist(),'value':value,'details':data,'gradient_errors':errors})
        records.append({'case':case,'exact_parameters':exact_parameters.tolist(),
                        'independent_complete_orders':exact,'balanced_support':rate,
                        'balanced_expected':math.sqrt(2)*ordinary,'gradient_controls':controls})
    names = ['tools/scaled_asymmetric_crossing_pressure.py','tools/check_scaled_asymmetric_pressure.py',
             'tools/asymmetric_crossing_pressure.py','tools/native_rectangular_smith_f4.py',
             'tools/native_rectangular_smith_f4.cpp','tools/two_place_f4_census.py',
             'tools/weighted_smith_oracle_f4.py','tools/check_joint_energy_conjugations.py',
             'evidence/two-place-free-gauge-asymmetric-manifest.json']
    result = {'status':'pass','pencil_cases_reconstructed':len(cases),
              'independent_complete_orders':orders_checked,'finite_derivative_comparisons':derivatives,
              'controls':records,'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Coordinate/gradient and input controls for a numerical search, '
              'not an exact improving pressure certificate or a global exclusion.'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','pencil_cases_reconstructed',
        'independent_complete_orders','finite_derivative_comparisons')}))


if __name__ == '__main__':
    main()
