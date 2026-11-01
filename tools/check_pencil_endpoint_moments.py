#!/usr/bin/env python3
"""Independent rational moments, complete orders and moment derivatives."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from check_joint_energy_conjugations import cut_matrix
from pencil_endpoint_moments import EndpointMoments
from search_f4_recentered_crossings import five_letter_source
from weighted_smith_oracle_f4 import WeightedSmithOracle


def variance(values):
    mean = sum(values)/len(values)
    return sum((v-mean)**2 for v in values)/len(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    model = EndpointMoments()
    source,_,_,_ = five_letter_source()
    letters = [0,2,3,4]
    rng,records,derivatives,orders = np.random.default_rng(160816),[],[],0
    for strong in (0,1):
        for trial in range(3):
            point = np.asarray(model.natural[12*strong:12*(strong+1)],dtype=float)
            if trial:
                point += rng.integers(-2,3,12)
            value,_,data = model.evaluate(point,strong)
            parameters = data['parameters']
            exact = []
            for place in range(2):
                h = [[0]+list(map(int,parameters[place*12+3*leg:place*12+3*(leg+1)])) for leg in range(4)]
                full = cut_matrix(source,place==0)
                ri = [x+5*u for u in letters for x in letters]
                ci = [y+5*v for v in letters for y in letters]
                if place:
                    ci = [v+5*u for u in letters for v in letters]
                matrix = [[full[i][j] for j in ci] for i in ri]
                row = [h[0][i%4]+h[2 if place==0 else 1][i//4] for i in range(16)]
                col = [h[1 if place==0 else 3][i%4]+h[3 if place==0 else 2][i//4] for i in range(16)]
                profile,_ = WeightedSmithOracle(matrix,infinity=place==0).profile(row,col)
                assert list(profile)==data['profiles'][place]
                exact.append(list(profile))
                orders += len(profile)
            r,s = [[Fraction(v) for v in p] for p in exact][strong],exact[1-strong]
            gx,gy = [[Fraction(v) for v in g] for g in data['physical_energies']]
            gap = sum(gx)/4+sum(gy)/4-sum(r)/16-min(s)
            delta = variance(r)-variance(gx)-variance(gy)
            assert gap==0 and sum(v==min(s) for v in s)==1 and float(delta)==value
            if trial==0:
                assert delta==1
            records.append({'strong_place':strong,'trial':trial,'data':data,
                            'exact_variance_difference':str(delta),'exact_first_order_gap':str(gap)})
        probe = np.asarray(model.natural[12*strong:12*(strong+1)])+rng.normal(size=12)*.03
        _,gradient,_ = model.evaluate(probe,strong)
        for j in range(12):
            shift = np.zeros(12)
            shift[j]=1e-6
            finite = (model.evaluate(probe+shift,strong)[0]-model.evaluate(probe-shift,strong)[0])/2e-6
            error = abs(float(finite-gradient[j]))
            assert error<1e-5,(strong,j,error)
            derivatives.append({'strong_place':strong,'coordinate':j,'error':error})
    root = Path(__file__).resolve().parents[1]
    names = ['tools/pencil_endpoint_moments.py','tools/check_pencil_endpoint_moments.py',
             'tools/pencil_principal_alphabets.py','tools/search_f4_recentered_crossings.py',
             'tools/native_rectangular_smith_f4.py','tools/native_rectangular_smith_f4.cpp',
             'tools/weighted_smith_oracle_f4.py','tools/check_joint_energy_conjugations.py']
    report = {'status':'pass','independent_complete_orders':orders,'rational_moment_checks':6,
              'finite_derivatives':len(derivatives),'controls':records,'derivatives':derivatives,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Checks the saturated simple-minimum locus and second moments/gradients. '
                      'The positive natural delta1 is not a global minimum or a factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','independent_complete_orders',
                                         'rational_moment_checks','finite_derivatives')}))


if __name__=='__main__':
    main()
