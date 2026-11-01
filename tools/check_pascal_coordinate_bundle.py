#!/usr/bin/env python3
"""Check complete-profile bundle gradients and held-out supporting planes."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from f4_pascal_flags import change_flag
from permuted_pascal_kernels import permuted_pascal_kernel
from rectangular_pascal_sectors import RectangularPascalPressure
from search_pascal_coordinate_bundle import log_profile_query


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(710710)
    seed = permuted_pascal_kernel([0,1,2,3])
    flag = [[int(i==j) for j in range(4)] for i in range(4)]
    flag[0][1] = 2
    model = RectangularPascalPressure(change_flag(seed,flag),[0,1,3],[0,1,2,3,7])
    points, records, derivatives = [], [], []
    for trial in range(12):
        theta = model.natural(.03)
        theta[:-1] += rng.normal(size=len(theta)-1)*.002
        rho,f,g,detail,physical = log_profile_query(model,theta)
        for j in (0,3,model.size,2*model.size):
            delta = np.zeros_like(theta)
            delta[j] = 1e-7
            finite = (log_profile_query(model,theta+delta)[1]-log_profile_query(model,theta-delta)[1])/2e-7
            error = abs(float(finite-g[j]))
            assert error < 2e-6, (trial,j,error)
            derivatives.append({'trial':trial,'coordinate':j,'error':error})
        assert abs((f+physical@theta)/(2*detail['entropy'])-rho) < 1e-12
        for oldx,oldf,oldg in points:
            assert f >= oldf+oldg@(theta-oldx)-1e-9
            assert oldf >= f+g@(oldx-theta)-1e-9
        points.append((theta,f,g))
        records.append({'parameters':theta.tolist(),'log_profile_product':f,'subgradient':g.tolist()})
    root = Path(__file__).resolve().parents[1]
    names = ['tools/search_pascal_coordinate_bundle.py','tools/check_pascal_coordinate_bundle.py',
             'tools/rectangular_pascal_sectors.py','tools/search_rectangular_core_bundle.py']
    report = {'status':'pass','nonbinary_query_points':12,'finite_derivatives':len(derivatives),
              'heldout_supporting_plane_inequalities':12*11,'records':records,'derivatives':derivatives,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Checks the log-profile subgradient adapter, supporting-plane inequalities '
                      'and fixed-type objective identity. No global optimum or stronger bound.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','nonbinary_query_points','finite_derivatives',
                                         'heldout_supporting_plane_inequalities')}))


if __name__ == '__main__':
    main()
