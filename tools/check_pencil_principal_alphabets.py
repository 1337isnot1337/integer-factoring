#!/usr/bin/env python3
"""Independent complete orders after unequal principal pencil restrictions."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_scaled_asymmetric_pressure import inputs
from native_rectangular_smith_f4 import RectangularOracle
from pencil_principal_alphabets import restrict
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import pack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest,cases = inputs(root)
    rng,records,orders,entries = np.random.default_rng(575026),[],0,0
    for index,ax,ay in ((0,[0,4],[0,1,2,3]),(33,[1,2,4],[0,2,3]),(66,[0,1,2,4],[1,3])):
        case,values = cases[index]
        seed,parameters,gx,gy = restrict(values,manifest['natural_heights_and_scale'],ax,ay)
        a,b = len(ax),len(ay)
        for v,u,y,x in product(range(a),range(b),range(b),range(a)):
            pair = seed['values'][x+a*(y+b*(u+b*v))]
            assert pair == values[ax[x]+5*(ay[y]+5*(ay[u]+5*ax[v]))]
            endpoint = pair[0].bit_count()%2 | (pair[1].bit_count()%2)<<1
            assert endpoint == int(x==v and y==u)
            entries += 1
        trial = np.asarray(parameters,dtype=int)
        trial[:-1] += rng.integers(-3,4,len(trial)-1)
        oracle = RectangularOracle(seed['values'],a,b)
        native = oracle.prefixes(trial,1)['orders']
        stride,offset = oracle.stride,oracle.offsets
        exact = []
        for place in range(2):
            h = [[0]+trial[place*stride+o:place*stride+o+q-1].tolist() for o,q in zip(offset,(a,b,b,a))]
            matrix = [[0]*(a*b) for _ in range(a*b)]
            for v,u,y,x in product(range(a),range(b),range(b),range(a)):
                r,c = (x+a*u,y+b*v) if place==0 else (x+a*y,v+a*u)
                matrix[r][c] = pack(seed['values'][x+a*(y+b*(u+b*v))])
            row = [int(h[0][i%a]+h[2 if place==0 else 1][i//a]) for i in range(a*b)]
            col = ([int(h[1][i%b]+h[3][i//b]) for i in range(a*b)] if place==0 else
                   [int(h[3][i%a]+h[2][i//a]) for i in range(a*b)])
            result,_ = WeightedSmithOracle(matrix,infinity=place==0).profile(row,col)
            exact += list(result)
        assert exact == native
        orders += len(exact)
        records.append({'case_index':index,'alphabet_x':ax,'alphabet_y':ay,
                        'parameters':trial.tolist(),'independent_orders':exact})
    names = ['tools/pencil_principal_alphabets.py','tools/check_pencil_principal_alphabets.py',
             'tools/check_scaled_asymmetric_pressure.py','tools/weighted_smith_oracle_f4.py',
             'tools/native_rectangular_smith_f4.py','tools/native_rectangular_smith_f4.cpp']
    report = {'status':'pass','selected_literal_entries_and_endpoints':entries,
              'independent_untruncated_orders':orders,'records':records,
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Exact unequal principal restrictions, endpoints and new-gauge full profiles. '
                      'No improving pressure or all-gauge exclusion.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','selected_literal_entries_and_endpoints',
                                         'independent_untruncated_orders')}))


if __name__ == '__main__':
    main()
