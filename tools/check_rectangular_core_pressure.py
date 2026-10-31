#!/usr/bin/env python3
"""Independently compare reduced Gibbs pressure with the full M4 tensor."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from check_rectangular_f4_repair import build
from native_rectangular_smith_f4 import RectangularOracle
from search_rectangular_core_pressure import CoreDefectPressure

def check(full_core=False):
    M=4;seed=build(tuple(range(M)),M,check_types=False)
    full=RectangularOracle(seed['values'],seed['q_x'],seed['q_u'],seed['denominator_degree'],endpoint_check=False)
    reduced=CoreDefectPressure(M,full_core=full_core);bulk=1<<(M-2);top=1<<(M-1)
    stride=10 if full_core else 7;dimension=5+2*stride
    mx=[((a%bulk)<<1)|(1,0,top,top|1)[a//bulk] for a in range(full.A)]
    mu=[((a%bulk)<<1)|(1,0,top)[a//bulk] for a in range(full.B)]
    rng=np.random.default_rng(9108 if full_core else 427);records=[]
    for case in range(1 if full_core else 2):
        theta=np.r_[rng.uniform(-.19,.19,dimension-1),.13+.04*case]
        parameters=[]
        for side in range(2):
            weights=([M]*M,[-M]*M,[2*a-M+1 for a in range(M)],[2*a-M+1 for a in range(M)]) if side==0 else ([-2*M+1+2*a for a in range(M)],[1+2*a for a in range(M)],[0]*M,[0]*M)
            extras=(theta[1],theta[0],theta[1],theta[0]) if side==0 else (theta[2],-theta[2],-theta[3],theta[3])
            o=4+stride*side;core=([0,*theta[o+7:o+10]] if full_core else [0]*4,[0,*theta[o+3:o+5],0],[0,*theta[o+5:o+7],0],[0,*theta[o:o+3]])
            for masks,ws,extra,cs in zip((mx,mu,mu,mx),weights,extras,core):
                values=[sum(ws[a]*theta[-1]+extra for a in range(M) if mask>>a&1)+cs[(mask&1)|((mask>>(M-2))&2)] for mask in masks]
                parameters.extend(value-values[0] for value in values[1:])
        parameters.append(theta[-1]);rho,gradient,detail=reduced.evaluate(theta);reference=full.types(parameters)[0]
        assert abs(rho-reference)<1e-10
        errors=[]
        for j in range(dimension):
            a=theta.copy();b=theta.copy();a[j]+=1e-6;b[j]-=1e-6
            errors.append(abs((reduced.evaluate(a)[0]-reduced.evaluate(b)[0])/2e-6-gradient[j]))
        assert max(errors)<1e-6
        records.append({'parameters':theta.tolist(),'reduced_rho':rho,'full_rho':reference,'pressure_error':abs(rho-reference),'derivative_coordinates':dimension,'maximum_derivative_error':max(errors),'entropy':detail['entropy']})
    sources=['rectangular_defect_profiles.py','rectangular_defect_coordinates.py','native_fock_f4.cpp','search_rectangular_defect_pressure.py','search_rectangular_core_pressure.py',Path(__file__).name]
    scope=f'{len(records)} generic M4 points; independent full-tensor Gibbs roots and all {len(records)*dimension} core/charge derivatives. Floating controls do not certify a new bound or a global optimum.'
    return {'status':'pass','scope':scope,'full_core':full_core,'cases':records,'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources}}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);parser.add_argument('--full-core',action='store_true');args=parser.parse_args()
    try:result=check(args.full_core)
    except Exception as error:
        if args.output:args.output.write_text(json.dumps({'status':'fail','error':repr(error)},indent=2)+'\n')
        raise
    rendered=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered)

if __name__=='__main__':main()
