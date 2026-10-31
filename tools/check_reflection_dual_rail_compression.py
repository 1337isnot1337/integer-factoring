#!/usr/bin/env python3
"""Independent minors and complete profiles for reflection compression."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import numpy as np
from reflection_dual_rail_gaussian import tensor,reflection_letters,gauged_kernel,logical_heights,profile_factors,physical_bit_energies
from selected_gaussian_coefficients import selected_tensor
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from prime_state_filter_pressure import gaussian_pressure
from search_reflection_gaussian_gauges import GaussianGaugePressure


def subset_sums(values):
    result=[0]
    for value in values:result += [x+value for x in result]
    return sorted(result)


def physical_parameters(interval,theta):
    m=interval//2;q=1<<m;size=4*m;result=[]
    for side in (1,0):
        for leg in range(4):
            weights=theta[side*size+leg*m:side*size+(leg+1)*m]
            result.extend(sum(weights[i] for i in range(m) if label>>i&1) for label in range(1,q))
    return result+[theta[-1]]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    coefficient_checks=endpoint_checks=0;physical_controls=[];kernel_controls=[];rate_controls=0;max_rate_error=0.;rng=random.Random(310753)
    for M in (4,8):
        actual=tensor(M);expected=selected_tensor(tuple(range(M)),M,reflection_letters(M))
        assert actual['values']==expected['values']
        coefficient_checks+=len(actual['values']);endpoint_checks+=len(actual['values'])
        q=actual['q_x'];oracle=RectangularOracle(actual['values'],q,q,endpoint_check=True)
        pressure=GaussianGaugePressure(M)
        natural=pressure.natural(1.)
        expected_orders=subset_sums(profile_factors(M));scalar=M*M//2
        _,_,orders=oracle.evaluate(physical_parameters(M,natural),[1/q]*q,[1/q]*q,1,True)
        assert orders[:q*q]==[x-scalar for x in expected_orders]
        assert orders[q*q:]==[x+scalar for x in expected_orders]
        physical_controls.append({'interval':M,'kind':'natural complete integer profiles','orders_per_cut':q*q})
        for trial in range(3):
            theta=pressure.natural(.13+trial*.07)
            theta[:-1]+=np.array([.04+rng.randrange(1,25)/1000 for _ in theta[:-1]])
            assert min(pressure.constraints@theta)>0
            small=pressure.evaluate(theta);full=oracle.types(physical_parameters(M,theta),1)[0]
            error=abs(small-full);assert error<1e-11,(M,trial,small,full)
            physical_controls.append({'interval':M,'kind':'independent physical mode gauges','trial':trial,'absolute_rate_error':error})
    # The native mask interface also covers its 32-bit endpoint exactly.
    # Both extreme isotropic minors equal det(A)=w^(M^2/2).
    for letter in (255,65280):
        extreme=selected_tensor(tuple(range(16)),16,(letter,))
        assert extreme['values']==[(1<<128,0)]
    for M in (4,8,16,32,64,128):
        m=M//2
        for place in ('zero','infinity'):
            kernel=gauged_kernel(M,place)
            rows=[leg*m+i for leg in ((0,1) if place=='zero' else (0,2)) for i in range(m)]
            columns=[i for i in range(4*m) if i not in rows]
            matrix=[[(kernel[i][j],0) for j in columns] for i in rows]
            precision=max(512,1<<(4*M-1).bit_length())
            orders=MatrixOracle(matrix,M,precision=precision).profile()['orders']
            assert orders==list(profile_factors(M)),(M,place,orders)
            kernel_controls.append({'interval':M,'place':place,'one_particle_rank':M,'orders':orders})
        for ell in (.007,.03,.13,.6):
            for rho in (.6,1/math.sqrt(2),.9):
                R=2*rho
                compressed=2*sum(np.logaddexp(0,-ell*x) for x in profile_factors(M))
                compressed-=2*R*sum(np.logaddexp(0,-ell*x/R) for x in physical_bit_energies(M))
                full=gaussian_pressure(2,M,0,2*ell,rho)/2
                error=abs(float(compressed)-full);max_rate_error=max(max_rate_error,error)
                assert error<1e-10,(M,ell,rho,compressed,full)
                rate_controls+=1
    result={'status':'pass','independent_original_minor_coefficients':coefficient_checks,'strict_endpoint_coefficients':endpoint_checks,'sixteen_mode_boundary_mask_controls':2,
            'complete_physical_controls':physical_controls,'uniform_formula_profile_controls':kernel_controls,
            'rate_identity_controls':rate_controls,'rate_identity_max_error':max_rate_error,
            'source_sha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
                ('reflection_dual_rail_gaussian.py','selected_gaussian_coefficients.py','selected_gaussian_coefficients.cpp',Path(__file__).name)},
            'scope':'Uniform polynomial compression and exact profiles; no stronger factoring exponent.'}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('source_sha256','uniform_formula_profile_controls')},indent=2))

if __name__=='__main__':main()
