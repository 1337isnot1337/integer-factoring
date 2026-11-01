#!/usr/bin/env python3
"""Exact rational-type improvement for the existing homogeneous cube kernel."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
from certify_natural_pressure_interval import log_interval


def certificate(bound, bits=192, terms=96):
    p=F(208173,10**6);r=F(434044514,10**9)
    Z=(1+r)*(1+r**2)*(1+r**4)*(1+r**5)
    lp=log_interval(p,bits,terms);lq=log_interval(1-p,bits,terms)
    lr=log_interval(r,bits,terms);lz=log_interval(Z,bits,terms)
    H=(-p*lp[1]-(1-p)*lq[1],-p*lp[0]-(1-p)*lq[0])
    N=(lz[0]-4*p*lr[1],lz[1]-4*p*lr[0])
    gap=(4*bound*H[0]-N[1],4*bound*H[1]-N[0])
    return {'bound':[bound.numerator,bound.denominator], 'physical_type':[p.numerator,p.denominator],
        'rational_r':[r.numerator,r.denominator], 'bits':bits,'terms':terms,
        'gap_interval_numerators':[[v.numerator,v.denominator] for v in gap],
        'gap_scale':2**bits,'passed':gap[0]>0,'rejected_at_this_point':gap[1]<0}


def main():
    target=F(624718442,10**9)
    a=certificate(target);b=certificate(target,224,112);negative=certificate(F(624718441,10**9))
    assert a['passed'] and b['passed'] and negative['rejected_at_this_point']
    # Literal determinant profile of the already constructed degree-six block.
    from cube_homogeneous_gaussian import augmentation_crossblock, exact_monomial_smith, kernel
    profiles=[exact_monomial_smith(augmentation_crossblock(kernel(),axis),axis)['orders'] for axis in range(3)]
    assert profiles==[[2,4,8,10]]*3,profiles
    names=['certify_cube_optimized_type.py','certify_natural_pressure_interval.py','cube_homogeneous_gaussian.py']
    report={'status':'passed','certificate':a,'higher_precision':b,'negative_control':negative,
        'exact_three_one_particle_profiles':profiles,
        'source_sha256':{n:sha256(Path('tools',n).read_bytes()).hexdigest() for n in names},
        'scope':'Stronger cube-interface upper bound for the existing family; not a factoring-exponent improvement or a global optimum certificate.'}
    Path('evidence/cube-optimized-type.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','strict_cube_rho':[target.numerator,target.denominator],'factoring_exponent_improved':False}))

if __name__=='__main__':main()
