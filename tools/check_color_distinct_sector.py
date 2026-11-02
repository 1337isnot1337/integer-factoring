#!/usr/bin/env python3
"""Controls and exact pressure certificate for the shared color-parity theorem."""
import argparse
from hashlib import sha256
from fractions import Fraction
from itertools import permutations,product
import json
from math import comb
from pathlib import Path
import cube_homogeneous_gaussian as cube


def colored_entry(left,right,axis,gaussian):
    first,second=cube.axis_parts(axis);masks=[0]*6
    for labels,vertices in ((left,first),(right,second)):
        for c,v in zip(labels,vertices):masks[c]|=1<<v
    value=cube.ONE
    for mask in masks[1:]:
        value=cube.multiply(value,gaussian[mask])
        if not value:break
    return value


def structural_control():
    kernel=cube.kernel();gaussian=cube.gaussian_tensor(kernel)
    assert gaussian[0]==cube.ONE
    assert all(not g for mask,g in enumerate(gaussian) if mask.bit_count()%2)
    assert all(sum(term)==3*mask.bit_count() for mask,g in enumerate(gaussian) for term in g)
    assert kernel[0][7]==cube.monomial(2,2,2)
    words=list(permutations((1,2,3,4)));records=[]
    for axis in range(3):
        first,second=cube.axis_parts(axis)
        assert all(cube.valuation(kernel[a][b],axis)==2 for a in first for b in second)
        endpoint_checks=0
        for left in words:
            for right in words:
                actual=colored_entry(left,right,axis,gaussian)
                expected=cube.ONE
                for i,color in enumerate(left):
                    expected=cube.multiply(expected,kernel[first[i]][second[right.index(color)]])
                assert actual==expected and cube.valuation(actual,axis)==8
                endpoint=all(right[second.index(v^7)]==left[i] for i,v in enumerate(first))
                assert len(actual)%2==endpoint;endpoint_checks+=1
        left=words[0]
        antipodal=tuple(left[first.index(v^7)] for v in second)
        assert colored_entry(left,antipodal,axis,gaussian)==cube.monomial(8,8,8)
        closed_checks=0
        for right in product(range(6),repeat=4):
            if sorted(right)!=[1,2,3,4]:
                assert not colored_entry(left,right,axis,gaussian)
                closed_checks+=1
        minor_checks=0
        for c,d in product(range(6),repeat=2):
            left=[0]*4;right=[0]*4;left[first.index(0)]=c;right[second.index(7)]=d
            actual=colored_entry(left,right,axis,gaussian)
            expected=cube.ONE if c==d==0 else kernel[0][7] if c==d else cube.ZERO
            assert actual==expected;minor_checks+=1
        records.append(dict(axis=axis,sector_order=24,sector_endpoint_checks=endpoint_checks,closed_sector_zero_checks=closed_checks,complement_minor_checks=minor_checks,smallest_sector_order=8))
    return records


def asymmetric_control():
    # The general theorem bounds the geometric mean; per-axis bounds need not hold.
    q=5;d=q-1;t=Fraction(1,2)
    orders=[6*sum(c!=0 for c in word) for word in product(range(q),repeat=4)]
    zx=sum(t**s for s in orders);zy=zz=q**4
    lower=sum(comb(d,r)*t**(2*r) for r in range(5))
    assert zx==(1+d*t**6)**4 and zx<lower and zx*zy*zz>=lower**3
    return dict(q=q,side_words_checked=len(orders),t=[1,2],single_axis_bound_fails=True,geometric_mean_bound_passes=True)


def main():
    controls=structural_control();asymmetric=asymmetric_control()
    certificate=dict(source='tools/check_cube_scalar_color_clone_floor.py',source_sha256=sha256(Path(__file__).with_name('check_cube_scalar_color_clone_floor.py').read_bytes()).hexdigest(),replayed_by='make verify-color-distinct-sector',newton_coefficients=50,minimum_q=34,rho_floor=[3,4],scope='Existing polynomial certificate reused by the ordinary parity-sector theorem; run the Make target for both certificate and structural controls.')
    result=dict(scope='Gaussian direct-sector controls and an asymmetric parity-preserving control supporting the ordinary parity-sector theorem, which reuses the separately replayed scalar polynomial certificate. No exhaustive tensor verification, specialized rank extrapolation, large-q matrix screen, stronger factoring exponent, or independent mathematical review.',original_color_checkpoint='1772970a7dbaab83eeb0607014a1111e46718d86',axes=controls,certificate=certificate,asymmetric_control=asymmetric,kernel_source_sha256=sha256(Path(cube.__file__).read_bytes()).hexdigest())
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(axes=len(controls),sector_endpoint_checks=sum(r['sector_endpoint_checks'] for r in controls),closed_sector_zero_checks=sum(r['closed_sector_zero_checks'] for r in controls),minor_checks=sum(r['complement_minor_checks'] for r in controls),reused_certificate_floor=[3,4],minimum_q=34)))


if __name__=='__main__':main()
