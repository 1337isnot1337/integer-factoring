#!/usr/bin/env python3
"""Exact all-permutation order controls and rational full-pressure comparisons."""
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import argparse,json,random
from pathlib import Path
from permuted_pascal_kernels import permutation_upset_orders


def exponential_residual(a, e, b, numerator, denominator, base):
    length = len(a)
    zero = infinity = physical = Fraction(1)
    for aa, ee, bb in zip(a, e, b):
        zero *= (1 + base ** (-numerator * aa)) ** 2
        infinity *= (1 + base ** (numerator * (2 * length - min(bb, ee))))
        infinity *= (1 + base ** (-numerator * max(bb, ee)))
        physical *= (1 + base ** (-denominator * (aa + bb - 2 * length)))
        physical *= (1 + base ** (-denominator * (aa + ee - 2 * length)))
    return (zero * infinity) ** denominator / physical ** numerator


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    basis_cases=0
    for image in permutations(range(8)):
        energy=sorted(15-2*v for v in permutation_upset_orders(image))
        assert all(value>=2*j+1 for j,value in enumerate(energy))
        basis_cases+=1
    affine_checks=0
    for length in (2,4,8):
        for aa in range(1,2*length,2):
            for ee in range(1,2*length,2):
                for numerator,denominator in ((3,2),(5,4),(7,5)):
                    assert denominator*(aa+ee-2*length)<numerator*aa
                    affine_checks+=1
                    for bb in range(-2*length,4*length+1):
                        if bb<=ee:assert denominator*(aa+bb-2*length)>numerator*(bb-2*length)
                        if bb>=ee:assert denominator*(aa+bb-2*length)<numerator*bb
                        affine_checks+=1
    rng=random.Random(202611010604);samples=list(permutations(range(4)))
    for _ in range(12):
        image=list(range(8));rng.shuffle(image);samples.append(tuple(image))
    comparisons=0;strict=0
    for image in samples:
        length=len(image);a=[2*length-1-2*v for v in permutation_upset_orders(image)];e=[2*v+1 for v in image];odd=list(range(1,2*length,2));order=sorted(range(length),key=lambda i:(a[i],e[i]));sorted_a=[a[i] for i in order];reordered_e=[e[i] for i in order]
        for numerator,denominator in ((3,2),(7,5)):
            base=Fraction(2*length+1,2*length)
            relaxed=exponential_residual(a,e,e,numerator,denominator,base)
            reduced=exponential_residual(odd,reordered_e,reordered_e,numerator,denominator,base)
            reference=exponential_residual(odd,odd,odd,numerator,denominator,base)
            assert relaxed>=reduced>=reference
            for beta in (odd,[2*i-length for i in range(length)],[4*i+1 for i in range(length)]):
                actual=exponential_residual(a,e,beta,numerator,denominator,base)
                assert actual>=relaxed
                comparisons+=1;strict+=int(actual>reference)
    root=Path(__file__).resolve().parents[1];sources=[Path(__file__),Path(__file__).with_name('permuted_pascal_kernels.py'),root/'research/permuted-pascal-cone-dominance.md']
    result={'status':'pass','all_length_eight_permutations':basis_cases,'exact_affine_derivative_sign_controls':affine_checks,'exact_rational_full_pressure_comparisons':comparisons,'strict_comparisons_to_ordered_kernel':strict,'scope':'All-scale comparison in the stated common-scale explicit cone. No assertion about arbitrary physical gauges or restricted alphabets.','source_sha256':{str(x.relative_to(root)):sha256(x.read_bytes()).hexdigest() for x in sources}}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))
if __name__=='__main__':main()
