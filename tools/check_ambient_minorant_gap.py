#!/usr/bin/env python3
"""Exact controls for the ambient positive-minorant norm-gap obstruction.

Finite compositions do not certify the asymptotic analytic bias estimate.
"""
import argparse
from fractions import Fraction
import json
from math import gcd,isqrt
from pathlib import Path
from binary_quadratic_forms import compose_forms,reduce_form
from prime_form_bias import prime_form_sample
from check_prime_translate_obstruction import recognition_failure


def is_prime(n):
    return n>1 and (n==2 or (n%2==1 and all(n%d for d in range(3,isqrt(n)+1,2))))


def walsh_control():
    n=16
    character=lambda x,c: (-1)**((x&c).bit_count()%2)
    mu=[Fraction(0)]+[Fraction(1,n-1)]*(n-1)
    bias=max(abs(sum(mu[x]*character(x,c) for x in range(n))) for c in range(1,n))
    records=[]
    for mix in (Fraction(1),Fraction(1,2)):
        # Convex combination of delta_0 and the zero-mean PD simplex kernel.
        f=[Fraction(1)]+[-(1-mix)/Fraction(n-1)]*(n-1)
        coefficients=[sum(f[x]*character(x,c) for x in range(n))/n for c in range(n)]
        assert all(v>=0 for v in coefficients) and sum(coefficients)==1
        mean=coefficients[0];expectation=sum(mu[x]*f[x] for x in range(n))
        assert all(f[x]<=0 for x in range(1,n))
        assert expectation<=0 and expectation>=mean-bias*(1-mean)
        assert mean<=bias/(1+bias)
        records.append(dict(mixture=str(mix),mean=str(mean),expectation=str(expectation),sharp=mean==bias/(1+bias)))
    assert records[0]['sharp']
    return dict(group='(Z/2)^4',bias=str(bias),sharp_mean=str(bias/(1+bias)),records=records)


def arithmetic_control(A,p,q):
    assert all(is_prime(v) for v in (A,p,q))
    assert p%(4*A*A)==1 and q%(4*A*A*p*p)==4*A*A*p*p-1
    assert p>=4*A*A+1 and q>=4*A*A*p*p-1 and q>p
    N=p*q;b=N.bit_length();C=(N+1)//(4*A);D=-N;H=isqrt(N//3)
    assert 4*A*C==N+1 and N%4==3 and gcd(A,N)==1
    assert p>(1<<(isqrt(b-1)+1))  # Actual stated small-prime preprocessing.
    assert N>4*A**4 and H>A*A and C>A**3
    target=(A,1,C);assert reduce_form(target)==target
    masses,records=prime_form_sample(D,A,12)
    assert masses
    controls=[]
    for source_prime,_,form,_ in records:
        translated=compose_forms(target,form)
        leading=translated[0]
        assert source_prime<A and leading==A*source_prime
        assert leading<=H and gcd(leading,N)==1
        assert recognition_failure(translated,D,H,A-1)=='leading_coefficient_not_smooth'
        controls.append(dict(prime=source_prime,source_form=form,translated_form=translated))
    return dict(A=A,factors=[p,q],N=N,bits=b,discriminant=D,target=target,full_height=H,smooth_bound=A-1,small_prime_preprocessing_cutoff=1<<(isqrt(b-1)+1),oriented_translates=len(controls),controls=controls)


def main():
    # Witnesses are checked by elementary trial primality, not trusted labels.
    witnesses=[(11,1453,2043650311),(13,677,4647456059),
               (17,3469,194757624823),(19,18773,6106829566511),
               (23,12697,682256839687),(29,10093,19533109708451)]
    controls=[arithmetic_control(*row) for row in witnesses]
    # Coarse Linnik exponent six: p=O(A^12), q=O(A^156).
    p_exponent=2*6;q_exponent=(2+2*p_exponent)*6;N_exponent=p_exponent+q_exponent
    assert (p_exponent,q_exponent,N_exponent)==(12,156,168)
    denominator=32*N_exponent;assert denominator==5376
    result=dict(scope='Exact finite controls for the ordinary ambient-minorant obstruction. Uniform prime bias and the infinite-family existence use the cited analytic estimates and least-prime theorem; they are not inferred from these examples.',walsh_control=walsh_control(),semiprime_controls=controls,oriented_translates_checked=sum(c['oriented_translates'] for c in controls),coarse_exponents=dict(p=p_exponent,q=q_exponent,N=N_exponent,minorant_binary_decay_denominator=denominator))
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(semiprimes=len(controls),oriented_translates_checked=result['oriented_translates_checked'],walsh_sharp_mean=result['walsh_control']['sharp_mean'],binary_decay_denominator=denominator)))


if __name__=='__main__':main()
