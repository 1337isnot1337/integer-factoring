#!/usr/bin/env python3
"""Independent divisor-sieve and streamed-arithmetic checks for the exact oracle."""
from dataclasses import asdict
import argparse
import hashlib
import json
from math import prod
from pathlib import Path
from time import monotonic

from exact_divisor_maximum import exact_divisor_maximum,first_primes,streamed_product_bit


def split_bound(n):
    low=n.bit_length()-1;high=(n-1).bit_length()
    return min([n-1,(n+1)//2+1]+[(1<<(s-1))*high+(n-1)//s for s in range(1,low+1)])


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--sieve-bits",type=int,default=18)
    parser.add_argument("--large",action="store_true");parser.add_argument("--output");args=parser.parse_args()
    if not 1<=args.sieve_bits<=22:parser.error("sieve bits must lie between1 and22")
    started=monotonic();limit=(1<<args.sieve_bits)-1
    counts=[0]*(limit+1)
    for divisor in range(1,limit+1):
        for multiple in range(divisor,limit+1,divisor):counts[multiple]+=1
    small=[];maximum=0;witness=0;next_bits=1
    for value in range(1,limit+1):
        if counts[value]>maximum:maximum=counts[value];witness=value
        if value==(1<<next_bits)-1:
            actual=exact_divisor_maximum(next_bits)
            assert (actual.maximum_divisors,actual.least_witness)==(maximum,witness)
            small.append(asdict(actual));next_bits+=1
    primes=first_primes(64)
    trial=[p for p in range(2,primes[-1]+1) if all(p%d for d in range(2,p))]
    assert primes==trial
    for k in range(1,65):assert primes[k-1]<=4*k*k
    # This independently checks the streamed-column carry recurrence.
    bit_cases=[[],[1],[2],[31],[2,3],[3,5,7],[2,3,5,7],
               [7,7,7,7],[1,2,1,3,1,2],[3,3,3,3,3,3]]
    bit_checks=0
    for factors in bit_cases:
        bound=max(1,max((x.bit_length() for x in factors),default=1))
        expected=prod(factors)
        for bit in range(max(1,len(factors)*bound)+2):
            observed=streamed_product_bit(len(factors),factors.__getitem__,bound,bit)
            assert observed==(expected>>bit)&1
            bit_checks+=1
    invalid=0
    for n in (0,-1,False,1.5):
        try:exact_divisor_maximum(n)
        except ValueError:invalid+=1
        else:raise AssertionError("invalid n accepted")
    table=[]
    for n in ([8,16,32,64,128] if args.large else [8,16,32,64]):
        actual=exact_divisor_maximum(n);record=asdict(actual)
        factors=first_primes(len(actual.exponents))
        assert actual.exponents==sorted(actual.exponents,reverse=True)
        assert prod(p**e for p,e in zip(factors,actual.exponents))==actual.least_witness<1<<n
        assert prod(e+1 for e in actual.exponents)==actual.maximum_divisors
        assert actual.optimal_count_bits<=split_bound(n)
        assert actual.maximum_divisors<=1<<actual.optimal_count_bits
        if actual.optimal_count_bits:assert actual.maximum_divisors>1<<(actual.optimal_count_bits-1)
        record["previous_split_count_bits"]=split_bound(n)
        record["least_strict_power_two_modulus_bits"]=actual.maximum_divisors.bit_length()
        table.append(record)
    result={"independently_sieved_integers":limit,"exact_bit_ranges_compared":len(small),
            "prime_bound_cases":64,"streamed_product_bits_compared":bit_checks,
            "invalid_inputs_rejected":invalid,"small_oracles":small,"comparison_table":table,
            "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(__file__).with_name("exact_divisor_maximum.py"))},
            "seconds":round(monotonic()-started,3)}
    text=json.dumps(result,indent=2)+"\n"
    if args.output:Path(args.output).write_text(text)
    print(json.dumps({k:v for k,v in result.items() if k!="small_oracles"},indent=2))


if __name__=="__main__":main()
