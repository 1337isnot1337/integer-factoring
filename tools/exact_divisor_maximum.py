#!/usr/bin/env python3
"""Exact maximum of tau(N) for 1 <= N < 2**n by canonical exponent patterns.

This fast implementation uses ordinary integer arithmetic and polynomial
workspace.  The companion proof supplies a lower-space streamed variant.
"""
from dataclasses import asdict,dataclass
import argparse
import json
from time import monotonic


def first_primes(count):
    if count<0:
        raise ValueError("prime count must be nonnegative")
    if not count:
        return []
    bound=4*count*count
    sieve=bytearray(b"\x01")*(bound+1)
    sieve[0:2]=b"\x00\x00"
    p=2
    while p*p<=bound:
        if sieve[p]:
            sieve[p*p:bound+1:p]=b"\x00"*((bound-p*p)//p+1)
        p+=1
    primes=[i for i,flag in enumerate(sieve) if flag]
    assert len(primes)>=count
    return primes[:count]


@dataclass
class DivisorMaximum:
    input_bits:int
    maximum_divisors:int
    optimal_count_bits:int
    least_witness:int
    exponents:list[int]
    canonical_patterns_checked:int


def exact_divisor_maximum(n):
    if not isinstance(n,int) or isinstance(n,bool) or n<1:
        raise ValueError("input bit bound must be a positive integer")
    primes=first_primes(n-1)
    limit=1<<n
    best_count=1;best_value=1;best_exponents=[];visited=0
    exponents=[]
    def visit(index,maximum_exponent,total_exponents,value,count):
        nonlocal best_count,best_value,best_exponents,visited
        visited+=1
        if count>best_count or (count==best_count and value<best_value):
            best_count=count;best_value=value;best_exponents=exponents.copy()
        if index>=len(primes):
            return
        power_value=value
        for exponent in range(1,min(maximum_exponent,n-1-total_exponents)+1):
            power_value*=primes[index]
            if power_value>=limit:
                break
            exponents.append(exponent)
            visit(index+1,exponent,total_exponents+exponent,power_value,count*(exponent+1))
            exponents.pop()
    visit(0,n-1,0,1,1)
    return DivisorMaximum(n,best_count,(best_count-1).bit_length(),best_value,best_exponents,visited)


def streamed_product_bit(factor_count,factor_at,leaf_bits,bit):
    """One product bit, with balanced recursion and no stored partial product.

    factor_at(i) must regenerate a nonnegative factor smaller than
    2**leaf_bits.  No memoization or factor array is used by this routine.
    The low-space theorem uses a trial-prime generator for that callback.
    """
    if factor_count<0 or leaf_bits<1 or bit<0:
        raise ValueError("invalid streamed-product parameters")
    def get(lo,hi,k):
        count=hi-lo
        if not count:
            return int(k==0)
        if k>=count*leaf_bits:
            return 0
        if count==1:
            value=factor_at(lo)
            assert value>=0 and value.bit_length()<=leaf_bits
            return (value>>k)&1
        mid=(lo+hi)//2
        left_bits=(mid-lo)*leaf_bits
        right_bits=(hi-mid)*leaf_bits
        carry=0
        digit=0
        for column in range(k+1):
            total=carry
            for i in range(max(0,column-right_bits+1),min(column,left_bits-1)+1):
                a=get(lo,mid,i)
                if a:
                    total+=get(mid,hi,column-i)
            digit=total&1
            carry=total>>1
        return digit
    return get(0,factor_count,bit)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("bits",type=int,nargs="+",default=[32])
    args=parser.parse_args();started=monotonic()
    results=[asdict(exact_divisor_maximum(n)) for n in args.bits]
    print(json.dumps({"results":results,"seconds":round(monotonic()-started,3)},indent=2))


if __name__=="__main__":main()
