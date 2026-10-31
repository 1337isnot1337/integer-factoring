#!/usr/bin/env python3
"""Check exact coefficient flags and complete finite frozen-row searches."""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from pathlib import Path
import random
import time
from frozen_core_row_profiles import compile_interval,combine_components
from search_frozen_core_complete_pressure import entry_flag,dot,screen,ProductComparison


def direct_minimum(matrix,shifts,metadata,gauge,M):
    values=[]
    for i,row in enumerate(matrix):
        for j,entry in enumerate(row):
            bits=entry[0]|entry[1]
            if bits:
                values.append(Fraction(shifts[i][j])+gauge*M*(metadata['row_charges'][i]-metadata['column_charges'][j])+(bits&-bits).bit_length()-1)
    assert values
    return min(values)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    rng=random.Random(310748);checks=0
    for M in (4,8):
        for beta in (0,1):
            _,compiled=compile_interval(M,beta)
            for gauge in (Fraction(-1),Fraction(0),Fraction(1,2)):
                flags=[entry_flag(comp,shifts,metadata,gauge,M)[0] for comp,shifts,metadata in compiled]
                for trial in range(20):
                    g=[rng.randrange(4) for _ in range(12)];g[10]=1^g[0]
                    for flag,(comp,shifts,metadata) in zip(flags,compiled):
                        predicted=next(entry['order'] for entry in flag if dot(entry['form'],g))
                        actual=direct_minimum(combine_components(comp,g),shifts,metadata,gauge,M)
                        assert actual==predicted,(M,beta,gauge,g,actual,predicted)
                        checks+=1
    scale=Fraction(40813173,25000000);runs=[]
    for beta in (0,1):
        _,compiled=compile_interval(4,beta)
        run=screen(compiled,4,(Fraction(0),Fraction(0)),scale,time.monotonic()+600,100000)
        assert run['complete'] and not run['improvements'] and run['total_affine_rows']==4**11
        runs.append({key:run[key] for key in ('total_affine_rows','distinct_tensor_rows','parameter_fiber_multiplicity','certifiably_pruned','profile_queries','profile_covered_rows','unresolved_rows','equal_to_gaussian','strictly_larger_than_gaussian','profile_digest_sha256')}|{'beta':beta})
    arithmetic_checks=0
    for _ in range(100):
        baseline=[[rng.randint(-10,10) for _ in range(6)] for side in range(2)]
        candidate=[[rng.randint(-10,10) for _ in range(6)] for side in range(2)]
        comparator=ProductComparison(Fraction(7,5),baseline)
        def exact(profiles):
            value=Fraction(1)
            for profile in profiles:value*=sum((Fraction(7,5)**(-order) for order in profile),Fraction())
            return value
        a,b=exact(candidate),exact(baseline)
        assert comparator.compare(candidate)==(a>b)-(a<b)
        arithmetic_checks+=1
    result={'status':'pass','independent_coefficient_minimum_checks':checks,'exact_rational_comparison_checks':arithmetic_checks,'complete_M4_runs':runs}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
