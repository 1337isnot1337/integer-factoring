#!/usr/bin/env python3
"""Independent all-minor controls for unrestricted binary polynomial synthesis."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random

from search_binary_two_cut_lattices import TwoCutLattice, BASE, HEIGHTS, flatten, labels
from search_binary_f4_flags import FlagSearch, flags as f4_flags, insert
from search_binary_flag_lattices import flags as binary_flags
from f4_lattice_linear import MUL, add, scale, coefficient


def pmul(first, second):
    result = {}
    for i, a in first.items():
        for j, b in second.items():
            result[i+j] = result.get(i+j,0) ^ MUL[a][b]
    return {power:value for power,value in result.items() if value}


def polynomial(pair):
    return {power:((pair[0] >> power) & 1) | (((pair[1] >> power) & 1) << 1)
            for power in range(max(pair).bit_length()) if ((pair[0] | pair[1]) >> power) & 1}


def literal_profiles(values, heights=HEIGHTS):
    result = []
    for side in range(2):
        matrix = [[{} for _ in range(4)] for _ in range(4)]
        for entry,pair in enumerate(values):
            row,column = flatten(entry,side)
            weight = sum(a*b for a,b in zip(labels(entry),heights[4*side:4*side+4]))
            matrix[row][column] = {(weight-power if side==0 else weight+power):value
                                   for power,value in polynomial(pair).items()}
        prefixes=[0]
        for size in range(1,5):
            minimum=None
            for rows in combinations(range(4),size):
                for columns in combinations(range(4),size):
                    determinant={}
                    for ordering in permutations(columns):
                        term={0:1}
                        for i,j in zip(rows,ordering):
                            term=pmul(term,matrix[i][j])
                        for power,value in term.items():
                            determinant[power]=determinant.get(power,0)^value
                    support=[power for power,value in determinant.items() if value]
                    if support:
                        order=min(support);minimum=order if minimum is None else min(minimum,order)
            assert minimum is not None
            prefixes.append(minimum)
        result.append(tuple(b-a for a,b in zip(prefixes,prefixes[1:])))
    return tuple(result)


def endpoint(values):
    for entry,pair in enumerate(values):
        x,y,u,v=labels(entry)
        assert ((pair[0].bit_count() & 1) | ((pair[1].bit_count() & 1) << 1)) == int(x==v and y==u)


def field_span(rows):
    space={(0,0)}
    for row in rows:
        encoded=(sum((value & 1) << i for i,value in enumerate(row)),
                 sum((value >> 1) << i for i,value in enumerate(row)))
        space={add(value,scale(encoded,multiplier)) for value in space for multiplier in range(4)}
    return tuple(sorted(space))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();rng=random.Random(310783)
    assert len(binary_flags())==315
    all_flags=f4_flags()
    keys=set()
    for frame in all_flags:
        key=tuple(field_span(frame[start:]) for start in (1,2,3))
        assert tuple(map(len,key))==(64,16,4)
        assert len(field_span(frame))==256
        assert key not in keys;keys.add(key)
    assert len(keys)==8925
    comparisons=literal_endpoint_coefficients=0
    controls=[]
    for degree in (5,8,12):
        lattice=TwoCutLattice(degree)
        solution=lattice.solve(lattice.orders);assert solution is not None
        free=[i for i in range(lattice.variables) if i not in solution]
        assert len(free)==1
        seeds=[lattice.seed(lattice.vector(solution,(color & 1) << free[0]),
                            lattice.vector(solution,(color >> 1) << free[0],True))
               for color in range(4)]
        seeds.extend(lattice.sample(solution,rng) for _ in range(4))
        for seed in seeds:
            endpoint(seed['values']);literal_endpoint_coefficients+=16
            assert lattice.check_witness(seed,lattice.orders)==literal_profiles(seed['values'])
            comparisons+=2
        controls.append({'degree':degree,'variables':lattice.variables,'base_lattice_dimension':lattice.variables-len(solution),
                         'all_F4_affine_line_assignments_replayed':4})
    laurent_controls=[]
    original=literal_profiles([(value,0) for value in BASE])
    for degree,clearing in ((9,1),(12,4),(16,7)):
        lattice=TwoCutLattice(degree,baseline_shift=clearing)
        shifted=tuple(tuple(value+(-clearing if side==0 else clearing) for value in original[side])
                      for side in range(2))
        assert lattice.profiles(lattice.seed(0))==shifted
        solution=lattice.solve(lattice.orders);assert solution is not None
        free=[i for i in range(lattice.variables) if i not in solution]
        assert len(free)==1
        seeds=[lattice.seed(lattice.vector(solution,(color & 1) << free[0]),
                            lattice.vector(solution,(color >> 1) << free[0],True))
               for color in range(4)]+[lattice.sample(solution,rng)]
        for seed in seeds:
            endpoint(seed['values']);literal_endpoint_coefficients+=16
            assert lattice.check_witness(seed,lattice.orders)==literal_profiles(seed['values'])
            assert lattice.profiles(seed)==shifted
            comparisons+=2
        laurent_controls.append({'degree_after_clearing':degree,'common_monomial_clearing':clearing,
                                 'base_lattice_dimension':lattice.variables-len(solution),
                                 'all_F4_affine_line_assignments_replayed':4,
                                 'complete_witness_replays':5})
    lattice=TwoCutLattice(5)
    # Independent non-Gaussian colored coefficients exercise both planes.
    for colors in product(range(4),repeat=3):
        low=high=0
        for variable,value in zip((1,18,63),colors):
            low|=(value&1)<<variable;high|=(value>>1)<<variable
        seed=lattice.seed(low,high)
        endpoint(seed['values']);literal_endpoint_coefficients+=16
        assert lattice.profiles(seed)==literal_profiles(seed['values'])
        comparisons+=2
    engine=FlagSearch(5)
    # A loose target admits arbitrary nonbinary flags and colored witnesses.
    nonbinary_controls=0
    for _ in range(12):
        frames=[rng.choice(all_flags) for _ in range(2)]
        targets=(0,)*8
        pivots={}
        for side,frame in enumerate(frames):
            system=engine.system(side,frame,targets[4*side:4*side+4]);assert system is not None
            for a,b,rhs in system.values():
                assert insert(pivots,(a|((rhs&1)<<engine.lattice.variables),b|((rhs>>1)<<engine.lattice.variables)),engine.lattice.variables)
        seed=engine.sample(pivots,rng)
        endpoint(seed['values']);literal_endpoint_coefficients+=16
        assert engine.verify(seed,frames,targets)==literal_profiles(seed['values'])
        comparisons+=2;nonbinary_controls+=1
    pruning_controls=0
    for clearing in (0,4):
        engine=FlagSearch(5 if clearing==0 else 12,clearing)
        for trial in range(12):
            frames=[rng.choice(all_flags) for _ in range(2)]
            targets=list(engine.lattice.orders)
            if trial%3==0:targets=[0]*8
            elif trial%3==2:targets[0]+=1;targets[3]-=1
            answers=[]
            for prune in (False,True):
                common=engine.common_system(targets) if prune else {}
                pivots=None if common is None else dict(common)
                if pivots is not None:
                    for side in range(2):
                        system=engine.system(side,frames[side],targets[4*side:4*side+4],common)
                        if system is None or not all(insert(pivots,(a|((rhs&1)<<engine.lattice.variables),
                                                                    b|((rhs>>1)<<engine.lattice.variables)),
                                                            engine.lattice.variables)
                                                     for a,b,rhs in system.values()):
                            pivots=None;break
                answers.append(pivots)
            assert (answers[0] is None)==(answers[1] is None)
            if answers[0] is not None:assert len(answers[0])==len(answers[1])
            if answers[1] is not None and trial%3==0:
                seed=engine.sample(answers[1],rng)
                endpoint(seed['values']);literal_endpoint_coefficients+=16
                assert engine.verify(seed,frames,targets)==literal_profiles(seed['values'])
                comparisons+=2
            pruning_controls+=1
    result={'status':'All independent exact controls passed; no pressure claim.',
            'binary_complete_flags':315,'F4_complete_flags':8925,'distinct_F4_tail_space_chains':len(keys),
            'complete_weighted_cut_profiles_by_all_minors':comparisons,
            'literal_endpoint_coefficients':literal_endpoint_coefficients,
            'base_lattice_controls':controls,'colored_free_coefficient_tables':64,
            'Laurent_clearing_controls':laurent_controls,
            'nonbinary_flag_witnesses':nonbinary_controls,
            'common_minimum_pruning_pair_comparisons':pruning_controls,
            'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                for name in ('search_binary_two_cut_lattices.py','search_binary_flag_lattices.py',
                             'search_binary_f4_flags.py',Path(__file__).name)}}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
