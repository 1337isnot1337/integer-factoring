#!/usr/bin/env python3
"""Independent full-tensor controls for complete-entry affine synthesis."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_block_smith_crossings import reverse_polynomial,smith_profile,toeplitz_rank
from check_physical_split_crossings import scalar_restriction
from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice,CompressedFullAlphabetFiber,cut_location
from reflection_dual_rail_gaussian import tensor
from search_binary_two_cut_lattices import TwoCutLattice,BASE,HEIGHTS


def independent_profiles(lattice,values):
    q,size=lattice.q,lattice.size;result=[];toeplitz=0
    for side in range(2):
        h=lattice.heights[4*side:4*side+4];scalar=lattice.scalar_shifts[side]
        matrix=[[(0,0)]*size for _ in range(size)]
        for entry,pair in enumerate(values):
            bits=lattice.labels[entry];i,j=cut_location(bits,q,side)
            weight=sum(row[label] for row,label in zip(h,bits))
            if side==0:
                shift=scalar-lattice.degree+weight
                matrix[i][j]=tuple(reverse_polynomial(part,lattice.degree)<<shift for part in pair)
            else:matrix[i][j]=tuple(part<<(scalar+weight) for part in pair)
        binary=scalar_restriction(matrix)
        bound=max(value.bit_length()-1 for row in binary for value in row)
        orders,precision=smith_profile(binary,bound)
        assert all(orders[i]==orders[i+1] for i in range(0,len(orders),2))
        for length in sorted({1,orders[0]+1,orders[len(orders)//2]+1,orders[-1]+1}):
            assert toeplitz_rank(binary,length)==sum(max(0,length-value) for value in orders)
            toeplitz+=1
        result.append(tuple(orders[i]-scalar for i in range(0,len(orders),2)))
    return tuple(result),toeplitz


def restricted_reference(alphabet):
    seed=tensor(4);q=len(alphabet)
    return tuple(seed['values'][alphabet[index%q]+4*alphabet[(index//q)%q]
                 +16*alphabet[(index//q**2)%q]+64*alphabet[(index//q**3)%q]]
                 for index in range(q**4))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--search-records',nargs='*',type=Path,default=[])
    args=parser.parse_args();rng=random.Random(311758)
    reference=TwoCutLattice(5)
    generalized=FullAlphabetTwoCutLattice([(value,0) for value in BASE],2,5,
                                         [(0,height) for height in HEIGHTS],raise_cap=4)
    assert generalized.orders==reference.orders
    binary_comparisons=0
    for targets in [reference.orders]+[tuple(value+int(i==j) for i,value in enumerate(reference.orders))
                                      for j in range(8)]:
        assert reference.solve(targets)==generalized.solve(targets);binary_comparisons+=1
    cases=[];endpoints=complete_profiles=toeps=feasibility=0
    for alphabet,energies in [((0,1,2),(0,12,4)),((0,1,2),(0,0,8)),
                              ((0,1,2,3),(0,12,4,16)),((0,1,2,3),(0,0,12,16))]:
        values=restricted_reference(alphabet);q=len(alphabet);degree=max(a.bit_length()-1 for a,b in values)
        zero=[0]*q;half=[value//2 for value in energies]
        model=FullAlphabetTwoCutLattice(values,q,degree,[zero,zero,half,half,half,half,zero,zero],raise_cap=4)
        floors=tuple(max(0,value-3) for value in model.orders)
        fiber=CompressedFullAlphabetFiber(model,floors)
        candidates=[model.orders,floors]
        for _ in range(4):
            candidate=[max(floors[i],value-rng.randrange(4)) for i,value in enumerate(model.orders)]
            for i in rng.sample(range(len(candidate)),3):candidate[i]=min(model.caps[i],candidate[i]+rng.randrange(1,4))
            candidates.append(tuple(candidate))
        samples=0
        for target in candidates:
            full=model.solve(target);small=fiber.solve(target);feasibility+=1
            assert (full is None)==(small is None)
            if full is None:continue
            assert len(full)==len(fiber.base_pivots)+len(small)
            witness=fiber.sample(small,rng)
            actual=model.check_witness(witness,target);endpoints+=q**4;complete_profiles+=2
            if samples<2:
                independent,count=independent_profiles(model,witness)
                assert actual==independent;toeps+=count;complete_profiles+=2
            samples+=1
        cases.append({'alphabet':alphabet,'energies':energies,'variables':model.variables,
                      'compressed_variables':fiber.variables,'feasible_full_tensor_witnesses':samples})
    negative_controls=0
    damaged=[(value,0) for value in BASE];damaged[0]=(damaged[0][0]^1,0)
    try:FullAlphabetTwoCutLattice(damaged,2,5,[(0,height) for height in HEIGHTS])
    except ValueError:negative_controls+=1
    else:raise AssertionError('Invalid crossing endpoint accepted.')
    for invalid in [tuple(value+1 for value in generalized.caps),(-1,)+generalized.orders[1:]]:
        try:generalized.solve(invalid)
        except ValueError:negative_controls+=1
        else:raise AssertionError('Incomplete coefficient window accepted.')
    search_replays=[]
    for path in args.search_records:
        record=json.loads(path.read_text());assert record['completed_requested_trials']
        for name,digest in record['source_sha256'].items():
            assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()==digest
        values=tensor(4)['values'];model=FullAlphabetTwoCutLattice(values,4,record['degree'],
            record['heights_infinity_then_zero'],raise_cap=record['arguments']['upper_window'])
        for row in record['retained_witnesses']:
            actual=model.check_witness(row['values'],row['targets'])
            assert actual==tuple(map(tuple,row['profiles']))
        search_replays.append({'record':path.name,'completed_trials':record['completed_trials'],
            'distinct_profile_classes':record['distinct_profile_classes'],
            'retained_full_tensor_replays':len(record['retained_witnesses'])})
    result={'status':'Complete-entry equations, compressed fibers, endpoints, and independent profiles passed.',
        'scope':'Exact sufficient row-lattice construction at arbitrary physical alphabet size; finite controls at two, three, and four symbols.',
        'binary_reference_system_comparisons':binary_comparisons,
        'full_and_compressed_feasibility_comparisons':feasibility,'literal_endpoint_entries':endpoints,
        'complete_native_and_independent_profiles':complete_profiles,'independent_toeplitz_checks':toeps,
        'invalid_endpoint_or_target_controls':negative_controls,'cases':cases,'search_record_replays':search_replays,
        'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            (Path(__file__).name,'full_alphabet_two_cut_lattices.py','weyl_lattice_filter_synthesis.py',
             'check_block_smith_crossings.py','check_physical_split_crossings.py')}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
