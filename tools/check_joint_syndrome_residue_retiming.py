#!/usr/bin/env python3
"""Exact residue/color, endpoint-cover, retiming, and transport controls."""

import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import random

from joint_syndrome_filter_synthesis import transform
from joint_syndrome_residue_retiming import (
    decompose_binary_components, endpoint_cover, retime_literal,
)


def literal_physical_beta(f, types, spacing):
    originals = []
    for state in range(1 << len(types)):
        syndrome = 0
        for bit, column in enumerate(types):
            if state & (1 << bit):
                syndrome ^= column
        originals.append(f[syndrome])
    result = []
    for state in range(len(originals)):
        out, subset = {}, state
        while True:
            power = state.bit_count()-subset.bit_count()
            shift = power
            while True:
                for degree, coefficient in originals[subset].items():
                    at = degree+spacing*shift
                    out[at] = out.get(at,0)^coefficient
                if shift == 0:
                    break
                shift = (shift-1)&power
            if subset == 0:
                break
            subset = (subset-1)&state
        result.append({degree: coefficient for degree,coefficient in out.items() if coefficient})
    return result


def exact_matching(zero_sets, infinity_sets, alpha, gamma):
    left, right = len(zero_sets),len(infinity_sets)
    source,sink=left+right,left+right+1
    size=sink+1
    capacity=[[F(0) for _ in range(size)] for _ in range(size)]
    for i,value in enumerate(alpha):
        capacity[source][i]=value
    for j,value in enumerate(gamma):
        capacity[left+j][sink]=value
    for i,s in enumerate(zero_sets):
        for j,t in enumerate(infinity_sets):
            if s&t:
                capacity[i][left+j]=F(1)
    flow=F(0)
    while True:
        previous=[None]*size;previous[source]=source;queue=deque([source])
        while queue and previous[sink] is None:
            node=queue.popleft()
            for other,value in enumerate(capacity[node]):
                if value>0 and previous[other] is None:
                    previous[other]=node;queue.append(other)
        if previous[sink] is None:
            break
        amount,node=F(1),sink
        while node!=source:
            amount=min(amount,capacity[previous[node]][node]);node=previous[node]
        node=sink
        while node!=source:
            parent=previous[node];capacity[parent][node]-=amount;capacity[node][parent]+=amount;node=parent
        flow+=amount
    return flow


def check():
    generator=random.Random(817371)
    counts={'original_profile_controls':0,'endpoint_cover_controls':0,
            'strict_component_reductions':0,'retimed_literal_vectors':0,
            'literal_physical_state_controls':0,'degree_cap_controls':0,
            'exact_transport_networks':0,'exact_transport_inequalities':0}
    records=[]
    for case in range(48):
        rank=2 if case<12 else 3
        types=tuple(range(1,1<<rank)) if case in (0,1,12,13) else tuple(1<<i for i in range(rank))
        spacing=generator.choice((1,3,5,127))
        f=[]
        for _ in range(1<<rank):
            row={degree:generator.randrange(1,8) for degree in generator.sample(range(4*spacing+1),generator.randrange(1,5))}
            endpoint=0
            for value in row.values():endpoint^=value
            if not endpoint:
                row[0]=row.get(0,0)^1
                row={degree:value for degree,value in row.items() if value}
            f.append(row)
        delta=transform(f,rank,types,spacing)
        original_orders,original_degrees=tuple(map(min,delta)),tuple(map(max,f))
        decomposition=decompose_binary_components(f,rank,types,spacing)
        positions=tuple(row['source_residue'] for row in decomposition['components'])
        lifted=retime_literal(decomposition,positions,spacing)
        assert lifted['actual_orders']==original_orders
        assert lifted['actual_degrees']==original_degrees
        counts['original_profile_controls']+=1
        cover=endpoint_cover(decomposition)
        assert len(cover['components'])<=1<<rank
        cover_positions=tuple(row['source_residue'] for row in cover['components'])
        reduced=retime_literal(cover,cover_positions,spacing)
        assert all(a>=b for a,b in zip(reduced['actual_orders'],original_orders))
        assert all(a<=b for a,b in zip(reduced['actual_degrees'],original_degrees))
        counts['endpoint_cover_controls']+=1
        counts['strict_component_reductions']+=int(reduced['actual_orders']!=original_orders or reduced['actual_degrees']!=original_degrees)
        length=len(cover['components'])
        retimings=[(0,)*length,(7,)*length,
                   tuple(7*i//max(1,length-1) for i in range(length)),
                   tuple(7-7*i//max(1,length-1) for i in range(length)),
                   tuple(generator.randrange(8) for _ in range(length))]
        for values in retimings:
            witness=retime_literal(cover,values,7)
            beta=literal_physical_beta(witness['f'],types,7)
            for state,row in enumerate(beta):
                shift=7*(rank-state.bit_count())
                assert {degree+shift:coefficient for degree,coefficient in row.items()}==witness['delta'][state]
            counts['retimed_literal_vectors']+=1
            counts['literal_physical_state_controls']+=len(beta)
            caps=[F(max(row),7) for row in witness['f']]
            for j,component in enumerate(cover['components']):
                upper=min([F(1)]+[caps[b]-degree for b,degree in enumerate(component['degrees']) if degree is not None])
                assert F(values[j],7)<=upper
            counts['degree_cap_controls']+=1
        records.append([rank,types,spacing,len(decomposition['components']),length,
                        original_orders,original_degrees,reduced['actual_orders'],reduced['actual_degrees']])

    # If the colors are discarded, colliding identical components cancel.
    control=decompose_binary_components(({0:1,1:2},{0:1,1:2}),1,(1,),3)
    collision=retime_literal(control,(0,0),7)
    assert collision['f']==({0:3},{0:3}) and all(collision['endpoint_values'])
    assert 1^1==0

    transport_records=[]
    for _ in range(64):
        components=generator.randrange(1,5);full=(1<<components)-1
        zero_sets=[generator.randrange(1,full+1) for _ in range(generator.randrange(2,7))]
        infinity_sets=[generator.randrange(1,full+1) for _ in range(generator.randrange(2,7))]
        left=[generator.randrange(1,8) for _ in zero_sets]
        right=[generator.randrange(1,8) for _ in infinity_sets]
        alpha=[F(value,sum(left)) for value in left];gamma=[F(value,sum(right)) for value in right]
        deficit=max(sum((weight for mask,weight in zip(zero_sets,alpha) if mask&~subset==0),F(0))
                    +sum((weight for mask,weight in zip(infinity_sets,gamma) if mask&subset==0),F(0))-1
                    for subset in range(1<<components))
        assert exact_matching(zero_sets,infinity_sets,alpha,gamma)==1-deficit
        counts['exact_transport_networks']+=1
        for _ in range(4):
            positions=[generator.randrange(8) for _ in range(components)]
            base=F(3,4);q=base**7
            zero=sum((weight*base**min(positions[j] for j in range(components) if mask>>j&1)
                      for mask,weight in zip(zero_sets,alpha)),F(0))
            infinity=sum((weight*base**(-max(positions[j] for j in range(components) if mask>>j&1))
                          for mask,weight in zip(infinity_sets,gamma)),F(0))
            ratio=zero*infinity
            assert ratio**deficit.denominator>=q**deficit.numerator
            counts['exact_transport_inequalities']+=1
        transport_records.append([components,zero_sets,infinity_sets,left,right,
                                  [deficit.numerator,deficit.denominator]])
    return {'status':'Exact independent-color reduction, endpoint-cover dominance, collision/reordering, degree caps, and transport controls passed; no factoring gain claimed.',
            'counts':counts,'profile_control_sha256':sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest(),
            'transport_control_sha256':sha256(json.dumps(transport_records,separators=(',',':')).encode()).hexdigest(),
            'independent_color_negative_control':collision,'seed':817371}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=check()
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['counts']))


if __name__=='__main__':main()
