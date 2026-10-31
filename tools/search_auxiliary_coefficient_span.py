#!/usr/bin/env python3
"""Exact F2 coefficient-span thresholds and rational rank-sum optima."""
from __future__ import annotations

import argparse
from hashlib import sha256
from heapq import heappop,heappush
from itertools import product
import json
from math import isqrt
from pathlib import Path

import check_cancellation_module_automata as module
from check_cancellation_module_crossing import local_crossing

ROOT = Path(__file__).resolve().parents[1]


def coefficient_matrices(local,degree,length):
    labels,dimension = 1 << length,len(next(iter(local.values())))
    matrices = [[[0]*(labels*labels) for _ in range(degree*length+1)] for _ in range(2)]
    endpoints = 0
    def visit(depth,words,state):
        nonlocal endpoints
        if depth == length:
            x,y,u,v = words
            value = state[0]
            assert value.bit_count() % 2 == int(x == v and y == u)
            endpoints += 1
            while value:
                bit = value & -value
                k = bit.bit_length()-1
                matrices[0][k][x*labels+y] |= 1 << (u*labels+v)
                matrices[1][k][x*labels+u] |= 1 << (y*labels+v)
                value ^= bit
            return
        for physical,matrix in local.items():
            following = [0]*dimension
            for a,b in product(range(dimension),repeat=2):
                following[b] ^= module.pmul(state[a],matrix[a][b])
            visit(depth+1,tuple(2*x+y for x,y in zip(words,physical)),tuple(following))
    visit(0,(0,0,0,0),tuple(int(i == 0) for i in range(dimension)))
    return matrices,endpoints


def combine(matrices,mask):
    result = [0]*len(matrices[0])
    while mask:
        bit = mask & -mask
        row = matrices[bit.bit_length()-1]
        result = [a^b for a,b in zip(result,row)]
        mask ^= bit
    return result


def reverse_rank(rows):
    basis = {}
    for row in reversed(rows):
        while row:
            pivot = row & -row
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def add_basis(basis,vector,tag):
    while vector:
        pivot = vector.bit_length()-1
        if pivot not in basis:
            basis[pivot] = vector,tag
            return
        row,other = basis[pivot]
        vector ^= row
        tag ^= other


def reduce_basis(basis,vector):
    tag = 0
    while vector:
        pivot = vector.bit_length()-1
        if pivot not in basis:
            return vector,tag
        row,other = basis[pivot]
        vector ^= row
        tag ^= other
    return 0,tag


def minimum_threshold(pairs):
    target = len(pairs)-1
    basis = {}
    for mask in sorted(range(len(pairs)),key=lambda w:(pairs[w][0]*pairs[w][1],w)):
        add_basis(basis,mask,1 << mask)
        remainder,tag = reduce_basis(basis,target)
        if remainder == 0:
            threshold = pairs[mask][0]*pairs[mask][1]
            witness = [w for w in range(len(pairs)) if tag >> w & 1]
            break
    lower_basis = {}
    for mask,(a,b) in enumerate(pairs):
        if a*b < threshold:
            add_basis(lower_basis,mask,0)
    dual = None
    for j in range(target.bit_length()):
        if j in lower_basis:
            continue
        candidate = 1 << j
        for pivot,(row,_) in sorted(lower_basis.items()):
            if (row & candidate).bit_count() % 2:
                candidate ^= 1 << pivot
        assert all((row & candidate).bit_count() % 2 == 0 for row,_ in lower_basis.values())
        if (target & candidate).bit_count() % 2:
            dual = candidate
            break
    assert dual is not None
    assert all((mask & dual).bit_count() % 2 == 0 for mask,(a,b) in enumerate(pairs) if a*b < threshold)
    value = 0
    for mask in witness:
        value ^= mask
        assert pairs[mask][0]*pairs[mask][1] <= threshold
    assert value == target
    return {'minimum_maximum_rank_product':threshold,'decomposition_masks':witness,
            'strict_lower_threshold_separating_covector':dual}


def shortest_sum(costs):
    target = len(costs)-1
    distances = [None]*len(costs)
    distances[0] = 0
    previous = {}
    queue = [(0,0)]
    while queue:
        cost,vertex = heappop(queue)
        if distances[vertex] != cost:
            continue
        if vertex == target:
            path = []
            while vertex:
                vertex,mask = previous[vertex]
                path.append(mask)
            potential = [cost if value is None else min(cost,value) for value in distances]
            assert potential[0] == 0 and potential[target] == cost
            assert all(potential[v^w] <= potential[v]+weight
                       for v in range(len(costs)) for w,weight in enumerate(costs))
            return cost,path,potential
        for mask,weight in enumerate(costs):
            other = vertex^mask
            candidate = cost+weight
            if distances[other] is None or candidate < distances[other]:
                distances[other] = candidate
                previous[other] = vertex,mask
                heappush(queue,(candidate,other))
    raise AssertionError('The full coefficient space must contain its endpoint.')


def rank_sum_bounds(pairs,precision):
    scale = 1 << precision
    squares = [a*b*scale**2 for a,b in pairs]
    floors = [isqrt(q) for q in squares]
    ceilings = [a+int(a*a != q) for a,q in zip(floors,squares)]
    lower,_,potential = shortest_sum(floors)
    upper,path,_ = shortest_sum(ceilings)
    assert lower <= upper
    target = 0
    for mask in path:
        target ^= mask
    assert target == len(pairs)-1 and sum(ceilings[w] for w in path) == upper
    return {'denominator':scale,'optimal_rank_sum_lower_numerator':lower,
            'optimal_rank_sum_upper_numerator':upper,'upper_decomposition_masks':path,
            'rank_sum_lower_potential':potential}


def verify_record(record,local,degree):
    matrices,endpoints = coefficient_matrices(local,degree,record['copies'])
    pairs = record['rank_pairs_by_functional_mask']
    assert len(pairs) == 1 << len(matrices[0])
    assert endpoints == record['endpoint_entries']
    for mask,pair in enumerate(pairs):
        assert [reverse_rank(combine(cut,mask)) for cut in matrices] == pair
    target,threshold = len(pairs)-1,record['minimum_maximum_rank_product']
    dual = record['strict_lower_threshold_separating_covector']
    assert (target & dual).bit_count() % 2 == 1
    assert all((w & dual).bit_count() % 2 == 0 for w,(a,b) in enumerate(pairs) if a*b < threshold)
    value = 0
    for mask in record['decomposition_masks']:
        assert pairs[mask][0]*pairs[mask][1] <= threshold
        value ^= mask
    assert value == target
    scale = record['denominator']
    squares = [a*b*scale**2 for a,b in pairs]
    floors = [isqrt(q) for q in squares]
    potential = record['rank_sum_lower_potential']
    assert len(potential) == len(pairs) and potential[0] == 0
    assert potential[target] == record['optimal_rank_sum_lower_numerator']
    assert all(potential[v^w] <= potential[v]+weight
               for v in range(len(pairs)) for w,weight in enumerate(floors))
    value,upper = 0,0
    for mask in record['upper_decomposition_masks']:
        value ^= mask
        upper += floors[mask]+int(floors[mask]**2 != squares[mask])
    assert value == target and upper == record['optimal_rank_sum_upper_numerator']
    assert record['upper_rank_sum_below_rho_seven_tenths'] == (upper**5 < scale**5*2**(7*record['copies']))
    assert potential[target]**2 > scale**2*2**(3*record['copies'])
    return len(pairs)


def analyze(name,local,degree,length,precision):
    matrices,endpoints = coefficient_matrices(local,degree,length)
    count = 1 << len(matrices[0])
    pairs = [[len(module.binary_basis(combine(cut,mask))) for cut in matrices]
             for mask in range(count)]
    threshold = minimum_threshold(pairs)
    sums = rank_sum_bounds(pairs,precision)
    checks = set(threshold['decomposition_masks']+sums['upper_decomposition_masks']+[count-1])
    for mask in checks:
        assert [reverse_rank(combine(cut,mask)) for cut in matrices] == pairs[mask]
    if length <= 2:
        independent,_ = module.coefficient_ranks(local,degree,length)
        assert [[pairs[1 << k][cut] for k in range(len(matrices[0]))] for cut in range(2)] == independent
    numerator,denominator = sums['optimal_rank_sum_upper_numerator'],sums['denominator']
    assert sums['optimal_rank_sum_lower_numerator']**2 > denominator**2*2**(3*length)
    return {'family':name,'copies':length,'degree_bound':degree*length,
            'coefficient_functionals_enumerated':count,'endpoint_entries':endpoints,
            'rank_pairs_by_functional_mask':pairs,**threshold,**sums,
            'upper_rank_sum_below_rho_seven_tenths':numerator**5 < denominator**5*2**(7*length),
            'independent_pivot_order_controls':len(checks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-copies',type=int,choices=(1,2,3,4),default=4)
    parser.add_argument('--precision',type=int,default=20)
    parser.add_argument('--families',nargs='+',choices=('elementary','module','pairing','mixed-loop'),
                        default=['elementary','module','pairing','mixed-loop'])
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    args = parser.parse_args()
    assert 1 <= args.precision <= 30
    constructors = {'elementary':lambda:module.elementary(1),'module':lambda:(local_crossing(),2),
                    'pairing':module.pairing,'mixed-loop':module.mixed_loop}
    if args.verify:
        saved = json.loads(args.verify.read_text())
        assert saved['prime'] == 2 and saved['status'] == 'pass'
        for name,digest in saved['source_sha256'].items():
            assert digest == sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
        checked = 0
        for record in saved['records']:
            local,degree = constructors[record['family']]()
            checked += verify_record(record,local,degree)
        print(json.dumps({'status':'pass','records':len(saved['records']),
                          'opposite_pivot_order_rank_pairs':checked,
                          'scope':'Every saved threshold and rank-sum bound verified from explicit witnesses, separating covectors, and shortest-path dual potentials; no search rerun needed.'}))
        return
    assert args.output
    records = []
    for name in args.families:
        local,degree = constructors[name]()
        for length in range(1,min(args.max_copies,3 if degree == 3 else 4)+1):
            record = analyze(name,local,degree,length,args.precision)
            records.append(record)
            print(json.dumps({k:v for k,v in record.items() if k not in
                              ('rank_pairs_by_functional_mask','rank_sum_lower_potential')}),flush=True)
    result = {'status':'pass','prime':2,'records':records,
              'scope':'Complete F2 coefficient recombinations at the stated finite lengths; exact minmax thresholds and rational optimal rank-sum bounds, without a uniform improved crossing or factoring claim.',
              'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                                for name in ('search_auxiliary_coefficient_span.py',
                                             'check_cancellation_module_automata.py',
                                             'check_cancellation_module_crossing.py')}}
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n')


if __name__ == '__main__':
    main()
