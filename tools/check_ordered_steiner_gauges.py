#!/usr/bin/env python3
"""Exact finite rank-sum certificates for ordered Steiner physical gauges."""
from __future__ import annotations

import argparse
from hashlib import sha256
import heapq
from itertools import product
import json
from math import isqrt
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from ordered_steiner_gauges import BOUNDARIES, FIXED_FIRST, GL2, IDENTITY, coefficient_rows

DENOMINATOR = 1 << 20
ROOT = Path(__file__).resolve().parents[1]


def rank_pairs(sides, low=False):
    size, dimension = len(sides[0][0]), len(sides[0])
    current, previous = [[0]*size for _ in range(2)], 0
    result = [(0, 0)]*(1 << dimension)
    for index in range(1, len(result)):
        code = index ^ (index >> 1)
        changed = (code ^ previous).bit_length()-1
        previous = code
        current = [[a ^ b for a, b in zip(rows, side[changed])] for rows, side in zip(current, sides)]
        result[code] = tuple(binary_rank(rows, low=low) for rows in current)
    return result


def weights(pairs, upper):
    result = []
    for r, c in pairs:
        square = r*c*DENOMINATOR**2
        weight = isqrt(square)
        result.append(weight+int(upper and weight*weight != square))
    return result


def shortest_paths(costs):
    size = len(costs)
    distance, parent = [10**30]*size, [None]*size
    distance[0], queue = 0, [(0, 0)]
    while queue:
        value, vertex = heapq.heappop(queue)
        if value != distance[vertex]:
            continue
        for code, weight in enumerate(costs[1:], 1):
            target = vertex ^ code
            candidate = value+weight
            if candidate < distance[target]:
                distance[target], parent[target] = candidate, (vertex, code)
                heapq.heappush(queue, (candidate, target))
    path, vertex = [], size-1
    while vertex:
        vertex, code = parent[vertex]
        path.append(code)
    return distance, path


def receipt(second, e, f, length):
    pairs = rank_pairs(coefficient_rows(FIXED_FIRST, second, e, f, length))
    potential, _ = shortest_paths(weights(pairs, False))
    upper_distances, path = shortest_paths(weights(pairs, True))
    return dict(second=list(second), e=e, f=f, length=length, ranks=pairs,
                lower_potential=potential, upper_path=path,
                lower_numerator=potential[-1], upper_numerator=upper_distances[-1])


def generate():
    records = [receipt((IDENTITY, a, b), e, f, 2)
               for a, b in product(GL2, repeat=2) for e, f in BOUNDARIES]
    best = min(records, key=lambda record: record['upper_numerator'])
    selected = receipt(best['second'], best['e'], best['f'], 3)
    hashes = {name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
              for name in ('ordered_steiner_gauges.py', 'check_ordered_steiner_gauges.py',
                           'check_steiner_cubic_crossings.py')}
    return dict(denominator=DENOMINATOR, first=list(FIXED_FIRST), two_position_cases=records,
                selected_three_position=selected, source_sha256=hashes,
                scope='Fixed first GL2(F2) table, all 36 normalized second tables and six valid boundaries; one selected three-position continuation.',
                factoring_exponent='1/2; unchanged')


def verify(data):
    assert data['denominator'] == DENOMINATOR and tuple(data['first']) == FIXED_FIRST
    expected = {(a, b, e, f) for a, b in product(GL2, repeat=2) for e, f in BOUNDARIES}
    assert len(data['two_position_cases']) == len(expected) == 216
    records = data['two_position_cases']+[data['selected_three_position']]
    rank_count = inequalities = 0
    for record in records:
        second, e, f, length = record['second'], record['e'], record['f'], record['length']
        if length == 2:
            assert second[0] == IDENTITY
            expected.remove((second[1], second[2], e, f))
        else:
            assert length == 3
        pairs = rank_pairs(coefficient_rows(FIXED_FIRST, second, e, f, length), low=True)
        assert pairs == list(map(tuple, record['ranks']))
        rank_count += 2*len(pairs)
        lower, upper = weights(pairs, False), weights(pairs, True)
        potential = record['lower_potential']
        assert len(potential) == len(pairs) and potential[0] == 0
        for vertex in range(len(pairs)):
            for code, cost in enumerate(lower):
                assert potential[vertex ^ code]-potential[vertex] <= cost
                inequalities += 1
        target = value = 0
        for code in record['upper_path']:
            target ^= code
            value += upper[code]
        assert target == len(pairs)-1 and value == record['upper_numerator']
        assert potential[-1] == record['lower_numerator'] <= value
        assert potential[-1]**2 > 3**(3*length)*DENOMINATOR**2
    assert not expected
    for name, digest in data['source_sha256'].items():
        assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest
    return dict(complete_two_position_cases=216, selected_three_position_cases=1,
                independently_replayed_cut_ranks=rank_count, exact_dual_inequalities=inequalities,
                minimum_two_position_lower_numerator=min(r['lower_numerator'] for r in data['two_position_cases']),
                denominator=DENOMINATOR, exact_comparison='Every retained lower rank-sum bound exceeds 3^(3m/2).',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/ordered-steiner-gauges.json')
    args = parser.parse_args()
    if args.generate:
        data = generate()
        args.evidence.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    else:
        data = json.loads(args.evidence.read_text())
    print(json.dumps(verify(data)))


if __name__ == '__main__':
    main()
