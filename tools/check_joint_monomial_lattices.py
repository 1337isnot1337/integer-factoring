#!/usr/bin/env python3
"""Complete coordinate-lattice modifications under the q5 pole bound."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import accumulate, product
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_joint_energy_conjugations import cut_matrix, weighted_coefficient
from check_joint_leading_tangent import realign
from joint_energy_conjugations import N, model
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def certificate():
    data = model()
    edges = []
    for r, row in enumerate(data['operator']):
        for c, (polynomial, high) in enumerate(row):
            assert not high
            if polynomial:
                low = (polynomial & -polynomial).bit_length()-1
                top = polynomial.bit_length()-1
                edges.extend(((r, c, low), (c, r, data['energies'][c]-top+8)))
    # Each edge r,c,w means p_c-p_r<=w. Its weight is nonnegative.
    assert all(w >= 0 for _, _, w in edges)
    distance = [[0 if r == c else float('inf') for c in range(N)] for r in range(N)]
    for r, c, w in edges:
        distance[r][c] = min(distance[r][c], w)
    for k in range(N):
        for r in range(N):
            for c in range(N):
                distance[r][c] = min(distance[r][c], distance[r][k]+distance[k][c])
    groups, unseen = [], set(range(N))
    while unseen:
        first = min(unseen)
        group = [j for j in sorted(unseen) if distance[first][j] == distance[j][first] == 0]
        groups.append(group)
        unseen.difference_update(group)
    assert groups == [[0], [1, 2, 3, 5, 10, 15], [4, 20],
                      [6, 7, 8, 11, 12, 13, 16, 17, 18, 24], [9, 14, 19, 21, 22, 23]]
    group_of = {value: i for i, group in enumerate(groups) for value in group}
    quotient = {}
    for r, c, w in edges:
        a, b = group_of[r], group_of[c]
        if a != b:
            quotient[a, b] = min(quotient.get((a, b), w), w)
    assert quotient == {(0, 3): 1, (3, 0): 0, (2, 3): 1, (3, 2): 1,
                        (1, 4): 1, (4, 1): 0}
    records, tensors = [], set()
    half = data['half']
    weights = [half[i//5] for i in range(N)]
    baseline_infinity, _ = WeightedSmithOracle(cut_matrix(data['values'], True), infinity=True).profile(weights, weights)
    baseline_zero, _ = WeightedSmithOracle(cut_matrix(data['values'], False)).profile(data['energies'], [0]*N)
    for a, c in product((0, 1), repeat=2):
        for b in range(a-1, a+2):
            charges = [0]*N
            for group, value in zip(groups, (0, 0, b, a, c)):
                for i in group:
                    charges[i] = value
            assert all(charges[j]-charges[i] <= w for i, j, w in edges)
            values = []
            for index, (low, high) in enumerate(data['values']):
                x, y, u, v = ((index//5**i) % 5 for i in range(4))
                shift = charges[x+5*y]-charges[v+5*u]
                if shift < 0:
                    assert not ((low | high) & ((1 << -shift)-1))
                    values.append((low >> -shift, high >> -shift))
                else:
                    values.append((low << shift, high << shift))
                assert ((values[-1][0].bit_count() % 2)
                        | ((values[-1][1].bit_count() % 2) << 1)) == int(x == v and y == u)
            tensors.add(tuple(values))
            infinity, _ = WeightedSmithOracle(cut_matrix(values, True), infinity=True).profile(weights, weights)
            zero, _ = WeightedSmithOracle(cut_matrix(values, False)).profile(data['energies'], [0]*N)
            assert infinity[0] == -8
            assert all(x <= y for x, y in zip(accumulate(infinity), accumulate(baseline_infinity)))
            assert all(x <= y for x, y in zip(accumulate(zero), accumulate(baseline_zero)))
            leading, high = weighted_coefficient(values, data['energies'], -8)
            assert not high
            matrix = realign(leading)
            rank = binary_rank([(matrix >> (N*r)) & ((1 << N)-1) for r in range(N)])
            records.append(dict(parameters=[a, b, c], infinity=infinity, zero=zero,
                                leading_realignment_rank=rank))
    assert len(records) == len(tensors) == 12
    assert sum(r['leading_realignment_rank'] == 1 for r in records) == 1
    return dict(zero_constraint_groups=groups,
                quotient_edges=[[*pair, w] for pair, w in sorted(quotient.items())],
                neutral_components=[[0, 2, 3], [1, 4]],
                baseline_profiles=[baseline_infinity, baseline_zero], records=records,
                fraction_free_profiles=26,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_monomial_lattices.py', 'check_joint_energy_conjugations.py',
                                            'joint_energy_conjugations.py', 'weighted_smith_oracle_f4.py',
                                            'check_joint_leading_tangent.py', 'search_f4_recentered_crossings.py')},
                scope='Every integer joint coordinate charge retaining polynomiality and no infinity pole below -8, on this reference at its natural heights. Binary profiles persist under any characteristic-two field extension.',
                conclusion='Exactly 12 effective tensors modulo two neutral component shifts; both profiles are prefix-dominated by the reference, so no candidate improves its cost at any scale.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-monomial-lattices.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == json.loads(json.dumps(result))
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_sha256', 'records', 'baseline_profiles')}))


if __name__ == '__main__':
    main()
