#!/usr/bin/env python3
"""Exact controls for the saturated zero-energy parity-core obstruction."""
import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
from random import Random

from check_fixed_field_descent import Field
from joint_six_to_five_retractions import crossing, source
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
FIELD = Field(7)


def multiply(a, b):
    out = {}
    for i in range((a.bit_length()+1)//2):
        for j in range((b.bit_length()+1)//2):
            out[i+j] = out.get(i+j, 0) ^ FIELD.mul((a >> (2*i)) & 3, (b >> (2*j)) & 3)
    return sum(v << (2*k) for k, v in out.items())


def order(p):
    assert p
    return ((p & -p).bit_length()-1)//2


def reverse(p, degree):
    assert (p.bit_length()-1)//2 <= degree
    return sum(((p >> (2*k)) & 3) << (2*(degree-k)) for k in range(degree+1))


def block_profile(matrix):
    result = []
    determinants = []
    for positions in ((0, 3), (1, 2)):
        i, j = positions
        det = multiply(matrix[i][i], matrix[j][j]) ^ multiply(matrix[i][j], matrix[j][i])
        first = min(order(matrix[a][b]) for a in positions for b in positions if matrix[a][b])
        result += [first, order(det)-first]
        determinants.append(det)
    return sorted(result), determinants


def certificate():
    d = 6
    # Cell saturation forces total 4d at both cuts. The reciprocal support
    # bound total(r)+4*min(s)<=4d forces both core minima to vanish.
    profiles = []
    for a, b, c in product(range(25), repeat=3):
        last = 24-a-b-c
        if 0 <= a <= b <= c <= last and a == 0 and a+b == 6 and a+b+c >= 15:
            profiles.append([a, b, c, last])
    assert profiles == [[0, 6, 9, 9]]
    parity_profiles = [sorted((a, 12-a, b, 12-b)) for a, b in product(range(7), repeat=2)]
    assert [0, 6, 9, 9] not in parity_profiles
    six_letter = json.loads((ROOT/'evidence/hereditary-q6-integer-target.json').read_text())
    target = six_letter['target']
    assert target[:4] == [0, 6, 7, 11] and six_letter['deficit'] == 6
    possible_core = [[0, 6, c, 18-c] for c in range(7, 10)]
    assert all(p not in parity_profiles for p in possible_core)
    # Formal coefficients of the two independent zero-paired leg charges.
    # At entry (x,y)|(v,u) they are (x-v,y-u).
    identities = []
    for side, positions, permutation in product((0, 1), ((0, 3), (1, 2)), permutations(range(2))):
        a = b = 0
        for i, row in enumerate(positions):
            column = positions[permutation[i]]
            if side == 0:
                x, y = divmod(row, 2)
                v, u = divmod(column, 2)
            else:
                x, u = divmod(row, 2)
                v, y = divmod(column, 2)
            a += x-v
            b += y-u
        assert (a, b) == (0, 0)
        identities.append([side, list(positions), list(permutation), a, b])
    rng = Random(319066)
    controls = []
    for degree in range(1, 5):
        for case in range(6):
            matrix = [[0]*4 for _ in range(4)]
            for i, j in product(range(4), repeat=2):
                if i.bit_count() % 2 == j.bit_count() % 2:
                    q = sum(rng.randrange(4) << (2*k) for k in range(degree))
                    matrix[i][j] = q ^ (q << 2) ^ int(i == j)
            opposite = [[0]*4 for _ in range(4)]
            for x, y, u, v in product(range(2), repeat=4):
                opposite[2*x+u][2*v+y] = reverse(matrix[2*x+y][2*v+u], degree)
            for cut in (matrix, opposite):
                for i, j in product(range(4), repeat=2):
                    value = 0
                    for k in range(degree+1):
                        value ^= (cut[i][j] >> (2*k)) & 3
                    assert value == int(i == j)
                exact, dets = block_profile(cut)
                independent, _ = WeightedSmithOracle(cut).profile([0]*4, [0]*4)
                assert list(independent) == exact
                assert all((p.bit_length()-1)//2 <= 2*degree for p in dets)
                controls.append(dict(degree=degree, case=case, profiles=exact, block_determinants=dets))
    data = source()
    parity = [label.bit_count() % 2 for label in data['letters']]
    joint = [parity[x] ^ parity[y] for y in range(6) for x in range(6)]
    assert all(p == 0 for i, row in enumerate(data['matrix']) for j, p in enumerate(row) if joint[i] != joint[j])
    retained = [parity[x] ^ parity[y] for x in range(5) for y in range(5)]
    escaped = crossing([(0, 57, 2)])
    violations = [(i, j) for i, row in enumerate(escaped['matrix']) for j, pair in enumerate(row)
                  if retained[i] != retained[j] and any(pair)]
    assert len(violations) == 6
    core = (0, 1, 5, 6)
    core_violations = [p for p in violations if p[0] in core and p[1] in core]
    assert core_violations == [(0, 1)]
    names = (Path(__file__).name, 'check_fixed_field_descent.py', 'weighted_smith_oracle_f4.py',
             'joint_six_to_five_retractions.py', 'polynomial_retraction_spaces.py')
    return dict(scope='All-field proof applies to the stated opposite-parity zero-energy core and both target prefix sectors; controls do not exclude other favorable profiles or same-parity zero labels.',
                deficit=d, target_prefixes=[0, 6, 15, 24], forced_core_profiles=profiles,
                block_determinant_valuation=12, parity_block_profile_allocations=len(parity_profiles),
                new_q6_core_possibilities=possible_core,
                formal_charge_identities=identities, literal_two_letter_tensors=24,
                independent_four_order_profiles=len(controls), independent_orders=4*len(controls),
                block_determinant_controls=2*len(controls), endpoint_entries=16*len(controls),
                source_joint_parity_entries=1296, escape_assignments=[[0, 57, 2]],
                parity_breaking_escape_entries=[list(p) for p in violations], controls=controls,
                parity_breaking_core_entries=[list(p) for p in core_violations],
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/zero-energy-parity-core.json').read_text())
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('formal_charge_identities', 'controls', 'source_sha256')}))


if __name__ == '__main__':
    main()
