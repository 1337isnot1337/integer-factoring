#!/usr/bin/env python3
"""Exact target rejection and literal hereditary adjugate-cap controls."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from check_joint_six_to_five_retractions import python_profiles
from check_joint_energy_conjugations import polynomial
from joint_six_to_five_retractions import crossing
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
TARGET = [0,6,9,9,14,19,19,20,22,22,22,22,24,24,24,24,30,34,34,36,36,37,38,40,41]


def caps(d, r, s):
    return [(4*sum(a[:3])+3*sum(b[:2]), 12*d) for a, b in ((r, s), (s, r))]


def normalize(r, s):
    d = -min(r)-min(s)
    assert d >= 0
    return d, [x-min(r) for x in r], [x-min(s) for x in s]


def literal_controls():
    rng = Random(0xA6D005)
    records = []
    for degree in range(1, 5):
        for case in range(8):
            matrix = []
            for i in range(4):
                row = []
                for j in range(4):
                    q = sum(rng.randrange(4) << (2*k) for k in range(degree))
                    row.append(q ^ (q << 2) ^ int(i == j))
                matrix.append(row)
            opposite = [[0]*4 for _ in range(4)]
            for x, y, u, v in product(range(2), repeat=4):
                opposite[2*x+u][2*v+y] = matrix[2*x+y][2*v+u]
            r, _ = WeightedSmithOracle(matrix).profile([0]*4, [0]*4)
            s, _ = WeightedSmithOracle(opposite, infinity=True).profile([0]*4, [0]*4)
            d, r, s = normalize(r, s)
            bounds = caps(d, r, s)
            assert all(a <= b for a, b in bounds)
            records.append(dict(degree=degree, case=case, deficit=d, r=r, s=s,
                                caps=[list(v) for v in bounds]))
    return records


def certificate():
    violations = caps(6, TARGET, TARGET)
    assert violations == [(78, 72), (78, 72)]
    # Normalization algebra: a,b are restricted-core minima relative to
    # global minima. The old q2,k4 cap gives the nonnegative extra 6b.
    # Compare coefficients of d,S3(u),S2(v),a,b exactly.
    # 4(S3(u)-3a)+3(S2(v)-2b) <= 12(d-a-b).
    left = [0, 4, 3, -12, -6]
    right = [12, 0, 0, -12, -12]
    moved = [l-r for l, r in zip(left, right)]
    assert moved == [-12, 4, 3, 0, 6]
    dense = literal_controls()
    inputs = json.loads((ROOT/'evidence/joint-six-to-five-retractions.json').read_text())
    core_indices = (0, 1, 5, 6)
    inherited = []
    for witness in inputs['witnesses']:
        code = crossing(witness['assignments'])
        global_r, global_s = python_profiles(code)
        d, r, s = normalize(global_r, global_s)
        assert all(a <= b for a, b in caps(d, r, s))
        # These zero-energy rows have zero Laurent clearing; use literal
        # packed F4 entries of the two-letter physical restriction.
        assert all(code['row_clearing'][i] == 0 for i in core_indices)
        core = [[code['matrix'][i][j] for j in core_indices] for i in core_indices]
        packed = [[polynomial(p) for p in row] for row in core]
        opposite = [[0]*4 for _ in range(4)]
        for x, y, u, v in product(range(2), repeat=4):
            opposite[2*x+u][2*v+y] = packed[2*x+y][2*v+u]
        core_r, _ = WeightedSmithOracle(packed).profile([0]*4, [0]*4)
        core_s, _ = WeightedSmithOracle(opposite, infinity=True).profile([0]*4, [0]*4)
        u = [x-min(global_r) for x in core_r]
        v = [x-min(global_s) for x in core_s]
        assert all(sum(u[:k]) >= sum(r[:k]) and sum(v[:k]) >= sum(s[:k]) for k in range(1, 5))
        assert min(u) >= 0 and min(v) >= 0
        exact = [(4*sum(a[:3])+3*sum(b[:2])+6*min(b), 12*d) for a, b in ((u, v), (v, u))]
        assert all(a <= b for a, b in exact)
        inherited.append(dict(assignments=witness['assignments'], deficit=d, global_r=r, global_s=s,
                              core_r=u, core_s=v,
                              strengthened_core_caps=[list(p) for p in exact]))
    names = (Path(__file__).name, 'joint_six_to_five_retractions.py', 'polynomial_retraction_spaces.py',
             'check_joint_six_to_five_retractions.py', 'weighted_smith_oracle_f4.py')
    return dict(scope='All-field hereditary zero-energy two-letter adjugate cap; excludes the previously retained q5 target and its average-prefix sector, not all targets at these energies.',
                paired_energies=[0,0,12,16,18], target=TARGET, deficit=6,
                target_caps=[list(v) for v in violations], target_violation=6,
                inherited_inequalities=['4*S3(r)+3*S2(s)<=12*d', '4*S3(s)+3*S2(r)<=12*d'],
                normalization_coefficients=moved, dense_two_letter_tensors=len(dense),
                dense_independent_orders=8*len(dense), restricted_five_letter_tensors=len(inherited),
                full_and_core_independent_orders=58*len(inherited), dense_controls=dense,
                inherited_controls=inherited,
                upstream_theorem_source='3670b7df8db2bca84bc1dc1ea492c7044fa93e60:research/graded-adjugate-caps.md, equation (1), q=2,k=4',
                input_sha256=sha256((ROOT/'evidence/joint-six-to-five-retractions.json').read_bytes()).hexdigest(),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/zero-energy-hereditary-adjugate.json').read_text())
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('target', 'dense_controls', 'inherited_controls', 'source_sha256')}))


if __name__ == '__main__':
    main()
