#!/usr/bin/env python3
"""Check the exact off-face dual identities for all balanced gauge vectors."""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from check_rectangular_uniform_entry_bounds import verify_supported_duals


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    model_path = root/'evidence/rectangular-full-leading-face-lp.json'
    witness_path = root/'evidence/rectangular-uniform-entry-inputs.json'
    model = json.loads(model_path.read_text())
    witnesses = json.loads(witness_path.read_text())['witnesses']
    duals = verify_supported_duals(model, witnesses)
    table = duals['leading_C_lower_bounds']
    for x in range(4):
        for y in range(3):
            assert Fraction(table[x][y]) >= (0 if (x, y) == (1, 1) else 2)
    gaps = []
    for certificate in model['C0_physical_gaps']:
        equality = dict(certificate['dual_equalities'])
        assert Fraction(equality[2]) == -1
        mass = sum(abs(Fraction(value)) for i, value in certificate['dual_inequalities']
                   if i < len(model['constraints']))
        assert mass == 2
        assert Fraction(certificate['optimum']) >= 2
        gaps.append({'letter': certificate['letter'], 'state': certificate['state'],
                     'constant_gap': certificate['optimum'],
                     'off_face_C_multiplier': '-1', 'entry_dual_mass': '2'})
    infinity = [(-2, 1), (-2, 3), (0, -3), (0, -1), (0, 1),
                (0, 3), (0, 5), (2, -5), (2, -3), (2, -1)]
    zero = [(0, 1), (0, 3), (0, 5), (2, -5), (2, -3), (2, -1),
            (0, 1), (0, 3), (2, -3), (2, -1)]
    negatives = []
    for a, b in infinity+zero:
        at_eight = 8*a+b
        assert at_eight
        assert not a or (a > 0) == (at_eight > 0)
        if at_eight < 0:
            negatives.append((a, b))
    assert tuple(map(sum, zip(*negatives))) == (-4, 0)
    assert len(infinity)+len(zero) == 20
    result = {'status': 'Exact dual and Gaussian constants support the uniform balanced-gauge theorem.',
              'scope': 'The baseline mixed rectangular partial repair; tau >= c log M with c > 1/sqrt(2).',
              'leading_C_lower_bounds': table, 'off_face_gap_certificates': gaps,
              'active_Gaussian_factors': 20, 'sum_negative_Gaussian_orders': '-4M',
              'physical_core_cardinalities': [4, 3],
              'coarse_pressure_constant': 'K0 = 20 log 2 + sqrt(2) log 12',
              'entry_offset_constant': 'Any B >= 0 satisfying all 28 eventual bounds v(entry) <= a M + B.',
              'coarse_pressure_bound': 'F >= T - (2B+8) ell - K0',
              'nonpositive_pressure_gap_bound': 'gap >= 2 tau - K0 - (2B+12) ell',
              'model_sha256': sha256(model_path.read_bytes()).hexdigest(),
              'witness_sha256': sha256(witness_path.read_bytes()).hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
