#!/usr/bin/env python3
"""Two exact direct sums control every polynomial energy-gauge tangent."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_controlled_permutation_smith import binary_rank
from check_joint_leading_tangent import certificate as constant_certificate, realign
from joint_energy_conjugations import N, Space, binary_product, model

ROOT = Path(__file__).resolve().parents[1]


def leading_minor_derivative(matrix, leading):
    base, delta = realign(leading), realign(matrix)
    mask = (1 << N)-1
    a = [(base >> (N*r)) & 1 for r in range(N)]
    b = base & mask
    rows = [(delta >> (N*r)) & mask for r in range(N)]
    assert base & 1
    result = 0
    for r in range(N):
        for c in range(N):
            value = (((rows[r] >> c) & 1)
                     ^ (((b >> c) & 1)*(rows[r] & 1))
                     ^ (a[r]*((rows[0] >> c) & 1))
                     ^ (a[r]*((b >> c) & 1)*(rows[0] & 1)))
            result ^= value << (N*r+c)
    return result


def certificate():
    data = model()
    constant = constant_certificate()
    c0, c1 = data['weighted'][-8], data['weighted'][-7]
    commutator = lambda a, b: binary_product(a, b) ^ binary_product(b, a)
    first = [commutator(1 << (N*r+c), c0) for r, c in data['positions']]
    second = [commutator(a, c1) for a in data['leading'].basis]
    first_tangent = [leading_minor_derivative(a, c0) for a in first]
    second_tangent = [leading_minor_derivative(a, c0) for a in second]
    assert binary_rank(first) == 59 and binary_rank(second) == 6
    assert binary_rank(first+second) == 65
    assert binary_rank(first_tangent) == 53 and binary_rank(second_tangent) == 6
    assert binary_rank(first_tangent+second_tangent) == 59
    assert len(data['leading'].basis) == 10 and len(data['full'].basis) == 4
    # Kernel of ad(C_1) on L is exactly the complete reference commutant.
    assert all(commutator(a, c0) == commutator(a, c1) == 0 for a in data['full'].basis)
    assert binary_rank(data['leading'].basis+data['full'].basis) == 10
    assert constant['tangent_dimension'] == 16
    return dict(energy_dimension=69, leading_commutant_dimension=10,
                full_commutant_dimension=4, constant_tangent_dimension=16,
                ordinary_image_ranks=[59, 6, 65], projected_image_ranks=[53, 6, 59],
                both_image_intersections_zero=True,
                next_commutator_basis=Space(second).basis,
                full_commutant_basis=data['full'].basis,
                all_degree_tangent_dimension='16+4D for polynomials of degree at most D',
                tangent_conclusion='Every positive-degree coefficient commutes with the complete reference; the constant coefficient lies in the physical-diagonal plus leading-centralizer tangent space.',
                formal_conclusion='Every identity-based formal polynomial energy gauge retaining the pole bound and leading rank is a constant energy gauge times a formal polynomial full-commutant unit.',
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_polynomial_tangent.py', 'check_joint_leading_tangent.py',
                                            'joint_energy_conjugations.py', 'search_f4_recentered_crossings.py')},
                scope='All characteristic-two fields and all polynomial degrees. Formal lifting is proved in the note; separate components of the gauge constraints are not classified.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-polynomial-tangent.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_sha256', 'next_commutator_basis', 'full_commutant_basis')}))


if __name__ == '__main__':
    main()
