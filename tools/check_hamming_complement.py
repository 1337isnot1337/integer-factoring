#!/usr/bin/env python3
"""Full inverse and complementary-minor checks for weighted crossing profiles."""

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_modular_recombination import insert, represent
from hamming_pfaffian import folded_matrix, profile
from hamming_sectors import complement_scalar, sector_labels
from check_modular_particle_hole import complement, power
from modular_pfaffian import polynomial_exact_quotient, polynomial_pfaffian, polynomial_product

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(1729)
    inverse_entries = 0
    core_records = []
    for bits in range(0, 6 if args.full else 5):
        k = 1 << bits
        for weights in [(1,)*(bits+1), (1,)+tuple(1+rng.randrange(3) for _ in range(bits))]:
            core = tuple(range(1, 2*k, 2))
            a = folded_matrix(core, core, core, core, weights)
            size = 4*k
            scalar = complement_scalar(weights)
            assert scalar.bit_count()%2 == 1
            reflection = [(i//k)*k+k-1-i%k for i in range(size)]
            for i, j in product(range(size), repeat=2):
                value = 0
                for t in range(size):
                    value ^= polynomial_product(a[i][t], a[reflection[t]][reflection[j]])
                assert value == (scalar if i == j else 0)
                inverse_entries += 1
            pfaffian_checked = k <= 8
            if pfaffian_checked:
                assert polynomial_pfaffian(a) == power(scalar, k)
            core_records.append({'core_modes': k, 'weights': weights, 'scalar': scalar,
                                 'full_pfaffian_checked': pfaffian_checked})
    profiles = 0
    span_comparisons = 0
    omitted_factor_fails = False
    for free, count in [(1, 1), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)]:
        k = 1 << (free+1)
        domain = sector_labels(free, count, 4 if (free, count) != (1, 1) else 2)
        for weights in [(1,)*(free+2), (1,)+tuple(1+rng.randrange(3) for _ in range(free+1))]:
            transformed = tuple(complement(label, k) for label in domain)
            scalar = complement_scalar(weights)
            exponent = k-4*count
            factor = power(scalar, abs(exponent))
            indices = list(product(range(len(domain)), repeat=4))
            if k >= 16 and not args.full:
                rng.shuffle(indices)
                indices = indices[:24]
            old_profiles, new_profiles = [], []
            for ix in indices:
                x, y, u, v = (domain[i] for i in ix)
                cx, cy, cu, cv = (transformed[i] for i in ix)
                old, new = profile(x, v, u, y, weights), profile(cx, cv, cu, cy, weights)
                if exponent >= 0:
                    assert new == polynomial_product(factor, old)
                    assert polynomial_exact_quotient(new, factor) == old
                else:
                    assert old == polynomial_product(factor, new)
                    assert polynomial_exact_quotient(old, factor) == new
                assert old.bit_count()%2 == new.bit_count()%2 == int(ix[0] == ix[3] and ix[1] == ix[2])
                omitted_factor_fails |= exponent != 0 and old != new
                profiles += 1
                old_profiles.append(old)
                new_profiles.append(new)
            def coefficient_basis(values):
                basis = {}
                for degree in range(max(values).bit_length()):
                    vector = sum(((value >> degree) & 1) << i for i, value in enumerate(values))
                    insert(basis, vector, 0)
                return basis
            old_basis, new_basis = coefficient_basis(old_profiles), coefficient_basis(new_profiles)
            assert len(old_basis) == len(new_basis)
            assert all(represent(old_basis, v) is not None for v, m in new_basis.values())
            assert all(represent(new_basis, v) is not None for v, m in old_basis.values())
            span_comparisons += 1
    assert omitted_factor_fails
    paths = [Path(__file__), ROOT/'tools/hamming_sectors.py', ROOT/'tools/hamming_pfaffian.py']
    result = {'status': 'weighted complement identity; no rank monotonicity asserted',
              'full': args.full, 'inverse_entries': inverse_entries,
              'complete_profile_pairs': profiles, 'coefficient_span_comparisons': span_comparisons,
              'core_cases': core_records,
              'omitted_factor_negative_control': omitted_factor_fails,
              'sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('core_cases', 'sha256')}))


if __name__ == '__main__':
    main()
