#!/usr/bin/env python3
"""Complete closed-point profiles of paired subalphabet restrictions.

The finite factorization is performed over F16 and descended by Frobenius
to F4; every descended factor is independently checked irreducible over
F4. All cut determinants are reconstructed exactly from those factors.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

from affine_place_profile_oracle import polynomial
from check_fused_closed_profiles import (
    add, canonical, from_canonical, gcd, irreducible_factors, multiply,
    polynomial_modulo, power_modulo, pscale,
)
from f16_identity_operator_profiles import INV
from f16_weighted_local_oracle import F16WeightedLocalOracle
from fused_six_vertex_crossings import power
from physical_split_crossings import flatten
from reflection_laurent_physical_frames import digest
from search_finite_fiber_restrictions import restriction
from search_fused_mixed_weighted import atomic_write


def factors_over_f4(value):
    factors = Counter(canonical(factor) for factor in irreducible_factors(value))
    result = []
    while factors:
        key, multiplicity = next(iter(factors.items()))
        del factors[key]
        conjugate = tuple(power(coefficient, 4) for coefficient in key)
        factor = from_canonical(key)
        if conjugate != key:
            assert factors.pop(conjugate) == multiplicity
            factor = multiply(factor, from_canonical(conjugate))
        assert set(factor.values()) <= set(range(4))
        degree = max(factor)
        x = polynomial_modulo({1: 1}, factor)
        current = x
        for j in range(1, degree + 1):
            current = power_modulo(current, 4, factor)
            if j <= degree // 2:
                assert max(gcd(factor, add(current, x))) == 0
        assert current == x
        result.append((canonical(factor), multiplicity))
    reconstructed = {0: value[max(value)]}
    for key, multiplicity in result:
        for _ in range(multiplicity):
            reconstructed = multiply(reconstructed, from_canonical(key))
    assert reconstructed == value
    return sorted(result, key=lambda item: (len(item[0]), item[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('evidence/deflated-finite-fiber-input.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    inputs = json.loads(raw)
    subsets = [list(value) for size in (2, 3) for value in combinations(range(4), size)]
    names = (Path(__file__).name, 'search_finite_fiber_restrictions.py',
             'physical_split_crossings.py', 'affine_place_profile_oracle.py',
             'check_fused_closed_profiles.py', 'f16_weighted_local_oracle.py',
             'f16_identity_operator_profiles.py', 'fused_six_vertex_crossings.py')
    report = {'scope': __doc__, 'input_sha256': sha256(raw).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'records': [], 'completed': False}
    started = time.monotonic()
    for first in subsets:
        for second in subsets:
            for source_index, case in enumerate(inputs['cases']):
                seed, content = restriction(case['seed'], first, second)
                n = len(first) * len(second)
                matrices = [[[polynomial(v) for v in row] for row in flatten(seed, cut)[0]]
                            for cut in (0, 1)]
                oracles = [F16WeightedLocalOracle(matrix, seed['degree']) for matrix in matrices]
                zero = [oracle.profile([0] * n, [0] * n, (0, 1)) for oracle in oracles]
                factors = [factors_over_f4(profile['determinant']) for profile in zero]
                keys = sorted({key for factor_list in factors for key, _ in factor_list},
                              key=lambda key: (len(key), key))
                local = []
                for key in keys:
                    orders = [list(oracle.profile([0] * n, [0] * n, key)['orders'])
                              for oracle in oracles]
                    local.append({'factor': list(key), 'residue_degree': len(key) - 1,
                                  'orders': orders})
                infinity = [list(oracle.profile([0] * n, [0] * n, None)['orders'])
                            for oracle in oracles]
                for cut in (0, 1):
                    assert sum(infinity[cut]) == n * seed['degree'] - max(zero[cut]['determinant'])
                    multiplicities = dict(factors[cut])
                    for entry in local:
                        assert sum(entry['orders'][cut]) == multiplicities.get(tuple(entry['factor']), 0)
                record = {'source_index': source_index, 'source_label': case['label'],
                          'subset_x': first, 'subset_y': second, 'degree': seed['degree'],
                          'tensor_sha256': digest(seed['values']), 'scalar_content': content,
                          'determinants': [profile['determinant'] for profile in zero],
                          'factorizations': [[{'factor': list(key), 'multiplicity': multiplicity}
                                               for key, multiplicity in family] for family in factors],
                          'closed_points': local, 'homogeneous_infinity_orders': infinity}
                report['records'].append(record)
                report['elapsed_seconds'] = time.monotonic() - started
                atomic_write(args.output, report)
                if len(report['records']) % 20 == 0:
                    print(json.dumps({'completed_records': len(report['records']),
                                      'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic() - started)
    atomic_write(args.output, report)
    nonzero = [row for row in report['records']
               if any(p['factor'] != [0, 1] for p in row['closed_points'])]
    print(json.dumps({'completed_records': len(report['records']),
                      'with_nonzero_singular_places': len(nonzero),
                      'elapsed_seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()
