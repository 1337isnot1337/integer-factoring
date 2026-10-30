#!/usr/bin/env python3
"""Independent dense checks for the two-patch trace crossing."""
import argparse
import hashlib
import json
import os
from pathlib import Path

from check_trace_bilinear_kernel import Field, irreducibility_by_division, rank
from two_patch_trace_probe import profile


def dense_check(t, polynomial):
    field = Field(polynomial)
    d = 1 << t
    powers = [field.power(field.alpha, j) for j in range(t)]
    coordinates = [sum(field.trace(field.mul(a, z)) << j
                       for j, a in enumerate(powers)) for z in range(d)]
    rows = [[0] * (d*d) for _ in range(2*t+1)]
    columns = [[0] * (d*d) for _ in range(2*t+1)]
    for x in range(d):
        for y in range(d):
            for u in range(d):
                v = x ^ y ^ u
                s, w = x ^ y, x ^ u
                error = s ^ w
                b0 = field.mul(s, error)
                b1 = field.mul(s ^ 1, error)
                polynomial_bits = 1
                # Multiply the linear factors directly, independently of the
                # probe's Lucas-theorem coefficient enumeration.
                for value in (coordinates[b0], coordinates[b1]):
                    for j in range(t):
                        if value >> j & 1:
                            polynomial_bits ^= polynomial_bits << 1
                assert (polynomial_bits.bit_count() & 1) == int(s == w)
                for k in range(2*t+1):
                    if polynomial_bits >> k & 1:
                        rows[k][x*d+y] |= 1 << (u*d+v)
                        columns[k][x*d+u] |= 1 << (y*d+v)
    row_ranks = [rank(matrix) for matrix in rows]
    column_ranks = [rank(matrix) for matrix in columns]
    expected = profile(t, polynomial)
    assert row_ranks == expected['row_ranks']
    assert column_ranks == expected['column_ranks']
    return {'t': t, 'field_polynomial': polynomial,
            'coefficient_cuts_checked': 2*(2*t+1),
            'endpoint_entries_checked': d**3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    profiles = []
    dense_checks = []
    for t, polynomial in ((2, 7), (3, 11), (4, 19), (5, 37), (6, 67)):
        irreducibility_by_division(polynomial)
        result = profile(t, polynomial)
        assert result['rank_sum_floor'] > (1 << (2*t))
        profiles.append(result)
        if t <= 4:
            dense_checks.append(dense_check(t, polynomial))
    here = Path(__file__).resolve().parent
    sources = [Path(__file__), here/'two_patch_trace_probe.py',
               here/'check_trace_bilinear_kernel.py']
    evidence = {'status': 'passed', 'profiles': profiles,
                'dense_checks': dense_checks,
                'coefficient_cuts_checked': sum(x['coefficient_cuts_checked']
                                               for x in dense_checks),
                'endpoint_entries_checked': sum(x['endpoint_entries_checked']
                                                for x in dense_checks),
                'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sources},
                'scope': 'Exact endpoint and finite ranks; all five tested raw rank sums exceed d^2. No asymptotic barrier or score improvement.'}
    rendered = json.dumps(evidence, indent=2, sort_keys=True) + '\n'
    temporary = args.output.with_suffix(args.output.suffix + '.tmp')
    with temporary.open('w') as out:
        out.write(rendered)
        out.flush()
        os.fsync(out.fileno())
    temporary.replace(args.output)
    saved = json.loads(args.output.read_text())
    assert len(saved['profiles']) == 5 and len(saved['dense_checks']) == 3
    print(json.dumps({'status': 'passed', 'profiles': 5,
                      'dense_coefficient_cuts': evidence['coefficient_cuts_checked'],
                      'endpoint_entries': evidence['endpoint_entries_checked']}))


if __name__ == '__main__':
    main()
