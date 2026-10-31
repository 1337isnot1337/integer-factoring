#!/usr/bin/env python3
"""Reconstruct every retained joint-filter witness from its exact bounds.

This checks the portable catalogue, rather than rerunning the much larger
Pareto search.  Feasibility algorithms have separate independent controls.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from joint_binary_filter_synthesis import synthesize as binary_synthesize
from joint_prime_filter_synthesis import synthesize as prime_synthesize


CATALOGUES = (
    'joint-binary-polynomial-n32-r3.json',
    'joint-binary-polynomial-n32-r4.json',
    'joint-binary-polynomial-n32-p3-r3.json',
    'joint-binary-polynomial-n32-r3-offset64.json',
    'joint-binary-polynomial-n32-r4-offset32.json',
    'joint-binary-polynomial-n32-r3-individual-offset32.json',
    'joint-qutrit-polynomial-n27-r2.json',
    'joint-qutrit-polynomial-n27-r2-offset27.json',
)


def digest(value):
    return sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'joint-filter-catalogue-replay.json')
    args = parser.parse_args()
    records = []
    total = endpoints = 0
    for name in CATALOGUES:
        data = json.loads((root / 'evidence' / name).read_text())
        rows = data['catalogue']
        assert digest(rows) == data['catalogue_sha256']
        assert len(rows) == data['distinct_ordered_beta_degree_profiles']
        assert sum(row['multiplicity'] for row in rows) == data['distinct_feasible_bound_pairs']
        seen = set()
        prime, modes, m = data.get('prime', 2), data['modes'], data['m']
        assert m == data['n'] * data['ramification']
        for row in rows:
            if prime == 2:
                witness = binary_synthesize(modes, m, row['source_degree_bounds'],
                                            row['source_beta_bounds'])
            else:
                witness = prime_synthesize(prime, modes, m, row['source_degree_bounds'],
                                           row['source_beta_bounds'])
            assert witness
            assert list(witness['actual_orders']) == row['beta_orders']
            assert list(witness['actual_degrees']) == row['degrees']
            assert list(witness['endpoint_values']) == row['endpoint_values']
            assert witness['endpoint_image_dimension'] == row['endpoint_image_dimension']
            assert witness['coefficient_field_order'] == row['coefficient_field_order']
            literal = [{'state': state, 'polynomial': sorted(polynomial.items())}
                       for state, polynomial in witness['alpha'].items()]
            assert digest(literal) == row['witness_polynomial_sha256']
            key = tuple(row['beta_orders']), tuple(row['degrees'])
            assert key not in seen
            seen.add(key)
            total += 1
            endpoints += len(row['endpoint_values'])
        records.append({'catalogue': name, 'literal_witnesses': len(rows),
                        'catalogue_sha256': data['catalogue_sha256']})
    result = {'status': 'Every retained joint-filter witness and catalogue hash replayed exactly.',
              'catalogues': records, 'literal_polynomial_witnesses': total,
              'nonzero_endpoint_values': endpoints}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'catalogues'}))


if __name__ == '__main__':
    main()
