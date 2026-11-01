#!/usr/bin/env python3
"""Exact all-scale winner among six F4 five-pair height-20 orderings.

The saved zero profiles were computed directly. Their infinity profiles
were derived from the orthogonal identity. This checker only replays the
integer comparison certificates; it does not recompute Smith profiles.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from check_tensor_coupled_quadratic_dominance import positive_difference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles', type=Path,
                        default=Path('evidence/tensor-f4-five-pair-height20-catalogue.json'))
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-f4-five-pair-dominance.json'))
    args = parser.parse_args()
    source = json.loads(args.profiles.read_text())
    assert source['completed'] and source['in_progress'] is None
    assert source['height_cap'] == 20 and source['modes'] == 14
    records = source['records']
    assert len(records) == 6 and {r['index'] for r in records} == set(range(6))
    reference = next(record for record in records if record['index'] == 0)
    assert reference['coefficients'] == [[1, 0, 1, 2, 2], [0, 1, 2, 1, 2]]
    physical = Counter(reference['profile']['heights'])
    base = Counter(reference['profile']['zero_orders_t'])
    certificates = []
    for record in records:
        profile = record['profile']
        assert profile['q'] == 54
        assert len(profile['zero_orders_t']) == len(profile['infinity_orders_t']) == 54 * 54
        assert profile['zero_orders_t'] == profile['infinity_orders_t']
        assert profile['computed_cuts'] == ['zero']
        assert profile['computed_smith_orders'] == 54 * 54
        assert Counter(profile['heights']) == physical
        certificate = positive_difference(Counter(profile['zero_orders_t']), base)
        assert certificate['comparison'] == ('equal' if record['index'] == 0 else 'strictly_larger')
        if record['index']:
            assert certificate['nonnegative_basis'] == 'bernstein'
            assert certificate['one_minus_power'] == 1
        certificates.append({'index': record['index'], 'coefficients': record['coefficients'],
                             **certificate})
    report = {'passed': True, 'scope': __doc__, 'reference_index': 0,
              'reference_coefficients': reference['coefficients'], 'profile_classes': 6,
              'direct_smith_orders_in_source': 6 * 54 * 54,
              'derived_infinity_orders_in_source': 6 * 54 * 54,
              'smith_orders_recomputed_by_this_checker': 0,
              'input_path': args.profiles.as_posix(),
              'input_sha256': sha256(args.profiles.read_bytes()).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'check_tensor_coupled_quadratic_dominance.py')},
              'certificates': certificates}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'reference_index': 0, 'strict_comparisons': 5,
                      'bernstein_degrees': [c['quotient_degree'] for c in certificates if c['index']],
                      'output': str(args.output)}))


if __name__ == '__main__':
    main()
