#!/usr/bin/env python3
"""Check the conservative norm parameters against exact density integers."""

from hashlib import sha256
import json
from pathlib import Path

from powered_prime_norm_parameters import powered_prime_parameters
from quasipolynomial_norm_parameters import (
    divisor_table_exponent, quasipolynomial_norm_parameters,
)


def divisor_controls():
    checked = 0
    for two in range(30):
        for three in range(18):
            for five in range(9):
                n = 2**two*3**three*5**five
                b = max(1, n.bit_length())
                tau = (two+1)*(three+1)*(five+1)
                bound = divisor_table_exponent(b)
                assert tau <= 1 << bound['divisor_count_exponent_bound']
                assert (b+1)*tau <= 1 << (2*bound['table_size_exponent_bound'])
                checked += 1
    return checked


def main():
    moderate = []
    for b in [2**j for j in range(14, 26)]+[2**17+33, 2**23-11]:
        current = quasipolynomial_norm_parameters(b)
        assert current['eligible']
        exact = powered_prime_parameters(b, current['smoothness_exponent'])
        assert exact['eligible']
        upper = current['density_bias_bits_upper_bound']
        assert 4*exact['density_denominator'] <= 1 << upper
        assert current['walk_steps'] >= exact['walk_steps']
        assert current['walk_steps']*current['sample_bias_bits'] >= upper
        assert current['ordered_walk_enumeration_exponent_bound'] == (current['smoothness_exponent']+1)*current['walk_steps']
        assert current['accepted_norm_log_two_bound'] >= exact['accepted_norm_log_two_bound']
        moderate.append(current)

    large = []
    for exponent in [32, 64, 128, 256, 512, 1024, 2048, 4096]:
        current = quasipolynomial_norm_parameters(1 << exponent)
        assert current['eligible']
        assert current['logarithm_ceiling'] == exponent
        assert current['sample_bias_bits'] > 0
        assert current['walk_steps']*current['sample_bias_bits'] >= current['density_bias_bits_upper_bound']
        current['input_bit_length_parameter_log_two'] = exponent
        large.append(current)

    ineligible = [quasipolynomial_norm_parameters(b) for b in [2, 16, 256, 1024, 8192]]
    assert all(not row['eligible'] for row in ineligible)
    for bad in [0, 1, -1, '16', 1.5]:
        try:
            quasipolynomial_norm_parameters(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Malformed bit-length parameter was accepted.')
    root = Path(__file__).resolve().parents[1]
    paths = ['tools/quasipolynomial_norm_parameters.py',
             'tools/check_quasipolynomial_norm_parameters.py',
             'tools/powered_prime_norm_parameters.py']
    report = {
        'status': 'pass',
        'scope': 'Exact existence budgets and processing bounds only; no smooth-norm selector.',
        'materialized_density_comparisons': moderate,
        'large_symbolic_parameters': large,
        'ineligible_finite_controls': ineligible,
        'divisor_count_controls': divisor_controls(),
        'malformed_input_controls': 5,
        'source_sha256': {path: sha256((root/path).read_bytes()).hexdigest() for path in paths},
    }
    output = root/'evidence/quasipolynomial-norm-parameters.json'
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'pass', 'exact_density_controls': len(moderate),
                     'large_symbolic_controls': len(large),
                     'divisor_count_controls': report['divisor_count_controls'],
                     'output': str(output)}))


if __name__ == '__main__':
    main()
