#!/usr/bin/env python3
"""Replay the stored large relation matrix, determinant and torsion certificate."""
import json
from math import isqrt
from pathlib import Path

from arithmetic_norm_sieve import find_pair, prime_base, roots_x, verify
from arithmetic_two_group import power
from binary_quadratic_forms import compose_forms, reduced_forms


def check_roots():
    controls = 0
    for discriminant in (-23, -84, -348, -66548):
        base = prime_base(discriminant, 23)
        for a, b, c in reduced_forms(discriminant):
            for y in range(4):
                for prime, root, _ in base:
                    expected = {x for x in range(prime) if (a*x*x+b*x*y+c*y*y) % prime == 0}
                    assert set(roots_x((a, b, c), y, prime, root)) == expected
                    controls += 1
    return controls


if __name__ == '__main__':
    evidence = Path(__file__).resolve().parents[1]/'evidence'
    result = {'small_residue_controls': check_roots(), 'certificates': []}
    for bits in (60, 80):
        certificate = json.loads((evidence/f'arithmetic-norm-sieve-{bits}bit.json').read_text())
        assert certificate['status'] == 'factor_found'
        result['certificates'].append(verify(certificate))
    failures = json.loads((evidence/'arithmetic-norm-search-exhaustion.json').read_text())['controls']
    for control in failures:
        integer, bound = int(control['integer']), control['smooth_bound']
        base = prime_base(-4*integer, bound)
        index = control['failed_generator_index']
        assert base[index][0] == control['generator_prime']
        target = power(base[index][2], 1 << (2*integer.bit_length()), (1, 0, integer), compose_forms)
        assert [str(value) for value in target] == control['target_form']
        pair, _, _ = find_pair(target, base, bound, isqrt(4*integer-1)//2,
                               control['x_bound'], control['y_bound'])
        assert pair is None
    result['replayed_fixed_box_exhaustions'] = len(failures)
    print(json.dumps(result, indent=2))
