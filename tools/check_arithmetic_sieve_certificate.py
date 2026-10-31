#!/usr/bin/env python3
"""Replay the stored large relation matrix, determinant and torsion certificate."""
import json
from math import isqrt
from pathlib import Path

from arithmetic_norm_sieve import find_pair, prime_base, roots_x, verify
from arithmetic_two_group import power
from binary_quadratic_forms import compose_forms, inverse_form, reduced_forms
from smooth_norm_relation import norm_to_prime_word


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
        a, b, c = target
        x_bound, y_bound = control['x_bound'], control['y_bound']
        radius = isqrt(4*integer-1)//2
        maximum_norm = a*x_bound*x_bound+abs(b)*x_bound*y_bound+c*y_bound*y_bound
        # Here every smooth norm would pass greedy balancing. Exhaustion
        # therefore concerns norm search itself, not the pair output format.
        assert bound*maximum_norm <= radius*radius
        pair, _, _ = find_pair(target, base, bound, isqrt(4*integer-1)//2,
                               control['x_bound'], control['y_bound'])
        assert pair is None
        x, y = control['outside_box_witness_vector']
        assert y > y_bound
        word = norm_to_prime_word(target, x, y, bound, radius*radius)
        assert word is not None and str(word['norm']) == control['outside_box_primitive_norm']
        assert word['word_length'] == control['outside_box_word_length']
        value = (1, 0, integer)
        for entry in word['word']:
            form, exponent = entry['form'], entry['exponent']
            value = compose_forms(value, power(form if exponent >= 0 else inverse_form(form),
                                               abs(exponent), (1, 0, integer), compose_forms))
        assert value == target
    result['replayed_fixed_box_exhaustions'] = len(failures)
    print(json.dumps(result, indent=2))
