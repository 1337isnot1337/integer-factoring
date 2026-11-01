#!/usr/bin/env python3
"""Independent dense controls for batched proper rational contractions."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random
from time import perf_counter

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_pole_polynomials import inverse_pole_polynomial_contraction, inverse_residue_block_batched
from pole_cauchy_residue_blocks import inverse_residue_block


ROOT = Path(__file__).resolve().parents[1]


def quad(value):
    return Quad(Fraction(value)) if isinstance(value, int) else Quad(*(Fraction(x) for x in value))


def rational_weight(a, pole, coefficients):
    inverse_value, power, result = 1/(a-quad(pole)), Quad(Fraction(1)), Quad()
    for coefficient in coefficients:
        power *= inverse_value
        result += quad(coefficient)*power
    return result


def check_scaled(actual, expected):
    if actual['working_precision'] == 0:
        assert expected.valuation() is None or expected.valuation() >= actual['absolute_precision']
        return
    scale = actual['scale_exponent']
    normalized = expected*(1 << -scale) if scale <= 0 else expected/(1 << scale)
    assert tuple(actual['residue']) == mod_pair(normalized, actual['working_precision'])
    assert actual['working_precision']+scale == actual['absolute_precision']


def controls():
    rng = random.Random(2026110115)
    contractions = blocks = old_agreements = 0
    examples = []
    for index in range(72):
        size = 1+index % 4
        left = tuple((a, 1) for a in sorted(rng.sample(list(range(1, 20, 2)), size)))
        right = ((1+2*(index % 3), size),)
        aa, bb = expand(left), expand(right)
        n = rng.randrange(31, 220, 2)
        k = max(abs(n-a*b) for a in aa for b in bb).bit_length()+1
        dense = inverse(matrix_at(n, aa, bb, k))
        c, d = (index % 5, -(1 << (index % 4))), ((index+1) % 5, -(1 << ((index+2) % 4)))
        alpha = [(rng.randrange(-4, 5), rng.randrange(-4, 5)) for _ in range(1+index % 5)]
        beta = [(rng.randrange(-4, 5), rng.randrange(-4, 5)) for _ in range(1+(index//5) % 4)]
        alpha[0] = (1, index % 3)
        beta[0] = (-2, 1)
        expected = sum((dense[j][i]*rational_weight(a, c, alpha)*rational_weight(b, d, beta)
                        for i, a in enumerate(aa) for j, b in enumerate(bb)), Quad())
        result = inverse_pole_polynomial_contraction(n, left, right, k, c, d, alpha, beta, 5)
        check_scaled(result, expected)
        assert result['assembled_product_jets'] <= 2
        contractions += 1
        if index < 24:
            precision = 1+index % 4
            row_bits, column_bits = index % 5, (index//5) % 4
            row_residue, column_residue = rng.randrange(16), rng.randrange(16)
            actual = inverse_residue_block_batched(n, left, right, k, row_residue, row_bits,
                                                    column_residue, column_bits, precision)
            literal = sum((dense[j][i] for i, a in enumerate(aa) for j, b in enumerate(bb)
                           if (a-row_residue) % (1 << row_bits) == 0
                           and (b-column_residue) % (1 << column_bits) == 0), Quad())
            literal *= 1 << actual['inverse_normalization_exponent']
            assert tuple(actual['residue']) == mod_pair(literal, precision)
            old = inverse_residue_block(n, left, right, k, row_residue, row_bits,
                                        column_residue, column_bits, precision)
            assert old['residue'] == actual['residue']
            blocks += 1
            old_agreements += 1
            examples.append({'n': n, 'size': size, 'precision': precision,
                             'residue': actual['residue'],
                             'new_working_precision': actual['batched_contraction']['working_precision'],
                             'old_maximum_working_precision': max(v['working_precision']
                                                                  for v in old['contraction_budgets']),
                             'new_assembled_jets': actual['batched_contraction']['assembled_product_jets'],
                             'old_repeated_pole_calls': old['repeated_pole_contractions']})
    # A high positive Q valuation omits both unit jets exactly.
    n, left, right, k = 15, ((1, 1),), ((3, 1),), 5
    c, d = 0, (15+(1 << 100), 32)
    actual = inverse_pole_polynomial_contraction(n, left, right, k, c, d, [1], [1], 5)
    expected = inverse(matrix_at(n, (1,), (3,), k))[0][0]/((1-quad(c))*(3-quad(d)))
    check_scaled(actual, expected)
    assert actual['assembled_product_jets'] == 0 and actual['product_unit_precision'] == 0
    positive_q = actual['product_valuation']
    # Factoring coefficient valuations can make the entire answer vanish.
    zero = inverse_pole_polynomial_contraction(n, left, right, k, (0, -1), (0, -1),
                                                [1 << 100], [1], 5)
    assert zero['working_precision'] == 0 and zero['residue'] == [0, 0]
    rejected = 0
    for alpha, beta, precision in (([], [1], 5), ([0], [1], 5), ([1], [(0, 0)], 5),
                                    ([1], [1], 0), ([(1, 2, 3)], [1], 5)):
        try:
            inverse_pole_polynomial_contraction(n, left, right, k, -1, -3, alpha, beta, precision)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('Invalid polynomial contract accepted')
    high_bits = []
    for a, b in ((1, 1), (3, 5), (5, 3)):
        value = inverse_residue_block_batched(15, ((1, 3),), ((1, 3),), 5, a, 3, b, 3, 5)
        assert value['residue'] == ([25, 25] if a == b == 1 else [16, 0])
        high_bits.append({'a': a, 'b': b, 'residue_mod_32': value['residue']})
    return {'dense_quadratic_coefficient_contractions': contractions,
            'dense_residue_blocks': blocks, 'independent_termwise_agreements': old_agreements,
            'positive_product_valuation_skip_control': positive_q,
            'whole_answer_vanishing_control': True, 'invalid_contract_rejections': rejected,
            'high_bit_controls': high_bits, 'precision_comparisons': examples}


def large_controls():
    n, size = (1 << 401)+45, (1 << 200)+13
    left, right = ((1, size),), ((3, size),)
    k = max(abs(n-3), abs(n-(2*size-1)*(2*size+1))).bit_length()+12
    records = []
    for row in (5, 7):
        previous = None
        for precision in (1, 5, 9):
            started = perf_counter()
            current = inverse_residue_block_batched(n, left, right, k, row, 7, 11, 9, precision)
            elapsed = perf_counter()-started
            if previous is not None:
                assert previous['residue'] == [x % (1 << previous['precision']) for x in current['residue']]
            records.append({'row_residue': row, 'precision': precision, 'seconds': elapsed,
                            'result': current})
            previous = current
    return {'n': n, 'input_bits': n.bit_length(), 'domain_size': size, 'k': k,
            'records': records, 'enumerated_matrix_nodes': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--small-only', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'claim': 'Batched proper rational weights and tight coefficient precision',
              'controls': controls(), 'factoring_exponent_improved': False}
    if not args.small_only:
        report['implicit_large_controls'] = large_controls()
    names = ('tools/check_pole_cauchy_pole_polynomials.py', 'tools/pole_cauchy_pole_polynomials.py',
             'tools/pole_cauchy_quadratic_selectors.py', 'tools/pole_cauchy_residue_blocks.py',
             'tools/modular_interval_jets.py')
    report['source_sha256'] = {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report['controls'].items()
                      if k not in ('precision_comparisons', 'high_bit_controls')}, indent=2))
    if 'implicit_large_controls' in report:
        print(json.dumps([{'row': item['row_residue'], 'precision': item['precision'],
                           'residue': item['result']['residue'], 'seconds': item['seconds'],
                           'working_precision': item['result']['batched_contraction']['working_precision']}
                          for item in report['implicit_large_controls']['records']], indent=2))


if __name__ == '__main__':
    main()
