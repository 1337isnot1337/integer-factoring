#!/usr/bin/env python3
"""Independent rational entropy witnesses for bare alphabets seven to nine."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from check_cell_minor_caps import Dyadic, log_at_least_one, point
from subalphabet_adjugate_caps import linear_caps


def independent_constraint(q, label):
    """Construct the displayed bare formulas without the general cap builder."""
    total = q*q
    if label[0] == 'prefix':
        assert label == ['prefix', q]
        return [F(i < q) for i in range(total)], F(q, 2)
    family, p, b, k = label
    n = p*b
    values = [F(0)]*total
    if family == 'square cofactor cell':
        assert p == b == q and k == n
        values[:n-1] = [F(n)]*(n-1)
        values[1] += n-1
    else:
        assert family in ('square cofactor shared line', 'rectangular cofactor shared line')
        assert b == q and k == q
        values[:n] = [F(k-(i >= n-k)) for i in range(n)]
        if p == b:
            assert family == 'square cofactor shared line'
            values[1] += k-1
        else:
            assert family == 'rectangular cofactor shared line'
            gap = n-1-p*p
            assert gap > 0
            values[1] += k*gap
    return values, F(k*(n-1))


def replay(root, entry):
    q, rho = entry['alphabet'], F(entry['rho_floor'])
    n, coefficients, total = q*q, [F(0)]*(q*q), F(0)
    reference = {(cap['family'], *cap['alphabet'], cap['size']): cap
                 for cap in linear_caps(1, [0]*q, [0]*q)}
    controls = []
    for label, weight in entry['integer_constraint_weights']:
        assert type(weight) is int and weight > 0
        row, rhs = independent_constraint(q, label)
        if label[0] != 'prefix':
            cap = reference[tuple(label)]
            assert rhs == cap['rhs']
            assert row == [a+b for a, b in zip(cap['first'], cap['second'])]
        coefficients = [a+weight*b for a, b in zip(coefficients, row)]
        total += weight*rhs
        controls.append({'label': label, 'integer_weight': weight, 'rhs': str(rhs)})
    probabilities = [value/(2*total) for value in coefficients]
    probabilities[0] = 1-sum(probabilities[1:])
    assert min(probabilities) > 0 and sum(probabilities) == 1
    assert 2*rho.numerator**2 > rho.denominator**2
    arithmetic, entropy = Dyadic(100), (0, 0)
    for probability, multiplicity in Counter(probabilities).items():
        log = log_at_least_one(arithmetic, point(arithmetic, 1/probability))
        entropy = arithmetic.add(entropy, arithmetic.multiply(
            point(arithmetic, multiplicity*probability), log))
    required = arithmetic.multiply(point(arithmetic, 2*rho),
                                  log_at_least_one(arithmetic, point(arithmetic, q)))
    margin = arithmetic.subtract(entropy, required)
    assert margin[0] > 0
    return {'alphabet': q, 'rho_floor': str(rho), 'weighted_rhs_before_normalization': str(total),
            'normalized_weighted_rhs': '1/2', 'constraints': controls,
            'probability_runs': [[str(p), m] for p, m in Counter(probabilities).items()],
            'complete_probability': list(map(str, probabilities)), 'precision': 100,
            'entropy_interval': list(map(str, entropy)),
            'entropy_margin_interval': list(map(str, margin)),
            'entropy_margin_decimal': [float(F(v, arithmetic.scale)) for v in margin],
            'strict_half_rate_integer_comparison': [2*rho.numerator**2, rho.denominator**2]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    records = json.loads((root/'evidence/subalphabet-bare-entropy-witnesses.json').read_text())
    assert [r['alphabet'] for r in records] == [7, 8, 9]
    result = [replay(root, r) for r in records]
    names = ['tools/check_subalphabet_bare_entropy.py', 'tools/subalphabet_adjugate_caps.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py',
             'evidence/subalphabet-bare-entropy-witnesses.json']
    report = {'status': 'pass', 'all_scale_entropy_certificates': result,
              'source_sha256': {n: sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope': 'All fields, all degrees and all scales for bare scalar local Smith '
              'pressure. The ordinary auxiliary projection transfer also applies. '
              'Independent physical gauges and improved factoring algorithms are outside scope.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': 'pass', 'floors': [[r['alphabet'], r['rho_floor']]
          for r in result], 'positive_entropy_margins': [r['entropy_margin_decimal'] for r in result]}))


if __name__ == '__main__':
    main()
