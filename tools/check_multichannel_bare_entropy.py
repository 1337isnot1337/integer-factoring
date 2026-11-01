#!/usr/bin/env python3
"""Independent exact all-scale bare-pressure witnesses on ten/eleven letters."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from check_cell_minor_caps import Dyadic, all_caps, log_at_least_one, point
from hereditary_adjugate_caps import constraints
from multichannel_cofactor_caps import linear_caps as multichannel_caps
from rectangular_adjugate_jet_caps import linear_caps as jet_caps

ROOT = Path(__file__).resolve().parents[1]


def independent_row(q, label):
    values = [F(0)]*(q*q)
    if label[0] == 'prefix':
        assert label == ['prefix', q]
        values[:q] = [F(1)]*q
        return values, F(q, 2)
    family, p, b, k = label
    n = p*b
    if family == 'rectangular jet cell':
        assert p < b == q and k == n
        values[:n-1] = [F(n)]*(n-1)
        values[1] += (n-p*p)*(n-1)
    else:
        assert k == b == q
        values[:n] = [F(k-(i >= n-k)) for i in range(n)]
        if family == 'square cofactor shared line':
            assert p == b
            values[1] += k-1
        else:
            assert family == 'multichannel cofactor shared line' and p < b
            h = p*p
            for i in range(n):
                values[i] += k*min(h, max(0, n-1-h*i))
    return values, F(k*(n-1))


def replay(entry):
    q, rho = entry['alphabet'], F(entry['rho_floor'])
    caps = jet_caps(1, [0]*q, [0]*q)+multichannel_caps(1, [0]*q, [0]*q)
    reference = {(cap['family'], *cap['alphabet'], cap['size']):
                 ([a+b for a, b in zip(cap['first'], cap['second'])], cap['rhs']) for cap in caps}
    for cap in constraints(1, [0]*q, [0]*q):
        if cap['family'].endswith('shared'):
            row = [sum(c*int(i < k) for k, c in cap['first_prefixes']+cap['second_prefixes']) for i in range(q*q)]
            m = cap['physical_size']
            reference['square cofactor shared line', m, m, cap['size']] = (row, cap['rhs'])
    total, coefficients = F(0), [F(0)]*(q*q)
    for label, weight in entry['integer_constraint_weights']:
        assert type(weight) is int and weight > 0
        row, rhs = independent_row(q, label)
        if label[0] == 'prefix':
            assert rhs == all_caps(1, [0]*q, [0]*q)[q][0]
        else:
            expected, bound = reference[tuple(label)]
            # Square prefix restoration may add a multiple of z1=0.
            assert rhs == bound and row[1:] == expected[1:]
        total += weight*rhs
        coefficients = [a+weight*b for a, b in zip(coefficients, row)]
    probabilities = [v/(2*total) for v in coefficients]
    probabilities[0] = 1-sum(probabilities[1:])
    assert min(probabilities) > 0 and sum(probabilities) == 1
    arithmetic, entropy = Dyadic(100), (0, 0)
    for p, multiplicity in Counter(probabilities).items():
        entropy = arithmetic.add(entropy, arithmetic.multiply(point(arithmetic, multiplicity*p),
            log_at_least_one(arithmetic, point(arithmetic, 1/p))))
    required = arithmetic.multiply(point(arithmetic, 2*rho), log_at_least_one(arithmetic, point(arithmetic, q)))
    margin = arithmetic.subtract(entropy, required)
    assert margin[0] > 0 and 2*rho.numerator**2 > rho.denominator**2
    return dict(alphabet=q, rho_floor=str(rho), constraint_count=len(entry['integer_constraint_weights']),
                combined_rhs=str(total), complete_probability=list(map(str, probabilities)),
                probability_runs=[[str(p), m] for p, m in Counter(probabilities).items()],
                normalized_weighted_rhs='1/2', precision=100, entropy_interval=list(entropy),
                entropy_margin_interval=list(margin), entropy_margin_decimal=[float(F(v, arithmetic.scale)) for v in margin],
                strict_half_comparison=[2*rho.numerator**2, rho.denominator**2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    entries = json.loads((ROOT/'evidence/multichannel-bare-entropy-witnesses.json').read_text())
    assert [entry['alphabet'] for entry in entries] == [10, 11]
    names = ['tools/check_multichannel_bare_entropy.py', 'tools/multichannel_cofactor_caps.py',
             'tools/rectangular_adjugate_jet_caps.py', 'tools/hereditary_adjugate_caps.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py',
             'evidence/multichannel-bare-entropy-witnesses.json']
    result = dict(scope='All fields, degrees and scales for bare scalar Smith pressure on ten/eleven letters; no independently gauged exclusion or stronger factoring algorithm.',
                  certificates=[replay(entry) for entry in entries],
                  source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/multichannel-bare-entropy-controls.json').read_text())
    print(json.dumps(dict(floors=[[r['alphabet'], r['rho_floor']] for r in result['certificates']],
                         entropy_margins=[r['entropy_margin_decimal'] for r in result['certificates']])))


if __name__ == '__main__':
    main()
