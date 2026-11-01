#!/usr/bin/env python3
"""Verify an integer q12 profile point and a rational-scale rate certificate."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_cell_minor_caps import all_caps
from hereditary_adjugate_caps import inequalities as square_caps
from higher_jet_word_caps import inequalities as word_caps
from inherited_second_jet_caps import inequalities as second_jet_caps
from multichannel_cofactor_caps import linear_caps as multichannel_caps
from rectangular_adjugate_jet_caps import linear_caps as jet_caps

ROOT = Path(__file__).resolve().parents[1]


def integer_digest(value):
    return sha256(value.to_bytes((value.bit_length()+7)//8, 'big')).hexdigest()


def certificate():
    path = ROOT/'evidence/q12-integer-word-target.json'
    data = json.loads(path.read_text())
    q, degree, r, e = data['alphabet'], data['deficit'], data['orders'], data['paired_energies']
    assert q == 12 and degree > 0 and type(degree) is int and e == [0]*q
    assert len(r) == q*q and r == sorted(r) and not r[0] and all(type(v) is int for v in r)
    assert all(sum(r[:k]) <= cap for k, (cap, _) in all_caps(degree, e, e).items())
    checks = {}
    for label, builder in (('square', square_caps), ('second_jet', second_jet_caps), ('word_jet', word_caps)):
        caps = builder(degree, r, r, e, e)
        assert min(cap['slack'] for cap in caps) >= 0
        checks[label] = dict(oriented_checks=len(caps), minimum_slack=str(min(c['slack'] for c in caps)))
    for label, builder in (('rectangular_jet', jet_caps), ('multichannel', multichannel_caps)):
        caps = builder(degree, e, e)
        slacks = [cap['rhs']-sum(c*v for c, v in zip(cap['first'], r))-
                  sum(c*v for c, v in zip(cap['second'], r)) for cap in caps]
        assert min(slacks) >= 0
        checks[label] = dict(compiled_checks=len(caps), minimum_slack=str(min(slacks)))
    numerator, denominator = data['rational_thermal_scale']
    bound, divisor = data['rho_upper_bound']
    assert (numerator, denominator, bound, divisor) == (17, 14, 707, 1000)
    highest = max(r)
    assert 2*highest >= degree
    polynomial = sum(numerator**(highest-a)*denominator**a for a in r)
    # B(s)=s^degree (sum s^(-r_i))^2. Clearing B(17/14)^250<12^707:
    left = polynomial**500
    right = q**707*numerator**(250*(2*highest-degree))*denominator**(250*degree)
    assert left < right and 2*bound**2 < divisor**2
    exact = dict(max_order=highest, left_bits=left.bit_length(), right_bits=right.bit_length(),
                 left_sha256=integer_digest(left), right_sha256=integer_digest(right))
    assert exact == data['integer_certificate']
    names = ['tools/check_q12_integer_word_target.py', 'tools/hereditary_adjugate_caps.py',
             'tools/higher_jet_word_caps.py', 'tools/inherited_second_jet_caps.py',
             'tools/multichannel_cofactor_caps.py', 'tools/rectangular_adjugate_jet_caps.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py']
    return dict(scope='Exact integer point of declared necessary constraints, with formal rate below707/1000; no realizing tensor or factoring improvement.',
                alphabet=q, integer_degree=degree, identical_profiles=True, checks=checks,
                rational_thermal_scale=[numerator, denominator], formal_rho_upper_bound='707/1000',
                exact_integer_comparison=exact, strict_half_comparison=[2*bound**2, divisor**2],
                input_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/q12-integer-word-target-controls.json').read_text())
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
