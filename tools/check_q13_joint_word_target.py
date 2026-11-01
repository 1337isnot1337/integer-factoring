#!/usr/bin/env python3
"""Exactly replay every declared compatibility family for a new q13 profile."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from box_adjugate_word_caps import inequalities as box_adjugate_word_caps
from check_bare_seven_to_twelve_letters import adjugate_violations, constraints, ln_bounds
from check_cell_minor_caps import all_caps
from hereditary_adjugate_caps import inequalities as square_caps
from higher_jet_word_caps import inequalities as word_caps
from inherited_second_jet_caps import inequalities as second_jet_caps
from multichannel_cofactor_caps import linear_caps as multichannel_caps
from rectangular_adjugate_jet_caps import linear_caps as rectangular_caps

ROOT = Path(__file__).resolve().parents[1]


def run():
    path = ROOT/'evidence/q13-joint-word-target.json'
    data = json.loads(path.read_text())
    q, degree, profile = data['alphabet'], data['M'], data['mu']
    assert q == 13 and degree == 2000 and data['identical_profiles'] is True
    assert len(profile) == q*q and all(type(x) is int for x in profile)
    assert profile == sorted(profile) and profile[0] == 0
    assert not adjugate_violations(q, profile, profile, degree)
    checks = {}
    caps = all_caps(degree, [0]*q, [0]*q)
    slacks = [bound-sum(profile[:k]) for k, (bound, _) in caps.items()]
    assert min(slacks) >= 0
    checks['cell_shared'] = dict(checks=len(caps), minimum_slack=str(min(slacks)))
    slacks = [size*degree-sum((a+b)*v for a, b, v in zip(first, second, profile))
              for _, first, second, size in constraints(q)]
    assert len(slacks) == 182 and min(slacks) >= 0
    checks['box'] = dict(oriented_checks=len(slacks), minimum_slack=min(slacks))
    for label, builder in (('square', square_caps), ('second_jet', second_jet_caps), ('word_jet', word_caps)):
        caps = builder(degree, profile, profile, [0]*q, [0]*q)
        slacks = [cap['slack'] for cap in caps]
        assert len(caps) == 1816 and min(slacks) >= 0
        checks[label] = dict(oriented_checks=len(caps), minimum_slack=str(min(slacks)))
    for label, builder in (('rectangular_jet', rectangular_caps), ('multichannel', multichannel_caps)):
        caps = builder(degree, [0]*q, [0]*q)
        slacks = [cap['rhs']-sum((a+b)*v for a, b, v in zip(cap['first'], cap['second'], profile))
                  for cap in caps]
        assert len(caps) == 8918 and min(slacks) >= 0
        checks[label] = dict(compiled_checks=len(caps), minimum_slack=str(min(slacks)), identical_orientations=True)
    # The original seven families remain valid checks of an incomplete
    # relaxation. A later composition applies box budgets to the adjugate too.
    later = box_adjugate_word_caps(degree, profile, profile)
    bad = [cap for cap in later if cap['slack'] < 0]
    assert len(later) == 1816 and len(bad) == 52 and min(cap['slack'] for cap in later) == -12956
    checks['later_box_adjugate_word'] = dict(status='excluded',oriented_checks=len(later),
                                            violated_checks=len(bad),minimum_slack=-12956)
    # Compute B exactly; round it upward before taking small rational log bounds.
    s = F(data['s'])
    assert s == F(414548, 413547)
    Z = sum(s**(-order) for order in profile)
    B = s**degree*Z*Z
    precision = 2**96
    upper = F(-((-B.numerator*precision)//B.denominator), precision)
    assert upper >= B
    _, log_upper = ln_bounds(upper)
    log_q_lower, _ = ln_bounds(F(q))
    numerator, denominator = data['rho_upper_bound']
    assert (numerator, denominator) == (7061, 10000)
    assert log_upper*denominator < 4*numerator*log_q_lower
    assert 2*numerator*numerator < denominator*denominator
    exact_rate = dict(thermal_scale=str(s), upward_dyadic_B=str(upper), dyadic_bits=96,
                      log_B_upper=str(log_upper), log_q_lower=str(log_q_lower),
                      formal_rho_upper_bound=str(F(numerator, denominator)),
                      threshold_integer_comparison=[2*numerator*numerator, denominator*denominator])
    names = ['tools/check_q13_joint_word_target.py', 'tools/check_bare_seven_to_twelve_letters.py',
             'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py',
             'tools/hereditary_adjugate_caps.py', 'tools/higher_jet_word_caps.py',
             'tools/inherited_second_jet_caps.py', 'tools/multichannel_cofactor_caps.py',
             'tools/rectangular_adjugate_jet_caps.py', 'tools/box_adjugate_word_caps.py']
    return dict(scope='Historical q13 point passes its original seven families but is excluded by later box-adjugate word constraints; formal rate only.',
                alphabet=q, degree=degree, profile_length=len(profile), checks=checks,
                exact_rate_certificate=exact_rate, original_square_box_adjugates='pass',
                input_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/q13-joint-word-target-controls.json').read_text())
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','exact_rate_certificate')}))


if __name__ == '__main__':
    main()
