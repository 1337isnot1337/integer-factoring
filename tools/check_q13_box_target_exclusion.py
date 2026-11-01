#!/usr/bin/env python3
"""Check newer cofactor exclusions of the historical thirteen-letter box target."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_bare_seven_to_twelve_letters import adjugate_violations, satisfies
from hereditary_adjugate_caps import inequalities as square_caps
from higher_jet_word_caps import inequalities as word_caps
from inherited_second_jet_caps import inequalities as second_jet_caps
from multichannel_cofactor_caps import linear_caps as multichannel_caps
from rectangular_adjugate_jet_caps import linear_caps as rectangular_caps

ROOT = Path(__file__).resolve().parents[1]


def run():
    path = ROOT/'evidence/bare-seven-to-twelve-letters.json'
    target = json.loads(path.read_text())['thirteen_target']
    q, degree, profile = 13, target['M'], target['mu']
    assert degree == 2000 and len(profile) == q*q and profile == sorted(profile) and profile[0] == 0
    assert not satisfies(q, profile, profile, degree)
    assert not adjugate_violations(q, profile, profile, degree)
    families = {}
    for label, builder, expected in (
            ('square', square_caps, (1816, 0, 523)),
            ('second_jet', second_jet_caps, (1816, 348, -76060)),
            ('word_jet', word_caps, (1816, 348, -213861))):
        caps = builder(degree, profile, profile, [0]*q, [0]*q)
        bad = [c for c in caps if c['slack'] < 0]
        minimum = min(c['slack'] for c in caps)
        assert (len(caps), len(bad), minimum) == expected
        if bad:
            witness = bad[0]
            assert (witness['physical_size'], witness['size'], witness['lhs'], witness['rhs']) == (13, 5, 1680456, 1680000)
        families[label] = dict(checks=len(caps), violated_checks=len(bad),
                               minimum_slack=int(minimum), violations=bad)
    for label, builder, expected in (
            ('rectangular_jet', rectangular_caps, 655),
            ('multichannel', multichannel_caps, 588)):
        caps = builder(degree, [0]*q, [0]*q)
        slacks = [cap['rhs']-sum(c*v for c, v in zip(cap['first'], profile))-
                  sum(c*v for c, v in zip(cap['second'], profile)) for cap in caps]
        assert len(caps) == 8918 and min(slacks) == expected
        families[label] = dict(checks=len(caps), violated_checks=0, minimum_slack=expected)
    names = ['tools/check_q13_box_target_exclusion.py',
             'tools/check_bare_seven_to_twelve_letters.py', 'tools/hereditary_adjugate_caps.py',
             'tools/higher_jet_word_caps.py', 'tools/inherited_second_jet_caps.py',
             'tools/multichannel_cofactor_caps.py', 'tools/rectangular_adjugate_jet_caps.py']
    return dict(scope='The recorded q13 box-relaxation profile is excluded by newer jet caps; the q13 criterion itself remains open.',
                alphabet=q, degree=degree, original_box_and_adjugate_checks='pass',
                input_sha256=sha256(path.read_bytes()).hexdigest(), families=families,
                source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    # Convert exact integral Fraction fields into the receipt's stable strings.
    result = json.loads(json.dumps(result, default=str))
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/q13-box-target-exclusion.json').read_text())
    print(json.dumps({'scope':result['scope'], 'families':{
        name:{k:v for k,v in data.items() if k != 'violations'}
        for name,data in result['families'].items()}}))


if __name__ == '__main__':
    main()
