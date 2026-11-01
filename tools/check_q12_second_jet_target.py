#!/usr/bin/env python3
"""Exact negative q12 point of the declared linear and nonlinear relaxation."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from check_cell_minor_caps import Dyadic, all_caps, pressure
from hereditary_adjugate_caps import inequalities as square_caps
from inherited_second_jet_caps import inequalities as second_jet_caps
from multichannel_cofactor_caps import linear_caps as multichannel_caps
from rectangular_adjugate_jet_caps import linear_caps as jet_caps

ROOT = Path(__file__).resolve().parents[1]


def certificate():
    path = ROOT/'evidence/q12-second-jet-target.json'
    point = json.loads(path.read_text())
    q, denominator = point['alphabet'], point['denominator']
    d, r = F(point['deficit'], denominator), [F(v, denominator) for v in point['orders']]
    e = point['paired_energies']
    assert q == 12 and e == [0]*q and len(r) == q*q and r == sorted(r) and not r[0]
    assert all(sum(r[:k]) <= bound for k, (bound, _) in all_caps(d, e, e).items())
    checks = {}
    for label, builder in (('square', square_caps), ('second_jet', second_jet_caps)):
        caps = builder(d, r, r, e, e)
        assert min(cap['slack'] for cap in caps) >= 0
        checks[label] = dict(oriented_checks=len(caps), minimum_slack=str(min(c['slack'] for c in caps)))
    for label, builder in (('rectangular_jet', jet_caps), ('multichannel', multichannel_caps)):
        caps = builder(d, e, e)
        slacks = [cap['rhs']-sum(c*v for c, v in zip(cap['first'], r))-
                  sum(c*v for c, v in zip(cap['second'], r)) for cap in caps]
        assert min(slacks) >= 0
        # The equal pair satisfies the identical reverse-cut constraints.
        checks[label] = dict(compiled_checks=len(caps), minimum_slack=str(min(slacks)))
    arithmetic = Dyadic(100)
    interval = pressure(arithmetic, d, r, e, e)
    assert list(interval) == point['pressure_interval']
    assert interval[1]*100 < -arithmetic.scale
    names = ['tools/check_q12_second_jet_target.py', 'tools/hereditary_adjugate_caps.py',
             'tools/inherited_second_jet_caps.py', 'tools/multichannel_cofactor_caps.py',
             'tools/rectangular_adjugate_jet_caps.py', 'tools/check_cell_minor_caps.py',
             'tools/charged_gaussian_floor.py']
    return dict(scope='A negative point of the declared necessary linear and second-jet inequalities; no tensor or stronger factoring algorithm.',
                alphabet=q, integer_deficit=point['deficit'], scale='1/'+str(denominator),
                checks=checks, pressure_interval=list(interval), precision=100, pressure_upper_bound='-1/100',
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
        assert result == json.loads((ROOT/'evidence/q12-second-jet-target-controls.json').read_text())
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
