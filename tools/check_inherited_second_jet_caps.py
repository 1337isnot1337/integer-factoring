#!/usr/bin/env python3
"""Replay graded nonlinear caps and reject a linear q12 relaxation point."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from check_cell_minor_caps import Dyadic, all_caps, literal_case, pressure, random_crossing
from hereditary_adjugate_caps import inequalities as square_caps
from inherited_second_jet_caps import inequalities, truncated_prefix
from multichannel_cofactor_caps import linear_caps as multichannel_caps
from rectangular_adjugate_jet_caps import linear_caps as jet_caps

ROOT = Path(__file__).resolve().parents[1]


def literal_controls():
    rng, count, records = Random(938625), 0, []
    for q in (2, 3, 4):
        for _ in range(3):
            tensor = random_crossing(rng, q, 2, .65)
            heights = [[[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2)]
            data = literal_case(rng, q, tensor, heights, cell_samples=0)
            for pair in data['profiles'].values():
                caps = inequalities(data['deficit'], pair['r'], pair['s'], data['paired_energies_X'], data['paired_energies_Y'])
                assert min(cap['slack'] for cap in caps) >= 0
                count += len(caps)
            records.append(data)
    return dict(tensors=9, unique_complete_profiles=36, inequalities=count,
                profiles_sha256=sha256(json.dumps(records, sort_keys=True, default=str).encode()).hexdigest())


def normalization_controls():
    rng, count = Random(610212), 0
    for n in (4, 9, 16):
        for _ in range(60):
            b = F(rng.randrange(10), 3)
            raw = [b]+sorted(b+F(rng.randrange(30), 4) for _ in range(n-1))
            local = [v-b for v in raw]
            for k in range(1, n+1):
                delta = truncated_prefix(raw, k)-truncated_prefix(local, k)
                assert k*b <= delta <= (2*k-1)*b
                assert k*(n-1)*b-delta >= (k*(n-3)+1)*b >= 0
                count += 1
    return count


def q12_control():
    path = ROOT/'evidence/inherited-second-jet-q12-point.json'
    point = json.loads(path.read_text())
    q, denominator = point['alphabet'], point['denominator']
    d, r = F(point['deficit'], denominator), [F(v, denominator) for v in point['orders']]
    e = point['paired_energies']
    assert q == 12 and e == [0]*q and len(r) == q*q and r == sorted(r) and not r[0]
    assert all(sum(r[:k]) <= cap for k, (cap, _) in all_caps(d, e, e).items())
    square = square_caps(d, r, r, e, e)
    assert min(cap['slack'] for cap in square) >= 0
    counts = []
    for builder in (jet_caps, multichannel_caps):
        caps = builder(d, e, e)
        assert all(sum(c*v for c, v in zip(cap['first'], r))+
                   sum(c*v for c, v in zip(cap['second'], r)) <= cap['rhs'] for cap in caps)
        counts.append(len(caps))
    arithmetic = Dyadic(100)
    interval = pressure(arithmetic, d, r, e, e)
    assert interval[1]*100 < -arithmetic.scale
    rejected = [cap for cap in inequalities(d, r, r, e, e) if cap['slack'] < 0]
    full = [cap for cap in rejected if cap['physical_size'] == q and cap['size'] in (q, q*q)]
    assert any(cap['size'] == q*q for cap in full)
    return dict(scope='Negative point of preceding linear constraints, excluded by the existing second-jet theorem; no tensor.',
                oriented_square_caps=len(square), rectangular_jet_caps=counts[0], multichannel_caps=counts[1],
                pressure_interval=list(interval), precision=100, pressure_upper_bound='-1/100',
                nonlinear_violations=len(rejected), full_alphabet_rejections=json.loads(json.dumps(full, default=str)),
                input_sha256=sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    names = ['tools/check_inherited_second_jet_caps.py', 'tools/inherited_second_jet_caps.py',
             'tools/hereditary_adjugate_caps.py', 'tools/rectangular_adjugate_jet_caps.py',
             'tools/multichannel_cofactor_caps.py', 'tools/check_cell_minor_caps.py', 'tools/charged_gaussian_floor.py']
    result = dict(scope='Existing universal second-jet theorem lifted to independent gradings and square restrictions; nonlinear conditions on the two actual profiles.',
                  literal_controls=literal_controls(), normalization_identities=normalization_controls(), q12=q12_control(),
                  source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/inherited-second-jet-controls.json').read_text())
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
