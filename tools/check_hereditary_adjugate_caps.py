#!/usr/bin/env python3
"""Exact hereditary-cap controls and a surviving integer six-letter target."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
from random import Random

from check_cell_minor_caps import Dyadic, all_caps, literal_case, pressure, random_crossing
from hereditary_adjugate_caps import constraints, inequalities

ROOT = Path(__file__).resolve().parents[1]


def raw_minima(tensor, heights):
    return [min((e if side == 0 else -e)+sum(heights[side][l][a[l]] for l in range(4))
                for a, p in tensor.items() for e in p) for side in (0, 1)]


def literal_controls():
    rng = Random(0xA6D0C0)
    tensors, profiles, cap_count, restricted, normalization_checks = 0, 0, 0, 0, 0
    records = []
    for q in (2, 3, 4):
        for case in range(4):
            heights = [[[rng.randint(-2, 2) for _ in range(q)] for _ in range(4)] for _ in range(2)]
            tensor = random_crossing(rng, q, 2, .65)
            result = literal_case(rng, q, tensor, heights, cell_samples=0)
            ex, ey, d = result['paired_energies_X'], result['paired_energies_Y'], result['deficit']
            for pair in result['profiles'].values():
                caps = inequalities(d, pair['r'], pair['s'], ex, ey)
                assert all(c['slack'] >= 0 for c in caps)
                cap_count += len(caps)
                profiles += 2
            tensors += 1
            # Independent normalization/restriction checks on proper cores.
            if case == 0:
                full_min = raw_minima(tensor, heights)
                for m in range(2, q):
                    ax = sorted(range(q), key=lambda i: (ex[i], i))[:m]
                    by = sorted(range(q), key=lambda i: (ey[i], i))[:m]
                    local = {}
                    for x, y, u, v in product(range(m), repeat=4):
                        key = (ax[x], by[y], by[u], ax[v])
                        if key in tensor:
                            local[x, y, u, v] = tensor[key]
                    local_heights = [[[heights[s][l][i] for i in (ax if l in (0, 3) else by)]
                                      for l in range(4)] for s in (0, 1)]
                    local_result = literal_case(rng, m, local, local_heights, cell_samples=0)
                    local_min = raw_minima(local, local_heights)
                    a, b = [u-v for u, v in zip(local_min, full_min)]
                    assert min(a, b) >= 0
                    assert local_result['deficit'] == d+min(ex[i] for i in ax)+min(ey[i] for i in by)-a-b
                    full_caps = [c for c in constraints(d, ex, ey) if c['physical_size'] == m]
                    local_caps = [c for c in constraints(local_result['deficit'],
                                                       local_result['paired_energies_X'],
                                                       local_result['paired_energies_Y']) if c['physical_size'] == m]
                    for name in result['profiles']:
                        g = result['profiles'][name]
                        core = local_result['profiles'][name]
                        u, v = [x+a for x in core['r']], [x+b for x in core['s']]
                        assert all(sum(u[:k]) >= sum(g['r'][:k]) and sum(v[:k]) >= sum(g['s'][:k]) for k in range(1, m*m+1))
                        for cap, lc in zip(full_caps, local_caps):
                            k = cap['size']
                            for first, second, low, lr, ls in ((u, v, b, core['r'], core['s']),
                                                              (v, u, a, core['s'], core['r'])):
                                def value(bound, r, s):
                                    return (sum(c*sum(r[:i]) for i, c in bound['first_prefixes'])+
                                            sum(c*sum(s[:i]) for i, c in bound['second_prefixes']))
                                extra = k*(m*m-3)+2
                                raw_slack = cap['rhs']-value(cap, first, second)-extra*low
                                local_slack = lc['rhs']-value(lc, lr, ls)
                                assert raw_slack == local_slack >= 0
                                normalization_checks += 1
                    restricted += 1
            records.append(result)
    return dict(tensors=tensors, complete_profiles=profiles, cap_checks=cap_count,
                proper_core_restrictions=restricted, normalization_identities=normalization_checks,
                record_sha256=sha256(json.dumps(records, sort_keys=True, default=str).encode()).hexdigest())


def subset_controls():
    rng = Random(15091)
    comparisons = 0
    for _ in range(4):
        x, y = sorted([0]+[F(rng.randrange(15), 3) for _ in range(4)]), sorted([0]+[F(rng.randrange(15), 4) for _ in range(4)])
        for m in range(2, 6):
            for cap in [c for c in constraints(1, x, y) if c['physical_size'] == m]:
                k = cap['size']
                for ax, by in product(combinations(x, m), combinations(y, m)):
                    constant = m*m-1+m*(sum(ax)+sum(by))
                    if cap['family'].endswith('cell'):
                        rhs = k*constant-sum(sorted(a+b for a in ax for b in by)[-k:])
                    else:
                        rhs = k*constant-max(sum(sorted(ax)[-k:])+k*max(by), sum(sorted(by)[-k:])+k*max(ax))
                    assert cap['rhs'] <= rhs
                    comparisons += 1
    return comparisons


def certificate():
    path = ROOT/'evidence/hereditary-q6-integer-target.json'
    point = json.loads(path.read_text())
    d, target, energies = point['deficit'], point['target'], point['energies']
    assert d == 6 and energies == [0, 0, 12, 16, 18, 24]
    assert len(target) == 36 and target == sorted(target) and target[0] == 0
    for k, (cap, _) in all_caps(F(d), energies, energies).items():
        assert sum(target[:k]) <= cap
    caps = inequalities(d, target, target, energies, energies)
    assert all(c['slack'] >= 0 for c in caps)
    scale = F(*point['scale'])
    arithmetic = Dyadic(point['pressure_precision'])
    interval = pressure(arithmetic, d*scale, [v*scale for v in target],
                        [v*scale for v in energies], [v*scale for v in energies])
    assert list(interval) == point['pressure_interval']
    assert interval[1]*250 < -arithmetic.scale
    for name, digest in point['source_sha256'].items():
        assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest, name
    old = json.loads((ROOT/'evidence/zero-energy-hereditary-adjugate.json').read_text())
    rejected = inequalities(old['deficit'], old['target'], old['target'], old['paired_energies'], old['paired_energies'])
    assert any(c['physical_size'] == 2 and c['size'] == 4 and c['slack'] == -6 for c in rejected)
    literal = literal_controls()
    subset = subset_controls()
    names = (Path(__file__).name, 'hereditary_adjugate_caps.py', 'check_cell_minor_caps.py', 'charged_gaussian_floor.py')
    return dict(scope='All-field hereditary necessary constraints and an exact negative point of their relaxation; no realizing tensor or stronger factoring algorithm.',
                literal_controls=literal, subset_minimum_comparisons=subset,
                q6_deficit=d, q6_paired_energies=energies, q6_target=target,
                q6_cap_checks=len(caps), q6_minimum_cap_slack=str(min(c['slack'] for c in caps)),
                pressure_upper_bound='-1/250', pressure_interval=list(interval), precision=arithmetic.precision,
                old_q5_target_rejected=True,
                input_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/hereditary-adjugate-controls.json').read_text())
    print(json.dumps({k: v for k, v in result.items() if k not in ('q6_target', 'source_sha256')}))


if __name__ == '__main__':
    main()
