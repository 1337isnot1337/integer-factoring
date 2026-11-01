#!/usr/bin/env python3
"""Exact growing-field descent, charge-area identities and finite span controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from check_controlled_permutation_smith import binary_rank
from check_fixed_field_descent import factor_controls
from check_subgroup_return_crossings import rank_pairs, upper_path
from ordered_charge_area import area_zero_cuts, clock, local_rows, statistics

ROOT = Path(__file__).resolve().parents[1]


def word_controls():
    rng = Random(312098)
    local = []
    for row, entries in enumerate(local_rows()):
        x, y = divmod(row, 3)
        for column, _, _, c, _ in entries:
            u, v = divmod(column, 3)
            assert abs(c) <= 2
            local.append((x, y, u, v))
    cases = 0
    for length in range(1, 9):
        field, powers = clock(length)
        order, half = len(powers), (len(powers)+1)//2
        psi = [powers[-half*j*j % order] for j in range(order)]
        phi = [powers[half*j*j % order] for j in range(order)]
        for trial in range(64):
            physical = [rng.choice(local) for _ in range(length)]
            if trial < 16:
                physical = [(x, y, x, y) for x, y, _, _ in physical]
            elif trial < 32:
                physical = [(x, y, y, x) for x, y, _, _ in physical]
            vector = psi[:]
            for x, y, u, v in reversed(physical):
                vector = [vector[(j+v) % order] for j in range(order)]
                vector = [field.mul(powers[(y-u)*j % order], z) for j, z in enumerate(vector)]
                vector = [vector[(j-x) % order] for j in range(order)]
                scalar = powers[-half*(x*x-v*v) % order]
                vector = [field.mul(scalar, z) for z in vector]
            value = 0
            for a, b in zip(phi, vector):
                value ^= field.mul(a, b)
            charge, area = statistics(physical)
            assert abs(charge) <= 2*length < order
            assert abs(area) <= 2*length*(length+1)
            assert value == (powers[area % order] if charge == 0 else 0)
            if trial < 32:
                assert (charge, area, value) == (0, 0, 1)
            # Independent double sum for the area (no prefix-state update).
            c = [x+y-u-v for x, y, u, v in physical]
            a = [x-v for x, _, _, v in physical]
            literal = -sum(v*dc for (_, _, _, v), dc in zip(physical, c))
            literal += sum(c[i]*a[j] for i in range(length) for j in range(i+1, length))
            assert literal == area
            cases += 1
    # Area-zero selection is itself order-sensitive.
    first, second = (0, 1, 0, 2), (0, 2, 1, 0)
    assert statistics((first, second)) == (0, 2)
    assert statistics((second, first)) == (0, 0)
    return cases


def generate():
    records = []
    for length, threshold in ((2, 576), (3, 92338)):
        cuts, support = area_zero_cuts(length)
        ranks = rank_pairs(cuts)
        records.append(dict(length=length, support=support, ranks=ranks, threshold=threshold,
                            detector=1 << (3 if length == 2 else 5),
                            upper_path=upper_path(ranks, threshold)))
    names = ('ordered_charge_area.py', 'check_ordered_charge_area.py',
             'check_subgroup_return_crossings.py', 'check_fixed_field_descent.py',
             'check_steiner_cubic_crossings.py')
    return dict(records=records,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in names},
                scope='Uniform identities and cost; complete F2 area-zero spans only at lengths two and three.',
                factoring_exponent='1/2; unchanged')


def verify(data):
    assert len(data['records']) == 2
    rank_count = 0
    for record, (length, threshold) in zip(data['records'], ((2, 576), (3, 92338))):
        assert record['length'] == length and record['threshold'] == threshold
        cuts, support = area_zero_cuts(length)
        assert support == record['support']
        assert binary_rank(cuts[1][0]) == 1
        ranks = rank_pairs(cuts, low=True)
        assert ranks == list(map(tuple, record['ranks']))
        detector = record['detector']
        assert detector.bit_count() == 1
        assert min(r*c for mask, (r, c) in enumerate(ranks) if mask & detector) == threshold
        endpoint = 0
        for mask in record['upper_path']:
            endpoint ^= mask
            assert ranks[mask][0]*ranks[mask][1] <= threshold
        assert endpoint == len(ranks)-1
        assert threshold >= 24**length
        rank_count += 2*len(ranks)
    assert 24**7 > 3**20 and 2*5**2 > 7**2
    fields = [factor_controls(clock(1 << k)[0]) for k in range(7)]
    assert [item['sharp_rank_increase'][1] for item in fields] == list(range(2, 9))
    for name, digest in data['source_sha256'].items():
        assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest
    return dict(independent_charge_area_words=word_controls(),
                extension_descent_controls=fields, independently_replayed_cut_ranks=rank_count,
                complete_finite_area_zero_cases=2, strict_finite_rho_floor='5/7 > 1/sqrt(2)',
                scope=data['scope'], factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/ordered-charge-area.json')
    args = parser.parse_args()
    if args.generate:
        data = generate()
        args.evidence.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    else:
        data = json.loads(args.evidence.read_text())
    print(json.dumps(verify(data)))


if __name__ == '__main__':
    main()
