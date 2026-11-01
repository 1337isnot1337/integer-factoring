#!/usr/bin/env python3
"""Finite inherited-height search in the strict eleven-parameter two-place pencil."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_joint_energy_conjugations import cut_matrix
from check_joint_endpoint_loops import operator_coefficients
from check_joint_zero_prefix import two_place_pencil
from joint_energy_conjugations import model
from native_weighted_smith_f4 import Oracle
from search_f4_recentered_crossings import five_letter_source, root_for_profiles
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'evidence/joint-two-place-pencil-search.json')
    parser.add_argument('--random-trials', type=int, default=0)
    parser.add_argument('--seed', type=int, default=70407)
    args = parser.parse_args()
    data = model()
    pencil = two_place_pencil(data, operator_coefficients(data))
    first, second = pencil['first_matrix'], pencil['second_matrix']
    source, heights, ex, ey = five_letter_source()
    parameters = [v for group in heights for v in group[1:]]+[1.0]
    half = [v//2 for v in ex]
    weights = [half[i//5] for i in range(25)]

    def profile(values):
        try:
            _, _, raw = Oracle(values, 5).raw_evaluate(parameters, want_spectra=True)
            orders = tuple(map(round, raw))
            assert all(abs(a-b) < 1e-8 for a, b in zip(raw, orders))
            return orders
        except ArithmeticError:
            infinity, _ = WeightedSmithOracle(cut_matrix(values, True), infinity=True).profile(weights, weights)
            zero, _ = WeightedSmithOracle(cut_matrix(values, False)).profile(data['energies'], [0]*25)
            return infinity+zero

    baseline = profile(source)
    reference = root_for_profiles((baseline[:25], baseline[25:]), ex, ey)
    classes, records = {baseline: reference}, []
    candidates = []
    for a in range(11):
        for b in range(1, 12-a):
            q = sum(1 << i for i in range(b+1) if i & b == i) << a
            for scalar in (1, 2, 3):
                candidates.append((dict(a=a, b=b, scalar=scalar),
                                   tuple(q if scalar >> p & 1 else 0 for p in range(2))))
    rng = random.Random(args.seed)
    for trial in range(args.random_trials):
        coefficients = [rng.randrange(4) for _ in range(11)]
        planes = tuple(sum(((v >> p) & 1) << i for i, v in enumerate(coefficients)) for p in range(2))
        candidates.append((dict(random_trial=trial), tuple(v ^ (v << 1) for v in planes)))
    for label, q in candidates:
        assert all(value.bit_count() % 2 == 0 for value in q)
        values = []
        for index, original in enumerate(source):
            x, y, u, v = ((index//5**i) % 5 for i in range(4))
            position = 25*(x+5*y)+v+5*u
            delta = tuple((value if first >> position & 1 else 0)
                          ^ (value << 4 if second >> position & 1 else 0) for value in q)
            values.append(tuple(a ^ b for a, b in zip(original, delta)))
        orders = profile(values)
        if orders not in classes:
            classes[orders] = root_for_profiles((orders[:25], orders[25:]), ex, ey)
        records.append(dict(**label, q_planes=list(q), profiles=[orders[:25], orders[25:]], **classes[orders]))
    names = ('search_joint_two_place_pencil.py', 'check_joint_zero_prefix.py',
             'joint_energy_conjugations.py', 'search_f4_recentered_crossings.py')
    result = dict(reference=reference, records=records, distinct_profiles=len(classes),
                  best=min(records, key=lambda record: record['rho']),
                  random_trials=args.random_trials, seed=args.seed,
                  source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                  scope='Finite F4 parameter catalogue at the inherited physical heights; pressure optimization uses floating arithmetic.',
                  factoring_exponent='1/2; unchanged')
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(json.dumps(dict(candidates=len(records), profiles=len(classes),
                          best_rho=result['best']['rho'], reference_rho=reference['rho'])))


if __name__ == '__main__':
    main()
