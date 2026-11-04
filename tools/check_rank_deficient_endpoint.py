#!/usr/bin/env python3
"""Deterministic screen of rank-one endpoint perturbations of the q4 seed.

Subtract a rank-one F4 tensor at degree zero and restore it at degree d.
The endpoint at w=1 is unchanged, while the constant term loses rank.
The natural pressure is fitted at the recorded q4 gauge. This is a finite
numerical screen, not a global optimum or a pressure certificate.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np
from scipy.optimize import minimize_scalar

from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from search_dense_coupled_reflection_frames import profile_rho
from search_physical_code_gauges import start_params

DEGREES = (1, 2, 3, 4, 6, 8)
TRIALS_PER_DEGREE = 1500
ENERGY = np.asarray([0, 12, 4, 16], float)
BASE_RHO = 0.7310809359736006
TOLERANCE = 1e-12


def natural_rho(values):
    oracle = RectangularOracle(values, 4, 4)
    unit = start_params(ENERGY, 1.0)
    orders = oracle.evaluate(unit, [.25] * 4, [.25] * 4, want_spectra=True)[2]
    fit = minimize_scalar(lambda s: profile_rho(orders, ENERGY, s),
                          bounds=(.002, 1), method='bounded',
                          options={'xatol': 1e-8})
    return float(fit.fun), float(fit.x)


def rank_one(rng, size=4):
    vectors = [[rng.randrange(4) for _ in range(size)] for _ in range(4)]
    vx, vy, vu, vv = vectors
    return [F4_PRODUCT[F4_PRODUCT[F4_PRODUCT[vx[x]][vy[y]]][vu[u]]][vv[v]]
            for x in range(size) for y in range(size)
            for u in range(size) for v in range(size)]


def perturb(base, tensor_value, degree):
    values = []
    for (low, high), entry in zip(base['values'], tensor_value):
        low_bit, high_bit = entry & 1, entry >> 1
        values.append((low ^ low_bit ^ (low_bit << degree),
                       high ^ high_bit ^ (high_bit << degree)))
    return dict(base, values=values, degree=max(base['degree'], degree),
                monomial_content=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    base = tensor(4)
    base_rho, base_scale = natural_rho(base['values'])
    if abs(base_rho - BASE_RHO) >= 1e-9:
        raise SystemExit('The q4 base pressure did not replay.')
    rng = random.Random(20261102)
    records = []
    for degree in DEGREES:
        best, best_trial = base_rho, None
        for trial in range(TRIALS_PER_DEGREE):
            tensor_value = rank_one(rng)
            for _ in range(rng.randrange(0, 3)):
                tensor_value = [left ^ right for left, right
                                in zip(tensor_value, rank_one(rng))]
            seed = perturb(base, tensor_value, degree)
            value, scale = natural_rho(seed['values'])
            if value < best - TOLERANCE:
                best, best_trial = value, trial
        if best < base_rho - TOLERANCE:
            raise SystemExit(f'Degree {degree} found a value below the base.')
        records.append({'compensation_degree': degree,
                        'trials': TRIALS_PER_DEGREE,
                        'best_natural_rho': best,
                        'best_trial': best_trial,
                        'improved': best_trial is not None})

    names = ('check_rank_deficient_endpoint.py', 'native_rectangular_smith_f4.py',
             'physical_split_crossings.py', 'reflection_dual_rail_gaussian.py',
             'search_dense_coupled_reflection_frames.py',
             'search_physical_code_gauges.py')
    result = {'status': 'pass',
              'scope': 'finite deterministic endpoint-perturbation screen; no pressure certificate',
              'base_natural_rho': base_rho,
              'base_natural_scale': base_scale,
              'random_seed': 20261102,
              'trials_per_degree': TRIALS_PER_DEGREE,
              'degrees': list(DEGREES),
              'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'base_natural_rho': base_rho,
                      'best_by_degree': [[row['compensation_degree'],
                                          row['best_natural_rho']]
                                         for row in records]}))


if __name__ == '__main__':
    main()
