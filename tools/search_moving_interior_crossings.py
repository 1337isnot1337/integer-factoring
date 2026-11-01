#!/usr/bin/env python3
"""Literal complete-profile screens that move arbitrary physical pairs.

Unequal heights and occupation numbers are allowed. The assigned physical
weight is charged through the actual full tensor; no compound shortcut is
used. Every fitted pressure is numerical evidence only.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations, product
import json
import math
from pathlib import Path
import time

from check_coherent_hot_boundary import logical, natural_profiles
from hot_code_charge_profiles import hot_energy
from moving_boundary_crossings import compress, reorder
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_moving_boundary_crossings import polynomial

ROOT = Path(__file__).resolve().parents[1]


def pressure(sides, energies):
    """Normalize exact scalar shifts before evaluating small Gibbs sums.

    logsumexp([0,-large,...]) can round to zero. log1p retains those
    terms, while exact minimum shifts avoid cancellation between cuts.
    """
    from scipy.optimize import brentq, minimize_scalar

    minima = [min(side) for side in sides]
    gaps = [[value - low for value in side] for side, low in zip(sides, minima)]
    energy_minimum = min(energies)
    energy_gaps = [value - energy_minimum for value in energies]
    slope = sum(minima) - 2 * energy_minimum
    assert slope <= 0  # Apply both entry-order bounds to the minimum-energy diagonal entry.

    def logarithm(values, scale):
        units = values.count(0)
        return math.log(units) + math.log1p(sum(math.exp(-scale * value) for value in values if value) / units)

    def objective(log_scale):
        scale = math.exp(log_scale)
        target = sum(logarithm(side, scale) for side in gaps) - scale * slope
        equation = lambda root: target - 2 * root * logarithm(energy_gaps, scale / root)
        upper = 4.
        while equation(upper) > 0:
            upper *= 2
        return brentq(equation, 1e-10, upper, xtol=1e-13) / 2

    fit = minimize_scalar(objective, bounds=(math.log(.005), math.log(1.5)), method='bounded',
                          options={'xatol': 1e-9, 'maxiter': 70})
    return {'rho': float(fit.fun), 'scale': math.exp(float(fit.x)),
            'completed': bool(fit.success), 'queries': int(fit.nfev)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field', type=int, choices=(2, 4), default=2)
    parser.add_argument('--degree', type=int, choices=(0, 1, 2, 3), default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    letters = [0, 1, 2, 3, 4, 5, 8, 6, 9]
    raw = restricted_tensor(16, list(map(logical, letters)))
    energies = list(map(hot_energy, letters))
    count = args.field ** (args.degree + 1)
    sources = ['tools/search_moving_interior_crossings.py', 'tools/moving_boundary_crossings.py',
               'tools/search_moving_boundary_crossings.py', 'tools/check_coherent_hot_boundary.py',
               'tools/restrict_implicit_reflection.py', 'tools/native_matrix_smith_f4.py',
               'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'coefficient_field': args.field, 'parameter_degree_bound': args.degree,
              'ambient_letters': letters, 'ambient_energies': energies, 'records': [],
              'completed': False, 'completed_cells': 0,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}}
    began = time.monotonic()
    best = None
    for first, second in combinations(range(9), 2):
        order = [i for i in range(9) if i not in (first, second)] + [first, second]
        original = reorder(raw, order)
        fixed = [energies[i] for i in order[:-2]]
        assigned = sorted({energies[first], (energies[first] + energies[second]) // 2, energies[second]})
        for ai, bi in product(range(count), repeat=2):
            a, b = polynomial(ai, args.field, args.degree), polynomial(bi, args.field, args.degree)
            seed = compress(original, 7, a, b)
            tensor_hash = digest(seed['values'])
            for energy in assigned:
                physical = fixed + [energy]
                sides = natural_profiles(seed, physical)
                fit = pressure(sides, physical)
                record = {'physical_pair': [letters[first], letters[second]], 'a': a, 'b': b,
                          'assigned_energy': energy, 'energies': physical, 'degree': seed['degree'],
                          'tensor_sha256': tensor_hash, 'orders_infinity_zero': sides,
                          'numerical_pressure': fit}
                report['records'].append(record)
                if best is None or fit['rho'] < best['numerical_pressure']['rho'] - 1e-12:
                    best = record
                    report['best'] = best
                    print(json.dumps({'evaluations': len(report['records']),
                                      'physical_pair': record['physical_pair'], 'a': a, 'b': b,
                                      'assigned_energy': energy, 'numerical_rho': fit['rho']}), flush=True)
        report.update(completed_cells=report['completed_cells'] + 1,
                      elapsed_seconds=time.monotonic() - began)
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        temporary.replace(args.output)
    report.update(completed=True, completed_evaluations=len(report['records']),
                  elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in ('completed', 'completed_cells', 'completed_evaluations', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    main()
