#!/usr/bin/env python3
"""Complete finite rank-one boundary-retraction screen for hot reflection codes.

All lower-energy states are retained. Every projective right vector and
every left vector pairing to one is tested at the final energy. Complete
compound profiles compute both natural cuts; numerical costs are search
evidence. Each improving tensor is independently constructed from the
original restricted Pfaffian tensor before it is retained.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from coherent_hot_boundary import compress
from hot_code_charge_profiles import hot_energy
from hot_compound_crossing import compound_crossing
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import field_dot, flatten
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_fused_mixed_weighted import atomic_write
from search_physical_code_gauges import start_params


def profile(seed, charges=None):
    matrix, rows, _ = flatten(seed, 0)
    if charges is None:
        return sorted(MatrixOracle(matrix, len(matrix)).profile(scale=2)['orders'])
    result = []
    for charge in sorted({charges[x] + charges[y] for x, y in rows}):
        keep = [i for i, (x, y) in enumerate(rows) if charges[x] + charges[y] == charge]
        block = [[matrix[i][j] for j in keep] for i in keep]
        result.extend(MatrixOracle(block, len(keep)).profile(scale=2)['orders'])
    # The paired code uses the same state charges, so row and column
    # labels have the same charge in flatten's two enumeration orders.
    for i, (x, y) in enumerate(rows):
        for j, (u, v) in enumerate(rows):
            if charges[x] + charges[y] != charges[u] + charges[v]:
                assert matrix[i][j] == (0, 0)
    return sorted(result)


def pressure(sides, energies):
    sides = [np.asarray(side) for side in sides]
    energies = np.asarray(energies)
    def objective(log_scale):
        scale = math.exp(log_scale)
        target = sum(float(logsumexp(-scale * side)) for side in sides)
        equation = lambda root: target - 2 * root * float(logsumexp(-scale * energies / root))
        return brentq(equation, 1e-10, 4, xtol=1e-13) / 2
    fit = minimize_scalar(objective, bounds=(math.log(.005), math.log(1.5)), method='bounded',
                          options={'xatol': 1e-9, 'maxiter': 70})
    return {'rho': float(fit.fun), 'scale': math.exp(float(fit.x)),
            'completed': bool(fit.success), 'queries': int(fit.nfev)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes', type=int, nargs='+', default=[10, 16])
    parser.add_argument('--seconds', type=float, default=480)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    labels = sorted(range(256), key=lambda value: (hot_energy(value), value))
    names = (Path(__file__).name, 'coherent_hot_boundary.py', 'hot_compound_crossing.py',
             'hot_code_charge_profiles.py', 'hot_code_charge_profiles.cpp',
             'physical_split_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'restrict_implicit_reflection.py')
    report = {'scope': __doc__, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'arguments': {'sizes': args.sizes, 'seconds': args.seconds},
              'records': [], 'winners': {}, 'completed': False, 'unresolved': 0}
    began, deadline = time.monotonic(), time.monotonic() + args.seconds
    exhausted = True
    for q in args.sizes:
        level = hot_energy(labels[q - 1])
        lower = [x for x in labels if hot_energy(x) < level]
        boundary = [x for x in labels if hot_energy(x) == level]
        if q != len(lower) + 1 or level >= 64:
            raise ValueError('Use the first retained letter of a boundary below energy 64.')
        letters = lower + boundary
        base = compound_crossing(letters)
        logical = lambda mask: sum(1 << (7 - a) for a in range(8) if mask >> a & 1)
        raw = restricted_tensor(16, list(map(logical, letters)))
        energies = list(map(hot_energy, lower)) + [level]
        charges = ([x.bit_count() for x in lower] + [boundary[0].bit_count()]
                   if len({x.bit_count() for x in boundary}) == 1 else None)
        rights = [v for v in product(range(4), repeat=len(boundary))
                  if any(v) and next(a for a in v if a) == 1]
        for right in rights:
            for left in product(range(4), repeat=len(boundary)):
                if field_dot(left, right) != 1:
                    continue
                if time.monotonic() >= deadline:
                    exhausted = False
                    break
                started = time.monotonic()
                normal = compress(base, len(lower), [left], [right])
                swapped = compress(base, len(lower), [left], [right], swap_y=True)
                sides = [profile(swapped, charges), profile(normal, charges)]
                assert all(value == round(value) for side in sides for value in side)
                sides = [[int(v) for v in side] for side in sides]
                fit = pressure(sides, energies)
                record = {'q': q, 'energy': level, 'lower_letters': lower,
                          'boundary_letters': boundary, 'left': left, 'right': right,
                          'energies': energies, 'orders_infinity_zero': sides,
                          'pressure': fit, 'seconds': time.monotonic() - started}
                if str(q) not in report['winners'] or fit['rho'] < report['winners'][str(q)]['pressure']['rho']:
                    seed = compress(raw, len(lower), [left], [right])
                    content = min((part & -part).bit_length() - 1
                                  for value in seed['values'] for part in value if part)
                    seed['values'] = [(a >> content, b >> content) for a, b in seed['values']]
                    seed['degree'] = max(max(a.bit_length(), b.bit_length()) - 1 for a, b in seed['values'])
                    actual = RectangularOracle(seed['values'], q, q).prefixes(start_params(np.asarray(energies), 1))['orders']
                    a = 128 - content
                    assert sorted(actual[:q*q]) == [v - a for v in sides[0]]
                    assert sorted(actual[q*q:]) == [v + a for v in sides[1]]
                    record.update(primitive_degree=seed['degree'], monomial_content=content,
                                  tensor_sha256=digest(seed['values']), literal_source_replay=True)
                    report['winners'][str(q)] = record
                    print(json.dumps({'q': q, 'left': left, 'right': right,
                                      'rho': fit['rho'], 'degree': seed['degree']}), flush=True)
                report['records'].append(record)
                report['elapsed_seconds'] = time.monotonic() - began
                atomic_write(args.output, report)
            if not exhausted:
                break
        if not exhausted:
            break
    report.update(completed=True, finite_list_exhausted=exhausted,
                  elapsed_seconds=time.monotonic() - began,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
