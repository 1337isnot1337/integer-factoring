#!/usr/bin/env python3
"""Finite unequal coherent-code screen with fully charged asymmetric costs."""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp
from certify_hot_charge_pressure import homogeneous_sum, integer_digest
from check_coherent_hot_boundary import full_maps, logical
from hot_code_charge_profiles import hot_energy
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import retract
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_fused_mixed_weighted import atomic_write


def natural_parameters(ex, ey):
    energies = [ex, ey, ey, ex]
    values = []
    for side in (0, 1):
        for leg, energy in enumerate(energies):
            values += [x / 2 for x in energy[1:]] if leg in ((2, 3) if side == 0 else (0, 1)) else [0] * (len(energy) - 1)
    return values + [1]


def support(sides, ex, ey, scale, p, q, h):
    a, b = p / h, q / h
    raw = a * float(logsumexp(-scale * np.asarray(sides[0]) / a))
    raw += b * float(logsumexp(-scale * np.asarray(sides[1]) / b))
    first, second = np.asarray(ex), np.asarray(ey)
    equation = lambda R: raw - R * (float(logsumexp(-scale * first / R)) + float(logsumexp(-scale * second / R)))
    return float(brentq(equation, 1e-10, 4, xtol=1e-13))


def exact_point(sides, ex, ey, scale, p, q, h):
    assert p*p + q*q == h*h
    z = Fraction(math.exp(-scale / (p*q))).limit_denominator(10000)
    n, d = z.numerator, z.denominator
    assert 0 < n < d
    pi, di = homogeneous_sum([q*h*x for x in sides[0]], n, d)
    pz, dz = homogeneous_sum([p*h*x for x in sides[1]], n, d)
    qx, dx = homogeneous_sum([p*q*x for x in ex], n, d)
    qy, dy = homogeneous_sum([p*q*x for x in ey], n, d)
    left, right = pi**p * pz**q, (qx*qy)**h
    difference = p*di + q*dz - h*(dx+dy)
    if difference >= 0:
        right *= d**difference
    else:
        left *= d**(-difference)
    return {'weights': [p, q, h], 'base': [n, d], 'intrinsic_scale': p*q*math.log(d/n),
            'strict_support_below_one': left < right,
            'comparison': 'P_infinity(z^(qh))^p P_zero(z^(ph))^q < (Q_X(z^(pq))Q_Y(z^(pq)))^h',
            'integer_bits': [left.bit_length(), right.bit_length()],
            'left_sha256': integer_digest(left), 'right_sha256': integer_digest(right),
            'absolute_difference_sha256': integer_digest(abs(left-right))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    labels = sorted(range(256), key=lambda x: (hot_energy(x), x))
    branches = [(3,4,5), (4,3,5), (5,12,13), (12,5,13), (20,21,29), (21,20,29)]
    names = (Path(__file__).name, 'physical_split_crossings.py', 'restrict_implicit_reflection.py',
             'reflection_dual_rail_gaussian.py', 'check_schur_gaussian_reduction.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'check_coherent_hot_boundary.py', 'certify_hot_charge_pressure.py')
    report = {'scope': 'Finite unequal coherent tensors with complete integer natural profiles and fully charged asymmetric Gibbs costs. Exact comparisons apply only at retained rational points.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'records': [], 'completed': False, 'branch_weights': branches}
    began = time.monotonic()
    for size, left_x, right_x in ((8, [1,1], [0,1]), (16, [1,1,0], [0,1,1])):
        level = hot_energy(labels[size-1])
        lower = [x for x in labels if hot_energy(x) < level]
        boundary = [x for x in labels if hot_energy(x) == level]
        b = len(boundary)
        raw = restricted_tensor(16, list(map(logical, lower+boundary)))
        identity = [[int(i == j) for j in range(b)] for i in range(b)]
        y_choices = [('lower_only', [], []), ('full_boundary', identity, identity)]
        if b == 3:
            for chosen in combinations(range(b), 2):
                selected = [identity[i] for i in chosen]
                y_choices.append(('coordinate_' + ''.join(map(str, chosen)), selected, selected))
            y_choices += [('coherent_01', [[1,1,0],[0,0,1]], [[0,1,0],[0,0,1]]),
                          ('coherent_12', [[1,0,0],[0,1,1]], [[1,0,0],[0,1,0]])]
        for x_name, xl, xr in (('coherent', [left_x], [right_x]),
                              ('coordinate', [identity[-1]], [identity[-1]])):
            lx, rx = full_maps(len(lower), b, xl, xr)
            for y_name, yl, yr in y_choices:
                started = time.monotonic()
                ly, ry = full_maps(len(lower), b, yl, yr)
                seed = retract(raw, lx, rx, ly, ry)
                qx, qy = seed['q_x'], seed['q_u']
                assert qx != qy
                content = min((part & -part).bit_length()-1 for value in seed['values'] for part in value if part)
                seed['values'] = [(a>>content, b>>content) for a,b in seed['values']]
                seed['degree'] = max(max(a.bit_length(), b.bit_length())-1 for a,b in seed['values'])
                ex = list(map(hot_energy, lower)) + [level]*len(xl)
                ey = list(map(hot_energy, lower)) + [level]*len(yl)
                values = RectangularOracle(seed['values'], qx, qy).prefixes(natural_parameters(ex,ey))['orders']
                n = qx*qy
                full = [sorted(values[:n]), sorted(values[n:])]
                assert min(full[0])+min(full[1]) == 0 and all(v == round(v) for side in full for v in side)
                sides = [[int(v-min(side)) for v in side] for side in full]
                trials = []
                for p,q,h in branches:
                    fit = minimize_scalar(lambda u: support(sides,ex,ey,math.exp(u),p,q,h),
                                          bounds=(math.log(.005),math.log(1.5)), method='bounded',
                                          options={'xatol':1e-8,'maxiter':60})
                    trials.append({'weights':[p,q,h],'support':float(fit.fun),
                                   'scale':math.exp(float(fit.x)),'completed':bool(fit.success)})
                best = min(trials,key=lambda r:r['support'])
                exact = exact_point(sides,ex,ey,best['scale'],*best['weights'])
                row = {'source_size':size,'x_kind':x_name,'y_kind':y_name,
                       'lower_letters':lower,'boundary_letters':boundary,
                       'x_left':xl,'x_right':xr,'y_left':yl,'y_right':yr,
                       'q_x':qx,'q_y':qy,'degree':seed['degree'],'scalar_content':content,
                       'tensor_sha256':digest(seed['values']),'energies_x':ex,'energies_y':ey,
                       'regular_profiles_infinity_zero':sides,'opposite_shifts':[min(side) for side in full],
                       'branch_trials':trials,'best_numerical_support':best,'exact_point':exact,
                       'seconds':time.monotonic()-started}
                report['records'].append(row)
                report['elapsed_seconds'] = time.monotonic()-began
                atomic_write(args.output,report)
                print(json.dumps({'dimensions':[qx,qy],'x':x_name,'y':y_name,'best_support':best['support'],
                                  'exact_below_one':exact['strict_support_below_one']}),flush=True)
                if exact['strict_support_below_one']:
                    report.update(completed=True,stopping_reason='Exact rational support below one; independent theorem audit required.')
                    atomic_write(args.output,report)
                    return
    report.update(completed=True,elapsed_seconds=time.monotonic()-began,
                  stopping_reason='Finite tensor/branch list exhausted; no retained rational point passes one.')
    atomic_write(args.output,report)


if __name__ == '__main__':
    main()
