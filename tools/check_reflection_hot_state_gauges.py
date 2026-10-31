#!/usr/bin/env python3
"""Independent original-tensor checks for full reflection hot-state gauges."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from check_rectangular_defect_support import Rat, ZERO, ONE, Z, identity, mm
from check_rectangular_f4_repair import scale
from check_reflection_hot_core_profiles import polynomial_seed
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import gauged_kernel, logical_heights, physical_bit_energies
from reflection_hot_state_support import full_images, check_support
from reflection_hot_state_gauges import HotStatePressure, normalize_mask, physical_index
from relabelled_gaussian_sums import flatten, substitute_power


def transport_control(interval, ramification, core):
    """Build all transported CAR columns from the original alternating kernel."""
    m, size, p = interval // 2, interval, ramification
    c, s = ONE + Z ** (2 * p), Z ** (4 * p)
    R = [[s ** (i & ~j) for j in range(m)] for i in range(m)]
    L = list(map(list, zip(*R)))
    U, Ui, V = [[[ZERO] * size for _ in range(size)] for _ in range(3)]
    for i in range(m):
        U[i][m+i] = U[m+i][i] = ONE
        Ui[i][m+i] = Ui[m+i][i] = ONE
        V[i][i] = V[m+i][m+i] = ONE
        for j in range(m):
            U[m+i][m+j] = Ui[i][j] = c * L[i][j]
            V[i][m+j] = c * R[i][j]
    assert mm(U, Ui) == mm(Ui, U) == mm(V, V) == identity(size)
    sigma = [Z ** (p * (2 * interval - 2 - 4 * i)) for i in range(m)]
    sigma += [Z ** (2 * interval * p) * value for value in sigma]
    reverse = list(range(m-1, -1, -1)) + list(range(size-1, m-1, -1))
    expected = full_images(interval, p, core)
    identities = 0
    for place, legs in (('zero', ((0, 1), (2, 3))), ('infinity', ((2, 0), (3, 1)))):
        kernel = [[Rat(substitute_power(entry, p)) for entry in row]
                  for row in gauged_kernel(interval, place)]
        rows, columns = [[leg * m + i for leg in pair for i in range(m)] for pair in legs]
        cross = [[kernel[i][j] for j in columns] for i in rows]
        internal_left = [[kernel[i][j] for j in rows] for i in rows]
        internal_right = [[kernel[i][j] for j in columns] for i in columns]
        assert mm(mm(U, cross), V) == [[sigma[i] if i == j else ZERO for j in range(size)] for i in range(size)]
        us, sv = mm(U, internal_left), mm(internal_right, V)
        for j in range(core):
            for group, which in enumerate(('first', 'second')):
                original = group * m + m - 1 - j
                images = {
                    f'left_{which}{j}_create': ([U[i][original] for i in reverse], [ZERO] * size),
                    f'left_{which}{j}_annihilate': ([us[i][original] for i in reverse], [Ui[original][i] for i in reverse]),
                    f'right_{which}{j}_create': ([sigma[i] * V[i][original] for i in reverse],
                                                       [sv[original][i] / sigma[i] for i in reverse]),
                    f'right_{which}{j}_annihilate': ([ZERO] * size, [V[original][i] / sigma[i] for i in reverse])}
                for name, actual in images.items():
                    assert actual == expected[name], (interval, p, core, place, name)
                    identities += 2 * size
    return {'interval': interval, 'ramification': p, 'core': core,
            'original_kernel_operator_identities': identities}


def original_weighted(interval, ramification, factor, core, mask, offsets, log_s, requested_precision=512):
    values, clearing = polynomial_seed(interval, ramification, factor, core)
    m, q, cq = interval // 2, 1 << (interval // 2), 1 << core
    core_labels = [sum(((state >> (m - 1 - j)) & 1) << j for j in range(core)) for state in range(q)]
    if mask is not None:
        values = [scale(value, mask[physical_index(tuple(core_labels[(index >> (leg*m)) & (q-1)]
                                                         for leg in range(4)), cq)])
                  for index, value in enumerate(values)]
    for index, value in enumerate(values):
        x, y, u, v = tuple((index >> (leg*m)) & (q-1) for leg in range(4))
        endpoint = (value[0].bit_count() & 1) + 2 * (value[1].bit_count() & 1)
        assert endpoint == int(x == v and y == u)
    bits = max(1, max(value.bit_length() for pair in values for value in pair))
    precision = max(requested_precision, 1 << (bits - 1).bit_length())
    profiles = []
    for place, name in enumerate(('zero', 'infinity')):
        h = logical_heights(interval, name)
        heights = [[ramification * sum(h[leg*m+j] for j in range(m) if state >> j & 1)
                    + offsets[place, leg, core_labels[state]] for state in range(q)] for leg in range(4)]
        first, second = (((0, 1), (2, 3)) if place == 0 else ((0, 2), (1, 3)))
        rw = [heights[first[0]][a] + heights[first[1]][b] for a, b in product(range(q), repeat=2)]
        cw = [heights[second[0]][a] + heights[second[1]][b] for a, b in product(range(q), repeat=2)]
        profile = MatrixOracle(flatten(values, m, bool(place)), q*q, precision=precision).profile(
            np.asarray(rw)*log_s, np.asarray(cw)*log_s, infinity=bool(place), scale=log_s)['orders']
        profiles.append(profile)
    natural = ramification * np.asarray(physical_bit_energies(interval))
    gx, gy = [], []
    for state in range(q):
        energy = sum(natural[i] for i in range(m) if state >> i & 1)
        label = core_labels[state]
        gx.append(energy + offsets[:, (0, 3), label].sum())
        gy.append(energy + offsets[:, (1, 2), label].sum())
    log_b = sum(float(logsumexp(-np.asarray(profile))) for profile in profiles)
    def equation(rate):
        return log_b - rate * (logsumexp(-log_s * np.asarray(gx) / rate)
                               + logsumexp(-log_s * np.asarray(gy) / rate))
    upper = 4.
    while equation(upper) > 0:
        upper *= 2
    rho = brentq(equation, 1e-9, upper, xtol=1e-13) / 2
    return profiles, rho, clearing, len(values)


def literal_control(interval, ramification, core, factor, mask, seed, cases=3, internal_regular=False):
    rng = np.random.default_rng(seed)
    mask = normalize_mask(mask, 1 << core)
    model = HotStatePressure(interval, factor, ramification, core, mask, internal_regular=internal_regular)
    endpoints, maximum_error = 0, 0.
    for case in range(cases):
        offsets = rng.integers(-19, 20, size=(2, 4, 1 << core)) / 4.
        offsets[:, :, 0] = 0
        if internal_regular:
            bit_weights = rng.integers(0, 20, size=(2, 4, core)) / 4.
            offsets = np.stack([sum((bit_weights[:, :, j] for j in range(core) if state >> j & 1),
                                   np.zeros((2, 4))) for state in range(1 << core)], axis=2)
        log_s = (1., .375, .625)[case % 3]
        expected, rho, clearing, endpoints = original_weighted(interval, ramification, factor, core, mask, offsets, log_s)
        actual, detail = model.evaluate(log_s, offsets)
        for place in range(2):
            normalized = np.asarray(detail['active_orders'][place])
            data = model.records[place][0]['data']
            normalized += log_s * (data['scalar_order'] + (1 if place == 0 else -1) * clearing)
            for spectator in data['spectator_orders']:
                normalized = np.r_[normalized, normalized + log_s * spectator]
            assert np.array_equal(np.sort(normalized), np.sort(expected[place])), (interval, ramification, core, case, place)
        maximum_error = max(maximum_error, abs(actual-rho))
        assert abs(actual-rho) < 2e-12
    return {'interval': interval, 'ramification': ramification, 'core': core,
            'internal_gaussians_removed': internal_regular,
            'nonconstant_masks': sorted(key for key in factor if key), 'coefficient_mask': mask is not None,
            'weighted_complete_cut_profiles': 2 * cases, 'endpoint_coefficients': endpoints,
            'maximum_pressure_error': maximum_error, 'metadata': model.records[0][1]}


def derivative_control(factor, mask, seed):
    rng = np.random.default_rng(seed)
    model = HotStatePressure(4, factor, mask=mask)
    offsets = rng.normal(size=(2, 4, 2))
    offsets[:, :, 0] = 0
    log_s, epsilon = .317, 1e-6
    _, (gradient, ds), _ = model.evaluate(log_s, offsets, True)
    errors = []
    for place, leg in product(range(2), range(4)):
        direction = np.zeros_like(offsets);direction[place, leg, 1] = epsilon
        difference = (model.evaluate(log_s, offsets+direction)[0] - model.evaluate(log_s, offsets-direction)[0]) / (2*epsilon)
        errors.append(abs(difference-gradient[place, leg, 1]))
    difference = (model.evaluate(log_s+epsilon, offsets)[0] - model.evaluate(log_s-epsilon, offsets)[0]) / (2*epsilon)
    errors.append(abs(difference-ds))
    assert max(errors) < 3e-8
    return {'derivatives': len(errors), 'maximum_error': max(errors)}


def widening_control():
    factor = {0:{0:1},15:{0:2,1:2}}
    model = HotStatePressure(4,factor,internal_regular=True)
    offsets = np.zeros((2,4,2));offsets[:,:,1] = 1000
    expected, rho, clearing, _ = original_weighted(4,1,factor,1,None,offsets,.001,8192)
    actual, detail = model.evaluate(.001,offsets)
    for place in range(2):
        data = model.records[place][0]['data']
        spectrum = np.asarray(detail['active_orders'][place]) + .001 * (data['scalar_order'] + (1 if place == 0 else -1) * clearing)
        assert np.max(np.abs(np.sort(spectrum)-np.sort(expected[place]))) < 1e-12
    assert abs(actual-rho) < 1e-12
    windows = [oracle.precision for side in model.oracles for oracle in side]
    assert max(windows) > 512
    return {'initial_elimination_precision':512,'final_elimination_precisions':windows,
            'maximum_complete_source_bits':max(md['defect_source_maximum_bits'] for _,side in model.records for md in side),
            'literal_cut_profiles':2,'maximum_height':1000,'pressure_error':abs(actual-rho)}


def verify(extended=False):
    supports = [check_support(M, P, k) for M, P, k in
                ((4,1,1),(8,1,1),(16,1,1),(32,1,1),(64,1,1),(16,3,1),(32,1,2),(64,1,2))]
    supports += [check_support(M, P, k, True) for M, P, k in
                 ((4,1,1),(8,1,1),(16,1,1),(32,1,1),(64,1,1),(32,3,2))]
    transports = [transport_control(M, P, k) for M, P, k in
                  ((4,1,1),(8,1,1),(8,3,2),(16,1,2),(32,1,1))]
    gaussian = {0:{0:1}}
    quartic = {0:{0:1},15:{0:2,1:2}}
    mixed = {0:{0:1},3:{-1:3,1:3},10:{0:2,2:2}}
    mask = lambda x,y,u,v: 1 if (x,y)==(v,u) else (2 if x+y+u+v==2 else 3)
    literal = [literal_control(4,1,1,f,None,100+i) for i,f in enumerate((gaussian,quartic,mixed))]
    literal += [literal_control(4,3,1,mixed,mask,105),literal_control(8,1,1,mixed,None,106),
                literal_control(8,1,1,quartic,mask,107)]
    literal += [literal_control(M,1,1,mixed,None,120+M,internal_regular=True) for M in (4,8)]
    if extended:
        factor = {0:{0:1},5:{0:2,1:2},85:{-1:3,2:3},255:{1:1,3:1}}
        literal += [literal_control(8,1,2,factor,None,108,cases=2)]
        literal += [literal_control(8,1,2,factor,None,109,cases=2,internal_regular=True)]
    derivatives = [derivative_control(f, m, 200+i) for i,(f,m) in
                   enumerate(((gaussian,None),(mixed,None),(quartic,mask)))]
    fast = HotStatePressure(4,gaussian,internal_regular=True)
    invalid = np.zeros((2,4,2));invalid[0,0,1] = -1
    try:
        fast.evaluate(.3,invalid)
    except ValueError:
        pass
    else:
        raise AssertionError('The regular-internal branch accepted a negative bit height.')
    root = Path(__file__).resolve().parent.parent
    files = ['research/reflection-hot-state-gauges.md', 'tools/reflection_hot_state_support.py',
             'tools/reflection_hot_state_coordinates.py','tools/reflection_hot_state_gauges.py',
             'tools/native_fock_f4_parity.cpp','tools/check_reflection_hot_state_gauges.py']
    return {'status':'pass','extended':extended,'support_controls':supports,'transport_controls':transports,
            'original_polynomial_controls':literal,'derivative_controls':derivatives,
            'adaptive_elimination_control':widening_control(),
            'original_complete_weighted_cut_profiles':sum(r['weighted_complete_cut_profiles'] for r in literal),
            'endpoint_coefficient_equalities':sum(r['endpoint_coefficients'] for r in literal),
            'source_sha256':{name:sha256((root/name).read_bytes()).hexdigest() for name in files},
            'scope':'Exact regular support and original weighted profiles; floating pressure and derivative controls; no improved exponent claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended',action='store_true')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = verify(args.extended)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ('status','original_complete_weighted_cut_profiles','endpoint_coefficient_equalities')}))


if __name__=='__main__':
    main()
