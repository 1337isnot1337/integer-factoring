#!/usr/bin/env python3
"""Screen the complete q96 profiles under unequal cut scales.

For a complete two-cut profile this computes the Gibbs-optimized
support value of the two-place asymmetric criterion.  The balanced
calibration is checked against the certified natural pressures.  The
interior screen keeps both supporting weights away from zero; a short
boundary sequence records the degenerate one-cut layer separately.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

PROFILES = (
    ('f8_hyperoval_full', 'evidence/f8-hyperoval-q96-full-profile.json', None),
    ('coherent_four_boundary', 'evidence/coherent-q96-four-boundary-candidate.json', None),
    ('coherent_coordinate', 'evidence/coherent-q96-coordinate-reference.json', None),
    ('coherent_hot_f4_hyperoval', 'evidence/coherent-hot-f4-hyperoval-q96.json', None),
    ('f8_retraction_678', 'evidence/f8-hyperoval-retraction-678.json', None),
    ('f8_retraction_789', 'evidence/f8-hyperoval-retraction-789.json', None),
)

EXPECTED = {
    'f8_hyperoval_full': ('certificate', 'evidence/f8-hyperoval-q96-pressure.json'),
    'coherent_four_boundary': ('pressure', 'evidence/coherent-q96-four-boundary-candidate.json'),
    'coherent_coordinate': (None, None),
    'coherent_hot_f4_hyperoval': ('pressure', 'evidence/coherent-hot-f4-hyperoval-q96.json'),
    'f8_retraction_678': ('numerical_pressure', 'evidence/f8-hyperoval-retraction-678.json'),
    'f8_retraction_789': ('numerical_pressure', 'evidence/f8-hyperoval-retraction-789.json'),
}


def load_profile(path):
    data = json.loads(Path(path).read_text())
    if 'profile' in data and isinstance(data['profile'], dict):
        data = data['profile']
    return {
        'zero': [int(x) for x in data['zero_orders_t']],
        'infinity': [int(x) for x in data['infinity_orders_t']],
        'heights': [int(x) for x in data['heights']],
    }


def expected_rho(name):
    kind, path = EXPECTED[name]
    if kind is None:
        return None
    data = json.loads(Path(path).read_text())
    if kind == 'certificate':
        numerator, denominator = data['certificate']['rho_strict_bound']
        return numerator / denominator
    value = data[kind]['rho'] if kind != 'pressure' else data['pressure']['rho']
    return float(value)


def support_root(profile, log_scale, angle):
    scale = math.exp(log_scale)
    a, b = math.cos(angle), math.sin(angle)
    if min(a, b) <= 0:
        raise ValueError('both supporting coefficients must be positive')
    zero = np.asarray(profile['zero'], dtype=float)
    infinity = np.asarray(profile['infinity'], dtype=float)
    heights = 2.0 * scale * np.asarray(profile['heights'], dtype=float)
    raw = a * float(logsumexp(-scale * infinity / a))
    raw += b * float(logsumexp(-scale * zero / b))

    def equation(rate):
        return raw - rate * (
            float(logsumexp(-heights / rate))
            + float(logsumexp(-heights / rate))
        )

    if equation(1.0) > 0:
        low, high = 1.0, 2.0
        while equation(high) > 0:
            high *= 2.0
            if high > 1e9:
                raise ArithmeticError('support root was not bracketed above one')
    else:
        low, high = 1e-12, 1.0
        if equation(low) <= 0:
            raise ArithmeticError('support root was not bracketed below one')
    return float(brentq(equation, low, high, xtol=1e-13, rtol=1e-15))


def optimize_scale(profile, angle, bounds=(-7.0, 3.0)):
    result = minimize_scalar(
        lambda u: support_root(profile, u, angle),
        bounds=bounds,
        method='bounded',
        options={'xatol': 1e-12},
    )
    return float(result.fun), float(result.x)


def balanced(profile):
    support, log_scale = optimize_scale(profile, math.pi / 4)
    return {
        'support': support,
        'rho': support / math.sqrt(2.0),
        'log_scale': log_scale,
        'scale': math.exp(log_scale),
    }


def refine(profile, log_scale, angle):
    for _ in range(8):
        log_scale = optimize_scale(profile, angle)[1]
        result = minimize_scalar(
            lambda theta: support_root(profile, log_scale, theta),
            bounds=(0.05, math.pi / 2 - 0.05),
            method='bounded',
            options={'xatol': 1e-12},
        )
        angle = float(result.x)
    return support_root(profile, log_scale, angle), log_scale, angle


def interior(profile):
    angles = np.linspace(0.05, math.pi / 2 - 0.05, 121)
    best = None
    for angle in angles:
        support, log_scale = optimize_scale(profile, float(angle))
        if best is None or support < best[0]:
            best = (support, log_scale, float(angle))
    support, log_scale, angle = refine(profile, best[1], best[2])
    return {
        'support': support,
        'log_scale': log_scale,
        'scale': math.exp(log_scale),
        'angle': angle,
        'a': math.cos(angle),
        'b': math.sin(angle),
    }


def boundary(profile):
    rows = []
    for angle in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2):
        support, log_scale = optimize_scale(profile, angle, bounds=(-17.0, 3.0))
        rows.append({
            'angle': angle,
            'support': support,
            'log_scale': log_scale,
            'scale': math.exp(log_scale),
        })
    return rows


def screen():
    began = time.monotonic()
    records = []
    for name, path, _ in PROFILES:
        profile = load_profile(path)
        bal = balanced(profile)
        fit = interior(profile)
        edge = boundary(profile)
        expected = expected_rho(name)
        record = {
            'name': name,
            'profile_path': path,
            'profile_sha256': sha256(Path(path).read_bytes()).hexdigest(),
            'channel_counts': [len(profile['zero']), len(profile['infinity'])],
            'height_count': len(profile['heights']),
            'expected_rho': expected,
            'balanced': bal,
            'calibration_error': None if expected is None else abs(bal['rho'] - expected),
            'interior': fit,
            'boundary': edge,
        }
        assert abs(fit['a'] ** 2 + fit['b'] ** 2 - 1.0) < 1e-12
        assert fit['support'] > 1.0, (name, fit['support'])
        assert all(row['support'] > 1.0 for row in edge), name
        if expected is not None:
            assert record['calibration_error'] < 1e-9, (name, record['calibration_error'])
        records.append(record)
    best = min(records, key=lambda row: row['interior']['support'])
    return {
        'scope': (
            'Deterministic numerical screen of six complete q96 profiles under the '
            'two-place asymmetric supporting-line criterion. The balanced calibration '
            'matches certified or recorded natural pressures. Interior angles keep '
            'both weights at least 0.05. No global optimality or exact pressure '
            'certificate is claimed; the factoring exponent is unchanged.'
        ),
        'profiles': records,
        'interior_best': {
            'name': best['name'],
            'support': best['interior']['support'],
            'angle': best['interior']['angle'],
            'a': best['interior']['a'],
            'b': best['interior']['b'],
            'scale': best['interior']['scale'],
        },
        'grid': {
            'angle_start': 0.05,
            'angle_stop': math.pi / 2 - 0.05,
            'angle_count': 121,
            'scale_bounds_log': [-7.0, 3.0],
            'boundary_angles': [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2],
        },
        'elapsed_seconds': time.monotonic() - began,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = screen()
    report['completed'] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(
        'best interior support {:.12f} on {} at angle {:.9f}'.format(
            report['interior_best']['support'],
            report['interior_best']['name'],
            report['interior_best']['angle'],
        )
    )


if __name__ == '__main__':
    main()
