#!/usr/bin/env python3
"""Check the resolved-subset pressure lower bound and the q94 exclusion."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import exp, isfinite, sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps

INPUT = Path('evidence/f8-hyperoval-q94-mode11-pilot.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
DROP = [90, 91]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    bad = [i for i, vector in enumerate(embedding)
           if any(word >> 11 & 1 for word in vector)]
    assert bad == DROP
    keep = [i for i in range(len(embedding)) if i not in bad]
    embedding = [embedding[i] for i in keep]
    retraction = [retraction[i] for i in keep]
    heights, charges = validate_maps(embedding, retraction, 11)
    assert heights == saved['heights']
    assert charges == saved['charges']
    assert min(heights) == 0
    assert max(heights) == 25
    assert len(embedding) == 94

    zero = np.asarray(saved['zero_orders_t'], dtype=float)
    infinity = np.asarray(saved['infinity_orders_t'], dtype=float)
    assert zero.shape == infinity.shape == (932,)
    assert (zero >= 0).all() and (infinity >= 0).all()
    h = np.asarray(heights, dtype=float)

    def target(scale: float) -> float:
        return float(logsumexp(-scale * zero) + logsumexp(-scale * infinity))

    def root(scale: float) -> float:
        equation = lambda rho: target(scale) - 4 * rho * float(logsumexp(-scale * h / rho))
        return float(brentq(equation, 1e-7, 2.0, xtol=1e-13))

    for scale in (0.05, 0.1, 0.2, 0.25, 0.5, 1.0):
        rho = root(scale)
        assert isfinite(rho) and rho > 0
        # F_t'(rho) = 4 H_t(rho) + 4 E_t[h]/rho > 0.
        weights = np.exp(-scale * h / rho)
        weights /= weights.sum()
        hterm = float(logsumexp(-scale * h / rho))
        derivative = 4 * hterm + 4 * float(weights @ h) / rho
        assert derivative > 0
        # Pointwise monotonicity: raising the target raises the root.
        target_plus = target(scale) + 1.0
        rho_plus = float(brentq(
            lambda r: target_plus - 4 * r * float(logsumexp(-scale * h / r)), 1e-7, 2.0))
        assert rho_plus > rho

    saved_fit = saved['numerical_pressure']
    pilot_scale = saved_fit['scale_t']
    pilot_rho = saved_fit['rho']
    assert abs(root(pilot_scale) - pilot_rho) < 1e-12
    assert abs(target(pilot_scale)
               - 4 * pilot_rho * float(logsumexp(-pilot_scale * h / pilot_rho))) < 1e-11
    assert saved_fit['completed'] is True
    assert pilot_rho > 1 / sqrt(2)

    receipt = {
        'passed': True,
        'q': 94,
        'modes': 11,
        'height_zero_present': True,
        'orders_nonnegative': True,
        'pilot_orders_per_cut': 932,
        'tail_counts': saved['tail_counts'],
        'pilot_rho': pilot_rho,
        'pilot_scale_t': pilot_scale,
        'excess_over_inverse_sqrt_two': pilot_rho - 1 / sqrt(2),
        'complete_pressure_lower_bound': pilot_rho,
        'complete_cannot_cross': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'scope': (
            'Replay of the recorded q94 cutoff-47 orders and the analytic '
            'resolved-subset lower bound. The pilot is above 1/sqrt(2), so '
            'the complete q94 pressure is above it too; no full q94 build is '
            'needed for the crossing question. No score change.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'q': 94, 'pilot_rho': pilot_rho,
                      'excess': pilot_rho - 1 / sqrt(2), 'cannot_cross': True}))


if __name__ == '__main__':
    main()
