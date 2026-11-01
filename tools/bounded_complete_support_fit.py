#!/usr/bin/env python3
"""Bounded exploratory SLSQP fits with every rejected query retained."""
from __future__ import annotations

import math
import time
import numpy as np
from scipy.optimize import minimize

from normalized_support_gap import normalized_gap
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure


def fit(seed, start, *, seconds=15., maxiter=100, balanced=False, height_bound=128.):
    model = ScaledAsymmetricPressure(seed)
    started = time.monotonic()
    point = np.asarray(start, dtype=float).copy()
    if balanced:
        point[-1] = 0.
    best, last, queries, rejected = None, None, 0, []

    def objective(x):
        nonlocal best, last, queries
        if time.monotonic()-started >= seconds:
            raise TimeoutError('per-trial cap checked between complete queries')
        queries += 1
        try:
            rate, gradient, details = model.evaluate(x)
        except (ArithmeticError, ValueError, OverflowError) as error:
            rejected.append({'query': queries, 'parameters': x.tolist(), 'error': str(error)})
            raise
        value, derivative = normalized_gap(rate, x, gradient)
        last = {'rate': float(rate), 'objective': float(value), 'parameters': x.tolist(), 'details': details}
        if best is None or value < best['objective']:
            best = last
        return value, derivative

    result = None
    try:
        result = minimize(objective, point, jac=True, method='SLSQP',
            bounds=[(-height_bound, height_bound)]*(model.parameters-2)
                   +[(-9., math.log(4.)), (0., 0.) if balanced else (-2., 2.)],
            options={'maxiter': maxiter, 'ftol': 1e-10})
        status, message = 'optimizer_returned', str(result.message)
    except TimeoutError as error:
        status, message = 'budget_reached_partial', str(error)
    except (ArithmeticError, ValueError, OverflowError) as error:
        status, message = 'query_rejected_partial', str(error)
    return {'best': best, 'last': last, 'initial_parameters': point.tolist(),
            'queries': queries, 'rejected_queries': rejected, 'status': status, 'message': message,
            'optimizer_success': bool(result.success) if result is not None else False,
            'balanced_angle_fixed': balanced, 'seconds': time.monotonic()-started}
