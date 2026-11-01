#!/usr/bin/env python3
"""Fixed-type proximal bundles with Gibbs and supporting-angle updates.

Adapt search_rectangular_core_bundle.py to actual asymmetric profile
pressure on the mixed graph tensors. All query points and rejected
coordinates are retained. QP solutions and stops are numerical only.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize, minimize_scalar

from asymmetric_profile_bundle import AsymmetricBundleOracle, raw_from_scaled, scaled_from_raw
from mixed_graph_retractions import FixedProfileSupport, literal_digest, tensor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.output.exists() or args.output.with_suffix('.jsonl').exists():
        parser.error('Refusing to overwrite a previous run.')
    root = Path(__file__).resolve().parents[1]
    source = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True):
        raise ValueError('Run from a clean frozen source.')
    report = {'status': 'running_incomplete', 'source_commit': source, 'manifest': manifest,
              'manifest_sha256': sha256(args.manifest.read_bytes()).hexdigest(), 'cases': [],
              'scope': 'Numerical fixed-type convex bundle steps and Gibbs/angle updates. '
                       'No strict pressure, global optimum or factoring improvement.'}
    begun = time.monotonic()

    def save():
        report['seconds'] = time.monotonic()-begun
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(report, indent=2)+'\n')
        temporary.replace(args.output)

    with args.output.with_suffix('.jsonl').open('x', buffering=1) as trace:
        def record(data):
            trace.write(json.dumps(data, separators=(',', ':'))+'\n')

        for registered in manifest['cases']:
            start = time.monotonic()
            deadline = start+manifest['seconds_per_case']
            literal = tensor(registered['coefficients'])
            assert literal_digest(literal) == registered['values_sha256']
            best = dict(rate=registered['rate'], parameters=registered['parameters'])
            state = {'label': registered['label'], 'status': 'running_incomplete', 'best': best,
                     'queries': 0, 'rejections': 0, 'master_failures': 0,
                     'minimum_heldout_plane_slack': 0., 'heldout_directed_comparisons': 0,
                     'epochs_completed': 0}
            report['cases'].append(state)
            save()
            rng = np.random.default_rng(registered['seed'])
            raw = raw_from_scaled(best['parameters'])
            logtan = float(best['parameters'][-1])
            bounds = Bounds(np.r_[np.full(32, -256.), math.exp(-9), -1e7],
                            np.r_[np.full(32, 256.), 4., 1e7])
            normal = np.zeros((64, 34))
            for j in range(32):
                normal[2*j, j], normal[2*j, 32] = 1., 64.
                normal[2*j+1, j], normal[2*j+1, 32] = -1., 64.
            metric = np.r_[np.ones(32), 22.]
            try:
                for epoch in range(manifest['epochs_ceiling']):
                    if time.monotonic() >= deadline:
                        raise TimeoutError
                    if epoch and epoch % 3 == 0:
                        raw = raw_from_scaled(best['parameters'])
                        raw[:-1] += rng.normal(0, .015, 32)
                        raw[:-1] = np.clip(raw[:-1], -63*raw[-1], 63*raw[-1])
                        logtan = float(best['parameters'][-1])
                    oracle = AsymmetricBundleOracle(literal, logtan)
                    cuts = []

                    def query(point):
                        nonlocal best
                        if time.monotonic() >= deadline:
                            raise TimeoutError
                        try:
                            result = oracle.query(point)
                        except (ArithmeticError, ValueError, OverflowError) as error:
                            state['rejections'] += 1
                            record({'case': registered['label'], 'epoch': epoch, 'status': 'query_rejected',
                                    'raw_parameters': point.tolist(), 'log_tangent': logtan, 'error': str(error)})
                            raise
                        rate, f, g, details, physical = result
                        for oldx, oldf, oldg in cuts:
                            slack = min(f-oldf-oldg@(point-oldx), oldf-f-g@(oldx-point))
                            state['minimum_heldout_plane_slack'] = min(state['minimum_heldout_plane_slack'], float(slack))
                            state['heldout_directed_comparisons'] += 2
                            if slack < -2e-7:
                                raise RuntimeError(f'Convex supporting-plane control failed: {slack}.')
                        cuts.append((point.copy(), f, g.copy()))
                        state['queries'] += 1
                        record({'case': registered['label'], 'epoch': epoch, 'status': 'complete_query',
                                'raw_parameters': point.tolist(), 'log_tangent': logtan,
                                'log_profile_value': f, 'subgradient': g.tolist(),
                                'physical_vector': physical.tolist(), 'rate': rate, 'details': details})
                        if rate < best['rate']:
                            best = {'rate': rate, 'parameters': scaled_from_raw(point, logtan).tolist(),
                                    'details': details}
                            state['best'] = best
                            print(json.dumps({'case': registered['label'], 'rate': rate,
                                              'query': state['queries'], 'seconds': time.monotonic()-start}), flush=True)
                        return result

                    result = query(raw)
                    f, physical, entropy = result[1], result[4], result[3]['entropy_sum']
                    anchor, current, mu, serious = raw.copy(), (f+physical@raw)/entropy, .03, 0
                    for iteration in range(manifest['iterations_per_epoch']):
                        if time.monotonic() >= deadline:
                            raise TimeoutError
                        A = np.asarray([np.r_[-g, 1.] for _, _, g in cuts])
                        lower = np.asarray([f-g@x for x, f, g in cuts])
                        initial = np.r_[anchor, max(f+g@(anchor-x) for x, f, g in cuts)]

                        def objective(x):
                            delta = (x[:33]-anchor)*metric
                            return (x[33]+physical@x[:33])/entropy+.5*mu*(delta@delta)

                        def derivative(x):
                            return np.r_[physical/entropy+mu*metric**2*(x[:33]-anchor), 1/entropy]

                        opt = minimize(objective, initial, jac=derivative, method='SLSQP',
                            bounds=bounds, constraints=[LinearConstraint(A, lower, np.inf),
                                                       LinearConstraint(normal, 0., np.inf)],
                            options={'maxiter': 200, 'ftol': 1e-11})
                        if (not np.all(np.isfinite(opt.x)) or min(A@opt.x-lower) < -1e-6
                                or min(normal@opt.x) < -1e-6):
                            state['master_failures'] += 1
                            record({'case': registered['label'], 'epoch': epoch, 'iteration': iteration,
                                    'status': 'master_unresolved', 'message': str(opt.message)})
                            break
                        candidate = opt.x[:33]
                        predicted = current-(opt.x[33]+physical@candidate)/entropy
                        if predicted < 2e-9 or max(abs(candidate-anchor)) < 1e-9:
                            result = query(anchor)
                            change = max(abs(result[4]-physical))
                            physical, entropy = result[4], result[3]['entropy_sum']
                            current = (result[1]+physical@anchor)/entropy
                            if change < 1e-7:
                                break
                            mu = .03
                            continue
                        try:
                            result = query(candidate)
                        except (ArithmeticError, ValueError, OverflowError):
                            mu = min(100., mu*4)
                            continue
                        fixed = (result[1]+physical@candidate)/entropy
                        actual = current-fixed
                        take = actual >= .2*predicted
                        if take:
                            anchor, current, serious = candidate.copy(), fixed, serious+1
                            mu = max(.0001, mu*.9)
                            if serious % 6 == 0:
                                physical, entropy = result[4], result[3]['entropy_sum']
                                current = (result[1]+physical@anchor)/entropy
                        else:
                            mu = min(10., mu*1.15)
                        record({'case': registered['label'], 'epoch': epoch, 'iteration': iteration,
                                'status': 'master_step', 'predicted': float(predicted), 'actual': float(actual),
                                'accepted': bool(take), 'mu': mu, 'master_success': bool(opt.success)})
                    raw = anchor
                    details = query(raw)[3]
                    fixed_profiles = FixedProfileSupport(details['scaled_orders'],
                        details['scaled_energy_x'], details['scaled_energy_y'])
                    angle_fit = minimize_scalar(lambda v: fixed_profiles.evaluate((0., v), False),
                        bounds=(-2., 2.), method='bounded', options={'xatol': 1e-9})
                    choices = [(float(angle_fit.fun), float(angle_fit.x)),
                               (fixed_profiles.evaluate((0., -2.), False), -2.),
                               (fixed_profiles.evaluate((0., 2.), False), 2.)]
                    value, logtan = min(choices)
                    if value < best['rate']:
                        best = {'rate': value, 'parameters': scaled_from_raw(raw, logtan).tolist()}
                        state['best'] = best
                    state['epochs_completed'] = epoch+1
                    state['seconds'] = time.monotonic()-start
                    save()
            except TimeoutError:
                state['status'] = 'complete_timed_budget'
            else:
                state['status'] = 'epoch_ceiling_reached'
            state['seconds'] = time.monotonic()-start
            save()
    report['status'] = 'complete_registered_attempts'
    report['trace_sha256'] = sha256(args.output.with_suffix('.jsonl').read_bytes()).hexdigest()
    save()
    print(json.dumps({'status': report['status'], 'seconds': report['seconds']}), flush=True)


if __name__ == '__main__':
    main()
