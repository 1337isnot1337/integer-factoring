#!/usr/bin/env python3
"""Bounded fixed-scale q13 inner relaxation; outputs profiles, never tensors."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from heapq import heappop, heappush
import json
from math import gcd, log, sqrt
from pathlib import Path
import signal
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize

from check_bare_seven_to_twelve_letters import constraints as boxes
from check_cell_minor_caps import all_caps
from hereditary_adjugate_caps import constraints as hereditary
from higher_jet_word_caps import word_cost_profile
from multichannel_cofactor_caps import linear_caps as multichannel
from rectangular_adjugate_jet_caps import linear_caps as rectangular

ROOT = Path(__file__).resolve().parents[1]


def word_majorants(q, source):
    """Affine prefixes from actual distinct ordered words and repeated cap slots.

    Selecting any k candidates bounds the minimum k-candidate sum above.
    Distinct ordered words can have the same count vector and must remain
    separate candidates. The serial field retains that multiplicity.
    """
    n = q*q
    source = tuple(map(F, source))
    assert len(source) == n and source == tuple(sorted(source)) and source[0] == 0
    cap, serial = q*source[1], 0
    zero = (0,)*n
    queue, selected = [(F(0), 0, serial, zero)], []
    while queue and len(selected) < n:
        cost, length, _, counts = heappop(queue)
        selected.append(counts)
        if length < q-1:
            for i in range(1, n):
                if cost+source[i] < cap:
                    new = list(counts)
                    new[i] += 1
                    serial += 1
                    heappush(queue, (cost+source[i], length+1, serial, tuple(new)))
    cap_row = (0, q)+(0,)*(n-2)
    selected += [cap_row]*(n-len(selected))
    prefixes, current = [zero], [0]*n
    for row in selected:
        current = [a+b for a, b in zip(current, row)]
        prefixes.append(tuple(current))
    actual = word_cost_profile(q, source)
    assert [sum(a*b for a, b in zip(row, source)) for row in selected] == list(actual)
    return prefixes


def majorant_controls():
    cases, comparisons = 0, 0
    for q in (2, 3, 4, 5, 7):
        n = q*q
        source = (0,)+tuple(1+(i*i+3*i)//7 for i in range(1, n))
        affine = word_majorants(q, source)
        for seed in range(12):
            values = (0,)+tuple(sorted((17*i*i+11*seed*i+seed*seed) % 43 for i in range(1, n)))
            actual = word_cost_profile(q, values)
            for k in range(n+1):
                assert sum(a*b for a, b in zip(affine[k], values)) >= sum(actual[:k])
                comparisons += 1
            cases += 1
    return dict(scope='Abstract profiles, not crossing tensors.', cases=cases, prefix_comparisons=comparisons)


def compile_caps(q, source):
    n, rows, counts = q*q, {}, {}

    def add(row, rhs, family):
        rhs = F(rhs)
        row = [int(a*rhs.denominator) for a in row]
        rhs = rhs.numerator
        common = gcd(rhs, gcd(*row))
        if common:
            row, rhs = [a//common for a in row], rhs//common
        key = tuple(row)
        if not any(key):
            assert rhs >= 0
            return
        # The first order is already fixed to zero by a bound.
        if not rhs and not any(row[1:]):
            return
        rows[key] = min(rhs, rows.get(key, rhs))
        counts[family] = counts.get(family, 0)+1

    for _, first, second, size in boxes(q):
        add([a+b for a, b in zip(first, second)], size, 'box')
    for k, (rhs, _) in all_caps(1, [0]*q, [0]*q).items():
        add([int(i < k) for i in range(n)], rhs, 'cell_shared')
    for family, builder in (('rectangular_jet', rectangular), ('multichannel', multichannel)):
        for cap in builder(1, [0]*q, [0]*q):
            add([a+b for a, b in zip(cap['first'], cap['second'])], cap['rhs'], family)
    word = {m:word_majorants(m, source[:m*m]) for m in range(2, q+1)}
    for cap in hereditary(1, [0]*q, [0]*q):
        first = [sum(c for k, c in cap['first_prefixes'] if i < k) for i in range(n)]
        square = [a+sum(c for k, c in cap['second_prefixes'] if i < k) for i, a in enumerate(first)]
        add(square, cap['rhs'], 'square')
        upper = word[cap['physical_size']][cap['size']]
        add([a+(upper[i] if i < len(upper) else 0) for i, a in enumerate(first)], cap['rhs'], 'word_majorant')
    A, rhs = [], []
    for row, bound in rows.items():
        factor = max(1, abs(bound))
        A.append([float(a)/factor for a in row])
        rhs.append(float(bound)/factor)
    for i in range(n-1):
        row = [0.]*n
        row[i], row[i+1] = 1., -1.
        A.append(row)
        rhs.append(0.)
    return np.array(A), np.array(rhs), dict(compiled_families=counts, unique_rows=len(rows), sorted_rows=n-1)


def run():
    controls = majorant_controls()
    target = json.loads((ROOT/'evidence/bare-seven-to-twelve-letters.json').read_text())['thirteen_target']
    q = 13
    source = [F(x, target['M']) for x in target['mu']]
    t = target['M']*log(float(F(target['s'])))
    compiled_start = time.monotonic()
    A, rhs, compiler = compile_caps(q, source)
    compile_seconds = time.monotonic()-compiled_start
    source = np.array(list(map(float, source)))
    products = A@source
    ratios = rhs[products > 0]/products[products > 0]
    shrink = min(1., min(ratios))*(1-1e-10)
    initial = source*shrink
    assert np.max(A@initial-rhs) < 1e-8
    trace = []
    best = initial.copy()
    fit_start = time.monotonic()

    def value(x):
        weights = np.exp(-t*x)
        return log(sum(weights)), -t*weights/sum(weights)

    def callback(x):
        nonlocal best
        excess = float(np.max(A@x-rhs))
        rho = (t+2*value(x)[0])/(4*log(q))
        trace.append(dict(iteration=len(trace)+1, seconds=time.monotonic()-fit_start,
                          formal_rho=rho, maximum_normalized_excess=excess))
        if excess <= 1e-8 and value(x)[0] < value(best)[0]:
            best = x.copy()

    def timeout(signum, frame):
        raise TimeoutError('registered 120-second fit limit')

    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(120)
    try:
        result = minimize(value, initial, jac=True, method='SLSQP',
                          bounds=Bounds([0.]*(q*q), [0.]+[float(q*q)]*(q*q-1)),
                          constraints=LinearConstraint(A, -np.inf, rhs), callback=callback,
                          options=dict(maxiter=150, ftol=1e-11))
        callback(result.x)
        termination = dict(success=bool(result.success), message=str(result.message),
                           iterations=int(result.nit), function_queries=int(result.nfev))
    except TimeoutError as error:
        termination = dict(success=False, message=str(error), iterations=len(trace))
    finally:
        signal.alarm(0)
    rho = (t+2*value(best)[0])/(4*log(q))
    return dict(scope='One numerical inner-relaxation fit; profiles are not tensors and failures are not universal obstructions.',
                alphabet=q, normalized_degree=1, fixed_thermal_scale=t,
                compile_seconds=compile_seconds, fit_seconds=time.monotonic()-fit_start,
                compiler=compiler, majorant_controls=controls, initial_uniform_shrink=shrink,
                initial_formal_rho=(t+2*value(initial)[0])/(4*log(q)),
                termination=termination, formal_rho=rho, threshold_margin=1/sqrt(2)-rho,
                maximum_normalized_excess=float(np.max(A@best-rhs)), orders=best.tolist(), trace=trace,
                continuation_gate=bool(rho < 1/sqrt(2)-1e-5),
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('orders','trace')}))


if __name__ == '__main__':
    main()
