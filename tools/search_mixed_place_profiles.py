#!/usr/bin/env python3
"""Explore the opposite-flag criterion using complete local Smith profiles.

All reported minima are numerical search results. Exact polynomial profiles
are inputs; a claimed crossing improvement still needs a rational certificate.
The optional allocation at a place may be zero, meaning that the place is
omitted. Point sets used by the two cuts always remain disjoint.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import exp, isfinite, log
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp


class Profile:
    def __init__(self, orders):
        if not orders or min(orders) < 0:
            raise ValueError("this unweighted search requires complete nonnegative profiles")
        counts = Counter(orders)
        self.orders = np.array(sorted(counts), dtype=float)
        self.logs = np.array([log(counts[int(a)]) for a in self.orders])
        self.n = len(orders)
        self.zero = counts.get(0, 0)
        self.top = max(orders)

    def value(self, x):
        return float(logsumexp(self.logs - x * self.orders))

    def inverse(self, target):
        """x with log(sum exp(-order*x))=target, or infinity if absent."""
        if target >= log(self.n):
            return 0.0
        if self.zero and target <= log(self.zero):
            return float("inf")
        if self.top == 0:
            return float("inf")
        high = 1.0
        while self.value(high) > target:
            high *= 2
            if high > 1e12:
                return float("inf")
        return float(brentq(lambda x: self.value(x) - target, 0, high,
                            xtol=2e-11, rtol=2e-12))


def envelope(profiles, t):
    """Minimize max_j U_j(exp(t/alpha_j)) with sum alpha_j=1.

    For a proposed log upper bound b, a place permits allocation at most
    t/inverse_j(b). These independent capacities must sum to at least one.
    A place with no positive capacity can be omitted. The optimum is the
    unique equality point, except for the constant-profile boundary.
    """
    if not profiles or t <= 0:
        raise ValueError("nonempty place set and positive outer log scale required")
    n = profiles[0].n
    if any(p.n != n for p in profiles):
        raise ValueError("both cuts must use complete profiles of equal dimension")
    if len(profiles) == 1:
        return profiles[0].value(t), [1.0]
    useful = [p for p in profiles if p.top]
    if not useful:
        return log(n), [1.0 / len(profiles)] * len(profiles)

    def capacity(bound):
        total = 0.0
        for p in useful:
            x = p.inverse(bound)
            if x == 0:
                return float("inf")
            if isfinite(x):
                total += t / x
        return total

    low = log(n) - t * sum(p.top for p in useful) - 1
    high = log(n) - 1e-13
    if capacity(high) < 1:
        # A very small t can be too close to the constant boundary for the
        # inversion tolerance. Its exact limiting upper bound is n.
        return log(n), [1.0 / len(profiles)] * len(profiles)
    bound = float(brentq(lambda b: capacity(b) - 1.0, low, high,
                         xtol=3e-10, rtol=3e-11))
    allocations = []
    for p in profiles:
        x = p.inverse(bound)
        allocations.append(t / x if isfinite(x) and x > 0 else 0.0)
    total = sum(allocations)
    if not total:
        # At very large scales the optimum can round to log(min zero count),
        # where the inverse has no finite floating representation. A literal
        # equal allocation on the lowest-floor places is still a valid upper
        # bound. It is not evidence of a rigorously computed optimum.
        floor = min(p.zero for p in useful)
        active = [i for i, p in enumerate(profiles) if p.top and p.zero == floor]
        allocations = [1.0 / len(active) if i in active else 0.0
                       for i in range(len(profiles))]
        total = 1.0
    allocations = [a / total for a in allocations]
    # Report the actual bound attained by the normalized floating allocation.
    attained = max(p.value(t / a) for p, a in zip(profiles, allocations) if a > 1e-14)
    return attained, allocations


def partitions(record):
    """Enumerate only points that have a nontrivial singular profile.

    A point useful to one cut only can be given to that cut without loss,
    since zero allocation remains permitted. Only shared singular points
    require binary choices. Omitted or full-rank points cannot help here.
    """
    row_only, col_only, shared = [], [], []
    for item in record['profiles']:
        if item['point'] == 1:
            continue
        first, second = item['cuts']
        if max(first) and max(second):
            shared.append(item)
        elif max(first):
            row_only.append(item)
        elif max(second):
            col_only.append(item)
    for choices in product((0, 1), repeat=len(shared)):
        row = row_only + [p for p, side in zip(shared, choices) if side == 0]
        col = col_only + [p for p, side in zip(shared, choices) if side == 1]
        if row and col:
            yield row, col


def optimize_partition(record, row, col, grid=17):
    row_profiles = [Profile(p['cuts'][0]) for p in row]
    col_profiles = [Profile(p['cuts'][1]) for p in col]
    degree, alphabet = record['degree'], record['alphabet']

    @lru_cache(maxsize=1024)
    def evaluate(log_t):
        t = exp(log_t)
        r, alpha = envelope(row_profiles, t)
        c, beta = envelope(col_profiles, t)
        return (degree * t + r + c) / (4 * log(alphabet)), alpha, beta, r, c

    logs = np.linspace(log(1e-4), log(12.0), grid)
    values = [evaluate(float(x))[0] for x in logs]
    candidates = [(values[i], float(logs[i])) for i in range(len(logs))]
    for i in range(1, len(logs) - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            fit = minimize_scalar(lambda x: evaluate(float(x))[0],
                                  bounds=(float(logs[i - 1]), float(logs[i + 1])),
                                  method='bounded', options={'xatol': 5e-7})
            candidates.append((float(fit.fun), float(fit.x)))
    _, best = min(candidates)
    value, alpha, beta, r, c = evaluate(best)
    return {'rho': value, 'outer_log_scale': exp(best),
            'outer_scale': exp(exp(best)),
            'row_places': [p['point'] for p in row],
            'column_places': [p['point'] for p in col],
            'row_allocations': alpha, 'column_allocations': beta,
            'row_log_envelope': r, 'column_log_envelope': c,
            'evaluations': evaluate.cache_info().currsize,
            'log_scale_window': [1e-4, 12.0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profiles', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--grid', type=int, default=17)
    parser.add_argument('--limit-partitions', type=int)
    args = parser.parse_args()
    raw = args.profiles.read_bytes()
    records = json.loads(raw)
    if isinstance(records, dict):
        records = records.get('cases', records.get('seeds'))
        if records is None:
            raise ValueError('profile catalogue needs a cases or seeds list')
    result = {'scope': 'numerical opposite-flag profile search; no global optimum or new factoring exponent',
              'input_sha256': sha256(raw).hexdigest(), 'seeds': []}
    started = time.monotonic()
    for record in records:
        entry = {'seed': {k: record[k] for k in ('m', 'theta', 'alphabet', 'degree') if k in record},
                 'completed_partitions': [], 'complete_partition_enumeration': True}
        for index, (row, col) in enumerate(partitions(record)):
            if args.limit_partitions is not None and index >= args.limit_partitions:
                entry['complete_partition_enumeration'] = False
                break
            answer = optimize_partition(record, row, col, args.grid)
            entry['completed_partitions'].append(answer)
        if entry['completed_partitions']:
            entry['best'] = min(entry['completed_partitions'], key=lambda x: x['rho'])
        result['seeds'].append(entry)
        result['elapsed_seconds'] = time.monotonic() - started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(result, indent=2) + '\n')
        temporary.replace(args.output)
        print(json.dumps({'seed': entry['seed'], 'partitions': len(entry['completed_partitions']),
                          'best': entry.get('best')}), flush=True)


if __name__ == '__main__':
    main()
