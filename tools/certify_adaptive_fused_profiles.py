#!/usr/bin/env python3
"""Generate and check a rational interval certificate for adaptive jet flags.

Floating arithmetic proposes allocations. Every exported certificate must
pass the independent integer-power inequalities in verify_certificate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import ceil, exp, lcm, log
from pathlib import Path

from adaptive_mixed_place import capacity_curve


def pair(value):
    value = Fraction(value)
    return [value.numerator, value.denominator]


def rational(value):
    return Fraction(*value)


def verify_certificate(certificate):
    q = certificate['alphabet']
    degree = certificate['degree']
    rho = rational(certificate['rho'])
    if q < 2 or degree < 1 or not Fraction(1, 2) <= rho < 1:
        raise ValueError('invalid seed or contraction parameter')
    row, column = certificate['row_profiles'], certificate['column_profiles']
    if not row or not column:
        raise ValueError('empty place set')
    row_points = [json.dumps(p['point']) for p in row]
    column_points = [json.dumps(p['point']) for p in column]
    if (len(set(row_points)) != len(row_points) or
            len(set(column_points)) != len(column_points) or
            not set(row_points).isdisjoint(column_points)):
        raise ValueError('place sets must be distinct and disjoint')
    for profile in row + column:
        if len(profile['orders']) != q * q or min(profile['orders']) < 0:
            raise ValueError('complete nonnegative local profiles required')
    previous, comparisons, all_passed, failed = Fraction(0), 0, True, []
    for index, interval in enumerate(certificate['intervals']):
        left, right = rational(interval['left']), rational(interval['right'])
        if left != previous or not left < right <= degree:
            raise ValueError('intervals must form an exact ordered cover')
        previous = right
        common = rho.denominator
        for side in ('row', 'column'):
            if not interval[side]:
                raise ValueError('a positive flag prefix needs a place')
            used = set()
            for allocation in interval[side]:
                position = allocation['profile']
                profiles = row if side == 'row' else column
                if position in used or not 0 <= position < len(profiles):
                    raise ValueError('invalid or repeated allocation place')
                used.add(position)
                delta, scale = rational(allocation['amount']), rational(allocation['scale'])
                if delta < 0 or scale < 1:
                    raise ValueError('negative jet budget or scale below one')
                common = lcm(common, delta.denominator)
        row_budget = sum((rational(a['amount']) for a in interval['row']), Fraction())
        col_budget = sum((rational(a['amount']) for a in interval['column']), Fraction())
        if row_budget < right or col_budget < degree - left:
            raise ValueError('an interval is not covered by its jet allocations')
        maxima = []
        for side, profiles in (('row', row), ('column', column)):
            terms = []
            for allocation in interval[side]:
                delta, scale = rational(allocation['amount']), rational(allocation['scale'])
                orders = profiles[allocation['profile']]['orders']
                partition = sum((scale ** -order for order in orders), Fraction())
                terms.append(partition ** common * scale ** int(delta * common))
                comparisons += 1
            maxima.append(max(terms))
        target_power = 4 * rho * common
        assert target_power.denominator == 1
        passed = maxima[0] * maxima[1] < q ** int(target_power)
        if not passed:
            all_passed = False
            failed.append(index)
    if previous != degree:
        raise ValueError('the interval cover stops before the full degree')
    return {'passed': all_passed, 'intervals': len(certificate['intervals']),
            'powered_local_bounds': comparisons, 'failed_intervals': failed,
            'rho': pair(rho)}


def propose_allocation(profiles, budget, amount_denominator, scale_denominator):
    import numpy as np
    from scipy.optimize import brentq
    from scipy.special import logsumexp

    exact_budget = Fraction(budget)
    budget = float(exact_budget)
    n = len(profiles[0]['orders'])
    curves = [capacity_curve(tuple(p['orders'])) for p in profiles]
    entropies = curves[0][0]
    capacities = sum((c[1] for c in curves), np.zeros_like(entropies))
    if budget >= capacities[-1] - 1e-12:
        return [{'profile': 0, 'amount': pair(exact_budget), 'scale': [1, 1]}]
    if budget <= 1e-13:
        index = min(range(len(profiles)), key=lambda i: profiles[i]['orders'].count(0))
        return [{'profile': index, 'amount': [0, 1], 'scale': [64, 1]}]
    unique, indices = np.unique(capacities, return_index=True)
    target = float(np.interp(budget, unique, entropies[indices]))
    output = []
    for index, profile in enumerate(profiles):
        counts = Counter(profile['orders'])
        energies = np.array(sorted(counts), dtype=float)
        log_counts = np.array([log(counts[int(a)]) for a in energies])
        if target <= log(counts.get(0, 1)) + 1e-13:
            continue

        def statistics(temperature):
            terms = log_counts - temperature * energies
            partition = float(logsumexp(terms))
            mean = float(np.exp(terms - partition) @ energies)
            return partition + temperature * mean, mean

        high = 1.0
        while statistics(high)[0] > target and high < 128:
            high *= 2
        temperature = brentq(lambda z: statistics(z)[0] - target, 0, high)
        mean = statistics(temperature)[1]
        delta = Fraction(ceil(mean * amount_denominator), amount_denominator)
        scale = Fraction(exp(min(temperature, log(64)))).limit_denominator(scale_denominator)
        output.append({'profile': index, 'amount': pair(delta), 'scale': pair(max(scale, 1))})
    if not output:
        raise ArithmeticError('numerical proposal has no active place')
    total = sum((rational(a['amount']) for a in output), Fraction())
    if total < exact_budget:
        largest = max(output, key=lambda a: rational(a['amount']))
        largest['amount'] = pair(rational(largest['amount']) + exact_budget - total)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles', type=Path)
    parser.add_argument('--allocation', type=Path)
    parser.add_argument('--alphabet', type=int, default=4)
    parser.add_argument('--rho', default='31/32')
    parser.add_argument('--intervals', type=int, default=256)
    parser.add_argument('--amount-denominator', type=int, default=8192)
    parser.add_argument('--scale-denominator', type=int, default=128)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify', type=Path)
    args = parser.parse_args()
    if args.verify:
        result = verify_certificate(json.loads(args.verify.read_text()))
        print(json.dumps(result))
        if not result['passed']:
            raise SystemExit(1)
        return
    if not args.profiles or not args.allocation or not args.output:
        parser.error('generation needs --profiles, --allocation and --output')
    records = json.loads(args.profiles.read_text())
    if isinstance(records, dict):
        records = records['cases']
    record = next(p for p in records if p['alphabet'] == args.alphabet)
    allocation = next(p for p in json.loads(args.allocation.read_text())
                      if p['alphabet'] == args.alphabet)['best']
    row_points, col_points = allocation['row_places'], allocation['column_places']
    by_point = {json.dumps(p['point']): p for p in record['profiles']}
    row = [{'point': point, 'orders': by_point[json.dumps(point)]['cuts'][0]} for point in row_points]
    col = [{'point': point, 'orders': by_point[json.dumps(point)]['cuts'][1]} for point in col_points]
    degree = record['degree']
    certificate = {'scope': 'exact sufficient adaptive mixed-place certificate; no post-half claim',
                   'alphabet': args.alphabet, 'degree': degree, 'rho': pair(Fraction(args.rho)),
                   'source_profiles_sha256': sha256(args.profiles.read_bytes()).hexdigest(),
                   'seed': {k: record[k] for k in ('m', 'theta', 'alphabet', 'degree') if k in record},
                   'row_profiles': row, 'column_profiles': col, 'intervals': []}
    for i in range(args.intervals):
        left, right = Fraction(i * degree, args.intervals), Fraction((i + 1) * degree, args.intervals)
        certificate['intervals'].append({
            'left': pair(left), 'right': pair(right),
            'row': propose_allocation(row, right, args.amount_denominator, args.scale_denominator),
            'column': propose_allocation(col, degree - left, args.amount_denominator, args.scale_denominator)})
    result = verify_certificate(certificate)
    print(json.dumps(result), flush=True)
    if not result['passed']:
        args.output.with_suffix('.failed.json').write_text(json.dumps(certificate, indent=2) + '\n')
        raise SystemExit(1)
    certificate['verification'] = result
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2) + '\n')


if __name__ == '__main__':
    main()
