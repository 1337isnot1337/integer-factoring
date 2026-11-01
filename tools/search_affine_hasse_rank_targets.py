#!/usr/bin/env python3
"""Exact higher-jet rank targets in a full polynomial affine frame span.

Complementary row/column rectangles vanish through prescribed Hasse orders
at F4 points2 and3. Optional original-coefficient constraints protect low
natural weighted jets at zero and infinity. Complete resulting tensors,
endpoints, all imposed jets, and numerical diagnostics are retained.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import gzip
from itertools import combinations
import json
import math
from pathlib import Path
import time

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, degree, scale
from screen_affine_finite_places import bitplane_rank, evaluate, rows
from search_physical_code_gauges import start_params
from search_reflection4_coupled_frames import save
from solve_dense_affine_scalar_content import Span


FULL = 6400
JET_COUNT = 16
JET_WIDTH = JET_COUNT * 256


def coefficient_planes(value):
    output = [(0, 0) for _ in range(25)]
    for entry in range(256):
        shift = entry * 25
        a, b = (value[0] >> shift) & ((1 << 25) - 1), (value[1] >> shift) & ((1 << 25) - 1)
        union = a | b
        while union:
            bit = union & -union
            union ^= bit
            exponent = bit.bit_length() - 1
            first, second = output[exponent]
            output[exponent] = first | (((a >> exponent) & 1) << entry), second | (((b >> exponent) & 1) << entry)
    return output


def expanded(value):
    coefficients = coefficient_planes(value)
    output = list(value)
    for point in (2, 3):
        powers = [1]
        for _ in range(24):
            powers.append(F4_PRODUCT[powers[-1]][point])
        offset = FULL + (point - 2) * JET_WIDTH
        for jet in range(JET_COUNT):
            answer = (0, 0)
            for exponent in range(jet, 25):
                if (jet & exponent) == jet:
                    answer = add(answer, scale(coefficients[exponent], powers[exponent - jet]))
            output[0] |= answer[0] << (offset + 256 * jet)
            output[1] |= answer[1] << (offset + 256 * jet)
    return tuple(output)


def polynomial_values(value):
    mask = (1 << 25) - 1
    return [((value[0] >> (25 * i)) & mask, (value[1] >> (25 * i)) & mask) for i in range(256)]


def packed(values):
    return tuple(sum(v[part] << (25 * i) for i, v in enumerate(values)) for part in (0, 1))


def add_column(system, value, index):
    reduced, witness = system.reduce(value, (1 << index, 0))
    if not any(reduced):
        return
    pivot = ((reduced[0] | reduced[1]) & -(reduced[0] | reduced[1])).bit_length() - 1
    c = ((reduced[0] >> pivot) & 1) + 2 * ((reduced[1] >> pivot) & 1)
    inv = (0, 1, 3, 2)[c]
    system.pivots[pivot] = scale(reduced, inv), scale(witness, inv)


def rectangle(point, cut, first, second, order):
    selected_rows, selected_columns = {5 * x for x in first}, {5 * x for x in second}
    entry_mask = 0
    for index in range(256):
        x, y, u, v = ((index // 4 ** j) % 4 for j in range(4))
        row, col = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
        if row not in selected_rows and col not in selected_columns:
            entry_mask |= 1 << index
    return sum(entry_mask << (FULL + (point - 2) * JET_WIDTH + 256 * jet) for jet in range(order))


def protection(zero_cutoff, infinity_cutoff):
    h = (0, 6, 2, 8)
    mask = 0
    for index in range(256):
        x, y, u, v = ((index // 4 ** j) % 4 for j in range(4))
        for exponent in range(25):
            if exponent + h[x] + h[y] < zero_cutoff or 24 - exponent + h[u] + h[v] < infinity_cutoff:
                mask |= 1 << (25 * index + exponent)
    return mask


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--span', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=2400)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.span.read_bytes()
    source = json.loads(gzip.decompress(raw))
    assert source['ambient_field_coordinates'] == FULL
    basis = [tuple(int(v, 16) for v in row['bitplanes']) for row in source['basis']]
    base = tensor(4)
    base_value = packed(base['values'])
    assert polynomial_values(base_value) == base['values']
    transformed = [expanded(v) for v in basis]
    target = expanded(base_value)
    # Independent Hasse evaluations check the packed transform before solving.
    from math import comb
    controls = 0
    for value in [base_value] + basis[:5]:
        out = expanded(value)
        values = polynomial_values(value)
        for index in (0, 1, 7, 63, 129, 255):
            for point in (2, 3):
                for jet in range(JET_COUNT):
                    total = 0
                    a, b = values[index]
                    for exponent in range(jet, 25):
                        if comb(exponent, jet) % 2:
                            coefficient = ((a >> exponent) & 1) + 2 * ((b >> exponent) & 1)
                            power = 1
                            for _ in range(exponent - jet):
                                power = F4_PRODUCT[power][point]
                            total ^= F4_PRODUCT[coefficient][power]
                    bit = FULL + (point - 2) * JET_WIDTH + 256 * jet + index
                    actual = ((out[0] >> bit) & 1) + 2 * ((out[1] >> bit) & 1)
                    assert actual == total
                    controls += 1
    supports = [(first, tuple(i for i in range(4) if i not in first)) for first in combinations(range(4), 2)]
    tasks = []
    for cutoffs in ((12, 20), (8, 16), (0, 0), (16, 24)):
        for order in (4, 2, 6, 8, 12):
            for first, second in supports:
                for reverse in (False, True):
                    other = (second, first) if reverse else (first, second)
                    tasks.append((cutoffs, [(2, 0, first, second, order), (3, 1, *other, order)]))
        for order in (4, 8, 12, 16):
            for first, second in supports:
                for point, cut in ((2, 0), (2, 1), (3, 0), (3, 1)):
                    tasks.append((cutoffs, [(point, cut, first, second, order)]))
    names = (Path(__file__).name, 'asymmetric_crossing_pressure.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'solve_dense_affine_scalar_content.py')
    report = {'schema_version': 1, 'scope': __doc__, 'span_sha256': sha256(raw).hexdigest(),
              'span_dimension': len(basis), 'span_packing': source['packing'],
              'requested_seconds': args.seconds, 'planned_targets': len(tasks),
              'independent_hasse_controls': controls, 'cases': [], 'retained_seeds': [],
              'feasible_targets': 0, 'unique_seeds': 0, 'unresolved_pressure_calls': 0,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}}
    began, seen = time.monotonic(), set()
    natural = start_params(np.asarray([0, 12, 4, 16]) * .2326459924883792, .2326459924883792)
    for index, (cutoffs, conditions) in enumerate(tasks):
        if time.monotonic() - began >= args.seconds:
            break
        fixed = protection(*cutoffs)
        finite = 0
        for condition in conditions:
            finite |= rectangle(*condition)
        mask = fixed | finite
        system = object.__new__(Span)
        system.pivots, system.increments = {}, []
        for j, column in enumerate(transformed):
            add_column(system, (column[0] & mask, column[1] & mask), j)
        rhs = target[0] & finite, target[1] & finite
        residual, word = system.reduce(rhs)
        record = {'target_index': index, 'protected_cutoffs': cutoffs, 'conditions': conditions,
                  'constraint_count': mask.bit_count(), 'constraint_rank': len(system.pivots),
                  'feasible': not any(residual), 'residual_bitplanes': [hex(v) for v in residual]}
        if not any(residual):
            report['feasible_targets'] += 1
            value = base_value
            for j, vector in enumerate(basis):
                c = ((word[0] >> j) & 1) + 2 * ((word[1] >> j) & 1)
                if c:
                    value = add(value, scale(vector, c))
            observed = expanded(value)
            assert (observed[0] & finite, observed[1] & finite) == (0, 0)
            assert (observed[0] & fixed, observed[1] & fixed) == (target[0] & fixed, target[1] & fixed)
            values = polynomial_values(value)
            assert [evaluate(v, 1) for v in values] == [evaluate(v, 1) for v in base['values']]
            digest = sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()
            record.update(seed_sha256=digest, coefficient_bitplanes=[hex(v) for v in word])
            if digest not in seen:
                seen.add(digest)
                actual_degree = max(map(degree, values))
                changed = dict(base, values=values, degree=actual_degree)
                finite_ranks = [[bitplane_rank(rows(values, cut, point)) for cut in (0, 1)] for point in (2, 3)]
                oracle = AsymmetricPressure(changed, 1)
                scores = []
                for multiplier in (1, .5, .2):
                    try:
                        scores.append(oracle.evaluate(np.r_[natural * multiplier / math.sqrt(2), math.pi / 4])[0] / math.sqrt(2))
                    except (ArithmeticError, ValueError):
                        report['unresolved_pressure_calls'] += 1
                score = min(scores) if scores else None
                record.update(degree=actual_degree, finite_ranks=finite_ranks, natural_rho=score)
                report['retained_seeds'].append(dict(changed, seed_sha256=digest, target_index=index,
                                                     finite_ranks=finite_ranks, natural_rho=score,
                                                     conditions=conditions, protected_cutoffs=cutoffs,
                                                     coefficient_bitplanes=record['coefficient_bitplanes']))
                all_seeds = report['retained_seeds']
                good_cost = sorted(all_seeds, key=lambda r: r['natural_rho'] or 1e6)[:32]
                high_order = sorted(all_seeds, key=lambda r: (-sum(c[-1] for c in r['conditions']), r['natural_rho'] or 1e6))[:32]
                report['retained_seeds'] = list({r['seed_sha256']: r for r in good_cost + high_order}.values())
        report['cases'].append(record)
        report.update(completed=False, elapsed_seconds=time.monotonic() - began, completed_targets=len(report['cases']), unique_seeds=len(seen))
        save(args.output, report)
        print(json.dumps({k: report[k] for k in ('completed_targets', 'feasible_targets', 'unique_seeds', 'elapsed_seconds')}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic() - began, completed_targets=len(report['cases']),
                  unique_seeds=len(seen), target_queue_exhausted=len(report['cases']) == len(tasks))
    save(args.output, report)


if __name__ == '__main__':
    main()
