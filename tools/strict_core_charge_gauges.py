#!/usr/bin/env python3
"""Arbitrary physical core gauges for the bounded strict charge projection."""
from __future__ import annotations

import argparse
import ctypes
from functools import lru_cache
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp

from check_rectangular_defect_support import Rat, Z, car_checks
from check_rectangular_defect_core_basis import (
    residue_operator, left_car, right_car, independent_rows,
)
from native_matrix_smith_f4 import MatrixOracle
from rectangular_defect_profiles import (
    Program, Expr, library, precision_sizes, retained_known_below, MASK,
    lcm2, quotient2,
)
from strict_core_charge_profiles import core_data, projection_program


@lru_cache(None)
def projection_data(place, interval, core, family='strict_low'):
    if family == 'strict_low':
        data = dict(core_data(place, interval, core))
        data['family'] = family
        data['physical_core_indices'] = list(range(core))
        return data
    if family != 'raw_extremes' or core != 2 or interval < 2:
        raise ValueError('The raw extreme family requires two core modes in an interval of length at least two.')
    from check_rectangular_full_core_support import full_core_data
    data = dict(full_core_data(place, interval))
    data['family'], data['core'] = family, core
    data['physical_core_indices'] = [0, interval - 1]
    images = {}
    for name, value in data['active_images'].items():
        prefix, mode, suffix = name.split('_')
        side = {'out': 'left', 'in': 'right', 'left': 'left', 'right': 'right'}[prefix]
        letter, original_index = mode[0], int(mode[1:])
        local_index = 0 if original_index == 0 else 1
        assert original_index in (0, interval - 1)
        images[f'{side}_{letter}{local_index}_{suffix}'] = value
    data['active_images'] = images
    left = ('v', 'y') if place == 'infinity' else ('y', 'x')
    right = ('u', 'x') if place == 'infinity' else ('v', 'u')
    data['core_names'] = {side: [f'{side}_{letter}{i}' for letter in letters for i in range(2)]
                          for side, letters in (('left', left), ('right', right))}
    active = set(data['active_indices'])
    data['spectator_charges'] = [1 if place == 'infinity' or i < interval else -1
                                for i in range(2 * interval) if i not in active]
    return data


@lru_cache(None)
def coordinate_basis(place, interval, core, side, family='strict_low'):
    """Certify regular core coordinates by a unit vacuum minor modulo z."""
    data = projection_data(place, interval, core, family)
    names = data['core_names'][side]
    images = {}
    for name in names:
        for suffix in ('create', 'annihilate'):
            image = data['active_images'][name + '_' + suffix]
            creation = [Rat(*entry) for entry in image['creation']]
            annihilation = [Rat(*entry) for entry in image['annihilation']]
            if side == 'right':
                creation = [value * Z ** (-order) for value, order in zip(creation, data['active_orders'])]
                annihilation = [value * Z ** order for value, order in zip(annihilation, data['active_orders'])]
            assert all(value.order() >= 0 for value in creation + annihilation if value)
            images[name + '_' + suffix] = (creation, annihilation)
    car_checks(images, names)
    active, rank = len(data['active_indices']), len(names)
    size, multiplicity = 1 << active, 1 << (active - rank)
    operators = {name: residue_operator(pair) for name, pair in images.items()}
    projector = [1 << i for i in range(size)]
    for name in names:
        occupied = left_car(left_car(projector, operators[name + '_annihilate']),
                            operators[name + '_create'])
        projector = [first ^ second for first, second in zip(projector, occupied)]
    rows, columns = independent_rows(projector)
    assert len(rows) == len(columns) == multiplicity
    minor = [sum(((projector[i] >> j) & 1) << k for k, j in enumerate(columns)) for i in rows]
    assert len(independent_rows(minor)[0]) == multiplicity
    left, right = [], [0] * size
    for mask in range(1 << rank):
        row_map, column_map = projector, projector
        for i, name in enumerate(names):
            if mask >> i & 1:
                row_map = right_car(row_map, operators[name + '_annihilate'])
                column_map = left_car(column_map, operators[name + '_create'])
        left.extend(row_map[i] for i in rows)
        for i, value in enumerate(column_map):
            right[i] |= sum(((value >> j) & 1) << (mask * multiplicity + k)
                            for k, j in enumerate(columns))
    assert len(independent_rows(left)[0]) == len(independent_rows(right)[0]) == size
    mode_charges = [1 if place == 'infinity' or i < interval else -1
                    for i in data['active_indices']]

    def charge(mask):
        return sum(value for i, value in enumerate(mode_charges) if mask >> i & 1)

    core_charges = []
    for name in names:
        creation, annihilation = images[name + '_create']
        values = {mode_charges[i] for i, value in enumerate(creation) if value}
        values |= {-mode_charges[i] for i, value in enumerate(annihilation) if value}
        assert len(values) == 1
        core_charges.append(values.pop())
    labels = {letter: [(mask >> (group * core)) & ((1 << core) - 1)
                       for mask in range(1 << rank) for _ in range(multiplicity)]
              for group, letter in enumerate([names[0].split('_')[1][0],
                                               names[core].split('_')[1][0]])}
    charges = lambda vacuum: [charge(index) + sum(value for i, value in enumerate(core_charges) if mask >> i & 1)
                               for mask in range(1 << rank) for index in vacuum]
    return {'names': names, 'regular_images': {
                name: {'creation': [value.pair() for value in pair[0]],
                       'annihilation': [value.pair() for value in pair[1]]}
                for name, pair in images.items()},
            'multiplicity': multiplicity, 'vacuum_rows': rows, 'vacuum_columns': columns,
            'row_charges': charges(rows), 'column_charges': charges(columns),
            'labels': labels, 'unit_basis_rank': size}


def coordinate_matrix(data, charges, profile_precision=None):
    left, right = [coordinate_basis(data['place'], data['M'], data['core'], side, data.get('family', 'strict_low'))
                   for side in ('left', 'right')]
    extended = dict(data)
    extended['active_images'] = dict(data['active_images'])
    names = []
    for side, basis in (('left', left), ('right', right)):
        prefix = 'basis_' + side + '_'
        names.append([prefix + name for name in basis['names']])
        extended['active_images'].update({prefix + name: value for name, value in basis['regular_images'].items()})
    program = Program(extended)
    size = 1 << program.k
    if program.k > 10:
        raise ValueError('Dense core-coordinate evaluation currently supports at most ten active modes.')
    diagonal = [sum(order for j, order in enumerate(data['active_orders']) if mask >> j & 1)
                for mask in range(size)]

    def vacuum(value, pairs):
        for name in pairs:
            value = program.add(value, program.number(value, name))
        return value

    left_names, right_names = names
    right_vectors = {0: vacuum(Expr(0), right_names)}
    outputs = {}
    for state in range(1 << len(right_names)):
        if state:
            bit = state & -state
            right_vectors[state] = program.car(right_vectors[state ^ bit],
                                               right_names[bit.bit_length() - 1] + '_create')
        initial = program.diagonal(right_vectors[state], diagonal)
        _, value = projection_program(data, charges, program, initial)
        annihilated = {0: value}
        for mask in range(1 << len(left_names)):
            if mask:
                bit = mask & -mask
                annihilated[mask] = program.car(annihilated[mask ^ bit],
                                                left_names[bit.bit_length() - 1] + '_annihilate')
            outputs[mask, state] = vacuum(annihilated[mask], left_names)
    base = min(value.base for value in outputs.values())
    denominator = 1
    for value in outputs.values():
        denominator = lcm2(denominator, value.denominator)
    output_ids = []
    for mask in range(1 << len(left_names)):
        for state in range(1 << len(right_names)):
            value = outputs[mask, state]
            factor = quotient2(denominator, value.denominator) << (value.base - base)
            output_ids.append(program.polynomial(value, factor).register)
    window, builder = precision_sizes(program, base, data['M'], profile_precision)
    words, output_words = builder // 64, window // 64
    instructions = (ctypes.c_int * (5 * len(program.instructions)))(
        *(value for row in program.instructions for value in row))
    packed = (ctypes.c_uint64 * (2 * words * len(program.coefficients)))(
        *(value >> (64 * w) & MASK for pair in program.coefficients for value in pair for w in range(words)))
    terms, offsets = [], [0]
    for car in program.cars:
        terms.extend(value for term in car for value in term)
        offsets.append(len(terms) // 3)
    cc = (ctypes.c_int * len(terms))(*terms)
    co = (ctypes.c_int * len(offsets))(*offsets)
    dd = (ctypes.c_int * size)(*(value - min(diagonal) for value in diagonal))
    ii = (ctypes.c_int * len(left['vacuum_rows']))(*left['vacuum_rows'])
    jj = (ctypes.c_int * len(right['vacuum_columns']))(*right['vacuum_columns'])
    rr = (ctypes.c_int * len(output_ids))(*output_ids)
    output = (ctypes.c_uint64 * (2 * output_words * size * size))()
    shifts = (ctypes.c_int64 * (size * size))()
    meta = (ctypes.c_int64 * 4)()
    lib = library(builder)
    integer = ctypes.POINTER(ctypes.c_int)
    unsigned = ctypes.POINTER(ctypes.c_uint64)
    signed = ctypes.POINTER(ctypes.c_int64)
    lib.fock_f4_coordinates_v2.argtypes = [ctypes.c_int, ctypes.c_int, integer, ctypes.c_int, unsigned,
        ctypes.c_int, integer, integer, integer, ctypes.c_int, integer, ctypes.c_int, integer,
        integer, ctypes.c_int, ctypes.c_int, unsigned, signed, signed]
    lib.fock_f4_coordinates_v2.restype = ctypes.c_int
    started = time.monotonic()
    status = lib.fock_f4_coordinates_v2(program.k, len(program.instructions), instructions,
        len(program.coefficients), packed, len(program.cars), co, cc, dd,
        left['multiplicity'], ii, right['multiplicity'], jj, rr, base, output_words, output, shifts, meta)
    if status:
        raise ArithmeticError(f'Core-coordinate construction failed: {status}, metadata={list(meta)}')
    matrix, shift_matrix = [], []
    for i in range(size):
        row, shift_row = [], []
        for j in range(size):
            index = i * size + j
            pair = tuple(sum(int(output[(2 * index + plane) * output_words + w]) << (64 * w)
                             for w in range(output_words)) for plane in range(2))
            assert not pair[1]
            row.append(pair)
            shift_row.append(int(shifts[index]))
        matrix.append(row)
        shift_matrix.append(shift_row)
    degree_bound = max(program.degree_bounds[index] for index in output_ids)
    metadata = {'matrix_order': size, 'registers': len(program.instructions),
                'profile_precision': window, 'builder_precision': builder,
                'unit_denominator_degree': denominator.bit_length() - 1,
                'program_degree_bound': degree_bound,
                'exact_polynomial_matrix': degree_bound < builder and int(meta[3]) < window,
                'known_below': retained_known_below(meta, window),
                'matrix_seconds': time.monotonic() - started}
    return matrix, shift_matrix, metadata, left, right


@lru_cache(None)
def gaussian_coordinates(place, interval, core, family):
    """The coordinate bases intertwine every physical core number operator.

    Therefore every diagonal charge projection is an exact entry mask of
    this one Gaussian coordinate matrix. No phase program is rebuilt.
    """
    data = projection_data(place, interval, core, family)
    return coordinate_matrix(data, (0,) * (1 << core))


class CoreChargePressure:
    def __init__(self, interval, charges, core=2, family='strict_low'):
        self.interval, self.core, self.q = interval, core, 1 << core
        self.charges = tuple(charges)
        if len(self.charges) != self.q or any(value not in (0, 1, 2) for value in self.charges):
            raise ValueError('Supply an F3-valued charge on every core state.')
        self.family = family
        self.records, self.oracles = [], []
        self.queries = 0
        for place in ('zero', 'infinity'):
            data = projection_data(place, interval, core, family)
            gaussian, shifts, metadata, left, right = gaussian_coordinates(place, interval, core, family)
            signed = lambda basis: [sum((1 if letter in ('x', 'u') else -1) * self.charges[labels[i]]
                                        for letter, labels in basis['labels'].items()) % 3
                                    for i in range(len(gaussian))]
            row_charges, column_charges = signed(left), signed(right)
            matrix = [[value if (row_charges[i] + column_charges[j]) % 3 == 0 else (0, 0)
                       for j, value in enumerate(row)] for i, row in enumerate(gaussian)]
            self.records.append((data, metadata, left, right))
            self.oracles.append(MatrixOracle(matrix, len(matrix), shifts,
                                            precision=metadata['profile_precision']))
        self.core_indices = self.records[0][0]['physical_core_indices']
        self.physical_spectators = [i for i in range(interval) if i not in self.core_indices]

    def evaluate(self, log_s, offsets=None, bulk=(0.0, 0.0), gradient=False):
        offsets = np.zeros((2, 4, self.q)) if offsets is None else np.asarray(offsets, dtype=float)
        if offsets.shape != (2, 4, self.q):
            raise ValueError('Core offsets have shape (two places, four legs, core alphabet).')
        letters = {'x': 0, 'y': 1, 'u': 2, 'v': 3}
        log_b = 0.0
        offset_derivative = np.zeros_like(offsets)
        bulk_derivative = np.zeros(2)
        scale_derivative = 0.0
        profiles = []
        for place, (record, oracle) in enumerate(zip(self.records, self.oracles)):
            data, metadata, left, right = record
            rows = bulk[place] * np.asarray(left['row_charges'], dtype=float)
            columns = np.zeros(len(right['column_charges']))
            for letter, labels in left['labels'].items():
                rows += offsets[place, letters[letter]][labels]
            for letter, labels in right['labels'].items():
                columns += offsets[place, letters[letter]][labels]
            result = oracle.profile(log_s * rows, log_s * columns, scale=log_s)
            orders = np.asarray(result['orders'])
            if not metadata['exact_polynomial_matrix']:
                threshold = log_s * (metadata['known_below'] + float(rows.min()) + float(columns.min()))
                if float(orders.max()) >= threshold:
                    raise ArithmeticError('Weighted profile exceeds its certified coefficient window.')
            log_z = float(logsumexp(-orders))
            probabilities = np.exp(-orders - log_z)
            log_b += log_z
            scale_derivative -= float(probabilities @ orders) / log_s
            row_indices, column_indices = np.asarray(result['rows']), np.asarray(result['cols'])
            bulk_derivative[place] -= log_s * float(probabilities @ np.asarray(left['row_charges'])[row_indices])
            for basis, indices in ((left, row_indices), (right, column_indices)):
                for letter, labels in basis['labels'].items():
                    values = np.asarray(labels)[indices]
                    offset_derivative[place, letters[letter]] -= log_s * np.bincount(
                        values, weights=probabilities, minlength=self.q)
            charges = np.asarray(data['spectator_charges'], dtype=float)
            spectator = np.asarray(data['spectator_orders'], dtype=float) + bulk[place] * charges
            log_b += float(np.logaddexp(0.0, -log_s * spectator).sum())
            occupation = np.exp(-np.logaddexp(0.0, log_s * spectator))
            scale_derivative -= float(occupation @ spectator)
            bulk_derivative[place] -= log_s * float(occupation @ charges)
            profiles.append(result['orders'])
        natural = np.asarray([4 * i - 2 * self.interval + 2 for i in range(self.interval)], dtype=float)
        populations = np.asarray([state.bit_count() for state in range(self.q)], dtype=float)
        core_energy = np.asarray([sum(natural[index] for i, index in enumerate(self.core_indices) if state >> i & 1)
                                  for state in range(self.q)], dtype=float)
        gx = core_energy + (bulk[0] + bulk[1]) * populations + offsets[:, (0, 3), :].sum(axis=(0, 1))
        gy = core_energy + (-bulk[0] + bulk[1]) * populations + offsets[:, (1, 2), :].sum(axis=(0, 1))
        spectator_physical = np.stack((natural[self.physical_spectators] + bulk[0] + bulk[1],
                                       natural[self.physical_spectators] - bulk[0] + bulk[1]))

        def equation(rate):
            physical = float(logsumexp(-log_s * gx / rate) + logsumexp(-log_s * gy / rate)
                             + np.logaddexp(0.0, -log_s * spectator_physical / rate).sum())
            return log_b - rate * physical

        if equation(1e-9) <= 0:
            raise ArithmeticError('A nonpositive pressure root requires separate review.')
        upper = 4.0
        while equation(upper) > 0:
            upper *= 2
        rate = brentq(equation, 1e-9, upper, xtol=1e-13)
        px = np.exp(-log_s * gx / rate - logsumexp(-log_s * gx / rate))
        py = np.exp(-log_s * gy / rate - logsumexp(-log_s * gy / rate))
        ps = np.exp(-np.logaddexp(0.0, log_s * spectator_physical / rate))
        entropy = -sum(float(p[p > 0] @ np.log(p[p > 0])) for p in (px, py))
        entropy += float(np.logaddexp(0.0, -log_s * spectator_physical / rate).sum()
                         + (ps * log_s * spectator_physical / rate).sum())
        self.queries += 1
        detail = {'rho': rate / 2, 'log_s': log_s, 'bulk': list(bulk), 'entropy': entropy,
                  'type_x_core': px.tolist(), 'type_y_core': py.tolist(),
                  'active_orders': profiles, 'log_profile_product': log_b}
        if not gradient:
            return rate / 2, detail
        if entropy < 1e-10:
            raise ArithmeticError('Physical entropy is too small for a stable derivative.')
        for place in range(2):
            offset_derivative[place] += log_s * np.stack((px, py, py, px))
        number_x, number_y = float(px @ populations + ps[0].sum()), float(py @ populations + ps[1].sum())
        bulk_derivative += log_s * np.array([number_x - number_y, number_x + number_y])
        scale_derivative += float(px @ gx + py @ gy + (ps * spectator_physical).sum())
        return rate / 2, (offset_derivative / (2 * entropy),
                          bulk_derivative / (2 * entropy), scale_derivative / (2 * entropy)), detail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=8)
    parser.add_argument('--family', choices=('strict_low', 'raw_extremes'), default='strict_low')
    parser.add_argument('--charges', nargs='+', type=int, default=[0, 0, 0, 1])
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--maxfun', type=int, default=250)
    parser.add_argument('--restarts', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    model = CoreChargePressure(args.interval, tuple(args.charges), family=args.family)
    dimension = 8 * (model.q - 1) + 3
    best, history = None, []
    rng = np.random.default_rng(597213)

    def save():
        args.output.write_text(json.dumps({'status': 'exploratory physical-core gauge search',
            'interval': args.interval, 'family': args.family, 'charges': args.charges, 'queries': model.queries,
            'seconds': time.monotonic() - started, 'best': best, 'restarts': history,
            'parameter_layout': 'core offsets 2x4x(q-1), zero bulk, infinity bulk, log(s)'}, indent=2) + '\n')

    def objective(parameters):
        nonlocal best
        if time.monotonic() - started > args.seconds:
            raise TimeoutError
        offsets = np.zeros((2, 4, model.q))
        offsets[:, :, 1:] = parameters[:-3].reshape((2, 4, model.q - 1))
        try:
            value, (do, db, ds), detail = model.evaluate(parameters[-1], offsets, parameters[-3:-1], True)
        except ArithmeticError:
            return 100.0, np.zeros(dimension)
        if best is None or value < best['rho'] - 1e-10:
            best = {**detail, 'parameters': parameters.tolist()}
            save()
            print(json.dumps({'query': model.queries, 'rho': value,
                              'seconds': time.monotonic() - started}), flush=True)
        return value, np.r_[do[:, :, 1:].reshape(-1), db, ds]

    try:
        for ell in (0.1, 0.3, 0.6):
            objective(np.r_[np.zeros(dimension - 1), ell])
        for restart in range(args.restarts):
            parameters = np.asarray(best['parameters'])
            if restart:
                parameters[:-1] += rng.normal(0, 0.25, dimension - 1)
            run = minimize(objective, parameters, jac=True, method='L-BFGS-B',
                           bounds=[(-32, 32)] * (dimension - 3)
                           + [(-4 * args.interval, 4 * args.interval)] * 2 + [(0.002, 3.0)],
                           options={'maxfun': args.maxfun, 'maxiter': args.maxfun,
                                    'ftol': 1e-13, 'gtol': 1e-8, 'maxls': 20})
            history.append({'restart': restart, 'rho': float(run.fun),
                            'success': bool(run.success), 'message': str(run.message),
                            'evaluations': int(run.nfev)})
            save()
    except TimeoutError:
        history.append({'status': 'time budget reached'})
    finally:
        save()


if __name__ == '__main__':
    main()
