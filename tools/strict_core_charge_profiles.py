#!/usr/bin/env python3
"""Bounded active profiles for strict binary crossings with a fixed low core.

The interval and core lengths are powers of two. A core of k physical
modes requires at most 4*k active fermion modes at either place. This
module constructs both regular concentration maps and applies arbitrary
F3-valued core-charge projections by an exact F4 phase sum.
"""
from __future__ import annotations

import argparse
import ctypes
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import logsumexp

from check_rectangular_defect_support import (
    ZERO, ONE, Z, UNIT, T, Rat, vec, pairing, car_checks, apply_shear,
)
from native_matrix_smith_f4 import MatrixOracle
from rectangular_defect_profiles import (
    Program, Expr, library, precision_sizes, retained_known_below, MASK,
)


def _validate(interval, core):
    if interval < 1 or interval & (interval - 1):
        raise ValueError('The interval length must be a positive power of two.')
    if core < 1 or core & (core - 1) or core > interval:
        raise ValueError('The core length must be a power of two within the interval.')


def strict_images(place, interval, core):
    """Exact CAR images after the strict Gaussian's regular diagonalization."""
    _validate(interval, core)
    m, size = interval, 2 * interval
    images = {}
    if place == 'infinity':
        orders = [1 - 2 * m + 2 * a for a in range(size)]
        left_letters, right_letters = ('v', 'y'), ('u', 'x')
        for i in range(core):
            images[f'left_v{i}_create'] = (vec(size, [(m + i, ONE)]), vec(size))
            images[f'left_v{i}_annihilate'] = (
                vec(size), vec(size, [(a, UNIT * T ** (i & ~a)) for a in range(m)]
                               + [(m + i, ONE)]))
            images[f'left_y{i}_create'] = (
                vec(size, [(i, ONE)]
                    + [(m + a, UNIT * T ** (a & ~i)) for a in range(m)]), vec(size))
            images[f'left_y{i}_annihilate'] = (vec(size), vec(size, [(i, ONE)]))
            images[f'right_u{i}_create'] = (
                vec(size, [(i, Z ** orders[i])]), vec(size))
            images[f'right_u{i}_annihilate'] = (
                vec(size), vec(size, [(i, Z ** (-orders[i]))]
                    + [(m + a, UNIT / Z * T ** (-(a & i))) for a in range(m)]))
            images[f'right_x{i}_create'] = (
                vec(size, [(a, UNIT * Z ** (1 - 2 * m + 2 * (a | i)))
                           for a in range(m)] + [(m + i, Z ** (1 + 2 * i))]), vec(size))
            images[f'right_x{i}_annihilate'] = (
                vec(size), vec(size, [(m + i, Z ** (-1 - 2 * i))]))
    elif place == 'zero':
        orders = [2 * a + 1 for _ in range(2) for a in range(m)]
        left_letters, right_letters = ('x', 'y'), ('v', 'u')
        for i in range(core):
            for first, second, offset, other in (('x', 'v', 0, m), ('y', 'u', m, 0)):
                images[f'left_{first}{i}_create'] = (
                    vec(size, [(offset + i, ONE)]), vec(size))
                images[f'left_{first}{i}_annihilate'] = (
                    vec(size, [(other + a, UNIT * T ** (a & i)) for a in range(m)]),
                    vec(size, [(offset + i, ONE)]))
                images[f'right_{second}{i}_create'] = (
                    vec(size, [(offset + i, Z ** (2 * i + 1))]),
                    vec(size, [(other + a, UNIT * T ** (a & i) / Z ** (2 * a + 1))
                               for a in range(m)]))
                images[f'right_{second}{i}_annihilate'] = (
                    vec(size), vec(size, [(offset + i, Z ** (-2 * i - 1))]))
    else:
        raise ValueError('The place is zero or infinity.')
    names = {
        'left': [f'left_{letter}{i}' for letter in left_letters for i in range(core)],
        'right': [f'right_{letter}{i}' for letter in right_letters for i in range(core)],
    }
    return images, orders, names


@lru_cache(None)
def core_data(place, interval, core=2):
    """Concentrate exactly while retaining regular Gaussian companions."""
    images, orders, names = strict_images(place, interval, core)
    length = interval // core
    traces = []
    active = set()
    if length == 1:
        active.update(range(2 * interval))
    else:
        step = T ** core
        unit = sum((step ** h for h in range(length - 1)), ZERO)
        assert unit.order() == 0
        for offset in (0, interval):
            for residue in range(core):
                index = lambda j: offset + residue + core * j
                if place == 'infinity':
                    active.update((index(0), index(1)))
                    for j in range(2, length):
                        apply_shear(images, index(j), index(1), step ** (j - 1), orders, traces)
                    for j in range(2, length):
                        apply_shear(images, index(1), index(j), ONE / unit, orders, traces)
                else:
                    top = length - 1
                    active.update((index(0), index(top)))
                    for j in range(1, top):
                        apply_shear(images, index(j), index(top), ONE, orders, traces)
                    for j in range(1, top):
                        apply_shear(images, index(top), index(j), step ** (top - j) / unit,
                                    orders, traces)
    active = sorted(active)
    outside = set(range(2 * interval)) - set(active)
    assert all(not creation[i] and not annihilation[i]
               for creation, annihilation in images.values() for i in outside)
    # The right physical generators are regular before transport through Sigma.
    for side in ('left', 'right'):
        for name in names[side]:
            for suffix in ('create', 'annihilate'):
                creation, annihilation = images[name + '_' + suffix]
                for i in active:
                    first, second = creation[i], annihilation[i]
                    if side == 'right':
                        first *= Z ** (-orders[i])
                        second *= Z ** orders[i]
                    assert all(value.order() >= 0 for value in (first, second) if value)
    return {
        'place': place, 'M': interval, 'core': core,
        'core_names': names, 'active_indices': active,
        'active_orders': [orders[i] for i in active],
        'spectator_orders': [orders[i] for i in range(2 * interval) if i not in active],
        'spectator_charges': [1 if place == 'infinity' or i < interval else -1
                              for i in range(2 * interval) if i not in active],
        'active_images': {
            name: {'creation': [creation[i].pair() for i in active],
                   'annihilation': [annihilation[i].pair() for i in active]}
            for name, (creation, annihilation) in images.items()},
        'regular_shears': traces,
    }


def affine_charge_representatives(core=2):
    """Quotient by adding a constant and multiplying by a nonzero F3 scalar."""
    q = 1 << core
    for tail in product(range(3), repeat=q - 1):
        charges = (0,) + tail
        first = next((value for value in charges if value), None)
        if first is None or first == 1:
            yield charges


def projection_program(data, charges, program=None, initial=None):
    core = data['core']
    if len(charges) != 1 << core or any(value not in (0, 1, 2) for value in charges):
        raise ValueError('Supply one F3 charge for every core occupation state.')
    program = Program(data) if program is None else program
    initial = Expr(0) if initial is None else initial
    omega = (1, 2, 3)

    def phase(value, side, letter, power):
        names = [f'{side}_{letter}{i}' for i in range(core)]
        answer = value
        for state, charge in enumerate(charges):
            coefficient = omega[(power * charge) % 3] ^ 1
            if coefficient:
                projected = program.unit(value, names, state, state)
                answer = program.add(answer, program.scale(projected, coefficient))
        return answer

    infinity = data['place'] == 'infinity'
    right = ('u', 'x') if infinity else ('v', 'u')
    left = ('v', 'y') if infinity else ('x', 'y')
    summands = []
    for power in range(3):
        value = initial
        for side, letters in (('right', right), ('left', left)):
            for letter in letters:
                signed = power if letter in ('x', 'u') else -power
                value = phase(value, side, letter, signed)
        summands.append(value)
    last = program.add(program.add(summands[0], summands[1]), summands[2])
    return program, last


def evaluate_program(data, program, last, profile_precision=None):
    size = 1 << program.k
    if program.k > 10:
        raise ValueError('Dense evaluation is limited to ten active modes; use the support theorem for larger cores.')
    diagonal = [sum(order for j, order in enumerate(data['active_orders']) if mask >> j & 1)
                for mask in range(size)]
    window, builder = precision_sizes(program, last.base + min(diagonal), data['M'], profile_precision)
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
    dd = (ctypes.c_int * size)(*diagonal)
    output = (ctypes.c_uint64 * (2 * output_words * size * size))()
    shifts = (ctypes.c_int64 * (size * size))()
    meta = (ctypes.c_int64 * 4)()
    started = time.monotonic()
    status = library(builder).fock_f4_evaluate_v2(
        program.k, len(program.instructions), instructions, len(program.coefficients), packed,
        len(program.cars), co, cc, last.register, dd, last.base, output_words, output, shifts, meta)
    if status:
        raise ArithmeticError(f'Core-charge program failed: {status}; metadata={list(meta)}')
    matrix, shift_matrix = [], []
    for i in range(size):
        row, shift_row = [], []
        for j in range(size):
            index = i * size + j
            pair = tuple(sum(int(output[(2 * index + plane) * output_words + w]) << (64 * w)
                             for w in range(output_words)) for plane in range(2))
            assert not pair[1], 'The F4 phase projection must descend to F2.'
            row.append(pair)
            shift_row.append(int(shifts[index]))
        matrix.append(row)
        shift_matrix.append(shift_row)
    metadata = {
        'matrix_order': size, 'active_modes': program.k, 'registers': len(program.instructions),
        'coefficient_count': len(program.coefficients), 'profile_precision': window,
        'builder_precision': builder, 'unit_denominator_degree': last.denominator.bit_length() - 1,
        'program_degree_bound': program.degree_bounds[last.register],
        'exact_polynomial_matrix': program.degree_bounds[last.register] < builder and int(meta[3]) < window,
        'known_below': retained_known_below(meta, window),
        'matrix_seconds': time.monotonic() - started,
    }
    return matrix, shift_matrix, metadata


def profile(place, interval, charges, core=2, profile_precision=None):
    data = core_data(place, interval, core)
    program, last = projection_program(data, tuple(charges))
    matrix, shifts, metadata = evaluate_program(data, program, last, profile_precision)
    result = MatrixOracle(matrix, len(matrix), shifts, precision=metadata['profile_precision']).profile()
    if not metadata['exact_polynomial_matrix'] and max(result['orders']) >= metadata['known_below']:
        raise ArithmeticError('The Smith profile exceeds the certified coefficient window.')
    orders = []
    for value in result['orders']:
        if value != round(value):
            raise ArithmeticError('The fixed natural profile must have integral orders.')
        orders.append(round(value))
    mode_charges = [1 if place == 'infinity' or i < interval else -1
                    for i in data['active_indices']]
    state_charges = [sum(value for j, value in enumerate(mode_charges) if state >> j & 1)
                     for state in range(len(matrix))]
    assert all(not any(value) or state_charges[i] == state_charges[j]
               for i, row in enumerate(matrix) for j, value in enumerate(row))
    return {
        'place': place, 'interval': interval, 'core': core, 'charges': list(charges),
        'active_orders': orders,
        'active_charges': [state_charges[i] for i in result['rows']],
        'spectator_orders': data['spectator_orders'],
        'spectator_charges': data['spectator_charges'], **metadata,
    }


def pressure(interval, records, log_s, zero_charge=0.0, infinity_charge=0.0,
             gradient=False):
    log_b = 0.0
    derivatives = np.zeros(3)
    for place, (record, shift) in enumerate(zip(records, (zero_charge, infinity_charge))):
        charges = np.asarray(record['active_charges'], dtype=float)
        orders = np.asarray(record['active_orders'], dtype=float) + shift * charges
        log_z = float(logsumexp(-log_s * orders))
        probabilities = np.exp(-log_s * orders - log_z)
        log_b += log_z
        derivatives[0] -= float(probabilities @ orders)
        derivatives[place + 1] -= log_s * float(probabilities @ charges)
        charges = np.asarray(record['spectator_charges'], dtype=float)
        orders = np.asarray(record['spectator_orders'], dtype=float) + shift * charges
        log_b += float(np.logaddexp(0.0, -log_s * orders).sum())
        probabilities = np.exp(-np.logaddexp(0.0, log_s * orders))
        derivatives[0] -= float(probabilities @ orders)
        derivatives[place + 1] -= log_s * float(probabilities @ charges)
    energies = np.asarray([4 * i - 2 * interval + 2 for i in range(interval)], dtype=float)
    physical = np.stack((energies + zero_charge + infinity_charge,
                         energies - zero_charge + infinity_charge))

    def equation(rate):
        return log_b - rate * float(np.logaddexp(0.0, -log_s * physical / rate).sum())

    lower, upper = 1e-9, 4.0
    if equation(lower) <= 0:
        raise ArithmeticError('A nonpositive pressure root requires separate exact review.')
    while equation(upper) > 0:
        upper *= 2
    rate = brentq(equation, lower, upper, xtol=1e-13)
    if not gradient:
        return rate / 2
    probabilities = np.exp(-np.logaddexp(0.0, log_s * physical / rate))
    entropy = float(np.logaddexp(0.0, -log_s * physical / rate).sum()
                    + (probabilities * log_s * physical / rate).sum())
    if entropy < 1e-10:
        raise ArithmeticError('Physical type entropy is too small for a stable derivative.')
    derivatives[0] += float((probabilities * physical).sum())
    derivatives[1] += log_s * float(probabilities[0].sum() - probabilities[1].sum())
    derivatives[2] += log_s * float(probabilities.sum())
    return rate / 2, derivatives / (2 * entropy)


def optimize_scale(interval, records):
    run = minimize_scalar(lambda x: pressure(interval, records, math.exp(x)),
                          bounds=(math.log(0.001), math.log(4.0)), method='bounded',
                          options={'xatol': 1e-9})
    return {'rho': float(run.fun), 'log_s': math.exp(float(run.x)),
            'evaluations': int(run.nfev), 'success': bool(run.success)}


def optimize_bulk(interval, records):
    initial = optimize_scale(interval, records)

    def objective(parameters):
        ell = math.exp(parameters[0])
        try:
            value, derivative = pressure(interval, records, ell, parameters[1], parameters[2], True)
        except ArithmeticError:
            return 100.0, np.zeros(3)
        derivative[0] *= ell
        return value, derivative

    candidates = []
    for j, k in ((0, 0), (1, 0), (-1, 0)):
        run = minimize(objective, [math.log(initial['log_s']), j, k], jac=True,
                       method='L-BFGS-B', bounds=[(math.log(0.001), math.log(4.0))]
                       + [(-4 * interval, 4 * interval)] * 2,
                       options={'maxfun': 300, 'maxiter': 300, 'ftol': 1e-13, 'gtol': 1e-9})
        candidates.append(run)
    best = min(candidates, key=lambda run: run.fun)
    return {'rho': float(best.fun), 'log_s': math.exp(float(best.x[0])),
            'zero_charge': float(best.x[1]), 'infinity_charge': float(best.x[2]),
            'evaluations': sum(int(run.nfev) for run in candidates),
            'success': bool(best.success), 'message': str(best.message)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=8)
    parser.add_argument('--core', type=int, default=2)
    parser.add_argument('--charges', type=int, nargs='+')
    parser.add_argument('--all-charges', action='store_true')
    parser.add_argument('--optimize-bulk', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    charges = list(affine_charge_representatives(args.core)) if args.all_charges else [
        tuple(args.charges) if args.charges is not None else (0, 0, 0, 1)]
    records = []
    started = time.monotonic()
    for charge in charges:
        profiles = [profile(place, args.interval, charge, args.core) for place in ('zero', 'infinity')]
        optimize = optimize_bulk if args.optimize_bulk else optimize_scale
        record = {'charges': list(charge), 'profiles': profiles, **optimize(args.interval, profiles)}
        records.append(record)
        result = {'status': 'exact natural profiles; exploratory floating pressure roots',
                  'interval': args.interval, 'core': args.core,
                  'seconds': time.monotonic() - started, 'records': records}
        args.output.write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps({key: value for key, value in record.items() if key != 'profiles'}), flush=True)


if __name__ == '__main__':
    main()
