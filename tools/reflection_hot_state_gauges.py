#!/usr/bin/env python3
"""Arbitrary physical hot-core state gauges on exact reflection defects."""
from __future__ import annotations

from itertools import product
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_hot_block_profiles import normalize_factor
from reflection_hot_state_coordinates import gaussian_coordinates


def physical_index(labels, alphabet):
    value = 0
    for label in labels:
        value = alphabet * value + label
    return value


def normalize_mask(mask, alphabet):
    if mask is None:
        return None
    if callable(mask):
        mask = tuple(mask(*labels) for labels in product(range(alphabet), repeat=4))
    else:
        mask = tuple(mask)
    if len(mask) != alphabet ** 4 or any(value not in (0, 1, 2, 3) for value in mask):
        raise ValueError('Supply an F4 mask on every physical core tuple.')
    if any(mask[physical_index((x, y, y, x), alphabet)] != 1
           for x, y in product(range(alphabet), repeat=2)):
        raise ValueError('The coefficient mask must retain every crossing diagonal with coefficient one.')
    return mask


def defect_block(coordinates, factor, sector, mask=None):
    data, left, right = coordinates['data'], coordinates['left'], coordinates['right']
    k, p, place = data['core'], data['ramification'], data['place']
    sign = 1 if place == 'zero' else -1
    gauged = (0, 1) if place == 'zero' else (2, 3)
    left_letters = [name.split('_')[1][0] for name in left['names'][::k]]
    right_letters = [name.split('_')[1][0] for name in right['names'][::k]]
    terms = []
    for physical_mask, polynomial in factor.items():
        left_mask = sum(1 << (group * k + j) for group, letter in enumerate(left_letters)
                        for j in range(k) if physical_mask >> ('xyuv'.index(letter) * k + j) & 1)
        right_mask = sum(1 << (group * k + j) for group, letter in enumerate(right_letters)
                         for j in range(k) if physical_mask >> ('xyuv'.index(letter) * k + j) & 1)
        height = p * sum(2 + 4 * j for leg in gauged for j in range(k)
                         if physical_mask >> (leg * k + j) & 1)
        terms += [(left_mask, right_mask, sign * degree + height, coefficient)
                  for degree, coefficient in polynomial.items()]
    rows, columns = coordinates['row_groups'][sector], coordinates['column_groups'][sector]
    matrix, shifts = [], []
    maximum_bits = 1
    for row in rows:
        row_state = row // left['multiplicity']
        out, offsets = [], []
        for column in columns:
            column_state = column // right['multiplicity']
            selected = []
            for a, b, shift, coefficient in terms:
                if row_state & a != a or column_state & b != b:
                    continue
                i, j = row - a * left['multiplicity'], column - b * right['multiplicity']
                parity = left['row_parities'][i]
                assert parity == right['column_parities'][j]
                local_row, local_column = coordinates['row_positions'][i], coordinates['column_positions'][j]
                value, high = coordinates['matrices'][parity][local_row][local_column]
                assert not high
                if value:
                    selected.append((value, coordinates['shifts'][parity][local_row][local_column] + shift,
                                     coefficient))
            if not selected:
                out.append((0, 0));offsets.append(0);continue
            base = min(shift for value, shift, coefficient in selected)
            planes = [0, 0]
            for value, shift, coefficient in selected:
                value <<= shift - base
                for bit in range(2):
                    if coefficient >> bit & 1:
                        planes[bit] ^= value
            if mask is not None:
                labels = tuple(left['labels'][letter][row] if letter in left['labels']
                               else right['labels'][letter][column] for letter in 'xyuv')
                coefficient = mask[physical_index(labels, 1 << k)]
                a, b = planes
                planes = ([0, 0] if coefficient == 0 else ([a, b] if coefficient == 1
                          else ([b, a ^ b] if coefficient == 2 else [a ^ b, a])))
            out.append(tuple(planes));offsets.append(base)
            maximum_bits = max(maximum_bits, max(planes).bit_length())
        matrix.append(out);shifts.append(offsets)
    metadata = dict(coordinates['metadata'])
    metadata['known_below'] += min(shift for a, b, shift, coefficient in terms)
    metadata['profile_precision'] = max(metadata['profile_precision'], 1 << (maximum_bits - 1).bit_length())
    if metadata['profile_precision'] > 8192:
        raise ArithmeticError('The complete retained defect array exceeds the native precision limit.')
    metadata['defect_source_maximum_bits'] = maximum_bits
    return matrix, shifts, metadata


class HotStatePressure:
    def __init__(self, interval, factor, ramification=1, core=1, mask=None, profile_precision=512,
                 internal_regular=False):
        self.interval, self.ramification, self.core, self.q = interval, ramification, core, 1 << core
        self.internal_regular = internal_regular
        self.factor = normalize_factor(factor, core)
        self.mask = normalize_mask(mask, self.q)
        if internal_regular and self.mask is not None:
            raise ValueError('Internal-factor removal supports exterior defects, not arbitrary coefficient masks.')
        self.records, self.oracles = [], []
        self.queries = 0
        for place in ('zero', 'infinity'):
            coordinates = gaussian_coordinates(interval, place, ramification, core, profile_precision, internal_regular)
            sides, metadata = [], []
            for sector in range(2):
                matrix, shifts, md = defect_block(coordinates, self.factor, sector, self.mask)
                sides.append(MatrixOracle(matrix, len(matrix), shifts, md['profile_precision']))
                metadata.append(md)
            self.records.append((coordinates, metadata))
            self.oracles.append(sides)

    def evaluate(self, log_s, offsets=None, gradient=False):
        if not math.isfinite(log_s) or log_s <= 0:
            raise ValueError('Use a positive finite logarithmic scale.')
        offsets = np.zeros((2, 4, self.q)) if offsets is None else np.asarray(offsets, dtype=float)
        if offsets.shape != (2, 4, self.q):
            raise ValueError('State offsets have shape (two places, four legs, core alphabet).')
        if self.internal_regular:
            for state in range(self.q):
                for bit in range(self.core):
                    if not state >> bit & 1 and np.any(offsets[:, :, state | (1 << bit)] < offsets[:, :, state]):
                        raise ValueError('The smaller core requires state heights monotone under adding an occupied bit.')
        letters = dict(zip('xyuv', range(4)))
        log_b, ds = 0.0, 0.0
        derivative = np.zeros_like(offsets)
        profiles = []
        for place, ((coordinates, metadata), oracles) in enumerate(zip(self.records, self.oracles)):
            left, right, data = coordinates['left'], coordinates['right'], coordinates['data']
            rows = np.zeros(left['unit_basis_rank'])
            columns = np.zeros(right['unit_basis_rank'])
            for letter, labels in left['labels'].items():
                rows += offsets[place, letters[letter]][labels]
            for letter, labels in right['labels'].items():
                columns += offsets[place, letters[letter]][labels]
            orders, row_indices, column_indices = [], [], []
            for sector, (oracle, md) in enumerate(zip(oracles, metadata)):
                rg, cg = coordinates['row_groups'][sector], coordinates['column_groups'][sector]
                while True:
                    try:
                        result = oracle.profile(log_s * rows[rg], log_s * columns[cg], scale=log_s)
                        break
                    except ArithmeticError as error:
                        if '(status 3)' not in str(error) or oracle.precision >= 8192:
                            raise
                        # The packed source is already complete at its reported
                        # width. Only the elimination workspace is widened.
                        # A truncated source still faces the independent
                        # absolute known-below check immediately afterwards.
                        oracle.precision *= 2
                md['smith_precision'] = oracle.precision
                if not md['exact_polynomial_matrix']:
                    known = log_s * (md['known_below'] + float(rows[rg].min()) + float(columns[cg].min()))
                    if max(result['orders']) >= known:
                        raise ArithmeticError('The weighted profile exceeds its certified retained coefficient window.')
                orders += result['orders']
                row_indices += [rg[i] for i in result['rows']]
                column_indices += [cg[i] for i in result['cols']]
            orders = np.asarray(orders)
            log_z = float(logsumexp(-orders))
            probabilities = np.exp(-orders - log_z)
            log_b += log_z
            ds -= float(probabilities @ orders) / log_s
            for basis, indices in ((left, row_indices), (right, column_indices)):
                for letter, labels in basis['labels'].items():
                    derivative[place, letters[letter]] -= log_s * np.bincount(
                        np.asarray(labels)[indices], weights=probabilities, minlength=self.q)
            spectators = np.asarray(data['spectator_orders'], dtype=float)
            log_b += float(np.logaddexp(0, -log_s * spectators).sum())
            ds -= float(spectators @ np.exp(-np.logaddexp(0, log_s * spectators)))
            profiles.append(orders.tolist())
        natural = self.ramification * np.asarray(physical_bit_energies(self.interval), dtype=float)
        hot = [self.interval // 2 - 1 - j for j in range(self.core)]
        cold = [i for i in range(len(natural)) if i not in hot]
        core_energy = np.asarray([sum(natural[hot[j]] for j in range(self.core) if state >> j & 1)
                                  for state in range(self.q)], dtype=float)
        gx = core_energy + offsets[:, (0, 3), :].sum(axis=(0, 1))
        gy = core_energy + offsets[:, (1, 2), :].sum(axis=(0, 1))
        cold_energy = natural[cold]
        def equation(rate):
            physical = float(logsumexp(-log_s * gx / rate) + logsumexp(-log_s * gy / rate)
                             + 2 * np.logaddexp(0, -log_s * cold_energy / rate).sum())
            return log_b - rate * physical
        if equation(1e-9) <= 0:
            raise ArithmeticError('A nonpositive pressure root requires separate review.')
        upper = 4.0
        while equation(upper) > 0:
            upper *= 2
        rate = brentq(equation, 1e-9, upper, xtol=1e-13)
        px = np.exp(-log_s * gx / rate - logsumexp(-log_s * gx / rate))
        py = np.exp(-log_s * gy / rate - logsumexp(-log_s * gy / rate))
        cold_occupations = np.exp(-np.logaddexp(0, log_s * cold_energy / rate))
        entropy = -sum(float(p[p > 0] @ np.log(p[p > 0])) for p in (px, py))
        entropy += 2 * float(np.logaddexp(0, -log_s * cold_energy / rate).sum()
                             + cold_occupations @ (log_s * cold_energy / rate))
        self.queries += 1
        detail = {'rho': rate / 2, 'log_s': log_s, 'active_orders': profiles,
                  'log_profile_product': log_b, 'type_x_core': px.tolist(),
                  'type_y_core': py.tolist(), 'entropy': entropy}
        if not gradient:
            return rate / 2, detail
        if entropy < 1e-10:
            raise ArithmeticError('The physical entropy is too small for stable derivatives.')
        for place in range(2):
            derivative[place] += log_s * np.stack((px, py, py, px))
        ds += float(px @ gx + py @ gy + 2 * cold_occupations @ cold_energy)
        return rate / 2, (derivative / (2 * entropy), ds / (2 * entropy)), detail
