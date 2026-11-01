#!/usr/bin/env python3
"""Original Fock minors with independently chosen B/Y and X/V alphabets.

Occupation blocks are square of total order |B| |X| in both cuts. Their
generic full rank follows from the principal Pascal endpoint; no Gaussian
profile substitution or deletion-closure assumption is used here.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import product
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp, softmax

from cube_f4_bivariate_smith import ZERO, add, multiply
from native_matrix_smith_f4 import MatrixOracle
from pascal_coordinate_sectors import SectorOracle
from weighted_pascal_flags import natural_heights


class CoordinateTensorBuilder:
    """Share a bounded cache of original minors across coordinate restrictions."""
    def __init__(self, seed):
        self.seed, self.n = seed, seed['n']
        if not 1 <= self.n <= 16:
            raise ValueError('Use one through sixteen ambient modes.')
        matrix = [[tuple(v) if isinstance(v, (tuple, list)) else (v, 0) for v in row]
                  for row in seed['matrix']]

        @lru_cache(1 << 20)
        def minor(rows, columns):
            if rows.bit_count() != columns.bit_count():
                return ZERO
            if not rows:
                return (1, 0)
            i = min((a for a in range(2*self.n) if rows >> a & 1),
                    key=lambda a: sum(matrix[a][b] != ZERO for b in range(2*self.n)
                                      if columns >> b & 1))
            remaining, value = columns, ZERO
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                j = bit.bit_length()-1
                if matrix[i][j] != ZERO:
                    value = add(value, multiply(matrix[i][j], minor(rows ^ (1 << i), columns ^ bit)))
            return value
        self.minor = minor

    def tensor(self, alphabet_b, alphabet_x):
        alphabets = [list(alphabet_b), list(alphabet_x)]
        for alphabet in alphabets:
            if (not 2 <= len(alphabet) <= 64 or len(set(alphabet)) != len(alphabet)
                    or any(type(w) is not int or not 0 <= w < 1 << self.n for w in alphabet)):
                raise ValueError('Supply two through 64 distinct words in each paired alphabet.')
        ab, ax = alphabets
        by_occupation = defaultdict(list)
        for v, y in product(ax, ab):
            by_occupation[v.bit_count()+y.bit_count()].append((v, y))
        values = {}
        before = self.minor.cache_info()
        for b, x in product(ab, ax):
            for v, y in by_occupation[b.bit_count()+x.bit_count()]:
                value = self.minor(v | (y << self.n), b | (x << self.n))
                if value != ZERO:
                    values[b, x, v, y] = value
        after = self.minor.cache_info()
        return {'alphabet_b':ab, 'alphabet_x':ax, 'values':values,
                'minor_cache_entries':after.currsize,
                'minor_cache_misses':after.misses-before.misses}


class RectangularSectorOracle(SectorOracle):
    def __init__(self, tensor, place, precision=512):
        self.place = place
        ab, ax = tensor['alphabet_b'], tensor['alphabet_x']
        a, b = len(ab), len(ax)
        ib, ix = {v:i for i,v in enumerate(ab)}, {v:i for i,v in enumerate(ax)}
        rows, columns = defaultdict(list), defaultdict(list)
        for i, word_b in enumerate(ab):
            for j, word_x in enumerate(ax):
                charge = word_b.bit_count() + (-1 if place == 0 else 1)*word_x.bit_count()
                rows[charge].append(i+a*j)
                columns[charge].append(j+b*i)
        blocks, locations = {}, {}
        for charge in sorted(rows):
            if len(rows[charge]) != len(columns[charge]):
                raise ArithmeticError('Occupation blocks must be square.')
            if len(rows[charge]) > 1024:
                raise ValueError('A selected occupation block exceeds the native order limit.')
            blocks[charge] = [[ZERO]*len(columns[charge]) for _ in rows[charge]]
            locations[charge] = ({v:i for i,v in enumerate(rows[charge])},
                                 {v:i for i,v in enumerate(columns[charge])})
        for (wb, wx, wv, wy), value in tensor['values'].items():
            charge = wb.bit_count() + (-1 if place == 0 else 1)*(wv if place == 0 else wx).bit_count()
            row = ib[wb]+a*ix[wv if place == 0 else wx]
            column = ix[wx if place == 0 else wv]+b*ib[wy]
            ri, ci = locations[charge]
            blocks[charge][ri[row]][ci[column]] = value
        self.blocks = [(c, rows[c], columns[c],
                        MatrixOracle(blocks[c], len(rows[c]), precision=precision)) for c in sorted(rows)]
        assert sum(len(r) for _, r, _, _ in self.blocks) == a*b


class RectangularPascalPressure:
    def __init__(self, seed, alphabet_b, alphabet_x, precision=512, builder=None):
        self.seed, self.n = seed, seed['n']
        if builder is not None and builder.seed != seed:
            raise ValueError('The shared minor cache belongs to a different source tensor.')
        tensor = (builder or CoordinateTensorBuilder(seed)).tensor(alphabet_b, alphabet_x)
        self.alphabets = [tensor['alphabet_b'], tensor['alphabet_x']]
        self.a, self.b = map(len, self.alphabets)
        self.leg_alphabets = [self.alphabets[i] for i in (0,1,1,0)]
        self.offsets = [0, self.a, self.a+self.b, self.a+2*self.b]
        self.size, self.queries, self.last = 2*(self.a+self.b), 0, None
        self.construction = {k:tensor[k] for k in ('minor_cache_entries','minor_cache_misses')}
        self.construction.update(nonzero_original_entries=len(tensor['values']),
            nonbinary_original_entries=sum(bool(v[1]) for v in tensor['values'].values()))
        self.cuts = [(RectangularSectorOracle(tensor, place, precision), legs)
                     for place, legs in enumerate(((0,2,1,3), (0,1,2,3)))]

    def heights(self, parameters):
        return [[parameters[p*self.size+o:p*self.size+o+len(letters)]
                 for o, letters in zip(self.offsets, self.leg_alphabets)] for p in range(2)]

    def physical(self, parameters):
        h = self.heights(parameters)
        return [h[0][0]+h[0][3]+h[1][0]+h[1][3], h[0][1]+h[0][2]+h[1][1]+h[1][2]]

    def lift_mode_parameters(self, parameters):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (8*self.n+1,):
            raise ValueError('Supply eight mode-height vectors and a common scale.')
        modes = parameters[:-1].reshape(2,4,self.n)
        return np.r_[[sum(modes[p,leg,bit] for bit in range(self.n) if w >> bit & 1)
                       for p in range(2) for leg in range(4) for w in self.leg_alphabets[leg]],
                      parameters[-1]]

    def natural(self, scale):
        return self.lift_mode_parameters(np.r_[[scale*v for p in natural_heights(self.seed) for v in p], scale])

    def rate(self, log_b, physical):
        def equation(rate):
            return log_b-rate*sum(float(logsumexp(-v/rate)) for v in physical)
        if equation(1e-8) <= 0:
            raise ArithmeticError('The positive Gibbs root is unresolved near zero entropy.')
        high = 4.
        while equation(high) > 0:
            high *= 2
            if high > 1e8:
                raise ArithmeticError('The paired pressure root was not bracketed.')
        return brentq(equation, 1e-8, high, xtol=1e-13), equation

    def fitted_natural_start(self):
        unit = self.natural(1.)
        self.evaluate(unit)
        orders = [np.asarray(cut['orders']) for cut in self.last['complete_original_profiles']]
        physical = self.physical(unit)
        def value(log_scale):
            scale = math.exp(log_scale)
            return self.rate(sum(float(logsumexp(-scale*v)) for v in orders),
                             [scale*v for v in physical])[0]/2
        fit = minimize_scalar(value, bounds=(-8., math.log(2.)), method='bounded')
        return unit*math.exp(fit.x)

    def evaluate(self, parameters, gradient=False):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (2*self.size+1,) or not np.all(np.isfinite(parameters)) or parameters[-1] <= 0:
            raise ValueError('Supply eight symbol-height vectors and a positive common scale.')
        h, scale = self.heights(parameters), parameters[-1]
        derivative, log_b, records = np.zeros_like(parameters), 0., []
        for place, (oracle, legs) in enumerate(self.cuts):
            row = [h[place][legs[0]][i]+h[place][legs[1]][j]
                   for j in range(self.b) for i in range(self.a)]
            col = [h[place][legs[2]][i]+h[place][legs[3]][j]
                   for j in range(self.a) for i in range(self.b)]
            result = oracle.profile(row, col, scale, infinity=bool(place))
            orders = np.asarray(result['orders'])
            log_b += float(logsumexp(-orders))
            probabilities = softmax(-orders)
            for weight, r, c in zip(probabilities, result['rows'], result['cols']):
                for leg, symbol in ((legs[0],r % self.a), (legs[1],r // self.a),
                                    (legs[2],c % self.b), (legs[3],c // self.b)):
                    derivative[place*self.size+self.offsets[leg]+symbol] -= weight
            derivative[-1] -= float(probabilities @ np.diff([0]+result['base_prefix']))
            records.append(result)
        physical = self.physical(parameters)
        rate, equation = self.rate(log_b, physical)
        probabilities = [softmax(-v/rate) for v in physical]
        entropy = sum(float(logsumexp(-v/rate)+(v/rate*p).sum()) for v,p in zip(physical, probabilities))
        if entropy < 1e-8:
            raise ArithmeticError('The physical entropy is too small for a stable rate.')
        for place in range(2):
            for leg, offset in enumerate(self.offsets):
                prob = probabilities[0 if leg in (0,3) else 1]
                derivative[place*self.size+offset:place*self.size+offset+len(prob)] += prob
        self.queries += 1
        self.last = {'rho':float(rate/2), 'entropy':entropy, 'critical_pressure':equation(math.sqrt(2)),
                     'parameters':parameters.tolist(), 'pair_probabilities':[p.tolist() for p in probabilities],
                     'complete_original_profiles':records}
        return (rate/2, derivative/(2*entropy)) if gradient else rate/2
