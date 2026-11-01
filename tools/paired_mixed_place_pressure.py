#!/usr/bin/env python3
"""Complete mixed-place pressure for unequal paired physical alphabets.

The four dimensions are (q_x,q_y,q_y,q_x).  One pair of physical types is
shared by every local place.  Conjugate labels are charged separately even
when a common irreducible factor supplies their identical local profile.
Floating evaluation is exploratory, not a threshold certificate.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.special import logsumexp, softmax

from search_fused_mixed_weighted import SearchDeadline, grouped_places, point_key


def checked_place_groups(baseline):
    """Reject repeated conjugates, including aliases of a rational point."""
    canonical, result = [], []
    for side in ('row', 'column'):
        places, allocations = baseline[side+'_places'], baseline[side+'_allocations']
        if len(places) != len(allocations):
            raise ValueError('Each labelled place needs its own allocation.')
        seen = set()
        for label, allocation in zip(places, allocations):
            if allocation <= 1e-7:
                continue
            factor = point_key(label)
            if factor == (1, 1):
                raise ValueError('Selected reconstruction places must avoid the strict endpoint.')
            if factor is None:
                index = 0
            else:
                degree = len(factor)-1
                if degree < 1 or factor[-1] != 1 or any(value not in range(4) for value in factor):
                    raise ValueError('Closed-place factors must be monic F4 polynomials.')
                index = 0 if isinstance(label, int) else int(label.rsplit(':', 1)[1])
                if not 0 <= index < degree:
                    raise ValueError('A conjugate index must lie below the residue degree.')
            key = factor, index
            if key in seen:
                raise ValueError('A conjugate place cannot be counted twice.')
            seen.add(key)
        canonical.append(seen)
        result.append(grouped_places(places, allocations))
    if not canonical[0].isdisjoint(canonical[1]):
        raise ValueError('The row and column reconstructions need disjoint conjugate places.')
    return result


class PairedMixedPressure:
    def __init__(self, oracle, baseline, train_allocations=True, deadline=float('inf')):
        self.oracle, self.deadline = oracle, deadline
        self.q_x = int(oracle.q_x)
        self.q_y = int(oracle.q_y)
        if min(self.q_x, self.q_y) < 2:
            raise ValueError('Both paired physical alphabets must have positive entropy capacity.')
        self.dimensions = (self.q_x, self.q_y, self.q_y, self.q_x)
        self.leg_offsets = np.cumsum([0]+list(self.dimensions))
        self.leg_slices = [slice(int(self.leg_offsets[i]), int(self.leg_offsets[i+1])) for i in range(4)]
        self.active_heights = np.asarray([int(self.leg_offsets[leg])+letter
                                          for leg, size in enumerate(self.dimensions)
                                          for letter in range(1, size)], dtype=int)
        self.height_size = len(self.active_heights)
        self.groups = checked_place_groups(baseline)
        self.train_allocations = train_allocations
        offset, self.height_slices = 0, []
        for groups in self.groups:
            current = []
            for _ in groups:
                current.append(slice(offset, offset+self.height_size))
                offset += self.height_size
            self.height_slices.append(current)
        self.type_slice = slice(offset, offset+self.q_x+self.q_y-2)
        self.type_slices = [slice(offset, offset+self.q_x-1),
                            slice(offset+self.q_x-1, offset+self.q_x+self.q_y-2)]
        offset = self.type_slice.stop
        self.allocation_slices = []
        for groups in self.groups:
            size = len(groups)-1 if train_allocations else 0
            self.allocation_slices.append(slice(offset, offset+size))
            offset += size
        self.scale_index, self.dimension = offset, offset+1
        self.queries, self.best = 0, None

    def height_vector(self, parameters):
        result = np.zeros(int(self.leg_offsets[-1]))
        result[self.active_heights] = parameters
        return result

    def initial(self, baseline, rng, start=0):
        values = np.zeros(self.dimension)
        scale = float(baseline['outer_log_scale'])
        if start:
            scale *= (.6, 1., 1.7, 2.5)[start % 4]
        values[self.scale_index] = np.clip(scale, .004, 3.5)
        for cut, groups in enumerate(self.groups):
            if self.train_allocations and len(groups) > 1:
                allocations = np.asarray([group['allocation'] for group in groups])
                values[self.allocation_slices[cut]] = np.clip(np.log(allocations[:-1]/allocations[-1]), -6, 6)
            if start:
                for current in self.height_slices[cut]:
                    values[current] = rng.normal(0, .12*scale, self.height_size)
        if start:
            slopes = ((.5, .5), (1., -1.), (-.7, -.3), (1.5, .3), (0., 0.))[start % 5]
            logits = np.concatenate([slope*np.arange(1, size) for slope, size in zip(slopes, (self.q_x, self.q_y))])
            values[self.type_slice] = np.clip(logits+rng.normal(0, .05, len(logits)), -5, 5)
        return values

    def bounds(self):
        bounds = [(-12., 12.)]*self.dimension
        for i in range(self.type_slice.start, self.type_slice.stop):
            bounds[i] = (-5., 5.)
        for current in self.allocation_slices:
            for i in range(current.start, current.stop):
                bounds[i] = (-6., 6.)
        bounds[self.scale_index] = (.002, 4.)
        return bounds

    def evaluate(self, values, smoothing=.01, record=True):
        if time.monotonic() >= self.deadline:
            raise SearchDeadline()
        if smoothing <= 0:
            raise ValueError('Use positive smoothing; the exact maximum is recorded separately.')
        values = np.asarray(values, dtype=float)
        scale = float(values[self.scale_index])
        if scale <= 0:
            raise ValueError('The common outer log scale must be positive.')
        types = [softmax(np.r_[0., values[current]]) for current in self.type_slices]
        entropies = [-float(p @ np.log(p)) for p in types]
        entropy = sum(entropies)
        if entropy < 1e-8:
            raise ArithmeticError('Physical entropy is unresolved.')
        mean = np.concatenate((types[0], types[1], types[1], types[0]))
        exact_total = total = self.oracle.degree*scale
        derivative = np.zeros(self.dimension)
        derivative[self.scale_index] = self.oracle.degree
        type_derivative = [np.zeros(self.q_x), np.zeros(self.q_y)]
        place_records = []
        for cut, groups in enumerate(self.groups):
            if self.train_allocations and len(groups) > 1:
                allocations = softmax(np.r_[values[self.allocation_slices[cut]], 0.])
            else:
                allocations = np.asarray([group['allocation'] for group in groups])
            costs, local = [], []
            for index, (group, allocation) in enumerate(zip(groups, allocations)):
                local_scale = scale*len(group['labels'])/allocation
                heights = self.height_vector(values[self.height_slices[cut][index]])
                log_sum, gradient_h, gradient_scale, orders = self.oracle.log_profile_with_gradient(
                    cut, group['point'], heights, local_scale)
                shift = float(mean @ heights)
                costs.append(log_sum+shift)
                local.append((heights, gradient_h+mean, gradient_scale, local_scale, orders, shift, log_sum))
            costs = np.asarray(costs)
            weights = softmax(costs/smoothing)
            total += smoothing*float(logsumexp(costs/smoothing))
            exact_total += float(np.max(costs))
            allocation_gradient = np.zeros(len(groups))
            cut_records = []
            for index, (group, allocation, weight, data) in enumerate(zip(groups, allocations, weights, local)):
                heights, grad_h, grad_scale, local_scale, orders, shift, log_sum = data
                derivative[self.height_slices[cut][index]] += weight*grad_h[self.active_heights]
                derivative[self.scale_index] += weight*grad_scale*local_scale/scale
                allocation_gradient[index] = -weight*grad_scale*local_scale
                type_derivative[0] += weight*(heights[self.leg_slices[0]]+heights[self.leg_slices[3]])
                type_derivative[1] += weight*(heights[self.leg_slices[1]]+heights[self.leg_slices[2]])
                cut_records.append({'point': None if group['point'] is None else list(group['point']),
                                    'conjugate_labels': group['labels'], 'group_allocation': float(allocation),
                                    'allocation_per_root': float(allocation/len(group['labels'])),
                                    'local_log_scale': float(local_scale), 'heights': heights.tolist(),
                                    'restored_orders': np.asarray(orders).tolist(), 'type_shift': shift,
                                    'log_profile_sum': float(log_sum), 'complete_cost': float(costs[index])})
            if self.train_allocations and len(groups) > 1:
                transformed = allocation_gradient-allocations*np.sum(allocation_gradient)
                derivative[self.allocation_slices[cut]] += transformed[:-1]
            place_records.append(cut_records)
        rho, exact_rho = total/(2*entropy), exact_total/(2*entropy)
        derivative /= 2*entropy
        for current, probabilities, h, gradient in zip(self.type_slices, types, entropies, type_derivative):
            cost_gradient = probabilities*(gradient-float(probabilities @ gradient))/(2*entropy)
            entropy_gradient = -probabilities*(np.log(probabilities)+h)
            derivative[current] = (cost_gradient-rho*entropy_gradient/entropy)[1:]
        self.queries += 1
        if record and (self.best is None or exact_rho < self.best['rho']):
            self.best = {'rho': float(exact_rho), 'smoothed_rho': float(rho),
                         'critical_pressure': float(exact_total-np.sqrt(2)*entropy),
                         'entropy_sum': entropy, 'paired_types': [p.tolist() for p in types],
                         'paired_dimensions': [self.q_x, self.q_y], 'outer_log_scale': scale,
                         'parameters': values.tolist(), 'place_records': place_records,
                         'all_physical_costs_included': True}
        return float(rho), derivative


class ScaledPairedMixedPressure(PairedMixedPressure):
    """Use physical heights h=t*H without changing the evaluated criterion."""
    def initial(self, baseline, rng, start=0):
        result = super().initial(baseline, rng, start)
        scale = result[self.scale_index]
        for cut in self.height_slices:
            for current in cut:
                result[current] /= scale
        return result

    def evaluate(self, values, smoothing=.01, record=True):
        values = np.asarray(values, dtype=float)
        actual = values.copy()
        scale = float(values[self.scale_index])
        for cut in self.height_slices:
            for current in cut:
                actual[current] *= scale
        previous = self.best
        value, derivative = super().evaluate(actual, smoothing, record)
        scale_derivative = derivative[self.scale_index]
        for cut in self.height_slices:
            for current in cut:
                scale_derivative += float(values[current] @ derivative[current])
                derivative[current] *= scale
        derivative[self.scale_index] = scale_derivative
        if self.best is not previous:
            self.best['optimizer_parameters'] = values.tolist()
            self.best['height_parameterization'] = 'Actual physical heights h equal common outer log scale t times optimized H.'
        return value, derivative
