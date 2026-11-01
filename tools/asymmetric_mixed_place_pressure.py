#!/usr/bin/env python3
"""Complete supporting-line pressure for arbitrary disjoint mixed places.

For a²+b²=1 the row and column Chernoff slopes are t/a and t/b.
Opposite reconstruction indices then contribute exactly D*t.  Physical
heights are branch-weighted, so the common physical type mean is restored
without an angle-dependent convention change.
"""
from __future__ import annotations

import math
import time

import numpy as np
from scipy.special import logsumexp, softmax

from search_fused_mixed_weighted import MixedPressure, SearchDeadline


class MixedSupportingPressure(MixedPressure):
    def __init__(self, *arguments, **keywords):
        super().__init__(*arguments, **keywords)
        self.angle_index = self.dimension
        self.dimension += 1

    def initial(self, baseline, rng, start):
        result = super().initial(baseline, rng, start)
        result[self.angle_index] = math.pi / 4
        return result

    def bounds(self):
        result = super().bounds()
        result[self.angle_index] = (.015, math.pi / 2 - .015)
        return result

    def evaluate(self, values, smoothing=.003, record=True):
        if time.monotonic() >= self.deadline:
            raise SearchDeadline()
        q, t = self.q, float(values[self.scale_index])
        angle = float(values[self.angle_index])
        branch_weights = (math.cos(angle), math.sin(angle))
        branch_derivatives = (-math.sin(angle), math.cos(angle))
        if min(branch_weights) <= 0:
            raise ValueError('Both supporting-line weights must be positive.')
        type_values = values[self.type_slice].reshape(2, q - 1)
        types = [softmax(np.r_[0.0, row]) for row in type_values]
        entropies = [-float(p @ np.log(p)) for p in types]
        entropy = sum(entropies)
        if entropy < 1e-8:
            raise ArithmeticError('The common physical type entropy is unresolved.')
        mean = np.concatenate((types[0], types[1], types[1], types[0]))
        total = exact_total = self.oracle.degree * t
        derivative = np.zeros(self.dimension)
        derivative[self.scale_index] = self.oracle.degree
        type_derivative = [np.zeros(q), np.zeros(q)]
        place_records = []
        for cut, groups in enumerate(self.groups):
            a, da = branch_weights[cut], branch_derivatives[cut]
            if self.train_allocations and len(groups) > 1:
                allocations = softmax(np.r_[values[self.allocation_slices[cut]], 0.0])
            else:
                allocations = np.asarray([group['allocation'] for group in groups])
            costs, data = [], []
            for index, (group, allocation) in enumerate(zip(groups, allocations)):
                local_scale = t * len(group['labels']) / (a * allocation)
                weighted_heights = self.height_vector(values[self.height_slices[cut][index]])
                heights = weighted_heights / a
                log_sum, grad_h, grad_s, orders = self.oracle.log_profile_with_gradient(
                    cut, group['point'], heights, local_scale)
                shift = float(mean @ weighted_heights)
                cost = a * log_sum + shift
                # Homogeneity gives this derivative even when a pivot
                # chamber changes; the selected chamber supplies a subgradient.
                channel_entropy = log_sum - float(grad_h @ heights) - grad_s * local_scale
                costs.append(cost)
                data.append((weighted_heights, heights, grad_h, grad_s, local_scale,
                             orders, shift, log_sum, channel_entropy))
            costs = np.asarray(costs)
            selected = softmax(costs / smoothing)
            total += smoothing * float(logsumexp(costs / smoothing))
            exact_total += float(np.max(costs))
            allocation_gradient = np.zeros(len(groups))
            records = []
            for index, (group, allocation, weight, local) in enumerate(zip(groups, allocations, selected, data)):
                weighted_h, h, gh, gs, local_scale, orders, shift, log_sum, channel_entropy = local
                derivative[self.height_slices[cut][index]] += weight * self.height_gradient(gh + mean)
                derivative[self.scale_index] += weight * gs * local_scale * a / t
                derivative[self.angle_index] += weight * da * channel_entropy
                allocation_gradient[index] = -weight * gs * local_scale * a
                type_derivative[0] += weight * (weighted_h[:q] + weighted_h[3 * q:])
                type_derivative[1] += weight * (weighted_h[q:2 * q] + weighted_h[2 * q:3 * q])
                records.append({'point': None if group['point'] is None else list(group['point']),
                                'conjugate_labels': group['labels'], 'group_allocation': float(allocation),
                                'allocation_per_root': float(allocation / len(group['labels'])),
                                'branch_weight': a, 'local_log_scale': float(local_scale),
                                'branch_weighted_heights': weighted_h.tolist(), 'physical_heights': h.tolist(),
                                'complete_orders': orders.tolist(), 'weighted_type_shift': shift,
                                'log_profile_sum': float(log_sum), 'weighted_complete_cost': float(costs[index]),
                                'channel_entropy': float(channel_entropy)})
            if self.train_allocations and len(groups) > 1:
                derivative[self.allocation_slices[cut]] += (
                    allocation_gradient - allocations * np.sum(allocation_gradient))[:-1]
            place_records.append(records)
        support, exact_support = total / entropy, exact_total / entropy
        derivative /= entropy
        type_gradient = []
        for p, h, gradient in zip(types, entropies, type_derivative):
            cost_gradient = p * (gradient - float(p @ gradient)) / entropy
            entropy_gradient = -p * (np.log(p) + h)
            type_gradient.extend((cost_gradient - support * entropy_gradient / entropy)[1:])
        derivative[self.type_slice] = type_gradient
        self.queries += 1
        if record and (self.best is None or exact_support < self.best['support_value']):
            self.best = {'support_value': float(exact_support), 'smoothed_support_value': float(support),
                         'threshold_margin': float(exact_total - entropy),
                         'weights': list(branch_weights), 'angle': angle,
                         'paired_types': [p.tolist() for p in types], 'entropy_sum': entropy,
                         'outer_log_scale': t, 'parameters': np.asarray(values).tolist(),
                         'place_records': place_records, 'all_physical_costs_included': True}
        return float(support), derivative


class ScaledMixedSupportingPressure(MixedSupportingPressure):
    """Use branch-weighted physical heights b=t*B for stable native windows."""
    def initial(self, baseline, rng, start):
        values = super().initial(baseline, rng, start)
        for cut in self.height_slices:
            for current in cut:
                values[current] /= values[self.scale_index]
        return values

    def evaluate(self, values, smoothing=.003, record=True):
        actual = np.asarray(values, float).copy()
        t = float(values[self.scale_index])
        for cut in self.height_slices:
            for current in cut:
                actual[current] *= t
        previous = self.best
        support, derivative = super().evaluate(actual, smoothing, record)
        dt = derivative[self.scale_index]
        for cut in self.height_slices:
            for current in cut:
                dt += float(np.asarray(values[current]) @ derivative[current])
                derivative[current] *= t
        derivative[self.scale_index] = dt
        if self.best is not previous:
            self.best['optimizer_parameters'] = np.asarray(values).tolist()
            self.best['height_parameterization'] = 'Branch-weighted heights b=t*B; actual local heights are b/a_cut.'
        return support, derivative
