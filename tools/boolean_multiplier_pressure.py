#!/usr/bin/env python3
"""Optimize the complete one-support multiplier envelope at fixed biases."""

from functools import lru_cache

import numpy as np

from boolean_multiplier_profiles import complete_multiplier_automaton


class MultiplierEnvelope:
    def __init__(self, maximum_boolean_degree):
        self.maximum_boolean_degree = maximum_boolean_degree
        self.automaton = complete_multiplier_automaton(maximum_boolean_degree)
        self.pattern_degrees = np.array(
            [0 if not pattern else pattern.bit_length() - 1
             for pattern in range(self.automaton['pattern_count'])])
        self.profiles = []
        for row in self.automaton['continuous_profiles']:
            self.profiles.append({'kind': 'continuous', **row})
        for tail_index, tail in enumerate(self.automaton['singleton_tails']):
            for kind in ('startup', 'cycle'):
                for row in tail[kind]:
                    self.profiles.append({
                        'kind': kind, 'pattern': tail['pattern'],
                        'active_patterns': (tail['pattern'],),
                        'outside_exponents': tail['outside_exponents'],
                        'divergence': row['divergence'],
                        'degree_offset': row['degree_offset'],
                        'period': tail['period'] if kind == 'cycle' else 0,
                        'tail_index': tail_index,
                    })
        self.active = np.array([
            [pattern in row['active_patterns'] for pattern in range(len(self.pattern_degrees))]
            for row in self.profiles])
        self.exponents = np.array([row['outside_exponents'] for row in self.profiles])
        self.offsets = np.array([row.get('degree_offset', 0) for row in self.profiles])
        self.continuous = np.array([row['kind'] == 'continuous' for row in self.profiles])
        self.lower = np.array([row.get('lower_degree', row.get('divergence', 0))
                               for row in self.profiles])
        self.upper = np.array([
            float('inf') if row.get('upper_degree') is None else row['upper_degree']
            for row in self.profiles])
        self.periods = np.array([row.get('period', 0) for row in self.profiles])
        self.fixed = self.automaton['fixed_profiles']

    @lru_cache(maxsize=2)
    def exponential_arrays(self, tau):
        outside = np.where(self.active, 0.0, np.exp(2 * tau * self.exponents))
        inside = np.where(self.active, np.exp(2 * tau * self.pattern_degrees[None, :]), 0.0)
        return outside, inside

    def evaluate(self, layer_probabilities, support_probabilities, tau):
        """Return minimum log cost, selected profile, and zero divergence.

        Rows are independent choices of the Boolean support/physical measure.
        The baseline is always available. The optimum along each infinite
        arithmetic progression needs only its two nearest allowed integers.
        """
        q = np.asarray(layer_probabilities, dtype=float)
        p = np.asarray(support_probabilities, dtype=float)
        assert q.ndim == 2 and q.shape[1] == len(self.pattern_degrees)
        assert p.shape == (len(q),) and tau > 0
        outside, inside = self.exponential_arrays(float(tau))
        outside_values, inside_values = q @ outside.T, q @ inside.T
        with np.errstate(divide='ignore', invalid='ignore'):
            log_a, log_b = np.log(outside_values), np.log(inside_values)
            log_p, log_not_p = np.log(p)[:, None], np.log1p(-p)[:, None]
            optimum = (log_b + log_not_p - log_a - log_p
                       - 2 * tau * self.offsets[None, :]) / (4 * tau)
        positions = np.maximum(optimum, self.lower[None, :])
        positions[:, self.continuous] = np.minimum(
            positions[:, self.continuous], self.upper[None, self.continuous])
        periodic = (~self.continuous) & (self.periods > 0)
        with np.errstate(invalid='ignore'):
            positions[:, periodic] = self.lower[None, periodic] + self.periods[None, periodic] * np.maximum(
                0, np.floor((positions[:, periodic] - self.lower[None, periodic])
                            / self.periods[None, periodic]))
        startup = (~self.continuous) & (self.periods == 0)
        positions[:, startup] = self.lower[None, startup]
        # Underflowed absent groups cannot improve a strictly positive measure.
        positions = np.where(np.isfinite(positions), positions, self.lower[None, :])

        def costs(divergences):
            theta = 2 * tau * divergences
            return (np.logaddexp(log_a, log_b - theta)
                    + np.logaddexp(log_not_p, log_p + theta + 2 * tau * self.offsets[None, :]))

        values = costs(positions)
        next_positions = positions.copy()
        next_positions[:, periodic] += self.periods[None, periodic]
        next_values = costs(next_positions)
        improve = next_values < values
        values = np.where(improve, next_values, values)
        positions = np.where(improve, next_positions, positions)
        indices = np.argmin(values, axis=1)
        minima = values[np.arange(len(q)), indices]
        chosen_positions = positions[np.arange(len(q)), indices]
        chosen = indices.copy()
        no_gain = minima >= -1e-12
        minima[no_gain] = 0
        chosen[no_gain] = -1
        chosen_positions[no_gain] = 0

        for fixed_index, row in enumerate(self.fixed):
            zero = q @ np.exp(-2 * tau * np.array(row['excesses']))
            infinity = 1 - p + p * np.exp(2 * tau * row['degree'])
            relative = np.log(zero) + np.log(infinity)
            better = relative < minima - 1e-12
            minima[better] = relative[better]
            chosen[better] = len(self.profiles) + fixed_index
            chosen_positions[better] = row['degree']
        return minima, chosen, chosen_positions

    def describe(self, index, divergence):
        if index < 0:
            return {'kind': 'baseline'}
        if index >= len(self.profiles):
            return {'kind': 'fixed', **self.fixed[index - len(self.profiles)]}
        row = self.profiles[index]
        return {**row, 'chosen_divergence': float(divergence),
                'chosen_degree': float(divergence + row.get('degree_offset', 0))}
