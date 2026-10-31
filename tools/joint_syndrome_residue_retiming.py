#!/usr/bin/env python3
"""Independent residue colors and convex continuous exponent retiming."""

from functools import lru_cache
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize
from scipy.special import logsumexp

from joint_syndrome_filter_synthesis import transform


def decompose_residues(polynomials, rank, column_types, spacing):
    polynomials = tuple({int(degree): int(coefficient) for degree, coefficient in row.items()
                         if coefficient} for row in polynomials)
    residues = sorted({degree % spacing for row in polynomials for degree in row})
    width = max(coefficient.bit_length() for row in polynomials for coefficient in row.values())
    components = []
    for residue in residues:
        f = tuple({degree // spacing: coefficient for degree, coefficient in row.items()
                   if degree % spacing == residue} for row in polynomials)
        delta = transform(f, rank, tuple(column_types), 1)
        endpoints = []
        for row in f:
            value = 0
            for coefficient in row.values():
                value ^= coefficient
            endpoints.append(value)
        components.append({'source_residue': residue, 'f': f, 'delta': delta,
                           'endpoint_values': tuple(endpoints),
                           'delta_orders': tuple(min(row) if row else None for row in delta),
                           'degrees': tuple(max(row) if row else None for row in f)})
    assert all(any(component['endpoint_values'][b] for component in components)
               for b in range(1 << rank))
    return {'rank': rank, 'column_types': tuple(column_types), 'source_spacing': spacing,
            'coefficient_bit_width': width, 'components': tuple(components),
            'source_positions': tuple(residue / spacing for residue in residues),
            'individually_admissible_components': tuple(i for i, row in enumerate(components)
                                                       if all(row['endpoint_values']))}


def decompose_binary_components(polynomials, rank, column_types, spacing):
    original = decompose_residues(polynomials, rank, column_types, spacing)
    components = []
    for component in original['components']:
        for plane in range(original['coefficient_bit_width']):
            f = tuple({degree: 1 for degree, coefficient in row.items() if coefficient & (1 << plane)}
                      for row in component['f'])
            if not any(f):
                continue
            delta = transform(f, rank, tuple(column_types), 1)
            endpoints = tuple(len(row) % 2 for row in f)
            components.append({'source_residue': component['source_residue'],
                               'source_coefficient_plane': plane, 'f': f, 'delta': delta,
                               'endpoint_values': endpoints,
                               'delta_orders': tuple(min(row) if row else None for row in delta),
                               'degrees': tuple(max(row) if row else None for row in f)})
    return {'rank': rank, 'column_types': tuple(column_types), 'source_spacing': spacing,
            'coefficient_bit_width': 1, 'components': tuple(components),
            'source_positions': tuple(row['source_residue'] / spacing for row in components),
            'individually_admissible_components': tuple(i for i, row in enumerate(components)
                                                       if all(row['endpoint_values']))}


def endpoint_cover(decomposition):
    """Choose at most one binary component for each as-yet-uncovered endpoint."""
    components = decomposition['components']
    masks = [sum(1 << b for b, value in enumerate(row['endpoint_values']) if value)
             for row in components]
    full = (1 << (1 << decomposition['rank'])) - 1
    covered, selected = 0, []
    while covered != full:
        chosen = max(range(len(masks)), key=lambda i: ((masks[i] & ~covered).bit_count(), -i))
        assert masks[chosen] & ~covered
        selected.append(chosen)
        covered |= masks[chosen]
    # Delete redundant choices while retaining the complete endpoint cover.
    for index in selected[:]:
        union = 0
        for other in selected:
            if other != index:
                union |= masks[other]
        if union == full:
            selected.remove(index)
    assert len(selected) <= 1 << decomposition['rank']
    reduced = dict(decomposition)
    reduced['components'] = tuple(components[index] for index in selected)
    reduced['source_positions'] = tuple(decomposition['source_positions'][index] for index in selected)
    reduced['individually_admissible_components'] = tuple(i for i, row in enumerate(reduced['components'])
                                                        if all(row['endpoint_values']))
    reduced['source_component_indices'] = tuple(selected)
    return reduced


def retime_literal(decomposition, numerators, spacing):
    components = decomposition['components']
    if len(numerators) != len(components) or any(value < 0 or value > spacing for value in numerators):
        raise ValueError('Use one offset in [0, spacing] for each residue component.')
    rank, types, width = (decomposition['rank'], decomposition['column_types'],
                          decomposition['coefficient_bit_width'])
    f = [dict() for _ in range(1 << rank)]
    for index, (component, offset) in enumerate(zip(components, numerators)):
        for syndrome, row in enumerate(component['f']):
            for degree, coefficient in row.items():
                target = spacing * degree + offset
                f[syndrome][target] = f[syndrome].get(target, 0) ^ (coefficient << (width * index))
    f = tuple({degree: coefficient for degree, coefficient in row.items() if coefficient} for row in f)
    delta = transform(f, rank, types, spacing)
    orders, degrees = tuple(map(min, delta)), tuple(map(max, f))
    expected_orders = tuple(min(spacing * component['delta_orders'][state] + offset
                                for component, offset in zip(components, numerators)
                                if component['delta_orders'][state] is not None)
                            for state in range(len(delta)))
    expected_degrees = tuple(max(spacing * component['degrees'][syndrome] + offset
                                 for component, offset in zip(components, numerators)
                                 if component['degrees'][syndrome] is not None)
                             for syndrome in range(1 << rank))
    assert orders == expected_orders and degrees == expected_degrees
    endpoints = []
    for row in f:
        value = 0
        for coefficient in row.values():
            value ^= coefficient
        endpoints.append(value)
    assert all(endpoints)
    return {'rank': rank, 'column_types': types, 'spacing': spacing,
            'f': f, 'delta': delta, 'actual_orders': orders, 'actual_degrees': degrees,
            'endpoint_values': tuple(endpoints), 'coefficient_field_degree': width * len(components),
            'retimed_numerators': tuple(numerators)}


class ConvexResiduePressure:
    """Search the exact convex profile family; floating minima are proposals."""
    def __init__(self, decomposition, degree_caps=None, order_targets=None):
        self.decomposition = decomposition
        self.rank = decomposition['rank']
        components = decomposition['components']
        self.count = len(components)
        self.orders = np.array([[float('inf') if value is None else value
                                 for value in component['delta_orders']] for component in components]).T
        self.degrees = np.array([[float('-inf') if value is None else value
                                  for value in component['degrees']] for component in components]).T
        self.zero_count, self.infinity_count = self.orders.shape[0], self.degrees.shape[0]
        self.position_lower, self.position_upper = np.zeros(self.count), np.ones(self.count)
        if degree_caps is not None:
            degree_caps = np.asarray(degree_caps, dtype=float)
            if degree_caps.shape != (self.infinity_count,):
                raise ValueError('Use one normalized degree cap for each coefficient.')
            self.position_upper = np.minimum(self.position_upper,
                                               np.min(degree_caps[:,None]-self.degrees,axis=0))
        if order_targets is not None:
            order_targets = np.asarray(order_targets, dtype=float)
            if order_targets.shape != (self.zero_count,):
                raise ValueError('Use one normalized delta target for each output.')
            self.position_lower = np.maximum(self.position_lower,
                                               np.max(order_targets[:,None]-self.orders,axis=0))
        if np.any(self.position_lower>self.position_upper):
            raise ValueError('The requested component offset box is empty.')
        dimension = self.count + self.zero_count + self.infinity_count
        rows, lower, upper = [], [], []
        for state in range(self.zero_count):
            for component in range(self.count):
                if np.isfinite(self.orders[state, component]):
                    row = np.zeros(dimension)
                    row[self.count + state], row[component] = 1, -1
                    rows.append(row); lower.append(-np.inf); upper.append(self.orders[state, component])
        for state in range(self.infinity_count):
            for component in range(self.count):
                if np.isfinite(self.degrees[state, component]):
                    row = np.zeros(dimension)
                    row[self.count + self.zero_count + state], row[component] = 1, -1
                    rows.append(row); lower.append(self.degrees[state, component]); upper.append(np.inf)
        self.constraint = LinearConstraint(np.array(rows), np.array(lower), np.array(upper))
        self.bounds = Bounds(np.r_[self.position_lower, np.full(dimension-self.count,-np.inf)],
                             np.r_[self.position_upper, np.full(dimension-self.count,np.inf)])
        if self.count <= 16:
            zero_sets = np.array([sum(1 << j for j, value in enumerate(row) if value == row.min())
                                  for row in self.orders], dtype=np.int64)
            infinity_sets = np.array([sum(1 << j for j, value in enumerate(row) if value == row.max())
                                      for row in self.degrees], dtype=np.int64)
            subsets = np.arange(1 << self.count, dtype=np.int64)
            self.zero_contained = (zero_sets[None,:] & ~subsets[:,None]) == 0
            self.infinity_complement_contained = (infinity_sets[None,:] & subsets[:,None]) == 0
        else:
            self.zero_contained = self.infinity_complement_contained = None

    def profiles(self, positions):
        positions = np.asarray(positions)
        return np.min(self.orders + positions, axis=1), np.max(self.degrees + positions, axis=1)

    def cost(self, positions, zero_probabilities, infinity_probabilities, tau):
        orders, degrees = self.profiles(positions)
        with np.errstate(divide='ignore'):
            return float(logsumexp(np.log(zero_probabilities)-2*tau*(orders-self.rank))
                         +logsumexp(np.log(infinity_probabilities)+2*tau*degrees))

    def center_bound(self, zero_probabilities, infinity_probabilities, tau):
        """Convex transport lower bound over the complete offset cube (floating)."""
        with np.errstate(divide='ignore'):
            z = np.log(zero_probabilities)-2*tau*(self.orders.min(axis=1)-self.rank)
            i = np.log(infinity_probabilities)+2*tau*self.degrees.max(axis=1)
        lz, li = logsumexp(z), logsumexp(i)
        if self.zero_contained is None:
            deficit = 1.0
        else:
            deficit = max(0.0, float(np.max(
                self.zero_contained @ np.exp(z-lz)
                +self.infinity_complement_contained @ np.exp(i-li)-1)))
        value = float(lz+li)
        return {'center_cost': value, 'transport_deficit': deficit,
                'relative_log_cost_lower_bound': value-2*tau*deficit}

    def optimize(self, zero_probabilities, infinity_probabilities, tau, initial=None):
        q, p = np.asarray(zero_probabilities), np.asarray(infinity_probabilities)
        if self.count == 1:
            positions = self.position_lower.copy()
            return {'positions': positions.tolist(), 'relative_log_cost': self.cost(positions,q,p,tau),
                    'optimizer_success': True, 'iterations': 0}
        positions = np.array(self.decomposition['source_positions'] if initial is None else initial)
        positions = np.clip(positions,self.position_lower,self.position_upper)
        orders, degrees = self.profiles(positions)
        start = np.r_[positions, orders, degrees]
        with np.errstate(divide='ignore'):
            log_q, log_p = np.log(q), np.log(p)

        def objective(value):
            z = log_q-2*tau*(value[self.count:self.count+self.zero_count]-self.rank)
            i = log_p+2*tau*value[self.count+self.zero_count:]
            lz, li = logsumexp(z), logsumexp(i)
            gradient = np.r_[np.zeros(self.count),-2*tau*np.exp(z-lz),2*tau*np.exp(i-li)]
            return float(lz+li), gradient

        result = minimize(objective, start, method='SLSQP', jac=True, bounds=self.bounds,
                          constraints=[self.constraint], options={'ftol':1e-11,'maxiter':160})
        positions = np.clip(result.x[:self.count],self.position_lower,self.position_upper)
        value = self.cost(positions,q,p,tau)
        return {'positions': positions.tolist(), 'relative_log_cost': value,
                'optimizer_success': bool(result.success), 'iterations': int(result.nit),
                'message': str(result.message)}
