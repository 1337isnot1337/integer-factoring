#!/usr/bin/env python3
"""Literal affine Boolean relabelings of symmetric polynomial filters.

An affine relabeling changes the coefficient polynomials f_k.  This module
re-evaluates the complete finite-difference transform; it never permutes an
old beta profile. Binary-extension coefficients use their additive planes.
"""

from hashlib import sha256
import json

from symmetric_binary_filter_synthesis import submasks


def affine_relabelings(modes=4):
    identity = tuple(1 << bit for bit in range(modes))
    matrices = [('identity', identity)]
    for target in range(modes):
        for source in range(modes):
            if target != source:
                rows = list(identity)
                rows[target] ^= 1 << source
                matrices.append((f'shear_{target}_{source}', tuple(rows)))
    matrices.append(('gray', tuple((1 << bit) ^ (1 << (bit + 1))
                                  if bit + 1 < modes else 1 << bit for bit in range(modes))))
    matrices.append(('inverse_gray', tuple(sum(1 << other for other in range(bit, modes))
                                          for bit in range(modes))))
    result = []
    for name, rows in matrices:
        labels = tuple(sum(((state & row).bit_count() % 2) << bit
                           for bit, row in enumerate(rows)) for state in range(1 << modes))
        assert sorted(labels) == list(range(1 << modes))
        for translation in range(1 << modes):
            mapping = tuple(label ^ translation for label in labels)
            result.append({'matrix': name, 'linear_rows': rows, 'translation': translation,
                           'labels': mapping, 'weights': tuple(label.bit_count() for label in mapping)})
    return tuple(result)


def polynomial_planes(polynomial, bits):
    return tuple(sum(1 << degree for degree, coefficient in polynomial.items()
                     if (coefficient >> bit) & 1) for bit in range(bits))


def planes_polynomial(planes):
    result = {}
    for bit, plane in enumerate(planes):
        while plane:
            low = plane & -plane
            degree = low.bit_length() - 1
            result[degree] = result.get(degree, 0) ^ (1 << bit)
            plane ^= low
    return result


class LiteralRelabeler:
    def __init__(self, witness):
        self.modes, self.m = witness['modes'], witness['m']
        self.bits = witness['coefficient_field_order'].bit_length() - 1
        self.values = tuple(polynomial_planes(value, self.bits) for value in witness['f'])
        self.degrees = witness['actual_degrees']
        self.endpoints = witness['endpoint_values']
        self.multiplied = {}
        for weight, planes in enumerate(self.values):
            for power in range(self.modes + 1):
                result = [0] * self.bits
                for shift in submasks(power):
                    for bit, plane in enumerate(planes):
                        result[bit] ^= plane << (self.m * shift)
                self.multiplied[weight, power] = tuple(result)

    def evaluate(self, weights, literal=False):
        assert len(weights) == 1 << self.modes
        orders = []
        beta = []
        for state in range(1 << self.modes):
            planes = [0] * self.bits
            for source in submasks(state):
                value = self.multiplied[weights[source], state.bit_count() - source.bit_count()]
                for bit in range(self.bits):
                    planes[bit] ^= value[bit]
            union = 0
            for plane in planes:
                union |= plane
            assert union
            orders.append((union & -union).bit_length() - 1)
            if literal:
                beta.append(tuple(planes))
        degrees = tuple(self.degrees[weight] for weight in weights)
        common = min(orders)
        assert all(degree >= common for degree in degrees)
        result = {'beta_orders': tuple(order - common for order in orders),
                  'degrees': tuple(degree - common for degree in degrees),
                  'endpoint_values': tuple(self.endpoints[weight] for weight in weights),
                  'removed_common_monomial_order': common}
        if literal:
            values = tuple(tuple(plane >> common for plane in self.values[weight]) for weight in weights)
            beta = tuple(tuple(plane >> common for plane in value) for value in beta)
            # A common x-power divides all f coordinates because the
            # transform and inverse are regular and unimodular at zero.
            assert all(not (plane & ((1 << common) - 1))
                       for weight in weights for plane in self.values[weight])
            result['f'] = tuple(planes_polynomial(value) for value in values)
            result['beta'] = tuple(planes_polynomial(value) for value in beta)
            payload = [[sorted(value.items()) for value in result[name]] for name in ('f', 'beta')]
            result['literal_polynomial_sha256'] = sha256(json.dumps(payload, separators=(',', ':')).encode()).hexdigest()
        return result
