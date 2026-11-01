#!/usr/bin/env python3
"""Shared letter gauges for homogeneous coherent physical retractions.

The matrix coefficients and the prescribed generic ranks are exact. Real
gauge weights are exploratory. A claimed rate improvement still needs an
independent rational gauge and pressure certificate.

For every pair of induced deletion edges j->k and l->m, on modes a,b,
the sufficient regularity cone is
    f[j]-f[k]+f[l]-f[m]+2*u*(a&b) >= 0.
The edge union covers both the embedding and the retraction spaces.
"""
from __future__ import annotations

from itertools import combinations_with_replacement

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from coherent_hot_compounds import support
from coherent_hot_compounds_fast import profiles as charge_profiles
from coherent_hot_mixed import profiles as parity_profiles, validate_maps
from native_matrix_smith_f4 import MatrixOracle


def deletion_edges(embedding, retraction, modes):
    """Return all nonzero induced bit-deletion entries in either map space."""
    embedding, retraction, _, _ = validate_maps(embedding, retraction, modes)
    edges = set()
    for basis, dual in ((embedding, retraction), (retraction, embedding)):
        for source, vector in enumerate(basis):
            for bit in range(modes):
                deleted = 0
                for word in support(vector):
                    if word >> bit & 1:
                        deleted ^= 1 << (word ^ (1 << bit))
                for target, functional in enumerate(dual):
                    if (deleted & functional).bit_count() % 2:
                        edges.add((bit, source, target))
    return sorted(edges)


def gauge_cone(q, vacuum, edges):
    """Rows act on (u, all non-vacuum f_j); nonnegative is sufficient."""
    indices = [j for j in range(q) if j != vacuum]
    position = {j: index + 1 for index, j in enumerate(indices)}
    strongest = {}
    for (a, source, target), (b, other, lower) in combinations_with_replacement(edges, 2):
        row = [0] * (q - 1)
        for label, sign in ((source, 1), (target, -1), (other, 1), (lower, -1)):
            if label != vacuum:
                row[position[label] - 1] += sign
        key = tuple(row)
        value = 2 * (a & b)
        strongest[key] = min(value, strongest.get(key, value))
    return np.asarray([(value,) + row for row, value in sorted(strongest.items())], dtype=float)


class SharedGaugeOracle:
    """Complete weighted cut profiles with derivatives in a pivot region."""

    def __init__(self, embedding, retraction, modes):
        charge_homogeneous = all(len({word.bit_count() for word in support(p | l)}) == 1
                                 for p, l in zip(embedding, retraction))
        constructor = charge_profiles if charge_homogeneous else parity_profiles
        self.natural = constructor(embedding, retraction, modes, keep_blocks=True)
        self.natural['weighted_block_decomposition'] = 'charge' if charge_homogeneous else 'parity'
        self.q = self.natural['q']
        self.heights = np.asarray(self.natural['heights'], dtype=float)
        vacuum = np.flatnonzero(self.heights == 0)
        if len(vacuum) != 1:
            raise ValueError('The gauge normalization needs a unique vacuum.')
        self.vacuum = int(vacuum[0])
        self.free = [j for j in range(self.q) if j != self.vacuum]
        self.edges = deletion_edges(embedding, retraction, modes)
        self.cone = gauge_cone(self.q, self.vacuum, self.edges)
        self.blocks = []
        for cut in self.natural['cuts']:
            blocks = []
            for block in cut['blocks']:
                matrix = [[(value, 0) for value in row] for row in block.pop('matrix')]
                blocks.append({'matrix': matrix, 'pairs': np.asarray(block['pairs'], dtype=int),
                               'precision': block['precision_t'],
                               'oracle': MatrixOracle(matrix, block['size'],
                                                      precision=block['precision_t'])})
            self.blocks.append(blocks)

    def _cut(self, index, u, f):
        orders, bases, labels = [], [], []
        for block in self.blocks[index]:
            pairs = block['pairs']
            weights = f[pairs].sum(axis=1)
            while True:
                try:
                    result = block['oracle'].profile(weights, weights, scale=u)
                    break
                except ArithmeticError:
                    block['precision'] *= 2
                    if block['precision'] > 8192:
                        raise
                    block['oracle'] = MatrixOracle(block['matrix'], len(pairs),
                                                   precision=block['precision'])
            orders.extend(result['orders'])
            bases.extend(np.diff([0] + result['base_prefix']))
            labels.extend(np.concatenate((pairs[result['rows']], pairs[result['cols']]), axis=1))
        orders = np.asarray(orders)
        bases = np.asarray(bases)
        labels = np.asarray(labels)
        if len(orders) != self.q * self.q:
            raise ArithmeticError('A complete cut did not contain q squared orders.')
        expected = u * bases + f[labels].sum(axis=1)
        if not np.allclose(orders, expected, rtol=1e-10, atol=1e-9):
            raise ArithmeticError('Weighted profile derivative reconstruction failed.')
        log_partition = float(logsumexp(-orders))
        probability = np.exp(-orders - log_partition)
        gradient = -np.bincount(labels.ravel(), np.repeat(probability, 4), minlength=self.q)
        return {'orders': orders.tolist(), 'base_orders': bases.tolist(),
                'labels': labels.tolist(), 'log_partition': log_partition,
                'gradient_u': -float(probability @ bases), 'gradient_f': gradient}

    def evaluate(self, parameters, *, enforce=True, details=False):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (self.q,) or not np.all(np.isfinite(parameters)):
            raise ValueError('Supply u and q minus one finite letter potentials.')
        u = float(parameters[0])
        if u <= 0:
            raise ValueError('The profile scale must be positive.')
        slacks = self.cone @ parameters
        minimum_slack = float(slacks.min()) if len(slacks) else 0.
        if enforce and minimum_slack < -1e-9:
            raise ValueError('The shared regularity cone failed.')
        f = np.zeros(self.q)
        f[self.free] = parameters[1:]
        cuts = [self._cut(index, u, f) for index in (0, 1)]
        target = sum(cut['log_partition'] for cut in cuts)
        energies = u * self.heights + 2 * f

        def equation(rho):
            return target - 4 * rho * float(logsumexp(-energies / rho))

        rho = float(brentq(equation, 1e-7, 2., xtol=1e-13))
        log_physical = float(logsumexp(-energies / rho))
        physical = np.exp(-energies / rho - log_physical)
        entropy = float(log_physical + physical @ energies / rho)
        if entropy <= 0:
            raise ArithmeticError('The pressure root has zero physical entropy.')
        gradient_f = (sum(cut['gradient_f'] for cut in cuts) + 8 * physical) / (4 * entropy)
        gradient_u = (sum(cut['gradient_u'] for cut in cuts)
                      + 4 * float(physical @ self.heights)) / (4 * entropy)
        gradient = np.asarray([gradient_u] + gradient_f[self.free].tolist())
        result = {'rho': rho, 'gradient': gradient, 'minimum_cone_slack': minimum_slack}
        if details:
            result.update(scale_t=u, potentials=f.tolist(), physical_energies=energies.tolist(),
                          physical_type=physical.tolist(), physical_entropy=entropy,
                          zero_orders=cuts[0]['orders'], infinity_orders=cuts[1]['orders'],
                          block_precisions=[[block['precision'] for block in blocks]
                                            for blocks in self.blocks])
        return result
