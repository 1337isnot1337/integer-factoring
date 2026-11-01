#!/usr/bin/env python3
"""Independent zero/infinity letter gauges with both physical costs charged.

This is a discovery oracle. Its finite real-weight evaluations do not by
themselves certify a rate improvement or a global optimum.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from coherent_weighted_gauges import SharedGaugeOracle


class AsymmetricGaugeOracle(SharedGaugeOracle):
    def __init__(self, embedding, retraction, modes):
        super().__init__(embedding, retraction, modes)
        rows, width = self.cone.shape
        self.asymmetric_cone = np.zeros((2 * rows, 2 * width - 1))
        self.asymmetric_cone[:rows, :width] = self.cone
        self.asymmetric_cone[rows:, 0] = self.cone[:, 0]
        self.asymmetric_cone[rows:, width:] = self.cone[:, 1:]

    def evaluate(self, parameters, *, enforce=True, details=False):
        parameters = np.asarray(parameters, dtype=float)
        if parameters.shape != (2 * self.q - 1,) or not np.all(np.isfinite(parameters)):
            raise ValueError('Supply u and two sets of q minus one finite potentials.')
        u = float(parameters[0])
        if u <= 0:
            raise ValueError('The common profile scale must be positive.')
        slacks = self.asymmetric_cone @ parameters
        minimum_slack = float(slacks.min()) if len(slacks) else 0.
        if enforce and minimum_slack < -1e-9:
            raise ValueError('At least one complete-cut regularity cone failed.')
        potentials = [np.zeros(self.q), np.zeros(self.q)]
        potentials[0][self.free] = parameters[1:self.q]
        potentials[1][self.free] = parameters[self.q:]
        cuts = [self._cut(i, u, potentials[i]) for i in (0, 1)]
        target = sum(cut['log_partition'] for cut in cuts)
        energies = u * self.heights + potentials[0] + potentials[1]

        def equation(rho):
            return target - 4 * rho * float(logsumexp(-energies / rho))

        rho = float(brentq(equation, 1e-7, 2., xtol=1e-13))
        log_physical = float(logsumexp(-energies / rho))
        physical = np.exp(-energies / rho - log_physical)
        entropy = float(log_physical + physical @ energies / rho)
        if entropy <= 0:
            raise ArithmeticError('The pressure root has zero physical entropy.')
        gradient_u = (sum(cut['gradient_u'] for cut in cuts)
                      + 4 * float(physical @ self.heights)) / (4 * entropy)
        gradients_f = [(cut['gradient_f'] + 4 * physical) / (4 * entropy) for cut in cuts]
        gradient = np.concatenate(([gradient_u], gradients_f[0][self.free], gradients_f[1][self.free]))
        result = {'rho': rho, 'gradient': gradient, 'minimum_cone_slack': minimum_slack}
        if details:
            result.update(scale_t=u, potentials=[f.tolist() for f in potentials],
                          physical_energies=energies.tolist(), physical_type=physical.tolist(),
                          physical_entropy=entropy, zero_orders=cuts[0]['orders'],
                          infinity_orders=cuts[1]['orders'],
                          block_precisions=[[block['precision'] for block in blocks]
                                            for blocks in self.blocks])
        return result
