#!/usr/bin/env python3
"""Complete weighted block profiles for literal cyclic Pauli tensors.

The blocks follow from exact polynomial zeros in both full cuts. Arbitrary
whole-symbol heights retain those zeros; no Gaussian profile formula or
charge-homogeneous height restriction is used. Global nested minors merge
the complete block profiles, including their attaining rows and columns.
"""
from __future__ import annotations

import heapq
import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from cyclic_pauli_crossings import charge_labels
from native_matrix_smith_f4 import MatrixOracle
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure


class PauliBlockOracle:
    def __init__(self, planes, metadata, precision=512):
        self.q, self.n = metadata['q'], metadata['q']**2
        self.stride, self.blocks = 4*(self.q-1), []
        self.queries = 0
        q, n = self.q, self.n
        for infinity in (True, False):
            matrices = ([p.reshape(q,q,q,q).transpose(0,3,2,1).reshape(n,n) for p in planes]
                        if infinity else planes)
            labels = charge_labels(metadata['k'], metadata['central_rows'], 1 if infinity else -1)
            if any(np.any(p[labels[:,None] != labels[None,:]]) for p in matrices):
                raise ArithmeticError('An exact off-charge entry is nonzero.')
            blocks = []
            for charge in sorted(set(map(int, labels))):
                indices = np.flatnonzero(labels == charge)
                matrix = [[tuple(int(p[i,j]) for p in matrices) for j in indices] for i in indices]
                blocks.append({'charge': charge, 'indices': indices, 'matrix': matrix,
                               'oracle': MatrixOracle(matrix, len(indices), precision=precision)})
            self.blocks.append(blocks)

    def prefixes(self, parameters, cut=1):
        parameters = np.asarray(parameters, dtype=float)
        if cut != 1 or parameters.shape != (2*self.stride+1,) or not np.all(np.isfinite(parameters)) or parameters[-1] <= 0:
            raise ValueError('Supply complete finite infinity/zero heights and positive scale, cut1.')
        self.queries += 1
        q, n, result = self.q, self.n, {'orders': [], 'base_prefix': [], 'rows': [], 'cols': []}
        pair_a, pair_b = np.arange(n)//q, np.arange(n)%q
        for side, blocks in enumerate(self.blocks):
            h = [np.r_[0., parameters[side*self.stride+leg*(q-1):side*self.stride+(leg+1)*(q-1)]]
                 for leg in range(4)]
            rw = h[0][pair_a]+h[2 if side == 0 else 1][pair_b]
            cw = h[3][pair_a]+h[1 if side == 0 else 2][pair_b]
            steps = []
            for block in blocks:
                indices = block['indices']
                answer = block['oracle'].profile(rw[indices], cw[indices], scale=parameters[-1], infinity=side == 0)
                previous, sequence = 0, []
                for order, base, row, col in zip(answer['orders'], answer['base_prefix'], answer['rows'], answer['cols']):
                    i, j = int(indices[row]), int(indices[col])
                    # Internal block columns use (V,Y) or (V,U); the common
                    # asymmetric feature oracle expects (Y,V) or (U,V).
                    native_col = (j % q)*q+j//q
                    sequence.append((order, base-previous, i, native_col))
                    previous = base
                steps.append(sequence)
            heap = [(sequence[0][0], i, 0) for i, sequence in enumerate(steps)]
            heapq.heapify(heap)
            cumulative, count = 0, 0
            while heap:
                _, i, j = heapq.heappop(heap)
                order, base, row, col = steps[i][j]
                cumulative += base
                result['orders'].append(order)
                result['base_prefix'].append(cumulative)
                result['rows'].append(row)
                result['cols'].append(col)
                count += 1
                if j+1 < len(steps[i]):
                    heapq.heappush(heap, (steps[i][j+1][0], i, j+1))
            if count != n:
                raise ArithmeticError('A complete promised block cut was not resolved.')
        return result


class PauliAsymmetricPressure(AsymmetricPressure):
    def __init__(self, planes, metadata, precision=512):
        self.seed, self.cut = {'q_x': metadata['q'], 'q_u': metadata['q']}, 1
        self.qx = self.qy = metadata['q']
        self.n = self.qx*self.qy
        self.stride = 4*(self.qx-1)
        self.alphabets = (self.qx,)*4
        self.offsets = np.arange(4)*(self.qx-1)
        self.parameters = 2*self.stride+2
        self.oracle = PauliBlockOracle(planes, metadata, precision)


class ScaledPauliPressure(ScaledAsymmetricPressure):
    def __init__(self, planes, metadata, precision=512):
        self.base = PauliAsymmetricPressure(planes, metadata, precision)
        self.parameters = self.base.parameters
