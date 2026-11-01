#!/usr/bin/env python3
"""Complete natural profiles for a physical constraint whose support grows.

Forbid joint occupation of physical modes 0 and M-1, and retain every
middle mode.  The alphabet has 3*2**(M-2) letters.  Regular concentration
leaves an eight-mode core, of ambient order 256 and generic rank 144,
with Gaussian spectators.  All polynomial coefficients are constructed
over F2 before the native Smith calculation.  Pressure searches are
numerical and are not strict certificates or global minima proofs.
"""

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from check_cancellation_module_automata import pmul
from coherent_full_tail_core import binary_rank, projector
from coherent_full_tail_core import pressure as gaussian_comparison
from native_matrix_smith_f4 import MatrixOracle


def images(M):
    """Conjugated physical operators in x0,x1,xp,xe,y0,y1,yp,ye order.

    Here p=M-2 and e=M-1.  Each group lists the creator/annihilator
    pair for physical mode 0 followed by that for physical mode M-1.
    The within-row and within-column Gaussian factors descend through
    this deletion-closed coordinate code, so only the compound is used.
    """
    assert isinstance(M, int) and M >= 4 and not M & (M - 1)
    c = 3  # 1+t
    first = sum(1 << (2 * j) for j in range(M - 1))
    second = sum(1 << (2 * j) for j in range(M - 2))

    def vector(terms):
        result = [0] * 8
        for i, value in terms:
            result[i] ^= value
        return result

    groups = [
        [
            (vector([(0, 1), (4, c), (5, c << 2)]), vector([(0, 1)])),
            (vector([(3, 1), (4, c), (6, c), (7, c)]),
             vector([(1, 1 << (2 * M - 4)), (3, 1)])),
        ],
        [
            (vector([(4, 1)]),
             vector([(4, 1), (0, c), (1, pmul(c, first)), (3, c)])),
            (vector([(7, 1)]),
             vector([(5, 1 << (2 * M - 4)), (7, 1),
                     (0, c << (2 * M - 2)), (1, c << (2 * M - 4)),
                     (2, pmul(c, second) << 2), (3, c)])),
        ],
        [
            (vector([(0, 1)]),
             vector([(0, 1), (4, c), (5, pmul(c, first) << 2),
                     (7, c << (2 * M - 2))])),
            (vector([(3, 1)]),
             vector([(1, 1), (3, 1), (4, c), (5, c),
                     (6, pmul(c, second)), (7, c)])),
        ],
        [
            (vector([(4, 1), (0, c), (1, c)]), vector([(4, 1)])),
            (vector([(7, 1), (0, c << (2 * M - 2)), (2, c << 2), (3, c)]),
             vector([(5, 1), (7, 1)])),
        ],
    ]
    for side in (groups[:2], groups[2:]):
        operators = side[0] + side[1]
        for i, (creator, _) in enumerate(operators):
            for j, (_, annihilator) in enumerate(operators):
                value = 0
                for a, b in zip(creator, annihilator):
                    value ^= pmul(a, b)
                assert value == int(i == j)
    retained = [0, 1, M - 2, M - 1]
    diagonal = [2 * j + 1 for j in retained]
    diagonal += [2 * M + 2 * j + 1 for j in retained]
    return groups, diagonal


def core_matrix(M):
    groups, diagonal = images(M)
    matrix = [[0] * 256 for _ in range(256)]
    for state in range(256):
        exponent = sum(value for j, value in enumerate(diagonal) if state >> j & 1)
        matrix[state][state] = 1 << exponent
    for group, operators in enumerate(groups):
        matrix = projector(matrix, operators, group >= 2)
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    # Sum the maximum degrees of the diagonal and all sixteen operators.
    assert degree <= 28 * M
    return matrix, degree


def core(M):
    matrix, degree = core_matrix(M)
    # (1+2z)^2 (1+z)^4; the generic-rank proof uses the original endpoint,
    # not evaluation of the concentration matrix at t=1.
    ranks = [1, 8, 26, 44, 41, 20, 4, 0, 0]
    result, blocks = [], []
    for charge, rank in enumerate(ranks):
        indices = [i for i in range(256) if i.bit_count() == charge]
        assert all(not matrix[i][j] for i in indices
                   for j in range(256) if j.bit_count() != charge)
        endpoint = [sum((matrix[i][j].bit_count() % 2) << column
                        for column, j in enumerate(indices)) for i in indices]
        endpoint_rank = binary_rank(endpoint)
        assert endpoint_rank == rank  # Additional finite control only.
        if not rank:
            assert all(not matrix[i][j] for i in indices for j in indices)
            continue
        block = [[(matrix[i][j], 0) for j in indices] for i in indices]
        coefficient_bits = max(value[0].bit_length() for row in block for value in row)
        precision = 128
        while precision < coefficient_bits:
            precision *= 2
        attempts = []
        while precision <= 8192:
            try:
                orders = MatrixOracle(block, rank, precision=precision).profile()['orders']
                attempts.append([precision, True])
                break
            except ArithmeticError:
                attempts.append([precision, False])
                precision *= 2
        else:
            raise ArithmeticError('Complete active profile exceeded the precision budget.')
        assert len(orders) == rank and all(value == round(value) and value >= 0 for value in orders)
        orders = list(map(round, orders))
        result.extend(orders)
        blocks.append(dict(charge=charge, size=len(indices), rank=rank,
                           endpoint_rank=endpoint_rank,
                           coefficient_degree=coefficient_bits - 1,
                           precision=precision, precision_attempts=attempts, orders=orders))
    assert len(result) == 144
    return sorted(result), blocks, degree


def pressure(M, orders):
    """Charge the actual hot heights 0,1,2M-1 and all middle heights."""
    import numpy as np
    from scipy.optimize import brentq, minimize_scalar
    from scipy.special import logsumexp

    spectators = [2 * j + 1 for j in range(2, M - 2)]
    spectators += [2 * M + 2 * j + 1 for j in range(2, M - 2)]
    physical = np.asarray([2 * j + 1 for j in range(1, M - 1)])
    orders, spectators = np.asarray(orders), np.asarray(spectators)

    def rho(log_u):
        u = math.exp(log_u)
        crossing = float(logsumexp(-u * orders) + np.logaddexp(0, -u * spectators).sum())

        def equation(rate):
            alphabet = float(logsumexp([0, -u / rate, -(2 * M - 1) * u / rate])
                             + np.logaddexp(0, -u * physical / rate).sum())
            return 2 * crossing - 4 * rate * alphabet

        return brentq(equation, 0.2, 1.5, xtol=1e-13)

    grid = np.linspace(math.log(1e-6), math.log(2), 101)
    values = [rho(x) for x in grid]
    candidates = [(values[0], float(grid[0])), (values[-1], float(grid[-1]))]
    brackets = []
    for i in range(1, len(grid) - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            fit = minimize_scalar(rho, bounds=(grid[i - 1], grid[i + 1]), method='bounded',
                                  options={'xatol': 1e-10})
            candidates.append((float(fit.fun), float(fit.x)))
            brackets.append(dict(rho=float(fit.fun), u=math.exp(float(fit.x)),
                                 success=bool(fit.success)))
    rate, log_u = min(candidates)
    return dict(rho=rate, u=math.exp(log_u), u_range=[1e-6, 2], grid_points=101,
                local_minima=brackets,
                scope='Numerical grid and local refinements; no global or strict certificate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', default='4,8,16,32,64,128,256,512')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sources = (Path(__file__).name, 'coherent_full_tail_core.py',
               'check_cancellation_module_automata.py', 'native_matrix_smith_f4.py',
               'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    report = dict(schema=1, complete=False, factoring_exponent_improved=False,
                  scope=__doc__, records=[],
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in sources})
    for M in map(int, args.modes.split(',')):
        start = time.monotonic()
        orders, blocks, degree = core(M)
        record = dict(M=M, orders=orders, blocks=blocks,
                      physical_alphabet_dimension=3 * 2 ** (M - 2),
                      physical_forbidden_pair=[0, M - 1],
                      active_ambient_order=256, active_generic_rank=144,
                      actual_coefficient_degree=degree, coefficient_degree_bound=28 * M,
                      raw_tensor_degree_w_bound=8 * M * M,
                      physical_hot_heights=[0, 1, 2 * M - 1],
                      physical_tail_heights=[2 * j + 1 for j in range(1, M - 1)],
                      spectator_orders=[2 * j + 1 for j in range(2, M - 2)]
                      + [2 * M + 2 * j + 1 for j in range(2, M - 2)],
                      pressure=pressure(M, orders),
                      gaussian_pressure=gaussian_comparison(M, orders, True),
                      seconds=time.monotonic() - start)
        report['records'].append(record)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
        print(dict(M=M, complete_active_orders=len(orders), rho=record['pressure']['rho'],
                   gaussian_rho=record['gaussian_pressure']['rho'], seconds=record['seconds']), flush=True)
    report['complete'] = True
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
