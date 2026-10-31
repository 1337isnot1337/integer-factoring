#!/usr/bin/env python3
"""Reconstruct every saved symmetric tile and replay complete root costs."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import fsum, lcm
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from prime_state_filter_pressure import gaussian_pressure
from search_symmetric_tiled_filters import polynomial_hash
from symmetric_binary_filter_synthesis import synthesize as binary_synthesize
from symmetric_prime_filter_synthesis import synthesize as prime_synthesize


def block_log_cost(prime, n, energies, row, tau, charge):
    log_coefficients = np.zeros(1)
    for energy in energies:
        expanded = []
        width = len(log_coefficients) + prime - 1
        for digit in range(prime):
            values = np.full(width, -np.inf)
            values[digit:digit + len(log_coefficients)] = log_coefficients - tau * energy * digit / n
            expanded.append(values)
        log_coefficients = logsumexp(np.asarray(expanded), axis=0)
    weights = np.arange(len(log_coefficients))
    physical = log_coefficients - logsumexp(log_coefficients)
    shifted = log_coefficients - tau * charge * weights
    shifted -= logsumexp(shifted)
    orders = np.asarray(row['beta_orders']) / row['source_m']
    degrees = np.asarray(row['degrees']) / row['source_m']
    return float(logsumexp(shifted + 2 * tau * (weights - orders))
                 + logsumexp(physical + 2 * tau * degrees))


def physical_entropy(prime, n, charge, scale, rho):
    centers = 4 * np.arange(n) - 2 * (n - 1)
    energies = np.r_[centers - charge, centers + charge]
    exponents = -energies[:, None] * scale / (2 * rho) * np.arange(prime)
    normalized = exponents - logsumexp(exponents, axis=1)[:, None]
    return float(-np.sum(np.exp(normalized) * normalized))


def check_catalogue(path):
    source = json.loads(path.read_text())
    catalogue = source['catalogue']
    prime, modes = source.get('prime', 2), source['block_size']
    assert len(catalogue) == source['catalogue_profiles']
    assert sha256(json.dumps(catalogue, separators=(',', ':')).encode()).hexdigest() == source['catalogue_sha256']
    normalized = set()
    for row in catalogue:
        arguments = (modes, row['source_m'], row['source_degree_bounds'], row['source_order_bounds'])
        witness = binary_synthesize(*arguments) if prime == 2 else prime_synthesize(prime, *arguments)
        assert witness
        assert tuple(row['beta_orders']) == witness['actual_orders']
        assert tuple(row['degrees']) == witness['actual_degrees']
        assert tuple(row['endpoint_values']) == witness['endpoint_values']
        assert row['coefficient_field_order'] == witness['coefficient_field_order']
        assert row['witness_polynomial_sha256'] == polynomial_hash(witness)
        normalized.add((tuple(Fraction(value, row['source_m']) for value in row['beta_orders']),
                        tuple(Fraction(value, row['source_m']) for value in row['degrees'])))
    assert len(normalized) == len(catalogue)
    common_ramification = lcm(*(row['source_ramification'] for row in catalogue))
    summaries = []
    for run in source['optimized_runs']:
        n = run['n']
        blocks = run['actual_mode_energy_blocks']
        untouched = run.get('untouched_spectator_mode_energies', [])
        assert all(len(block) == modes for block in blocks)
        assert sorted([value for block in blocks for value in block] + untouched) == list(range(1, 2 * n, 2))
        indices = run['selected_catalogue_indices']
        assert len(indices) == len(blocks)
        assert sum(index >= 0 for index in indices) == run['nonbaseline_blocks']
        values = []
        for block, index in zip(blocks, indices):
            if index < 0:
                assert index == -1
                values.append(0.0)
            else:
                assert index < len(catalogue)
                row = catalogue[index]
                assert n * common_ramification % row['source_m'] == 0
                values.append(block_log_cost(prime, n, block, row, run['tau'], run['scaled_charge']))
        relative = fsum(values)
        assert abs(relative - run['relative_log_cost']) < 1e-8
        assert all(abs(a - b) < 1e-8 for a, b in zip(values, run['block_relative_log_costs']))
        root_residual = gaussian_pressure(prime, n, run['charge'], run['scale'], run['rho']) + relative
        critical = gaussian_pressure(prime, n, run['charge'], run['scale']) + relative
        assert abs(root_residual) < 1e-7
        assert abs(critical - run['full_critical_pressure']) < 1e-8
        entropy = physical_entropy(prime, n, run['charge'], run['scale'], run['rho'])
        assert entropy > 1e-8
        summaries.append({'n': n, 'partition': run['partition'], 'rho': run['rho'],
                          'root_residual': root_residual, 'physical_entropy': entropy,
                          'full_critical_pressure': critical,
                          'nonbaseline_blocks': run['nonbaseline_blocks']})
    return {'file': path.name, 'prime': prime, 'block_size': modes,
            'exact_reconstructed_polynomial_profiles': len(catalogue),
            'common_odd_ramification': common_ramification,
            'complete_root_replays': summaries}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalogues', nargs='*', type=Path)
    parser.add_argument('--output', type=Path, default=root / 'evidence' / 'symmetric-filter-catalogue-check.json')
    args = parser.parse_args()
    names = ('symmetric-tiled-filters-r4.json', 'symmetric-tiled-filters-r8.json',
             'symmetric-qutrit-tiled-filters-r2.json', 'symmetric-qutrit-tiled-filters-r4.json')
    paths = args.catalogues or [root / 'evidence' / name for name in names]
    records = [check_catalogue(path) for path in paths]
    result = {'status': 'Every catalogue polynomial and each reported complete positive-root evaluation replayed.',
              'scope': 'Finite exact witness reconstruction and floating search-value replay; no all-family pressure claim.',
              'exact_reconstructed_polynomial_profiles': sum(row['exact_reconstructed_polynomial_profiles'] for row in records),
              'completed_root_replays': sum(len(row['complete_root_replays']) for row in records),
              'catalogues': records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
