#!/usr/bin/env python3
"""T-primary cut dimensions of the twelve-site coefficient algebra.

The class sum T satisfies x(x-1)^2(x+1)^2.  Its trace powers on the row
and column cut modules therefore determine the dimensions of the three
generalized eigenspaces modulo five, together with the nilpotent trace
corrections.  The result is a checksum and a structural invariant, not
an exact rank computation.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRIME = 5


def load(name):
    return json.loads((ROOT / 'evidence' / name).read_text())


def product(regular, left, right):
    return (np.tensordot(left, regular, axes=(0, 0)) @ right) % PRIME


def t_power_traces(regular, unit, class_sum, trace_vector, count=6):
    powers = [unit.copy()]
    for _ in range(count - 1):
        powers.append(product(regular, powers[-1], class_sum))
    return [int(power @ trace_vector % PRIME) for power in powers]


def primary_dimensions(traces):
    rows = [
        [1, pow(4, k, PRIME), k % PRIME, k * pow(4, k, PRIME) % PRIME]
        for k in range(1, 5)
    ]
    target = traces[1:5]
    solutions = []
    for candidate in itertools.product(range(PRIME), repeat=4):
        if all(sum(a * b for a, b in zip(row, candidate)) % PRIME == value
               for row, value in zip(rows, target)):
            solutions.append(candidate)
    if len(solutions) != 1:
        raise ValueError(f'expected one T-primary solution, got {solutions}')
    d1, d4, m1, m4 = solutions[0]
    d0 = (traces[0] - d1 - d4) % PRIME
    predicted = (d1 + d4 * pow(4, 5, PRIME) + 5 * m1
                 + 5 * pow(4, 5, PRIME) * m4) % PRIME
    if predicted != traces[5]:
        raise ValueError('T-primary trace recurrence does not reproduce T^5')
    return {
        'dimensions_by_eigenvalue': {'0': d0, '1': d1, '4': d4},
        'nilpotent_trace_corrections': {'1': m1, '4': m4},
        't5_check': predicted,
    }


def idempotent_group_sums(traces, joint_characters):
    groups = {'0': [], '1': [], '4': []}
    for index, (t_value, _) in joint_characters.items():
        value = str(t_value)
        if value not in groups:
            raise ValueError(f'unexpected T character {t_value}')
        groups[value].append(index)
    return {
        value: sum(int(traces[int(index)]) for index in indices) % PRIME
        for value, indices in groups.items()
    }


def check():
    table = load('bounded-spin-twelve-coefficient-table.json')
    cut = load('bounded-spin-twelve-cut-trace-form.json')
    rt = load('bounded-spin-twelve-RT-coefficient.json')
    regular = np.array(table['regular'], dtype=np.int64)
    unit = np.array(table['unit'], dtype=np.int64)
    class_sum = np.array(rt['class_sum_coordinates'], dtype=np.int64)
    joint = rt['class_sum_separator_joint_characters']
    result = {}
    for cut_name in ('row', 'column'):
        trace_vector = np.array(cut['trace_vectors'][cut_name], dtype=np.int64)
        traces = t_power_traces(regular, unit, class_sum, trace_vector)
        primary = primary_dimensions(traces)
        group_sums = idempotent_group_sums(cut['idempotent_trace_residues'][cut_name], joint)
        if group_sums != primary['dimensions_by_eigenvalue']:
            raise ValueError(f'{cut_name}: idempotent group sums do not match T-primary dimensions')
        result[cut_name] = {
            'power_traces_0_to_5': traces,
            **primary,
            'idempotent_group_sums': group_sums,
        }
    return {
        'status': 'pass',
        'length': 12,
        'prime': PRIME,
        'basis_indices': table['basis_indices'],
        'class_sum_coordinates': rt['class_sum_coordinates'],
        'cuts': result,
        'scope': ('T-primary dimensions and nilpotent trace corrections of the '
                  'twelve-site row and column cut modules, modulo five.  They '
                  'are congruences, not exact dimensions; individual central '
                  'idempotent ranks, the two-cut pairing, and the score remain open.'),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({
        'status': result['status'],
        'row': result['cuts']['row']['dimensions_by_eigenvalue'],
        'column': result['cuts']['column']['dimensions_by_eigenvalue'],
        'scope': result['scope'],
    }, indent=2))


if __name__ == '__main__':
    main()
