#!/usr/bin/env python3
"""Bounded complete-profile experiments for homogeneous F_2 crossings.

The C++ path uses exact carryless polynomial arithmetic and fraction-free
Pfaffian elimination.  Small profiles are independently compared with the
Python determinant algorithm and the original normal-order recursion.
"""

import argparse
from functools import lru_cache
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess
import tempfile

from check_angle_pfaffian import pfaffian_polynomial
from check_homogeneous_central_selection import independent_rows
from check_homogeneous_central_symbols import binomial_two
from check_modular_crossing import column_bits, determinant_operator_two, rank_two, regular_partitions
from check_partition_channels import operator


def selected_domains(j: int, h: int, wanted: int | None = None):
    parameters = regular_partitions(2, j)
    partners = regular_partitions(2, h)
    d = len(parameters) if wanted is None else wanted
    if not 1 <= d <= 64:
        raise ValueError("the bounded full-tensor checker supports 1 to 64 labels")
    rows = [sum(determinant_operator_two(a, (), b) << k
                for k, b in enumerate(partners)) for a in parameters]
    selected = independent_rows(rows, d)
    if selected is None:
        raise ValueError(f"vacuum block ({j},{h}) does not have rank {d}")
    indices, columns = selected
    a = tuple(parameters[i] for i in indices)
    b = tuple(partners[i] for i in columns)
    return a, b, a


def encode_input(j, h, energy, a, b, states):
    lines = [f"{j} {h} {energy} {len(a)}"]
    for domain in (a, b, states):
        lines.extend(" ".join(map(str, (len(parts), *parts))) for parts in domain)
    return "\n".join(lines) + "\n"


@lru_cache(None)
def python_profile(a, u, b, y):
    j, h, energy = sum(a), sum(b), sum(y)
    assert sum(u) == energy and j <= h
    vertices = ([('U', t) for t in u] + [('A', n) for n in a]
                + [('B', m) for m in b] + [('Y', t) for t in y])
    matrix = [[0] * len(vertices) for _ in vertices]
    for row, (first, n) in enumerate(vertices):
        for col in range(row + 1, len(vertices)):
            second, m = vertices[col]
            value = 0
            if (first, second) == ('U', 'A') and (m - 1) & ~(n - 1) == 0:
                value = binomial_two(m) << n
            elif (first, second) == ('A', 'Y') and n & m == 1:
                value = binomial_two(n)
            elif (first, second) == ('A', 'B') and n & m == 1:
                value = 1 << n
            elif (first, second) == ('U', 'B') and (m - 1) & ~(n - 1) == 0:
                value = binomial_two(m)
            elif (first, second) == ('B', 'Y') and n & m == 1:
                value = binomial_two(n) << m
            elif (first, second) == ('U', 'Y') and n == m:
                value = 1 << m
            matrix[row][col] = matrix[col][row] = value
    return pfaffian_polynomial(matrix)


def original_profile(a, u, b, y):
    j, h, energy = sum(a), sum(b), sum(y)
    result = 0
    indices = [i for i in range(h + j + 1) if i & ~(h + j) == 0]
    for shift in range(-energy - j, energy + h + 1):
        value = 0
        for index in indices:
            middle = energy + h - shift - index
            for state in regular_partitions(2, middle):
                value ^= (operator(a, u, state) & 1) & (operator(b, state, y) & 1)
        result |= value << (shift + energy + j)
    return result


def compile_checker():
    source = Path(__file__).with_suffix('.cpp')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    cpu_info = Path('/proc/cpuinfo')
    accelerated = cpu_info.exists() and 'pclmulqdq' in cpu_info.read_text()
    flags = ['-mpclmul', '-msse2'] if accelerated else []
    binary = Path(tempfile.gettempdir()) / (
        'homogeneous_profiles_' + digest[:16] + ('_clmul' if accelerated else '_scalar'))
    if not binary.exists():
        subprocess.run(['g++', '-std=c++17', '-O3', *flags,
                        str(source), '-o', str(binary)], check=True)
    return binary, digest


def run_case(binary, j, h, wanted=None, emit=False):
    a, b, states = selected_domains(j, h, wanted)
    command = [str(binary)] + (['--emit-profiles'] if emit else [])
    result = subprocess.run(command, input=encode_input(j, h, j, a, b, states),
                            text=True, stdout=subprocess.PIPE, check=True)
    report = json.loads(result.stdout)
    report['operator_partitions'] = [list(x) for x in a]
    report['partner_partitions'] = [list(x) for x in b]
    report['state_partitions'] = [list(x) for x in states]
    if emit:
        profiles = report.pop('profiles_words')
        index = 0
        degree_bound = 3 * j + h
        minimum_valuation = degree_bound + 1
        jet_indices = [1 << k for k in range(degree_bound.bit_length())]
        jet_rows = {k: [0] * (len(a) ** 2) for k in jet_indices}
        for x, y, v, u in product(range(len(a)), repeat=4):
            expected = python_profile(a[x], states[u], b[v], states[y])
            actual = sum(word << (64 * k) for k, word in enumerate(profiles[index]))
            assert actual == expected, (j, h, x, y, u, v)
            residual = actual
            if y == u and determinant_operator_two(a[x], (), b[v]):
                residual ^= 1 << (2 * j)
            translated = 0
            while residual:
                bit = residual & -residual
                translated ^= binomial_two(bit.bit_length() - 1)
                residual ^= bit
            if translated:
                minimum_valuation = min(minimum_valuation,
                                        (translated & -translated).bit_length() - 1)
            for k in jet_indices:
                if (translated >> k) & 1:
                    jet_rows[k][x * len(a) + y] ^= 1 << (v * len(a) + u)
            index += 1
        report['independent_python_profiles'] = index
        assert minimum_valuation == report['minimum_residual_taylor_valuation']
        assert [[k, rank_two(jet_rows[k]), rank_two(column_bits(jet_rows[k], len(a)))]
                for k in jet_indices] == report['power_two_residual_jet_ranks']
    return report


def verify(binary):
    original = 0
    for j, h in ((1, 1), (3, 5), (4, 4), (7, 9), (8, 8)):
        a, b, states = selected_domains(j, h)
        for x, y, v, u in product(range(len(a)), repeat=4):
            assert python_profile(a[x], states[u], b[v], states[y]) == original_profile(
                a[x], states[u], b[v], states[y])
            original += 1
    records = [run_case(binary, j, h, emit=True)
               for j, h in ((1, 1), (3, 5), (8, 8), (15, 17))]
    return {'original_operator_profiles': original,
            'cpp_python_profiles': sum(r['independent_python_profiles'] for r in records),
            'small_cases': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', default='', help='comma-separated operator:partner weights')
    parser.add_argument('--limit-labels', type=int)
    parser.add_argument('--skip-verification', action='store_true')
    parser.add_argument('--output')
    args = parser.parse_args()
    binary, digest = compile_checker()
    report = {'status': 'pass', 'scope': 'finite complete-piece ranks; no improved factoring exponent',
              'cpp_sha256': digest,
              'rank_invariance': 'the nonsingular vacuum minor is inverted on the x leg; both cut ranks are unchanged'}
    if not args.skip_verification:
        report['verification'] = verify(binary)
    records = []
    for case in filter(None, args.cases.split(',')):
        j, h = map(int, case.split(':'))
        records.append(run_case(binary, j, h, args.limit_labels))
        report['complete_profiles'] = records
        if args.output:
            Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        print(json.dumps({'case': case, 'labels': records[-1]['labels'],
                          'rank_sum': records[-1]['rank_sum'],
                          'threshold': records[-1]['finite_seed_threshold'],
                          'maximum_rank_product': records[-1]['maximum_rank_product'],
                          'seconds': records[-1]['seconds']}), flush=True)
    report['complete_profiles'] = records
    encoded = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        Path(args.output).write_text(encoded)
    print(json.dumps({'status': 'pass', 'cases': len(records),
                      'output': args.output}), flush=True)


if __name__ == '__main__':
    main()
