#!/usr/bin/env python3
"""Audit complete central ranks via the small kernels of their two factors.

Default: weights 8,16,32, compared to the existing complete-matrix evidence.
--full adds weight 64. These are finite exact ranks, not an asymptotic bound.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
import tempfile

from check_modular_crossing import regular_partitions, modular_operator

ROOT = Path(__file__).resolve().parents[1]


def write_input(path, weight):
    labels = regular_partitions(2, weight)
    states = regular_partitions(2, 2*weight)
    with path.open('w') as out:
        for partitions in (labels, states):
            out.write(f'{len(partitions)}\n')
            for partition in partitions:
                out.write(' '.join(map(str, (len(partition), *partition)))+'\n')
    return labels, states


def audit_coefficients(binary, path, labels, states):
    rng = random.Random(20261030)
    triples = [(rng.randrange(len(labels)), rng.randrange(len(labels)),
                rng.randrange(len(states))) for _ in range(300)]
    # Include all diagonal vacuum pairings, even when random nonzero tests are rare.
    triples += [(a, a, 0) for a in range(len(labels))]
    text = ''.join(f'{a} {y} {j}\n' for a, y, j in triples)
    result = subprocess.run([str(binary), str(path), '--coefficients'],
                            input=text, capture_output=True, text=True, check=True)
    actual = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
    expected = [(modular_operator(2, labels[a], states[j], labels[y]),
                 modular_operator(2, labels[a], labels[y], states[j]),
                 modular_operator(2, labels[a], (), labels[y]), int(a == y))
                for a, y, j in triples]
    assert actual == expected
    return len(expected)*4


def run(weight, binary, temporary, log_directory):
    path = temporary/f'input-{weight}.txt'
    labels, states = write_input(path, weight)
    checked = audit_coefficients(binary, path, labels, states)
    result = subprocess.run([str(binary), str(path)],
                            capture_output=True, text=True, check=True)
    records = [json.loads(line) for line in result.stdout.splitlines()]
    if log_directory:
        log_directory.mkdir(parents=True, exist_ok=True)
        (log_directory/f'central-rank-{weight}.jsonl').write_text(result.stdout)
    last, summary = records[-2:]
    assert last['pairs'] == len(labels)**2 or min(
        last['left_rank'], last['right_rank']) == len(states)+1
    columns = len(states)+1
    assert last['left_rank']+summary['kernel_left'] == columns
    assert last['right_rank']+summary['kernel_right'] == columns
    assert summary['central_rank'] == (last['left_rank']+last['right_rank']-columns
                                      + summary['kernel_pair_rank'])
    return {
        'weight': weight, 'label_count': len(labels), 'intermediate_channels': len(states),
        'augmented_channels': columns, 'processed_pairs': last['pairs'],
        'left_rank': last['left_rank'], 'right_rank': last['right_rank'],
        'kernel_left': summary['kernel_left'], 'kernel_right': summary['kernel_right'],
        'kernel_pair_rank': summary['kernel_pair_rank'], 'central_rank': summary['central_rank'],
        'independent_source_coefficient_checks': checked,
        'input_sha256': sha256(path.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--log-directory', type=Path)
    args = parser.parse_args()
    source = ROOT/'tools/modular_central_rank.c'
    reference_path = ROOT/'evidence/modular-crossing.json'
    reference = {row['weight']: row['central_row_rank'] for row in
                 json.loads(reference_path.read_text())['homogeneous_characteristic_two']}
    with tempfile.TemporaryDirectory(prefix='factoring-central-') as name:
        temporary = Path(name)
        binary = temporary/'rank'
        subprocess.run(['cc', '-O3', '-std=c11', '-Wall', '-Wextra', '-Werror',
                        '-o', str(binary), str(source)], check=True)
        weights = (8, 16, 32, 64) if args.full else (8, 16, 32)
        rows = []
        for weight in weights:
            row = run(weight, binary, temporary, args.log_directory)
            if weight in reference:
                assert row['central_rank'] == reference[weight]
            rows.append(row)
            print(json.dumps(row), flush=True)
    result = {
        'status': 'finite exact complete-piece ranks; no asymptotic factoring improvement',
        'full': args.full, 'cases': rows,
        'sha256': {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                   for path in (source, Path(__file__), reference_path)},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
