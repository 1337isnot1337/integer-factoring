#!/usr/bin/env python3
"""Complete homogeneous central ranks on explicit binary label prefixes.

The first 2^floor(log2 D) ascending regular partitions are selected. Their
operator duals are rows of the full weight-h inverse, not a selected inverse.
These finite data do not establish an asymptotic bound.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

from check_modular_central_rank import ROOT, write_input, audit_coefficients
from check_modular_crossing import (
    homogeneous_family_two, homogeneous_piece_two, rank_two, column_bits,
)


def run(weight, binary, temporary):
    path = temporary/f'binary-{weight}.txt'
    labels, states = write_input(path, weight)
    selected = 1 << (len(labels).bit_length()-1)
    all_labels, duals = homogeneous_family_two(weight)
    assert labels == all_labels
    indices = {label: i for i, label in enumerate(labels)}
    with path.open('a') as out:
        out.write(f'{selected}\n')
        for terms in duals[:selected]:
            row = [indices[term] for term in terms]
            out.write(' '.join(map(str, (len(row), *row)))+'\n')
    coefficient_checks = audit_coefficients(binary, path, labels, states)
    result = subprocess.run([str(binary), str(path), '--binary-prefix'],
                            check=True, capture_output=True, text=True)
    records = [json.loads(line) for line in result.stdout.splitlines()]
    progress, row, column = records[-3:]
    assert progress['states'] == progress['selected'] == selected
    assert row['orientation'] == 0 and column['orientation'] == 1
    channels = len(states)+1
    for orientation, entry in enumerate((row, column)):
        left, right = progress['ranks'][2*orientation:2*orientation+2]
        assert entry['kernel_dimensions'] == [channels-left, channels-right]
        assert entry['complete_rank'] == left+right-channels+entry['kernel_pair_rank']
    direct_checked = weight <= 32
    if direct_checked:
        direct = homogeneous_piece_two(weight, 0, labels[:selected], duals[:selected])
        assert rank_two(direct) == row['complete_rank']
        assert rank_two(column_bits(direct, selected)) == column['complete_rank']
    return {
        'weight': weight, 'full_label_count': len(labels), 'binary_label_count': selected,
        'binary_width': selected.bit_length()-1, 'augmented_channels': channels,
        'factor_ranks_row_then_column': progress['ranks'],
        'row_rank': row['complete_rank'], 'column_rank': column['complete_rank'],
        'row_kernel_dimensions': row['kernel_dimensions'],
        'column_kernel_dimensions': column['kernel_dimensions'],
        'row_kernel_pair_rank': row['kernel_pair_rank'],
        'column_kernel_pair_rank': column['kernel_pair_rank'],
        'independent_source_coefficient_checks': coefficient_checks,
        'direct_complete_matrix_checked': direct_checked,
        'input_sha256': sha256(path.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = ROOT/'tools/modular_central_rank.c'
    with tempfile.TemporaryDirectory(prefix='factoring-binary-central-') as name:
        temporary = Path(name)
        binary = temporary/'rank'
        subprocess.run(['cc', '-O3', '-std=c11', '-Wall', '-Wextra', '-Werror',
                        '-o', str(binary), str(source)], check=True)
        weights = (8, 16, 32, 64) if args.full else (8, 16, 32)
        rows = []
        for weight in weights:
            row = run(weight, binary, temporary)
            rows.append(row)
            print(json.dumps(row), flush=True)
    result = {
        'status': 'finite complete ranks on explicit binary prefixes; no asymptotic bound',
        'full': args.full, 'cases': rows,
        'sha256': {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                   for path in (source, Path(__file__),
                                ROOT/'tools/check_modular_central_rank.py')},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
