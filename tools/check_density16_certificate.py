#!/usr/bin/env python3
"""Replay the complete binary product-density certificate with a C++20 compiler."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cxx', default='g++')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = root / 'tools' / 'check_density16_native.cpp'
    matrix = json.loads((root / 'evidence' / 'universal-entangler16.json').read_text())
    with tempfile.TemporaryDirectory(prefix='factoring-density-') as directory:
        binary = Path(directory) / 'check-density'
        result = Path(directory) / 'result.json'
        command = [args.cxx, '-O3', '-std=c++20', str(source), '-o', str(binary)]
        subprocess.run(command, check=True)
        subprocess.run([str(binary), '0', '32768', str(result)], check=True,
                       stdout=subprocess.PIPE, text=True)
        data = json.loads(result.read_text())
        assert data['field'] == 2
        assert data['start'] == 0 and data['end'] == 32768
        assert data['checked'] == 2**30
        assert data['matrix_rows'] == matrix['matrix_rows']
        assert data['product_densities'] == 0 and data['witnesses'] == []
        assert data['ordinary_conjugation_checks'] == 96
        assert data['realignment_rank_controls'] == 67
        data['source_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
        data['binary_sha256'] = hashlib.sha256(binary.read_bytes()).hexdigest()
        data['wrapper_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    data['compile_command'] = f'{args.cxx} -O3 -std=c++20 tools/check_density16_native.cpp -o <temporary-binary>'
    data['scope'] = ('No normalized common product density over F2 for the recorded W; '
                     'does not exclude extension-field densities.')
    print(json.dumps(data, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
