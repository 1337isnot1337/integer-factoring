#!/usr/bin/env python3
"""Compile and independently replay the optional polynomial certificates."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_json(command):
    completed = subprocess.run(command, check=True, text=True,
                               stdout=subprocess.PIPE)
    return json.loads(completed.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cxx', default='g++')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    tools = root / 'tools'
    evidence = root / 'evidence'
    names = ['smith_mpo_search.cpp', 'smith_walk_search.cpp',
             'smith_walk_probe.cpp', 'ternary_resolvent_probe.cpp',
             'schur_span_recombination.cpp']
    manifest = json.loads((evidence / 'smith-native-source-manifest.json').read_text())
    for name, expected in manifest.items():
        assert digest(tools / name) == expected
    with tempfile.TemporaryDirectory(prefix='factoring-polynomial-') as directory:
        temporary = Path(directory)
        programs = {}
        for name, standard in [('smith_walk_probe', 'c++17'),
                               ('ternary_resolvent_probe', 'c++20'),
                               ('schur_span_recombination', 'c++20')]:
            program = temporary / name
            subprocess.run([args.cxx, '-O3', '-std=' + standard,
                            str(tools / (name + '.cpp')), '-o', str(program)],
                           check=True)
            programs[name] = program
        smith = run_json([sys.executable, str(tools / 'check_degree_sensitive_smith_paths.py'),
                          '--native-probe', str(programs['smith_walk_probe'])])
        ternary_path = temporary / 'ternary.json'
        subprocess.run([str(programs['ternary_resolvent_probe']), str(ternary_path)],
                       check=True, stdout=subprocess.PIPE, text=True)
        ternary = json.loads(ternary_path.read_text())
        recorded = json.loads((evidence / 'ternary-resolvent-probe.json').read_text())
        assert ternary['records'] == recorded['records']
        audit = run_json([sys.executable, str(tools / 'check_ternary_resolvent_probe.py'),
                          '--evidence', str(ternary_path)])
        schur = run_json([sys.executable, str(tools / 'check_schur_span_recombination.py'),
                          '--native', str(programs['schur_span_recombination']),
                          '--pieces', str(evidence / 'homogeneous-h4-pieces.json'),
                          str(evidence / 'full-h4-pieces.json')])
    report = {
        'status': 'passed',
        'source_sha256': {name: digest(tools / name) for name in names},
        'wrapper_sha256': digest(Path(__file__)),
        'smith_crosscheck': smith['independent_native_crosscheck'],
        'ternary_toeplitz_checks': audit['count'],
        'ternary_records_match': True,
        'schur_span_records': schur['records'],
        'scope': 'Exact finite native replay and independent rational or dense-matrix controls; no improved exponent.'
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
