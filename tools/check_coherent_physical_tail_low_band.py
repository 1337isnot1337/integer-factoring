#!/usr/bin/env python3
"""Rebuild every leading-band witness and check its minor in independent C++."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np

from coherent_physical_tail_low_band import decomposed, F4, LINES, MUL, power, rowbasis


def verify(path):
    saved = json.loads(path.read_text())
    assert saved['complete'] and not saved['factoring_exponent_improved']
    assert saved['field_polynomial'] == 0x11d and saved['F4_embedding'] == list(F4)
    assert saved['specialization_t'] == 2 and saved['physical_hot_heights'] == [0, 1, 3, 4]
    assert saved['projective_functionals'] == [list(row) for row in LINES]
    assert len(saved['records']) == 8 and saved['exact_pair_checks'] == 57800
    for name, digest in saved['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest
    values = np.array([[F4[c] for c in row] for row in LINES], dtype=np.uint8)
    height = np.array([power(2, -h) for h in (0, 1, 3, 4)], dtype=np.uint8)
    witnesses = bytearray()
    records = iter(saved['records'])
    for side in ('row', 'column'):
        W, coordinates = decomposed(2, side)
        assert W.shape == (16, 256, 16)
        for weights in ((False, False), (False, True), (True, False), (True, True)):
            record = next(records)
            assert record['side'] == side and record['weights'] == list(weights)
            assert record['coordinate_rows'] == coordinates.tolist()
            assert record['shape'] == list(W.shape) and not record['failures']
            assert record['rank_histogram'] == {'16': 7225}
            first = MUL[values, height] if weights[0] else values
            second = MUL[values, height] if weights[1] else values
            for left in range(85):
                partial = np.zeros((4, 256, 16), dtype=np.uint8)
                for a in range(4):
                    for b in range(4):
                        partial[b] ^= MUL[first[left, a], W[a + 4 * b]]
                for right in range(85):
                    matrix = np.bitwise_xor.reduce(MUL[second[right, :, None, None], partial], axis=0)
                    # Select a minor through column elimination of the transpose.
                    # Its determinant is then checked by a separate scalar engine.
                    pivots = rowbasis(matrix.T)
                    assert len(pivots) == 16
                    witnesses.extend(matrix[pivots].tobytes())
    assert len(witnesses) == 57800 * 16 * 16
    source = Path(__file__).with_name('coherent_physical_tail_low_band_check.cpp')
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    with tempfile.TemporaryDirectory(prefix='coherent-physical-tail-witnesses-') as temporary:
        binary = Path(temporary) / 'check'
        subprocess.run(compiler + ['-O2', '-std=c++17', str(source), '-o', str(binary)],
                       check=True, capture_output=True)
        run = subprocess.run([str(binary), '57800'], input=witnesses,
                             check=True, capture_output=True)
    native = json.loads(run.stdout)
    assert native['independent_nonsingular_minors'] == 57800
    return dict(status='passed', factoring_exponent_improved=False,
                exact_field_pair_checks=57800, native_check=native,
                witness_bytes_sha256=sha256(witnesses).hexdigest(),
                receipt_sha256=sha256(path.read_bytes()).hexdigest(),
                source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, source.name)},
                scope='All 85 by 85 nonzero F4 functional pairs, both sides and all four ambient-height weighting patterns; generic leading-band rank16, with no complete high-band profile claim.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.evidence)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
