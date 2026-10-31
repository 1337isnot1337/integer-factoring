#!/usr/bin/env python3
"""Adversarial controls for exported Fock coefficient precision metadata.

The ordinary native window rejects this example when its window is kept.
Padding the exported arrays cannot restore discarded coefficients; the
retained common precision bound must reject the resulting false profile.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import rectangular_defect_profiles as profiles
from native_matrix_smith_f4 import MatrixOracle


def synthetic_data():
    images = {}
    for mode in range(2):
        for creation, label in ((True, 'create'), (False, 'annihilate')):
            images[f'm{mode}_{label}'] = {
                'creation': [(int(creation and j == mode), 1) for j in range(2)],
                'annihilation': [(int(not creation and j == mode), 1) for j in range(2)],
            }
    return {'M': 4, 'place': 'zero', 'active_indices': [0, 1],
            'active_orders': [0, 0], 'active_images': images}


class SyntheticProgram(profiles.Program):
    def build(self, *args, **kwargs):
        # [[1,1,0,0], [1,1+z^150,z^200,0],
        #  [0,z^200,z^300,0], [0,0,0,1]].
        terms = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0),
                 (1, 1, 150), (1, 2, 200), (2, 1, 200), (2, 2, 300), (3, 3, 0)]
        expressions = [self.scale(self.unit(profiles.Expr(0), ['m0', 'm1'], i, j), z=degree)
                       for i, j, degree in terms]
        value = expressions[0]
        for term in expressions[1:]:
            value = self.add(value, term)
        return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    with patch.object(profiles, 'Program', SyntheticProgram):
        matrix, shifts, short = profiles.active_matrix(synthetic_data(), profile_precision=128)
        complete_matrix, complete_shifts, complete = profiles.active_matrix(
            synthetic_data(), profile_precision=512)

    assert not short['exact_polynomial_matrix']
    assert complete['exact_polynomial_matrix']
    assert short['known_below'] == 128
    assert short['builder_known_below'] == 1024
    exact = MatrixOracle(complete_matrix, 4, complete_shifts, precision=512).profile()['orders']
    assert exact == [0., 0., 150., 250.]

    try:
        MatrixOracle(matrix, 4, shifts, precision=128).profile()
    except ArithmeticError:
        native_window_rejected = True
    else:
        raise AssertionError('The original window unexpectedly resolved this profile.')

    padded = MatrixOracle(matrix, 4, shifts, precision=512).profile()['orders']
    assert padded == [0., 0., 200., 200.]
    assert padded != exact
    assert max(padded) < short['builder_known_below']
    assert max(padded) >= short['known_below']

    # The retained cutoff transforms by the common scale and a lower
    # bound on row/column weights. An exact longer-window replay is valid.
    row_weights, column_weights, scale = [0.25]*4, [-0.5]*4, 0.5
    weighted = MatrixOracle(matrix, 4, shifts, precision=512).profile(
        row_weights, column_weights, scale=scale)['orders']
    cutoff = short['known_below']*scale+min(row_weights)+min(column_weights)
    assert max(weighted) >= cutoff
    replay = MatrixOracle(complete_matrix, 4, complete_shifts, precision=512).profile(
        row_weights, column_weights, scale=scale)['orders']
    assert replay == [-0.25, -0.25, 74.75, 124.75]

    result = {'status': 'retained coefficient precision is separate from builder precision',
              'source_polynomial_degree': 300, 'short_output_window': 128,
              'short_builder_precision': short['builder_precision'],
              'retained_known_below': short['known_below'],
              'builder_known_below': short['builder_known_below'],
              'exact_profile': exact, 'padded_truncated_profile': padded,
              'same_window_native_rejection': native_window_rejected,
              'old_common_bound_would_accept_padding': True,
              'new_common_bound_rejects_padding': True,
              'weighted_padded_profile': weighted, 'weighted_retained_cutoff': cutoff,
              'weighted_exact_replay': replay,
              'scope': 'metadata hardening for truncated arrays; no failure reproduced in the standard same-window native path'}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
