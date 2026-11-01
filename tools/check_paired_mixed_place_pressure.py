#!/usr/bin/env python3
"""Cross-check unequal paired pressure, gradients, and tensor symmetrization."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.special import softmax

from affine_place_profile_oracle import AffinePlaceOracle
from f4_closed_place_oracle import F4ClosedPlaceOracle
from f4_regular_smith_frames import polynomial_product
from paired_mixed_place_pressure import PairedMixedPressure, ScaledPairedMixedPressure
from search_finite_fiber_restrictions import restriction
from search_fused_mixed_weighted import MixedPressure
from search_fused_physical_bases import ScaledHeightPressure


def pair_swapped_tensor(seed):
    """K tensor the copy obtained by swapping its two paired species."""
    a, b = seed['q_x'], seed['q_u']
    q = a*b
    def entry(x, y, u, v):
        return seed['values'][x+a*y+a*b*u+a*b*b*v]
    values = []
    for v in range(q):
        for u in range(q):
            for y in range(q):
                for x in range(q):
                    first = entry(x % a, y % b, u % b, v % a)
                    second = entry(y // b, x // a, v // a, u // b)
                    values.append(polynomial_product(first, second))
    degree = max(plane.bit_length()-1 for value in values for plane in value)
    assert degree == 2*seed['degree']
    return {'values': values, 'q_x': q, 'q_u': q, 'dimensions': [q]*4,
            'degree': degree, 'strict_endpoint': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/deflated-finite-fiber-input.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    seed = json.loads(raw)['cases'][0]['seed']
    rng = np.random.default_rng(61219)
    ordinary = {'row_places': [0, 2], 'column_places': [None], 'row_allocations': [.7, .3],
                'column_allocations': [1.], 'outer_log_scale': .6}
    equal_old, equal_new = AffinePlaceOracle(seed), F4ClosedPlaceOracle(seed)
    equal_errors = []
    for first, second in ((PairedMixedPressure(equal_new, ordinary),
                           MixedPressure(equal_old, ordinary, 'full', False, True, float('inf'))),
                          (ScaledPairedMixedPressure(equal_new, ordinary),
                           ScaledHeightPressure(equal_old, ordinary, 'full', False, True, float('inf')))):
        for case in range(6):
            values = first.initial(ordinary, rng, case)
            values += rng.normal(0, .15, len(values))
            values[-1] = .6+case/10
            a, da = first.evaluate(values)
            b, db = second.evaluate(values)
            error = max(abs(a-b), float(np.max(abs(da-db))))
            assert error < 2e-10
            equal_errors.append(error)
    factor = (1, 0, 2, 1)
    labels = [f'root:{factor}:{i}' for i in range(3)]
    baseline = {'row_places': [0]+labels[:2], 'column_places': [None, labels[2]],
                'row_allocations': [.5, .25, .25], 'column_allocations': [.7, .3],
                'outer_log_scale': .6}
    controls, gradient_coordinates, gradient_errors = [], 0, []
    for first, second in (([0, 1], [0, 1, 3]), ([0, 1, 3], [0, 1])):
        restricted, _ = restriction(seed, first, second)
        oracle = F4ClosedPlaceOracle(restricted)
        for pressure_class in (PairedMixedPressure, ScaledPairedMixedPressure):
            pressure = pressure_class(oracle, baseline)
            for trial in range(3):
                values = pressure.initial(baseline, rng, trial+1)
                values += rng.normal(0, .12, len(values))
                values[-1] = .6+trial/7
                value, analytic = pressure.evaluate(values, .07, False)
                numeric = np.zeros(len(values))
                for index in range(len(values)):
                    offset = np.zeros(len(values))
                    offset[index] = 1e-6
                    numeric[index] = (pressure.evaluate(values+offset, .07, False)[0]
                                      -pressure.evaluate(values-offset, .07, False)[0])/2e-6
                error = float(np.max(abs(numeric-analytic)))
                assert error < 5e-5, error
                gradient_coordinates += len(values)
                gradient_errors.append(error)
            controls.append({'paired_dimensions': [len(first), len(second)],
                             'class': pressure_class.__name__, 'parameter_count': pressure.dimension,
                             'gradient_vectors': 3})
    restricted, _ = restriction(seed, [0, 1], [0, 1, 3])
    product_seed = pair_swapped_tensor(restricted)
    first = PairedMixedPressure(F4ClosedPlaceOracle(restricted), baseline)
    second = PairedMixedPressure(F4ClosedPlaceOracle(product_seed), baseline)
    values = first.initial(baseline, rng, 1)
    values += rng.normal(0, .18, len(values))
    values[-1] = .7
    combined = second.initial(baseline, rng, 0)
    for cut in range(2):
        for index in range(len(first.groups[cut])):
            heights = first.height_vector(values[first.height_slices[cut][index]])
            blocks = [heights[current] for current in first.leg_slices]
            expanded = []
            for leg, other in enumerate((1, 0, 3, 2)):
                expanded.extend(blocks[leg][a]+blocks[other][b]
                                for b in range(len(blocks[other])) for a in range(len(blocks[leg])))
            combined[second.height_slices[cut][index]] = np.asarray(expanded)[second.active_heights]
        combined[second.allocation_slices[cut]] = values[first.allocation_slices[cut]]
    px, py = [softmax(np.r_[0., values[current]]) for current in first.type_slices]
    product_types = [np.asarray([a*b for b in py for a in px]),
                     np.asarray([a*b for b in px for a in py])]
    for current, probabilities in zip(second.type_slices, product_types):
        combined[current] = np.log(probabilities[1:]/probabilities[0])
    combined[-1] = values[-1]
    first.evaluate(values)
    second.evaluate(combined)
    error = abs(first.best['rho']-second.best['rho'])
    assert error < 2e-10
    channel_errors = []
    for old_cut, new_cut in zip(first.best['place_records'], second.best['place_records']):
        for old, new in zip(old_cut, new_cut):
            expected = sorted(a+b for a, b in product(old['restored_orders'], repeat=2))
            mismatch = max(abs(a-b) for a, b in zip(expected, new['restored_orders']))
            assert mismatch < 2e-8
            channel_errors.append(mismatch)
    rejections = 0
    bad = [ordinary | {'column_places': [0]},
           ordinary | {'column_places': ['root:(0, 1):0']},
           ordinary | {'row_places': [0, 0]},
           ordinary | {'row_places': [0, 'root:(1, 0, 2, 1):3']},
           ordinary | {'row_places': [0, 1]}]
    for proposal in bad:
        try:
            PairedMixedPressure(equal_new, proposal)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError('Invalid or repeated place was accepted.')
    names = (Path(__file__).name, 'paired_mixed_place_pressure.py', 'f4_closed_place_oracle.py',
             'search_fused_mixed_weighted.py', 'search_fused_physical_bases.py',
             'search_finite_fiber_restrictions.py', 'f4_regular_smith_frames.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'scope': __doc__, 'passed': True, 'input_sha256': sha256(raw).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'equal_alphabet_controls': len(equal_errors), 'largest_equal_error': max(equal_errors),
              'unequal_gradient_coordinates': gradient_coordinates, 'largest_gradient_error': max(gradient_errors),
              'unequal_controls': controls, 'pair_swapped_pressure_error': error,
              'pair_swapped_profile_count': len(channel_errors), 'largest_convolution_error': max(channel_errors),
              'pair_swapped_endpoint_entries': len(product_seed['values']), 'invalid_place_rejections': rejections}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
