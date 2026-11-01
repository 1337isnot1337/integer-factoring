#!/usr/bin/env python3
"""Combine the two completed W18 threshold censuses for every piece count.

Checks the full affine-one coverage ledger and every retained full matrix.
The existing affine-zero checker supplies its independent count controls.
Neither checker repeats all billions of native threshold rejections.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from math import isqrt
import json
import lzma
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from check_fused_f4_envelope import dense_rank
from check_fused_ordered_projection import certificate as binary_certificate
from check_fused_scalar_projection import DUAL, TARGET
from f4_lattice_linear import add, scale
from fused_f4_envelope import Envelope
from joint_energy_conjugations import Space
from run_fused_f4_subspace import candidate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    archive = root/'evidence/fused-f4-one-census.json.xz'
    packet = json.loads(lzma.decompress(archive.read_bytes()))
    registration, summary = packet['registration'], packet['summary']
    manifest = registration['manifest']
    zero_archive = root/'evidence/fused-f4-zero-census.json.xz'
    zero = json.loads(lzma.decompress(zero_archive.read_bytes()))
    assert registration['manifest'] == zero['registration']['manifest']
    assert registration['affine_value'] == summary['affine_value'] == 1
    assert registration['threshold'] == summary['threshold'] == 501
    assert summary['status'] == zero['summary']['status'] == 'complete'
    directions, total = manifest['directions'], manifest['registered']
    assert total == 4**17 and len(directions) == 17
    space = Space([manifest['offset']]+directions)
    assert len(space.basis) == 18
    assert (manifest['offset'] & DUAL).bit_count() % 2 == 1
    assert all((v & DUAL).bit_count() % 2 == 0 for v in directions)
    for name, expected in manifest['source_sha256'].items():
        assert sha256(registration['as_run_sources'][name].encode()).hexdigest() == expected
    assert sha256((root/'evidence/fused-f4-subspace-input.bin').read_bytes()).hexdigest() == manifest['block_input_sha256']
    envelope = Envelope()
    counts, histogram, position, retained_controls, two_piece_controls = Counter(), Counter(), 0, [], []
    def independent_ranks(point):
        exact = []
        for low,high in envelope.matrices(point):
            matrix = [[((low >> (81*r+c))&1) | (((high >> (81*r+c))&1)<<1)
                       for c in range(81)] for r in range(81)]
            exact.append(dense_rank(matrix))
        assert tuple(exact) == envelope.ranks(point)
        return exact
    for chunk in packet['chunks']:
        data, retained = chunk['summary'], chunk['retained_text']
        assert sha256((json.dumps(data,indent=2,sort_keys=True)+'\n').encode()).hexdigest() == chunk['summary_sha256']
        assert sha256(retained.encode()).hexdigest() == data['retained_sha256']
        assert data['status'] == 'complete' and data['start'] == position
        assert data['stop'] > position and data['visited'] == data['stop']-position
        assert data['zeros'] == data['scalar_orbit_skips'] == 0
        assert data['visited'] == data['retained']+data['rejected']
        records, previous = [json.loads(line) for line in retained.splitlines()], position-1
        assert len(records) == data['retained']
        for record in records:
            p, point, ranks = record['position'], tuple(record['coordinates']), tuple(record['ranks'])
            assert previous < p < data['stop'] and p >= position
            previous = p
            assert point == candidate(p,manifest,1)
            assert (point[0]&DUAL).bit_count() % 2 == 1
            assert (point[1]&DUAL).bit_count() % 2 == 0
            for plane in point:
                assert space.binary_coordinates(plane) is not None
            exact = independent_ranks(point)
            assert tuple(exact) == ranks
            assert ranks[0]*ranks[1] <= 501
            retained_controls.append({'position':p,'coordinates':point,'independent_ranks':exact})
            histogram[ranks] += 1
            product = ranks[0]*ranks[1]
            for scalar in (1,2,3):
                complement = add((TARGET,0),scale(point,scalar))
                assert complement != (0,0)
                complement_ranks = independent_ranks(complement)
                integer_cost_lower_bound = isqrt(product)+isqrt(complement_ranks[0]*complement_ranks[1])
                assert integer_cost_lower_bound >= 28
                two_piece_controls.append({'coordinates':point,'scalar':scalar,
                    'complement_coordinates':complement,'complement_ranks':complement_ranks,
                    'integer_cost_lower_bound':integer_cost_lower_bound})
        for key in ('visited','retained','rejected','zeros','scalar_orbit_skips'):
            counts[key] += data[key]
        position = data['stop']
    assert position == total and len(packet['chunks']) == summary['chunks'] == 512
    assert summary['registered'] == total
    assert all(summary[k] == v for k,v in counts.items())
    assert summary['histogram'] == [{'ranks':list(k),'count':v} for k,v in sorted(histogram.items())]
    assert counts['retained'] == len(retained_controls) > 0
    minimum = min(r*c for r,c in histogram)
    assert minimum == 162
    assert (TARGET & DUAL).bit_count() % 2 == 1
    assert space.binary_coordinates(TARGET) is not None
    assert envelope.ranks((TARGET,0)) == (81,81)
    # The zero gate independently proves the same W and reconstructs its
    # source block input and every chunk's normalization/zero counts.
    with TemporaryDirectory(prefix='fused-f4-both-censuses-') as temporary:
        receipt = Path(temporary)/'zero-controls.json'
        subprocess.run([sys.executable,str(root/'tools/check_fused_f4_zero_census.py'),
                        '--output',str(receipt)],check=True,cwd=root,capture_output=True,text=True)
        zero_controls = json.loads(receipt.read_text())
    assert zero_controls['nonzero_rank_product_floor_in_declared_affine_zero_space'] == 56
    assert total+zero_controls['exact_native_threshold_tests'] == (4**18-1)//3
    # Three pieces: at least one has nonzero f, product>=162; the rest
    # have product>=56. Their lower bound lies strictly between27and28.
    assert 4*162 > 25**2 and 16*56 > 29**2
    assert 162 < 13**2 and 4*56 < 15**2
    # Two pieces outside the catalogue: sqrt502+sqrt56>22+7>28.
    assert 502 > 22**2 and 56 > 7**2 and 22+7 > 28
    # The nine-letter threshold is below27 because sqrt2<3/2.
    assert 2*4 < 3**2 and 9**3 == 27**2 and 28 < 81
    binary = binary_certificate()
    assert binary['minimum_nonzero_product_lower_bound']==56
    assert binary['smallest_listed_product']==162
    assert binary['minimum_listed_complement_product']==3969
    assert isqrt(162)+isqrt(3969)>28
    report = {'status':'pass','affine_one_positions':total,'affine_zero_positions':total,
              'affine_one_native_threshold_tests':counts['visited'],
              'affine_zero_native_threshold_tests':zero_controls['exact_native_threshold_tests'],
              'complete_projective_W18_points':total+zero_controls['exact_native_threshold_tests'],
              'affine_one_exact_rank_product_minimum':minimum,
              'affine_zero_nonzero_rank_product_lower_bound':56,
              'retained_independent_controls':retained_controls,
              'two_piece_complement_controls':two_piece_controls,
              'independent_f4_full_matrix_ranks':2*(len(retained_controls)+len(two_piece_controls)),
              'single_piece_crossing_ranks':[81,81],
              'all_piece_count_cost_lower_bound':'S >= sqrt(162)+2*sqrt(56) > 27 > 9^sqrt(2)',
              'balanced_finite_block_parameter_floor':'rho_star > 3/4',
              'minimum_catalogued_two_piece_integer_cost':min(c['integer_cost_lower_bound'] for c in two_piece_controls),
              'proof':'A nonzero functional value can be scaled to one without changing ranks. '
                      'Every decomposition of E in W has such a piece. A single piece is E, '
                      'of cost 81. For two nonzero pieces outside the catalogue, '
                      'cost >= sqrt(502)+sqrt(56)>29; inside it every scalar complement '
                      'gives cost >= 28. At least three nonzero pieces cost '
                      '>= sqrt(162)+2sqrt(56), a number strictly between 27 and 28.',
              'binary_full_33_dimensional_refinement':{'source_certificate':binary,
                  'rank_sum_lower_bound':'sqrt(162)+2*sqrt(56)>27',
                  'balanced_parameter_floor':'rho_star>3/4',
                  'scope':'Same piece-count refinement of the existing complete binary catalogue. No new binary exhaustive replay.'},
              'zero_controls':zero_controls,
              'archive_sha256':sha256(archive.read_bytes()).hexdigest(),
              'zero_archive_sha256':sha256(zero_archive.read_bytes()).hexdigest(),
              'scope':'All decompositions over F4 with every piece in this fixed18-dimensional '
                      'subspace of the33-dimensional fused envelope, for arbitrary piece count. '
                      'Coverage counts and retained ranks are independently checked; billions '
                      'of native threshold rejections are not independently repeated. '
                      'The bound concerns the balanced finite-block rank-sum criterion. '
                      'The full 33-dimensional F4 space, asymmetric criteria, other fields/generators '
                      'remain open. No factoring improvement.',
              'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in
                  ['tools/check_fused_f4_decomposition_obstruction.py','tools/run_fused_f4_subspace.py',
                   'tools/check_fused_f4_zero_census.py','tools/check_fused_ordered_projection.py','tools/fused_f4_envelope.py',
                   'tools/check_fused_f4_envelope.py','tools/check_fused_scalar_projection.py']}}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','affine_one_positions',
          'affine_zero_positions','affine_one_exact_rank_product_minimum','all_piece_count_cost_lower_bound')}))


if __name__ == '__main__':
    main()
