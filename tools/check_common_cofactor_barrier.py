#!/usr/bin/env python3
"""Complete common-cofactor witnesses, verified by projected determinants."""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from check_joint_determinant_barrier import determinant,power,rectangle_multiplicity
from check_multipoint_hermite_crossings import multiply,rank,interpolation_weights
from check_multipoint_rank_two_seeds import orbit_coverage
from check_transpose_intertwiner_barrier import independent_jet_ranks
from common_cofactor import find_certificate,matrices
from search_multipoint_rank_one_seeds import seed_from_values
from search_multipoint_rank_two_seeds import low_rank_representatives,low_rank_tensors

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def exterior_lift(vector):
    # A 3x4 matrix with the prescribed three-row Pluecker coordinates.
    pivot = next(i for i,a in enumerate(vector) if a)
    normal = [(-1)**i*a % 3 for i,a in enumerate(vector)]
    inverse = pow(normal[pivot],-1,3)
    rows = []
    for j in range(4):
        if j == pivot:
            continue
        row = [int(i == j) for i in range(4)]
        row[pivot] = -normal[j]*inverse % 3
        rows.append(row)
    rows[0] = [a*vector[pivot] % 3 for a in rows[0]]
    for j in range(4):
        matrix = [[[row[k]] for k in range(4) if k != j] for row in rows]
        assert determinant(matrix,3)[0] == vector[j]
    return tuple(map(tuple,rows))


def projected_determinant(matrix,left,right):
    a,b = exterior_lift(tuple(left)),exterior_lift(tuple(right))
    projected = [[[sum(a[i][r]*matrix[r][s][k]*b[j][s]
                       for r in range(4) for s in range(4)) % 3 for k in range(3)]
                  for j in range(3)] for i in range(3)]
    return determinant(projected,3)


def generate(path):
    reps,values = low_rank_representatives(),low_rank_tensors()
    records = []
    for first,a in enumerate(reps):
        for second,b in enumerate(values):
            seed = seed_from_values(a,b)
            certificate = find_certificate(seed)
            assert certificate is not None,(first,second)
            records.append({'first':first,'second':second,**certificate})
    result = {'total':len(records),'records':records,
              'scope':'All normalized ternary binary quadratics with two-cut rank product at most two at 0 and 2.'}
    path.write_text(json.dumps(result,separators=(',',':'))+'\n')


def representation_controls():
    random = Random(314159)
    checks = 0
    for prime in (2,3,5):
        for _ in range(4):
            matrix = [[[a,(int(i == j)-a) % prime]
                       for j in range(3) for a in [random.randrange(prime)]]
                      for i in range(3)]
            p = determinant([[matrix[i][j] for j in range(1,3)] for i in range(1,3)],prime)
            p = power(p,2,prime)
            weights = [random.randrange(prime) for _ in range(5)]
            if sum(a*b for a,b in zip(weights,p)) % prime == 0:
                k = next(i for i,a in enumerate(p) if a)
                weights[k] = (weights[k]+1) % prime
            assert sum(a*b for a,b in zip(weights,p)) % prime
            work = [[[1]]]
            for _ in range(4):
                n = len(work)
                work = [[multiply(work[i//3][j//3],matrix[i%3][j%3],prime)
                         for j in range(3*n)] for i in range(3*n)]
            operator = [[sum(a*b for a,b in zip(weights,value)) % prime
                         for value in row] for row in work]
            assert rank(operator,prime) >= rectangle_multiplicity(2,2)
            checks += 1
    return checks


def actual_jet_controls():
    path = ROOT/'evidence/asymmetric-quadratic-actual-jets.json'
    data = json.loads(path.read_text())
    catalogue = json.loads((ROOT/'evidence/multipoint-rank-two-seeds.json').read_text())['catalogue']
    assert data['tested_witnesses'] == len(data['records']) == 195
    assert data['copies'] == [1,2,3]
    assert len(data['excluded_equal_determinant_signatures']) == 126
    ordered = sorted(data['records'],key=lambda x:x['results'][-1]['best_rank_product'])
    selected = {record['signature_index'] for record in ordered[:6]}
    for record in reversed(data['records']):
        if len(selected) == 12:
            break
        selected.add(record['signature_index'])
    independent = 0
    minima = {1:100000,2:100000,3:100000}
    entries = 0
    for record in data['records']:
        seed = seed_from_values(*catalogue[record['signature_index']]['values'])
        for result in record['results']:
            copies = result['copies']
            counts = result['multiplicities']
            assert sum(counts) == 2*copies+1
            assert result['weights'] == interpolation_weights(2*copies,[0,2,None],counts,3)
            minima[copies] = min(minima[copies],result['best_rank_product'])
            entries += result['complete_reconstructed_entries']
            if record['signature_index'] not in selected:
                continue
            for j,point in enumerate((0,2,None)):
                for cut in (0,1):
                    # The independent oracle uses degree three. At infinity
                    # its extra homogeneous zero adds exactly copies orders.
                    padded = independent_jet_ranks(seed,copies,point,cut)
                    if point is None:
                        assert not any(padded[:copies])
                        actual = padded[copies:]
                    else:
                        actual = padded[:2*copies+1]
                        assert not any(padded[2*copies+1:])
                    assert actual == [pair[cut] for pair in result['complete_jet_cut_ranks'][j]]
                    independent += len(actual)
    assert minima == {1:4,2:64,3:816}
    assert entries == 851760
    return {'witnesses':195,'copies':[1,2,3],'independent_cut_rank_values':independent,
            'complete_reconstructed_entries':entries,'finite_best_rank_products':minima,
            'scope':'One retained witness per unequal-determinant signature, not every coefficient seed.',
            'evidence_sha256':sha256(path.read_bytes()).hexdigest()}


def check(path):
    saved = json.loads(path.read_text())
    reps,values = low_rank_representatives(),low_rank_tensors()
    assert saved['total'] == len(saved['records']) == len(reps)*len(values) == 32265
    general_left = []
    for number,record in enumerate(saved['records']):
        first,second = divmod(number,len(values))
        assert (record['first'],record['second']) == (first,second)
        seed = seed_from_values(reps[first],values[second])
        for key in ('row_left','row_right','column_left','column_right'):
            assert len(record[key]) == 4 and any(record[key])
            assert all(x in (0,1,2) for x in record[key])
        row,col = matrices(seed)
        p = projected_determinant(row,record['row_left'],record['row_right'])
        q = projected_determinant(col,record['column_left'],record['column_right'])
        assert p == q == record['polynomial']
        assert sum(p) % 3 == 1
        if sum(bool(x) for x in record['row_left']) > 1 or sum(bool(x) for x in record['column_left']) > 1:
            general_left.append([first,second])
    assert general_left == [[7,29],[7,56]]
    coverage = orbit_coverage()
    representation = representation_controls()
    jets = actual_jet_controls()
    return {'status':'pass','complete_normalized_seeds':len(saved['records']),
            'complete_actual_coefficient_seeds':len(values)**2,
            'finite_value_orbit_coverage':coverage,
            'direct_projected_determinant_identities':2*len(saved['records']),
            'nonzero_endpoint_identities':len(saved['records']),
            'general_left_projection_cases':general_left,
            'literal_positive_characteristic_cofactor_rank_checks':representation,
            'actual_jet_controls':jets,
            'leading_rho_floor':'log(3)/(2 log(2))',
            'scope':saved['scope'],
            'source_sha256':{name:sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                              for name in ('common_cofactor.py','check_common_cofactor_barrier.py',
                                           'check_joint_determinant_barrier.py',
                                           'search_multipoint_rank_one_seeds.py',
                                           'search_multipoint_rank_two_seeds.py',
                                           'check_multipoint_rank_two_seeds.py',
                                           'check_transpose_intertwiner_barrier.py',
                                           'check_multipoint_hermite_crossings.py',
                                           'search_asymmetric_quadratic_jets.py',
                                           'search_multipoint_actual_jets.py')},
            'catalogue_sha256':sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue',type=Path,default=ROOT/'evidence/common-cofactor-quadratics.json')
    parser.add_argument('--generate',type=Path)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    path = args.catalogue
    if args.generate:
        generate(args.generate)
        path = args.generate
    result = check(path)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
