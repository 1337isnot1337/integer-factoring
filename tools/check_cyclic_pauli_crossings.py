#!/usr/bin/env python3
"""Independent packed polynomial and complete weighted-profile Pauli controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np

from check_qutrit_weyl_crossing import field4_multiply
from cyclic_pauli_crossings import construct, literal_seed
from cyclic_pauli_pressure import ScaledPauliPressure
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply, scalar_multiply
from windowed_f4_smith import pack


def independent_monomial(a, b):
    length, permutation, phases = len(a), [], []
    for index in range(3**length):
        source = [(index//3**(length-1-i)) % 3 for i in range(length)]
        destination, phase = 0, 1
        for ai, bi, si in zip(a, b, source):
            destination = 3*destination+(si+ai) % 3
            phase = field4_multiply(phase, (1,2,3)[(2*ai*bi+bi*si) % 3])
        permutation.append(destination)
        phases.append(phase)
    for i in range(len(permutation)):
        j, coefficient = i, 1
        for _ in range(3):
            coefficient = field4_multiply(coefficient, phases[j])
            j = permutation[j]
        assert j == i and coefficient == 1
    return permutation, phases


def independent_columns(planes, metadata, columns):
    n, checked = metadata['q']**2, 0
    actions = [(independent_monomial(a,b), c, degree) for a,b,c,degree in reversed(metadata['word'])]
    for column in columns:
        vector = [int(i == column) for i in range(n)]
        for (permutation, phases), coefficient, degree in actions:
            def action(values):
                result = [0]*n
                for i, value in enumerate(values):
                    if value:
                        result[permutation[i]] ^= scalar_multiply(value, phases[i])
                return result
            first = action(vector)
            second = action(first)
            projected = [a^b^c for a,b,c in zip(vector, first, second)]
            correction = [scalar_multiply(value, coefficient) for value in projected]
            vector = [(value << (2*degree)) ^ c ^ (c << 2) for value,c in zip(vector,correction)]
        for row, value in enumerate(vector):
            assert value == pack((int(planes[0][row,column]), int(planes[1][row,column])))
            checked += 1
    return checked


def independent_weighted_orders(model, parameters):
    q, stride = model.base.qx, model.base.stride
    answer = []
    for side, blocks in enumerate(model.base.oracle.blocks):
        h = [[0]+list(map(int,parameters[side*stride+leg*(q-1):side*stride+(leg+1)*(q-1)]))
             for leg in range(4)]
        orders = []
        for block in blocks:
            labels = list(map(int, block['indices']))
            rw = [h[0][i//q]+h[2 if side == 0 else 1][i%q] for i in labels]
            cw = [h[3][j//q]+h[1 if side == 0 else 2][j%q] for j in labels]
            matrix = [[pack(value) for value in row] for row in block['matrix']]
            result,_ = WeightedSmithOracle(matrix, infinity=side == 0).profile(rw,cw,int(parameters[-1]))
            orders.extend(result)
        answer.extend(sorted(orders))
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, begun = np.random.default_rng(27320261101), time.monotonic()
    entries, orders, finite, records, gradient_errors = 0, 0, 0, [], []
    base1, meta1 = construct(1)
    base2, meta2 = construct(2)
    packed1 = [[pack((int(base1[0][i,j]),int(base1[1][i,j]))) for j in range(9)] for i in range(9)]
    tensor_entries = 0
    for i in range(81):
        x,y = divmod(i,9)
        first_row,second_row = (x//3)*3+y//3, (x%3)*3+y%3
        for j in range(81):
            v,u = divmod(j,9)
            first_col,second_col = (v//3)*3+u//3, (v%3)*3+u%3
            value = polynomial_multiply(packed1[first_row][first_col],packed1[second_row][second_col])
            assert value == pack((int(base2[0][i,j]),int(base2[1][i,j])))
            tensor_entries += 1
    fixtures = [(1, []),
        (2,[{'position':1,'a':[1,1,1,1],'b':[1,0,0,1],'coefficient':2}]),
        (3,[{'position':1,'a':[1,1,0,1,1,0],'b':[1,0,1,0,1,0],'coefficient':2}])]
    for k, loops in fixtures:
        planes, metadata = construct(k,loops)
        q, n = metadata['q'], metadata['q']**2
        columns = list(range(n)) if k == 1 else [0,n//2,n-1] if k == 2 else [0,n-1]
        entries += independent_columns(planes,metadata,columns)
        model = ScaledPauliPressure(planes,metadata)
        parameters = np.r_[rng.integers(-2,3,8*(q-1)),1]
        native = model.base.oracle.prefixes(parameters,1)['orders']
        exact = independent_weighted_orders(model,parameters)
        assert native == exact
        orders += len(exact)
        record = {'metadata': metadata, 'parameters': parameters.tolist(),
                  'independent_complete_weighted_orders': exact,
                  'exact_unweighted_bare_order_sum': int(sum(exact)-q*sum(parameters[:-1])),
                  'block_sizes': [[len(b['indices']) for b in side] for side in model.base.oracle.blocks]}
        if k == 2:
            full = ScaledAsymmetricPressure(literal_seed(planes,metadata))
            assert full.base.oracle.prefixes(parameters,1)['orders'] == native
            comparisons = []
            for logtan in (-1.2,.4,1.1):
                point = np.r_[rng.normal(0,.2,64),math.log(.08),logtan]
                a,ga,_ = model.evaluate(point)
                b,gb,_ = full.evaluate(point)
                assert abs(a-b) < 1e-10 and max(abs(ga-gb)) < 2e-7
                comparisons.append({'point':point.tolist(),'rate':a,'gradient_error':float(max(abs(ga-gb)))})
            record['full_native_comparisons'] = comparisons
        if k >= 2:
            point = np.r_[rng.normal(0,.17,8*(q-1)),math.log(.07),.3]
            _,gradient,_ = model.evaluate(point)
            coordinates = (list(range(len(point))) if k == 2 else
                sorted({(q-1)*leg+j for leg in range(8) for j in (0,q-2)} | {len(point)-2,len(point)-1}))
            errors = []
            for coordinate in coordinates:
                delta = np.zeros_like(point)
                delta[coordinate] = 1e-6
                finite_value = (model.evaluate(point+delta,False)[0]-model.evaluate(point-delta,False)[0])/2e-6
                error = abs(float(finite_value-gradient[coordinate]))
                assert error < 2e-5, (k,coordinate,error)
                errors.append(error)
                finite += 1
            gradient_errors.extend(errors)
            record['finite_gradient_control'] = {'point':point.tolist(),'coordinates':coordinates,'errors':errors}
        records.append(record)
        print(json.dumps({'fixture_k':k,'independent_orders':len(exact),'seconds':time.monotonic()-begun}),flush=True)
    names = ['cyclic_pauli_crossings.py','cyclic_pauli_pressure.py','check_cyclic_pauli_crossings.py',
             'regularized_transpose_words.py','check_qutrit_weyl_crossing.py','weighted_smith_oracle_f4.py',
             'native_matrix_smith_f4.py','native_matrix_smith_f4.cpp','native_rectangular_smith_f4.cpp',
             'scaled_asymmetric_crossing_pressure.py','asymmetric_crossing_pressure.py']
    report = {'status':'pass','independent_packed_polynomial_entries':entries,
              'full_two_cell_tensor_identity_entries':tensor_entries,
              'independent_untruncated_weighted_orders':orders,
              'finite_gradient_comparisons':finite,'maximum_gradient_error':max(gradient_errors),
              'records':records,'seconds':time.monotonic()-begun,
              'source_sha256':{n:sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope':'Exact literal coefficients, full direct-sum weighted profiles at new integer gauges, '
                       'and numerical derivative controls. All q27 weighted orders are independently '
                       'recomputed untruncated block by block; no favorable pressure or score is claimed.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','independent_packed_polynomial_entries',
            'independent_untruncated_weighted_orders','finite_gradient_comparisons','seconds')}))


if __name__ == '__main__':
    main()
