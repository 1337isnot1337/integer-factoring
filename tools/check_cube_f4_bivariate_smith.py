#!/usr/bin/env python3
"""Replay exact generic F4 row-frame certificates, without specializations."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from random import Random
import subprocess

import cube_f4_bivariate_smith as engine
from cube_bivariate_smith import native_program

# Literal flat rank-two seed from the concurrent search.  The fixture is
# retained here so replay never depends on temporary files or sampling.
MONOMIALS = tuple((x,y,6-x-y) for x in range(7) for y in range(7-x))
BINARY_SEED = (
    32768,0,0,0,0,0,0,0,0,294921,294912,0,32776,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
    0,2129992,2129920,0,64,0,0,0,0,0,0,65,0,2129920,2129992,0,
    32776,0,0,0,0,0,0,0,0,294912,294921,0,32768,0,0,0,
    0,136609792,134217728,0,2129920,0,0,0,0,0,0,2129920,0,134217729,136609792,0,
    294912,0,0,0,0,0,0,0,0,32776,32768,0,294921,0,0,0,
    0,0,0,294912,0,32776,32768,0,0,0,0,0,0,0,0,294921,
    0,134217792,136609792,0,2129992,0,0,0,0,0,0,2129992,0,136609792,134217793,0,
    0,134217793,136609792,0,2129992,0,0,0,0,0,0,2129992,0,136609792,134217792,0,
    294921,0,0,0,0,0,0,0,0,32768,32776,0,294912,0,0,0,
    0,0,0,294921,0,32768,32776,0,0,0,0,0,0,0,0,294912,
    0,136609792,134217729,0,2129920,0,0,0,0,0,0,2129920,0,134217728,136609792,0,
    0,0,0,32768,0,294921,294912,0,0,0,0,0,0,0,0,32776,
    0,2129992,2129920,0,65,0,0,0,0,0,0,64,0,2129920,2129992,0,
    0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
    0,0,0,32776,0,294912,294921,0,0,0,0,0,0,0,0,32768,
)
EXPECTED = [0,2,2,4,4,4,6,6,6,6,8,8,8,10,10,12]


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def certificate_digest(record):
    serial = dict(record)
    serial['frames'] = [[[format(a,'x'),format(b,'x')] for a,b in row]
                        for row in record['frames']]
    return digest(serial)


def seed_tensor():
    assert len(BINARY_SEED) == 256
    return [{e:1 for i,e in enumerate(MONOMIALS) if code>>i&1}
            for code in BINARY_SEED]


def binary_reference(matrix):
    """Existing F2 engine on two literal parity blocks, without edits."""
    assert all(b == 0 for row in matrix for a,b in row)
    orders = []
    for parity in range(2):
        labels = [i for i in range(16) if i.bit_count()%2 == parity]
        data = '8\n'+'\n'.join(format(matrix[i][j][0],'x')
                                 for i in labels for j in labels)+'\n'
        result = subprocess.run([str(native_program('cube_bivariate_smith'))],
                                input=data,text=True,capture_output=True)
        assert result.returncode == 0,result.stderr
        orders += json.loads(result.stdout)['orders']
    return sorted(orders)


def change_units(matrix, rng, rounds):
    """Known invertible constant row and column operations over F4."""
    result = [list(row) for row in matrix]
    n = len(result)
    operations = []
    for index in range(rounds):
        a,b = rng.sample(range(n),2)
        c = rng.randrange(1,4)
        if index%2:
            for row in result:
                row[a] = engine.add(row[a],engine.scale(row[b],c))
            operations.append(['column_add',a,b,c])
        else:
            result[a] = [engine.add(x,engine.scale(y,c))
                         for x,y in zip(result[a],result[b])]
            operations.append(['row_add',a,b,c])
    result[0] = [engine.scale(value,2) for value in result[0]]
    operations.append(['row_scale',0,2])
    for row in result:
        row[-1] = engine.scale(row[-1],3)
    operations.append(['column_scale',n-1,3])
    return result,operations


def literal_matrix(tensor,axis):
    first = [i for i in range(8) if not (i>>axis&1)]
    second = [i for i in range(8) if i>>axis&1]
    other = next(i for i in range(3) if i != axis)
    matrix = [[(0,0) for _ in range(16)] for _ in range(16)]
    for mask,polynomial in enumerate(tensor):
        row = sum(((mask>>p)&1)<<i for i,p in enumerate(first))
        column = sum(((mask>>p)&1)<<i for i,p in enumerate(second))
        value = [0,0]
        for exponent,c in polynomial.items():
            for bit in range(2):
                if c>>bit&1:
                    value[bit] ^= 1 << (exponent[axis]+256*exponent[other])
        matrix[row][column] = tuple(value)
    return matrix


def check():
    rng = Random(2026110117)
    seed = seed_tensor()
    binary_controls = []
    transformed_controls = []
    parity_breaking_controls = []
    known_controls = []
    products = determinants = 0
    for axis in range(3):
        matrix = literal_matrix(seed,axis)
        assert matrix == engine.tensor_matrix(seed,axis)
        record = engine.profile_tensor(seed,axis)
        verified = engine.verify_tensor(seed,axis,record)
        assert record['orders'] == EXPECTED == binary_reference(matrix)
        binary_controls.append(dict(axis=axis,orders=record['orders'],
                                    certificate_sha256=certificate_digest(record),verification=verified))
        products += verified['polynomial_product_entries'];determinants += 2
        changed,operations = change_units(matrix,rng,32)
        assert any(value[1] for row in changed for value in row)
        updated = engine.profile_matrix(changed)
        checked = engine.verify_matrix(changed,updated)
        assert updated['orders'] == EXPECTED
        transformed_controls.append(dict(axis=axis,operations=operations,orders=updated['orders'],
                                         certificate_sha256=certificate_digest(updated),verification=checked))
        products += checked['polynomial_product_entries'];determinants += 2

    # Literal nonbinary tensor directions in odd physical parity preserve
    # the crossing endpoint because each added polynomial evaluates to 0.
    mixed = copy.deepcopy(seed)
    mixed[1] = {(6,0,0):2,(0,6,0):2}
    mixed[128] = {(0,0,6):3,(3,3,0):3}
    assert any(mask.bit_count()%2 and polynomial for mask,polynomial in enumerate(mixed))
    for axis in range(3):
        record = engine.profile_tensor(mixed,axis)
        verified = engine.verify_tensor(mixed,axis,record)
        parity_breaking_controls.append(dict(axis=axis,orders=record['orders'],
                                             certificate_sha256=certificate_digest(record),verification=verified))
        products += verified['polynomial_product_entries'];determinants += 2

    # Profiles are known before elimination: constant F4 unit transforms
    # preserve the diagonal x orders and the nonzero units in F4[y].
    for index,n in enumerate((2,3,4,8,16,16)):
        orders = sorted(rng.randrange(7) for _ in range(n))
        diagonal = [[engine.ZERO for _ in range(n)] for _ in range(n)]
        for i,order in enumerate(orders):
            d = rng.randrange(7-order)
            term = engine.scale((1 << (order+256*d),0),rng.randrange(1,4))
            if d:
                term = engine.add(term,engine.scale((1 << order,0),rng.randrange(1,4)))
            diagonal[i][i] = term
        matrix,operations = change_units(diagonal,rng,4*n)
        record = engine.profile_matrix(matrix)
        verified = engine.verify_matrix(matrix,record)
        assert record['orders'] == orders
        known_controls.append(dict(index=index,dimension=n,orders=orders,operations=operations,
                                   certificate_sha256=certificate_digest(record),verification=verified))
        products += verified['polynomial_product_entries'];determinants += 2

    # All four F4 specializations are singular, but this generic matrix
    # has complete orders (0,2).  Sampling cannot certify this example.
    vanishes = (1 << (4*256)) ^ (1 << 256)
    trap = [[(vanishes,0),engine.ZERO],[engine.ZERO,(1 << 2,0)]]
    trap_record = engine.profile_matrix(trap)
    trap_verified = engine.verify_matrix(trap,trap_record)
    assert trap_record['orders'] == [0,2]
    for value in range(4):
        power = 1
        for _ in range(4):
            power = engine.F4_MUL[power][value]
        assert power ^ value == 0
    products += trap_verified['polynomial_product_entries'];determinants += 2

    rejected = 0
    bad = copy.deepcopy(trap_record)
    bad['frames'][0] = [engine.ZERO,engine.ZERO]
    try:
        engine.verify_matrix(trap,bad)
    except AssertionError:
        rejected += 1
    else:
        raise AssertionError('A singular row frame was accepted.')
    bad = copy.deepcopy(trap_record)
    bad['orders'][-1] += 1;bad['prefixes'][-1] += 1
    try:
        engine.verify_matrix(trap,bad)
    except AssertionError:
        rejected += 1
    else:
        raise AssertionError('An incorrect row order was accepted.')
    for matrix in ([[engine.ZERO]], [[(1 << 7,0)]]):
        try:
            engine.profile_matrix(matrix)
        except (ArithmeticError,ValueError):
            rejected += 1
        else:
            raise AssertionError('An invalid matrix contract was accepted.')
    invalid = copy.deepcopy(seed)
    invalid[0][(2,2,2)] = 2
    try:
        engine.profile_tensor(invalid,0)
    except ValueError:
        rejected += 1
    else:
        raise AssertionError('An invalid crossing endpoint was accepted.')
    paths = [Path(__file__),Path(engine.__file__),Path(engine.__file__).with_suffix('.cpp'),
             Path(engine.__file__).with_name('cube_bivariate_smith.cpp'),
             Path(engine.__file__).with_name('cube_bivariate_smith.py'),
             Path(engine.__file__).with_name('cube_homogeneous_gaussian.py')]
    return dict(schema=1,scope='Exact F4(y) Smith certificates; no new factoring exponent claim',
                field='F2[alpha]/(alpha^2+alpha+1)',stride=256,input_total_degree_bound=6,
                binary_seed_sha256=digest(BINARY_SEED),binary_seed_controls=binary_controls,
                nonbinary_unit_transform_controls=transformed_controls,
                parity_breaking_tensor_controls=parity_breaking_controls,
                independent_known_profile_controls=known_controls,
                all_F4_specializations_singular_control=dict(orders=trap_record['orders'],verification=trap_verified),
                exact_profile_count=16,verified_polynomial_product_entries=products,
                nonzero_generic_constant_determinants=determinants,rejected_controls=rejected,
                source_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('evidence/cube-f4-bivariate-smith.json'))
    args = parser.parse_args()
    result = check()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:result[key] for key in ('exact_profile_count','verified_polynomial_product_entries',
                                                'nonzero_generic_constant_determinants','rejected_controls')}))


if __name__ == '__main__':
    main()
