#!/usr/bin/env python3
"""Direct ring sums and field-reduction controls for Teichmuller carries."""
from hashlib import sha256
from itertools import product
import argparse
import json
from math import comb
from pathlib import Path
import random

from teichmuller_cyclic_carries import (UnramifiedRing, polynomial_sum,
                                      parity_lift_coefficients, residue_root_count, system_root_count)
from joint_frobenius_carries import joint_trace
from check_carry_fourier_elimination import BinaryField

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/teichmuller-cyclic-carry-controls.json'
POLYNOMIALS = {1: (0, 1), 2: (1, 1, 1), 3: (1, 1, 0, 1), 4: (1, 1, 0, 0, 1)}


def evaluate(ring, variables, t, factor, point):
    row, terms = factor
    value = 0
    for coefficient, h in terms:
        columns = len(h)//variables
        term = coefficient
        for v, x in enumerate(point):
            exponent = sum(h[v*columns+j] << (row+t*j) for j in range(columns))
            term = ring.mul(term, ring.power(x, exponent))
        value = ring.add(value, term)
    return value


def direct_sum(ring, width, variables, t, factors, mask):
    q = 1 << width
    domain = [ring.teichmuller(x) for x in range(ring.field_size)]
    domain = [x for x in domain if ring.power(x, q) == x]
    assert len(domain) == q
    answer, assignments = 0, 0
    for point in product(*([x for x in domain if x or not nz] for nz in mask)):
        value = 1
        for factor in factors:
            value = ring.mul(value, evaluate(ring, variables, t, factor, point))
        answer = ring.add(answer, value)
        assignments += 1
    return answer, assignments


def run():
    rng = random.Random(202611017)
    counts = dict(polynomial_queries=0, direct_assignments=0, field_reductions=0,
                  ring_homomorphism_checks=0, teichmuller_moments=0, norm_queries=0,
                  parity_polynomials=0, residue_root_queries=0, system_root_queries=0)
    digest = sha256()
    configurations = [(1, 1, 1, 1), (2, 2, 1, 2), (2, 2, 2, 2),
                      (3, 3, 1, 1), (3, 3, 3, 2), (4, 2, 1, 2), (4, 2, 2, 2)]
    models = []
    for precision in (1, 2, 3, 4):
        for degree in (1, 2, 3, 4):
            ring = UnramifiedRing(precision, POLYNOMIALS[degree])
            for _ in range(8):
                a, b = (rng.randrange(1 << (precision*degree)) for _ in range(2))
                assert ring.frobenius(ring.add(a, b)) == ring.add(ring.frobenius(a), ring.frobenius(b))
                assert ring.frobenius(ring.mul(a, b)) == ring.mul(ring.frobenius(a), ring.frobenius(b))
                orbit = a
                for _ in range(degree):
                    orbit = ring.frobenius(orbit)
                assert orbit == a
                counts['ring_homomorphism_checks'] += 3
            domain = [ring.teichmuller(x) for x in range(ring.field_size)]
            assert len(set(domain)) == ring.field_size
            assert all(ring.power(x, ring.field_size) == x and ring.frobenius(x) == ring.power(x, 2)
                       for x in domain)
            for exponent in range(2*ring.field_size+1):
                observed = 0
                for x in domain:
                    observed = ring.add(observed, ring.power(x, exponent))
                moment = ring.field_size if exponent == 0 else ring.field_size-1 if exponent % (ring.field_size-1) == 0 else 0
                assert observed == ring.integer(moment)
                counts['teichmuller_moments'] += 1
            models.append(dict(degree=degree, precision=precision, polynomial=list(POLYNOMIALS[degree]),
                               frobenius_generator=ring.coefficients(ring.frobenius_generator)))
        for degree, width, t, variables in configurations:
            ring = UnramifiedRing(precision, POLYNOMIALS[degree])
            dimension = variables*(width//t)
            field = BinaryField(sum(c << i for i, c in enumerate(POLYNOMIALS[degree])))
            for trial in range(4):
                if trial == 0:
                    factors = []
                elif trial == 1:
                    factors = [(0, [(ring.integer(3), (0,)*dimension)])]
                else:
                    factors = []
                    for _ in range(2):
                        terms = []
                        for _ in range(3):
                            h = [0]*dimension
                            for _ in range(rng.randrange(trial)):
                                h[rng.randrange(dimension)] += 1
                            terms.append((rng.randrange(1 << (precision*degree)), tuple(h)))
                        factors.append((rng.randrange(t), terms))
                for mask in product((False, True), repeat=variables):
                    expected, assignments = direct_sum(ring, width, variables, t, factors, mask)
                    observed = polynomial_sum(ring, width, t, variables, factors, mask)
                    assert observed == expected, (precision, degree, width, t, variables, trial, mask)
                    counts['polynomial_queries'] += 1
                    counts['direct_assignments'] += assignments
                    if precision == 1:
                        assert joint_trace(field, degree, width, t, variables, factors, [], mask) == observed
                        counts['field_reductions'] += 1
                    digest.update(json.dumps([precision, degree, width, t, variables, factors, mask, observed]).encode())
    for degree in (2, 3):
        for precision in (2, 3, 6):
            ring = UnramifiedRing(precision, POLYNOMIALS[degree])
            base = [(ring.integer(3), (0,)), (ring.generator, (1,)), (1, (2,))]
            factors = [(r, [(ring.frobenius(c, r), h) for c, h in base]) for r in range(degree)]
            expected, assignments = direct_sum(ring, degree, 1, degree, factors, (False,))
            observed = polynomial_sum(ring, degree, degree, 1, factors)
            assert observed == expected and not any(ring.coefficients(observed)[1:])
            counts['norm_queries'] += 1
            counts['direct_assignments'] += assignments
    ring = UnramifiedRing(6, POLYNOMIALS[2])
    full_constant = polynomial_sum(ring, 2, 2, 1, [])
    punctured_constant = polynomial_sum(ring, 2, 2, 1, [], True)
    assert full_constant == 4 and punctured_constant == 3
    naive_root_count = ring.sub(ring.integer(4), polynomial_sum(
        ring, 2, 2, 1, [(r, [(ring.integer(3), (0,))]) for r in range(2)]))
    assert naive_root_count == 32  # True root count of constant 3 modulo 2 is zero.
    assert ring.power(ring.integer(2), 2) != ring.add(ring.power(1, 2), ring.power(1, 2))
    corrected = residue_root_count(ring, 1, [(ring.integer(3), (0,))])
    assert corrected == 0
    for precision in range(1, 13):
        coefficients = parity_lift_coefficients(precision)
        expanded = [0]*(2*precision)
        n = 2*precision-1
        for k in range(precision, 2*precision):
            for j in range(n-k+1):
                expanded[k+j] += comb(n, k)*comb(n-k, j)*(-1)**j
        assert coefficients == expanded
        values = range(1 << precision) if precision <= 8 else (0, 1, 2, 3, (1 << precision)-1)
        for a in values:
            value = sum(c*a**j for j, c in enumerate(coefficients)) % (1 << precision)
            assert value == a % 2
        counts['parity_polynomials'] += 1
    for degree in (1, 2, 3):
        for precision in ((1, 2) if degree == 3 else (1, 2, 3, 4, 6)):
            ring = UnramifiedRing(precision, POLYNOMIALS[degree])
            terms = [(ring.integer(3), (0,)), (ring.generator, (1,)), (1, (2,))]
            expected = sum(ring.mod_two(evaluate(ring, 1, degree, (0, terms), (ring.teichmuller(x),))) == 0
                           for x in range(ring.field_size)) % ring.modulus
            assert residue_root_count(ring, 1, terms) == expected
            counts['residue_root_queries'] += 1
    ring = UnramifiedRing(3, POLYNOMIALS[2])
    linear = [(1, (1, 0)), (ring.generator, (0, 1)), (1, (0, 0))]
    assert residue_root_count(ring, 2, linear) == 4
    counts['residue_root_queries'] += 1
    for precision in (1, 2, 3):
        for second in ([(1, (1,))], [(1, (1,)), (1, (0,))]):
            equations = [[(1, (1,))], second]
            expected = 1 if len(second) == 1 else 0
            assert system_root_count(precision, POLYNOMIALS[1], 1, equations) == expected
            counts['system_root_queries'] += 1
    assert system_root_count(1, POLYNOMIALS[2], 1, [[(1, (1,))], [(2, (1,))]]) == 1
    assert system_root_count(8, POLYNOMIALS[3], 2, []) == 64
    counts['system_root_queries'] += 2
    base = UnramifiedRing(4, POLYNOMIALS[1])
    assert base.unit_inverse(3) == 11 and base.power(3, 1) != 11
    paths = ['tools/teichmuller_cyclic_carries.py', 'tools/check_teichmuller_cyclic_carries.py',
             'tools/joint_frobenius_carries.py', 'tools/check_carry_fourier_elimination.py']
    return dict(scope='Polynomial Teichmuller sums and precision-linear residue-root indicators; inverse-factor compilation and integer word encoding remain open.',
                seed=202611017, models=models, counts=counts, query_sha256=digest.hexdigest(),
                boundary_controls=dict(full_constant=full_constant, punctured_constant=punctured_constant),
                non_root_indicator_control=dict(precision=6, constant=3, true_roots=0, naive_lifted_sum=naive_root_count,
                                                corrected_root_count=corrected),
                source_sha256={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.generate:
        RECEIPT.write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert receipt == json.loads(RECEIPT.read_text())
    print(json.dumps({key: receipt[key] for key in ('scope', 'counts', 'boundary_controls', 'non_root_indicator_control')}))


if __name__ == '__main__':
    main()
