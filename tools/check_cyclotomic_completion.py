#!/usr/bin/env python3
"""Exact lifted-field, norm, lattice and supplied-divisor coefficient controls."""
import argparse
from hashlib import sha256
from itertools import product as cartesian_product
import json
from math import gcd, isqrt
from pathlib import Path
from random import Random

from cyclotomic_completion import (
    bits, completion_decoder, hensel_factor, local_remainder, reconstruct_coefficient,
)
from cyclotomic_character_decoder import (
    binary_cyclotomic, polynomial_divmod, polynomial_remainder,
)
from odd_character_residues import (
    multiply_cyclotomic, primitive_generator, root_of_unity_powers,
)
from sparse_faithful_character import decode_sparse_faithful, primitive_component
from sparse_divisor_residues import support_polynomial
from divisor_moment_recovery import smallest_integer_root
from check_factoring_primitive_recipients import determinant

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT/'evidence/one-cyclotomic-completion.json'


def all_components(order):
    remaining, factors = binary_cyclotomic(order), []
    while remaining != 1:
        field = primitive_component(remaining)
        factors.append(field)
        remaining, rest = polynomial_divmod(remaining, field)
        assert rest == 0
    return factors


def local_character_powers(order, factor, precision):
    """Direct powers of -T in the actual finite coefficient ring."""
    modulus = 1 << precision
    current = (1,)+(0,)*(len(factor)-2)
    powers = []
    for _ in range(2*order):
        powers.append(current)
        shifted = local_remainder((0,)+current, factor, precision)
        current = tuple(-x % modulus for x in shifted)
    assert current == powers[0]
    assert powers[order] == ((modulus-1),)+(0,)*(len(factor)-2)
    return powers


def supplied(integer, prime, field, bound):
    """Fixture provider lists divisors; this is explicitly not a fast engine."""
    order = (prime-1)//2
    factor, precision, _, *_ = completion_decoder(order, field, bound)
    local = local_character_powers(order, factor, precision)
    phi, global_powers = root_of_unity_powers(order)
    generator = primitive_generator(prime)
    logs = {pow(generator, j, prime): j for j in range(prime-1)}
    observed, truth, mask = [0]*len(local[0]), [0]*(len(phi)-1), 0
    divisors = [r for r in range(1, integer+1) if integer % r == 0]
    assert gcd(integer, prime) == 1 and len(divisors) <= bound
    for divisor in divisors:
        exponent = logs[divisor % prime]
        observed = [(x+y) % (1 << precision) for x, y in zip(observed, local[exponent])]
        truth = [x+(-1)**exponent*y for x, y in zip(truth, global_powers[exponent % order])]
        mask ^= 1 << (exponent % order)
    return tuple(observed), tuple(truth), mask, divisors


def compute():
    records, digest, actual, generic, norm_controls = [], sha256(), 0, 0, 0
    random = Random(20261103)
    for prime in (7, 11, 19, 31, 43, 47):
        order = (prime-1)//2
        phi, _ = root_of_unity_powers(order)
        for field in all_components(order):
            factor, precision, height, reduced, _, norms, stats, kernel = completion_decoder(order, field, 12)
            dimension, degree = len(phi)-1, len(factor)-1
            assert abs(determinant(reduced)) == 1 << (precision*degree)
            for row in kernel:
                assert not any(local_remainder(row, factor, precision))
            # Independently build the full integral multiplication matrix and
            # check the divisibility used in the separation proof.
            word = kernel[-1]
            multiplication = [list(multiply_cyclotomic(word, (0,)*i+(1,), phi))
                              for i in range(dimension)]
            norm = abs(determinant(multiplication))
            assert norm and norm % (1 << (precision*degree)) == 0
            assert norm**2 <= (dimension*sum(x*x for x in word))**dimension
            norm_controls += 1
            points = 0
            for integer in range(2, 81):
                if gcd(integer, prime) != 1:
                    continue
                observed, truth, mask, _ = supplied(integer, prime, field, 12)
                recovered = reconstruct_coefficient(order, field, 12, observed)
                assert recovered == truth
                binary = sum((x & 1) << i for i, x in enumerate(recovered))
                assert binary == polynomial_remainder(mask, binary_cyclotomic(order))
                digest.update(json.dumps([integer, prime, field, observed, recovered]).encode())
                actual += 1
                points += 1
            for _ in range(8):
                truth = tuple(random.randint(-height, height) for _ in range(dimension))
                observed = local_remainder(truth, factor, precision)
                assert reconstruct_coefficient(order, field, 12, observed) == truth
                generic += 1
            smallest_norm = min(norms)
            records.append(dict(prime=prime, order=order, field=field, dimension=dimension,
                                residue_degree=degree, precision=precision, coefficient_bound=height,
                                actual_points=points, lll=stats,
                                smallest_gram_schmidt_square=[smallest_norm.numerator, smallest_norm.denominator]))

    field = all_components(3)[0]
    factor, precision, *_ = completion_decoder(3, field, 2)
    for truth in cartesian_product(range(-2, 3), repeat=2):
        assert reconstruct_coefficient(3, field, 2, local_remainder(truth, factor, precision)) == truth
        generic += 1

    collisions = []
    symmetric = 0
    f503 = 1986479280909191
    for j in range(51):
        if f503 >> j & 1:
            symmetric ^= 1 << ((j-24) % 251)
            symmetric ^= 1 << ((j-26) % 251)
    pairs = [j for j in range(1, 126) if symmetric >> j & 1]
    left503 = sum((1 << j) | (1 << (251-j)) for j in pairs[:7])
    for order, field, left, right in ((23, 2787, 35, 2752),
                                     (251, f503, left503, symmetric ^ left503)):
        factor = hensel_factor(order, field, 2)
        a, b = (local_remainder(bits(word, order), factor, 2) for word in (left, right))
        assert a != b and all((x-y) % 2 == 0 for x, y in zip(a, b))
        collisions.append(dict(order=order, field=field, words=[left, right],
                               same_mod_two=True, different_mod_four=True))

    # The old complete supplied-data CRT control, now obtaining every full
    # faithful observation from just one lifted field. No cross-level roots
    # or divisor sets are supplied to the CRT consumer.
    crt_records, crt_models = [], []
    levels = (23, 47, 59, 83, 107, 167)
    for integer in (8, 9, 11, 15):
        coefficients, modulus, count = [], 1, None
        for prime in levels:
            order = (prime-1)//2
            field = primitive_component(binary_cyclotomic(order))
            if integer == 8:
                state = completion_decoder(order, field, 4)
                crt_models.append(dict(prime=prime, field=field, dimension=len(state[3]),
                                       residue_degree=len(state[0])-1, precision=state[1]))
            observed, _, _, _ = supplied(integer, prime, field, 4)
            recovered = reconstruct_coefficient(order, field, 4, observed)
            binary = sum((x & 1) << i for i, x in enumerate(recovered))
            mask = decode_sparse_faithful(prime-1, binary, 4)
            generator = primitive_generator(prime)
            positions = [pow(generator, j, prime) for j in range(prime-1)
                         if mask >> (j % order) & 1]
            nodes = {x*x % prime for x in positions}
            if count is None:
                count = len(nodes)
                coefficients = [0]*count
            assert len(nodes) == count and len(positions) == 2*count
            polynomial = support_polynomial(nodes, prime)
            coefficients = [x+modulus*((y-x)*pow(modulus, -1, prime) % prime)
                            for x, y in zip(coefficients, polynomial[:-1])]
            modulus *= prime
        assert modulus > 2*(integer*integer+1)**4
        polynomial = [x if 2*x <= modulus else x-modulus for x in coefficients]+[1]
        quotient = [0]*count
        quotient[-1] = 1
        for j in range(count-2, -1, -1):
            quotient[j] = polynomial[j+1]+quotient[j+1]
        assert polynomial[0]+quotient[0] == 0
        root = smallest_integer_root(quotient, integer*integer)
        factor = isqrt(root)
        assert factor*factor == root and 1 < factor <= integer and integer % factor == 0
        assert all(factor % d for d in range(2, isqrt(factor)+1))
        crt_records.append(dict(integer=integer, smallest_factor=factor,
                                polynomial=polynomial, levels=levels, modulus=modulus))

    _, powers105 = root_of_unity_powers(105)
    assert max(abs(x) for row in powers105 for x in row) > 1
    rejected = 0
    tests = [lambda: hensel_factor(2, 7, 2), lambda: hensel_factor(23, 15, 2),
             lambda: hensel_factor(23, 2787, 0), lambda: hensel_factor(23, 7, 2),
             lambda: completion_decoder(23, 2787, 0),
             lambda: reconstruct_coefficient(3, 7, 1, (0,)),
             lambda: reconstruct_coefficient(3, 7, 1, (-1, 0)),
             lambda: reconstruct_coefficient(3, 7, 1, (2, 0))]
    for test in tests:
        try:
            test()
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('Invalid input or out-of-box residue was accepted')
    names = ['tools/cyclotomic_completion.py', 'tools/check_cyclotomic_completion.py',
             'tools/check_factoring_primitive_recipients.py', 'tools/cyclotomic_character_decoder.py',
             'tools/odd_character_residues.py', 'tools/sparse_faithful_character.py',
             'tools/sparse_divisor_residues.py', 'tools/divisor_moment_recovery.py',
             'tools/check_carry_fourier_elimination.py']
    return json.loads(json.dumps(dict(
        scope='Exact Hensel, integral norm, LLL and one-completion reconstruction controls. Literal-divisor fixtures supply the data; no fast original-index coefficient engine, new Lean, independent mathematical review or stronger complete bound.',
        models=records, actual_divisor_coefficients=actual, generic_bounded_vectors=generic,
        integer_norm_controls=norm_controls, observations_sha256=digest.hexdigest(),
        lifted_collisions=collisions, supplied_crt_factorizations=crt_records,
        crt_completion_models=crt_models,
        rejected_inputs=rejected,
        composite_order_coefficient_bound=dict(order=105, monomial_coefficient_bound=max(abs(x) for row in powers105 for x in row), lattice_reconstruction_replayed=False),
        optional_interrupted_pilot=dict(order=105, dimension=48, residue_degree=12,
                                        exit_code=143, scope='Private exact Python LLL pilot deliberately stopped after about seven minutes before output; not a pass. No larger reconstruction run is included.'),
        source_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in names})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.write:
        EXPECTED.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(EXPECTED.read_text())
    print('One cyclotomic completion: Hensel, norm, exact lattice and supplied CRT controls PASS')
