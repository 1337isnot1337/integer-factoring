#!/usr/bin/env python3
"""Exact sparse-code, decoder and CRT controls; no fast coefficient engine."""
import argparse
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path
from random import Random

from check_carry_fourier_elimination import irreducible
from cyclotomic_character_decoder import (
    modular_compose, modular_inverse, modular_power, modular_product,
    parameters, polynomial_divmod, polynomial_product, polynomial_remainder,
)
from odd_character_residues import primitive_generator, root_of_unity_powers
from sparse_faithful_character import (
    acquire_rough_prime_levels, acquire_sieved_prime_levels, decode_sparse_faithful,
    least_prime_factor, primitive_component,
)
from sparse_divisor_residues import support_polynomial
from divisor_moment_recovery import smallest_integer_root

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/sparse-faithful-character.json'


def supplied_faithful(order, mask):
    """Independent integer cyclotomic root powers, reduced only at the end."""
    modulus, powers = root_of_unity_powers(order)
    values = [0]*(len(modulus)-1)
    for j in range(order//2):
        if mask >> j & 1:
            values = [x+y for x, y in zip(values, powers[j])]
    return sum((x & 1) << j for j, x in enumerate(values))


def run():
    kernel_records, kernel_words = [], 0
    for order in (6, 10, 12, 14, 18, 22, 28, 30, 42, 50, 54, 60, 70, 88):
        e, u, phi, modulus = parameters(order)
        width, degree = order//2, modulus.bit_length()-1
        weights = []
        for multiple in range(1, 1 << (width-degree)):
            word = polynomial_product(multiple, modulus)
            assert word.bit_length() <= width
            weights.append(word.bit_count())
            kernel_words += 1
        prime = least_prime_factor(u)
        witness = sum(1 << (e*j*(u//prime)) for j in range(prime))
        assert min(weights) == prime == witness.bit_count()
        assert polynomial_remainder(witness, modulus) == 0
        kernel_records.append(dict(order=order, e=e, odd_part=u,
                                   exhaustive_nonzero_kernel_words=len(weights),
                                   minimum_weight=min(weights), attaining_word=witness))

    rng, digest, decoded = Random(20261102), sha256(), 0
    decode_records = []
    for order in (2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 28, 30, 44, 46, 58, 88, 98, 210, 264):
        e, u, phi, modulus = parameters(order)
        field = primitive_component(phi)
        assert field == 3 or irreducible(field)
        assert polynomial_remainder(phi, field) == 0
        if u > 1:
            assert modular_power(2, u, field) == 1
            assert all(modular_power(2, u//p, field) != 1
                       for p in range(3, u+1, 2) if u % p == 0
                       and all(p % d for d in range(2, isqrt(p)+1)))
        bound = min(4, (least_prime_factor(u)-1)//2) if u > 1 else min(4, e)
        width, masks = order//2, {0}
        exhaustive = width <= 14
        if exhaustive:
            masks.update(mask for mask in range(1 << width) if mask.bit_count() <= bound)
        else:
            for size in range(1, bound+1):
                masks.update(sum(1 << j for j in rng.sample(range(width), size)) for _ in range(30))
        for mask in sorted(masks):
            observed = supplied_faithful(order, mask)
            assert observed == polynomial_remainder(mask, modulus)
            assert decode_sparse_faithful(order, observed, bound) == mask
            decoded += 1
            digest.update(json.dumps([order, bound, mask, observed]).encode())
        decode_records.append(dict(order=order, e=e, odd_part=u, sparsity=bound,
                                   auxiliary_field_polynomial=field,
                                   complete_sparse_ball=exhaustive, vectors=len(masks)))

    actual, skipped = 0, 0
    for prime in (7, 11, 17, 23, 29, 47, 59, 83):
        order = prime-1
        e, u, phi, modulus = parameters(order)
        bound = min(4, (least_prime_factor(u)-1)//2) if u > 1 else min(4, e)
        generator = primitive_generator(prime)
        logs = {pow(generator, j, prime): j for j in range(order)}
        _, powers = root_of_unity_powers(order)
        for integer in range(2, 97):
            divisors = [d for d in range(1, integer+1) if integer % d == 0]
            mask, coordinates = 0, [0]*(modulus.bit_length()-1)
            for divisor in divisors:
                if divisor % prime:
                    exponent = logs[divisor % prime]
                    mask ^= 1 << (exponent % (order//2))
                    coordinates = [x+y for x, y in zip(coordinates, powers[exponent])]
            if mask.bit_count() > bound:
                skipped += 1
                continue
            observation = sum((x & 1) << j for j, x in enumerate(coordinates))
            assert decode_sparse_faithful(order, observation, bound) == mask
            actual += 1

    # The sharp distance gives two colliding sparse words outside the strict
    # certificate. Phi_3=1+X+X^2: X^2 and 1+X have the same observation.
    _, _, _, modulus = parameters(6)
    assert polynomial_remainder(4, modulus) == polynomial_remainder(3, modulus)
    assert polynomial_remainder(1, 3) == polynomial_remainder(2, 3)
    assert decode_sparse_faithful(8, 1, 1) == 1
    assert decode_sparse_faithful(8, 2, 1) == 2
    # Even at a reduced eligible level, one field component can lose sparse
    # input data. The auxiliary field is used after retaining the full input.
    _, _, phi, full_modulus = parameters(46)
    component = primitive_component(phi)
    assert component == 2787 and component.bit_count() == 7
    left, right = 35, 2752
    assert left ^ right == component and (left.bit_count(), right.bit_count()) == (3, 4)
    assert polynomial_remainder(left, component) == polynomial_remainder(right, component) == 35
    assert polynomial_remainder(left, full_modulus) != polynomial_remainder(right, full_modulus)
    assert decode_sparse_faithful(46, left, 4) == left
    assert decode_sparse_faithful(46, right, 4) == right
    # Actual divisor data has a known twisted-reciprocal promise. At q47
    # inversion exchanges the only two fields, so one supplies both.
    companion, remainder = polynomial_divmod(phi, component)
    assert remainder == 0 and companion == 3189
    assert int(bin(component)[2:][::-1], 2) == companion
    generator = primitive_generator(47)
    logs = {pow(generator, j, 47): j % 23 for j in range(46)}
    reciprocal_digest, reciprocal_vectors = sha256(), 0
    for integer in range(2, 400):
        if integer % 47 == 0:
            continue
        shift, mask = logs[integer % 47], 0
        for divisor in range(1, integer+1):
            if integer % divisor == 0:
                mask ^= 1 << logs[divisor % 47]
        assert mask == sum(((mask >> j) & 1) << ((shift-j) % 23) for j in range(23))
        center = 12*shift % 23
        normalized = sum(((mask >> j) & 1) << ((j-center) % 23) for j in range(23))
        assert normalized & 1 == int(isqrt(integer)**2 == integer)
        observed = supplied_faithful(46, mask)
        first = polynomial_remainder(observed, component)
        second = modular_product(modular_power(2, shift, companion),
            modular_compose(first, modular_power(2, 22, companion), companion), companion)
        correction = modular_product(second ^ polynomial_remainder(first, companion),
                                     modular_inverse(component, companion), companion)
        recovered = first ^ polynomial_product(component, correction)
        assert recovered == observed
        reciprocal_digest.update(json.dumps([integer, shift, first, recovered]).encode())
        reciprocal_vectors += 1
    assert reciprocal_vectors == 390
    assert all(word != sum(((word >> j) & 1) << ((c-j) % 23) for j in range(23))
               for word in (left, right) for c in range(23))

    # A reciprocal sparse promise still does not always supply every field.
    # q503 is eligible even at sparsity14: 2*14 < leastprime(251).
    _, odd, _, phi503 = parameters(502)
    factor503 = primitive_component(phi503)
    assert odd == 251 and factor503 == 1986479280909191
    assert factor503 == int(bin(factor503)[2:][::-1], 2)
    symmetric_kernel = 0
    for j in range(51):
        if factor503 >> j & 1:
            symmetric_kernel ^= 1 << ((j-24) % 251)
            symmetric_kernel ^= 1 << ((j-26) % 251)
    pairs = [j for j in range(1, 126) if symmetric_kernel >> j & 1]
    real_left = sum((1 << j) | (1 << (251-j)) for j in pairs[:7])
    real_right = symmetric_kernel ^ real_left
    for word in (real_left, real_right):
        assert word.bit_count() == 14 and word & 1 == 0
        assert word == sum(((word >> j) & 1) << (-j % 251) for j in range(251))
    assert polynomial_remainder(real_left, factor503) == polynomial_remainder(real_right, factor503)
    assert polynomial_remainder(real_left, phi503) != polynomial_remainder(real_right, phi503)
    assert 2*14 < least_prime_factor(251)
    rejected = []
    for order, observation, bound in ((6, 3, 2), (22, 1 << 10, 2), (22, 3, 0)):
        try:
            decode_sparse_faithful(order, observation, bound)
        except ValueError:
            rejected.append([order, observation, bound])
        else:
            raise AssertionError('Invalid or uncertified sparse input accepted')

    acquisitions = [acquire_rough_prime_levels(s, 4) for s in (1, 2, 4)]
    for record in acquisitions:
        values = [x['prime'] for x in record['levels']]
        assert len(set(values)) == record['count']
        for row in record['levels']:
            q = row['prime']
            assert q % 4 == 3 and q % record['auxiliary_prime'] == row['auxiliary_residue']
            assert all(q % p == 2 for p in record['forbidden_odd_primes'])
            assert gcd(row['progression_residue'], record['progression_modulus']) == 1
            assert all(q % d for d in range(2, isqrt(q)+1))
            u = (q-1)//2
            assert u == 1 or least_prime_factor(u) > 2*record['sparsity']

    sieved = [acquire_sieved_prime_levels(s, 4) for s in (1, 2, 4)]
    for record in sieved:
        assert len(set(record['levels'])) == record['count']
        for q in record['levels']:
            assert q % 4 == 3 and all(q % d for d in range(2, isqrt(q)+1))
            assert (q-1)//2 == 1 or least_prime_factor((q-1)//2) > 2*record['sparsity']
    crt_density_controls = 0
    for limit in (100, 500):
        primes = [q for q in range(3, limit+1, 4)
                  if all(q % d for d in range(2, isqrt(q)+1))]
        for d in range(1, 40, 2):
            assert sum(((q-1)//2) % d == 0 for q in primes) == sum(
                q % (4*d) == 2*d+1 for q in primes)
            phi = sum(gcd(j, d) == 1 for j in range(1, d+1))
            assert sum(gcd(j, 4*d) == 1 for j in range(1, 4*d+1)) == 2*phi
            assert sum(gcd(j, 4*d) == 1 and j % 4 == 3
                       for j in range(1, 4*d+1)) == phi
            crt_density_controls += 1

    # A finite end-to-end supplied-data pipeline on six eligible levels.
    # These small controls certify their no-collision property independently;
    # they do not replace the uniform L-level count in the written theorem.
    crt_records = []
    for integer in (9, 11, 15):
        divisors = [d for d in range(1, integer+1) if integer % d == 0]
        coefficients, crt_modulus = [0]*len(divisors), 1
        for q in (23, 47, 59, 83, 107, 167):
            order = q-1
            generator = primitive_generator(q)
            logs = {pow(generator, j, q): j for j in range(order)}
            _, powers = root_of_unity_powers(order)
            coordinates = [0]*len(powers[0])
            for d in divisors:
                coordinates = [x+y for x, y in zip(coordinates, powers[logs[d % q]])]
            observed = sum((x & 1) << j for j, x in enumerate(coordinates))
            mask = decode_sparse_faithful(order, observed, 4)
            positions = [pow(generator, j, q) for j in range(order)
                         if mask >> (j % (order//2)) & 1]
            assert len(positions) == 2*len(divisors)
            nodes = {x*x % q for x in positions}
            assert len(nodes) == len(divisors)
            reduced = support_polynomial(nodes, q)
            inverse = pow(crt_modulus, -1, q)
            coefficients = [x+crt_modulus*((y-x)*inverse % q)
                            for x, y in zip(coefficients, reduced[:-1])]
            crt_modulus *= q
        assert crt_modulus > 2*(integer*integer+1)**4
        polynomial = [x if 2*x <= crt_modulus else x-crt_modulus for x in coefficients]+[1]
        exact = [1]
        for divisor in divisors:
            following = [0]*(len(exact)+1)
            for j, c in enumerate(exact):
                following[j] -= divisor*divisor*c
                following[j+1] += c
            exact = following
        assert polynomial == exact
        quotient = [0]*(len(polynomial)-1)
        quotient[-1] = 1
        for j in range(len(quotient)-2, -1, -1):
            quotient[j] = polynomial[j+1]+quotient[j+1]
        root = smallest_integer_root(quotient, integer*integer)
        factor = isqrt(root)
        assert factor*factor == root and factor == divisors[1]
        crt_records.append(dict(target=integer, supplied_divisor_count=len(divisors),
                                eligible_levels=[23, 47, 59, 83, 107, 167],
                                reconstructed_divisor_square_polynomial=polynomial,
                                crt_modulus=crt_modulus, smallest_factor=factor))
    names = ['tools/sparse_faithful_character.py', 'tools/check_sparse_faithful_character.py',
             'tools/cyclotomic_character_decoder.py', 'tools/odd_character_residues.py',
             'tools/check_carry_fourier_elimination.py', 'tools/sparse_divisor_residues.py',
             'tools/divisor_moment_recovery.py']
    return dict(scope='Exact sparse faithful-character decoding and fixed CRT acquisition controls. Divisor observations are supplied from explicit divisors; no fast composite-index coefficient engine, asymptotic analytic replay or stronger factoring bound.',
                kernel_controls=kernel_records, decoding_controls=decode_records,
                prime_acquisition_controls=acquisitions, rejected_inputs=rejected,
                sieved_prime_acquisition_controls=sieved,
                sieve_remainder_margin=dict(lower_sieve_exponent=[67, 200],
                                             bombieri_vinogradov_exponent=[1, 2]),
                supplied_crt_factorizations=crt_records,
                collision_outside_certificate=dict(order=6, words=[4, 3]),
                reduced_quotient_collision=dict(order=8, words=[1, 2], reduced_modulus=3),
                one_component_collision=dict(order=46, sparsity=4, words=[left, right],
                                             component_modulus=component,
                                             shared_component_observation=35),
                arithmetic_reciprocity=dict(prime=47, components=[component, companion],
                    supplied_divisor_vectors=reciprocal_vectors,
                    records_sha256=reciprocal_digest.hexdigest(),
                    generic_q47_words_have_no_reciprocal_center=True),
                reciprocal_one_component_collision=dict(prime=503, sparsity=14,
                    words=[real_left, real_right], component_modulus=factor503,
                    scope='Reciprocal sparse vectors, not asserted to be actual divisor data'),
                counts=dict(nonzero_kernel_words=kernel_words, sparse_roundtrips=decoded,
                            actual_divisor_vectors=actual, actual_vectors_outside_promise=skipped,
                            acquired_prime_levels=sum(len(r['levels']) for r in acquisitions),
                            sieved_prime_levels=sum(len(r['levels']) for r in sieved),
                            fixed_four_crt_density_controls=crt_density_controls,
                            complete_supplied_crt_factorizations=len(crt_records)),
                vector_sha256=digest.hexdigest(),
                source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in names},
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
    print(json.dumps(dict(scope=receipt['scope'], counts=receipt['counts'])))


if __name__ == '__main__':
    main()
