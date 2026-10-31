"""Deterministic finite norm search and dominant relation certificates.

Successful certificates are exact. The fixed search box has no proved
uniform success or improved factoring-time bound. Run the CLI explicitly;
the large certificate replay has its own Make target.
"""
import argparse
import hashlib
import json
from math import gcd, isqrt
import sys

from arithmetic_two_group import power, project_to_two_part, two_group_basis
from binary_quadratic_forms import ambiguous_factor, compose_forms, inverse_form, reduce_form
from bounded_smooth_norm import norm_to_pair, smooth_factors
from check_positive_relation_kernel import determinant


def prime_base(discriminant, bound):
    sieve = bytearray(b'\x01')*(bound+1)
    sieve[:2] = b'\x00\x00'
    for prime in range(2, isqrt(bound)+1):
        if sieve[prime]:
            sieve[prime*prime::prime] = b'\x00'*(((bound-prime*prime)//prime)+1)
    base = []
    for prime in range(2, bound+1):
        if not sieve[prime] or gcd(prime, discriminant) > 1:
            continue
        roots = [middle for middle in range(prime+1)
                 if (middle*middle-discriminant) % (4*prime) == 0]
        if roots:
            middle = roots[0]
            base.append((prime, middle, reduce_form((prime, middle, (middle*middle-discriminant)//(4*prime)))))
    return base


def roots_x(form, y, prime, root):
    a, b, c = form
    if prime == 2:
        return [x for x in range(2) if (a*x*x+b*x*y+c*y*y) % prime == 0]
    if a % prime:
        return {((-b+sign*root)*y*pow(2*a, -1, prime)) % prime for sign in (-1, 1)}
    if b*y % prime:
        return [(-c*y*y*pow(b*y, -1, prime)) % prime]
    return range(prime) if c*y*y % prime == 0 else []


def find_pair(form, base, bound, radius, x_bound, y_bound):
    a, b, c = form
    if gcd(a, b*b-4*a*c) == 1 and smooth_factors(a, bound) is not None:
        pair = norm_to_pair(form, 1, 0, bound, radius)
        if pair is not None:
            return pair, (1, 0), 1
    candidates = 1
    for y in range(1, y_bound+1):
        values = [a*x*x+b*x*y+c*y*y for x in range(x_bound+1)]
        residual = values[:]
        for prime, root, _ in base:
            for residue in roots_x(form, y, prime, root):
                for x in range(residue, len(values), prime):
                    while residual[x] > 1 and residual[x] % prime == 0:
                        residual[x] //= prime
        for x, remainder in enumerate(residual):
            candidates += 1
            if remainder != 1 or values[x] > radius*radius:
                continue
            pair = norm_to_pair(form, x, y, bound, radius)
            if pair is not None:
                return pair, (x, y), candidates
    return None, None, candidates


def word_vector(form, base):
    remaining, middle, _ = form
    result = []
    for prime, root, _ in base:
        exponent = 0
        while remaining % prime == 0:
            exponent += 1
            remaining //= prime
        sign = 1 if (middle-root) % (2*prime) == 0 else -1
        assert not exponent or (middle-sign*root) % (2*prime) == 0
        result.append(sign*exponent)
    assert remaining == 1
    return result


def relation_row(index, diagonal, pair, base):
    left, right = (word_vector(form, base) for form in pair)
    row = [diagonal*int(index == j)-left[j]+right[j] for j in range(len(base))]
    assert row[index] > sum(abs(value) for j, value in enumerate(row) if j != index)
    return row


def matrix_hash(rows):
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


def extraction(base, annihilator, integer):
    identity = (1, 0, integer)
    generators = [form for _, _, form in base]
    projected = project_to_two_part(generators, annihilator, identity, compose_forms)
    basis = two_group_basis(projected, identity, compose_forms, inverse_form, integer.bit_length()+3)
    tops = [power(form, 1 << (order-1), identity, compose_forms) for form, order in basis]
    factor = next((factor for form in tops for factor in [ambiguous_factor(form, integer)] if factor), None)
    return basis, factor


def generate(integer, bound, x_bound, y_bound):
    assert integer > 1 and integer % 2 == 1 and bound >= 2 and x_bound >= 1 and y_bound >= 1
    discriminant = -4*integer
    radius = isqrt(-discriminant-1)//2
    base = prime_base(discriminant, bound)
    identity = (1, 0, integer)
    diagonal = 1 << (2*integer.bit_length())
    rows, vectors = [], []
    for index, (_, _, generator) in enumerate(base):
        target = power(generator, diagonal, identity, compose_forms)
        pair, vector, candidates = find_pair(target, base, bound, radius, x_bound, y_bound)
        print(f'relation {index+1}/{len(base)}: {candidates} inspected candidates', file=sys.stderr, flush=True)
        if pair is None:
            return {'status': 'search_exhausted', 'integer': str(integer), 'failed_generator_index': index,
                    'completed_relations': len(rows), 'x_bound': x_bound, 'y_bound': y_bound}
        rows.append(relation_row(index, diagonal, pair, base))
        vectors.append(vector)
    annihilator = abs(determinant(rows))
    basis, factor = extraction(base, annihilator, integer)
    return {
        'scope': 'Finite exact relation/torsion certificate; fixed-box search has no uniform success or improved asymptotic guarantee',
        'format': 'dominant-norm-relations-v1', 'status': 'factor_found' if factor else 'relations_verified_no_split',
        'integer': str(integer), 'smooth_bound': bound, 'x_bound': x_bound, 'y_bound': y_bound,
        'factor_base_size': len(base), 'diagonal_hex': hex(diagonal), 'row_witness_vectors': vectors,
        'relation_matrix_sha256': matrix_hash(rows),
        'relation_determinant_hex': hex(annihilator), 'determinant_bits': annihilator.bit_length(),
        'determinant_two_valuation': (annihilator & -annihilator).bit_length()-1,
        'extracted_two_group_basis': [[[str(value) for value in form], order] for form, order in basis],
        'recovered_factor': str(factor) if factor else None,
    }


def verify(certificate):
    integer, bound = int(certificate['integer']), certificate['smooth_bound']
    radius = isqrt(4*integer-1)//2
    identity, diagonal = (1, 0, integer), int(certificate['diagonal_hex'], 16)
    base = prime_base(-4*integer, bound)
    vectors = certificate['row_witness_vectors']
    assert len(vectors) == len(base) == certificate['factor_base_size']
    rows = []
    for index, ((_, _, generator), (x, y)) in enumerate(zip(base, vectors)):
        assert 0 <= x <= certificate['x_bound'] and 0 <= y <= certificate['y_bound']
        target = power(generator, diagonal, identity, compose_forms)
        pair = norm_to_pair(target, x, y, bound, radius)
        assert pair is not None
        row = relation_row(index, diagonal, pair, base)
        value = identity
        for exponent, (_, _, form) in zip(row, base):
            value = compose_forms(value, power(form if exponent >= 0 else inverse_form(form), abs(exponent), identity, compose_forms))
        assert value == identity
        rows.append(row)
    assert matrix_hash(rows) == certificate['relation_matrix_sha256']
    annihilator = abs(determinant(rows))
    assert hex(annihilator) == certificate['relation_determinant_hex']
    assert annihilator.bit_length() == certificate['determinant_bits']
    assert (annihilator & -annihilator).bit_length()-1 == certificate['determinant_two_valuation']
    basis, factor = extraction(base, annihilator, integer)
    encoded = [[[str(value) for value in form], order] for form, order in basis]
    assert encoded == certificate['extracted_two_group_basis']
    assert factor is not None and str(factor) == certificate['recovered_factor'] and integer % factor == 0
    return {'integer': str(integer), 'verified_relations': len(rows), 'determinant_bits': annihilator.bit_length(),
            'determinant_two_valuation': certificate['determinant_two_valuation'],
            'two_group_basis_order_logs': [order for _, order in basis], 'proper_factor': str(factor)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('integer', type=int)
    parser.add_argument('--smooth-bound', type=int, default=1000)
    parser.add_argument('--x-bound', type=int, default=4096)
    parser.add_argument('--y-bound', type=int, default=64)
    arguments = parser.parse_args()
    print(json.dumps(generate(arguments.integer, arguments.smooth_bound, arguments.x_bound, arguments.y_bound), indent=2))
