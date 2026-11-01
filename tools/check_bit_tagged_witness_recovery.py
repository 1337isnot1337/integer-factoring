#!/usr/bin/env python3
"""Exact finite checks for bit-tagged deterministic witness recovery.

The finite evaluator is the explicit multilinear polynomial of a supplied
witness list. It checks the recovery transformations, not a new fast
factoring circuit. All production recovery loops are sequential.
"""

from itertools import combinations, product
from math import isqrt
import json


def ceil_div(first, second):
    return (first + second - 1) // second


def ceil_log2(value):
    return (value - 1).bit_length()


def tag(number, width, prime):
    return (1,) + tuple((number >> bit) & 1 for bit in range(width))


def add(first, second, prime):
    return tuple((a + b) % prime for a, b in zip(first, second))


def scale(value, scalar, prime):
    return tuple(scalar * coordinate % prime for coordinate in value)


def multiply(first, second, prime):
    return (first[0] * second[0] % prime,) + tuple(
        (first[0] * b + second[0] * a) % prime
        for a, b in zip(first[1:], second[1:]))


def evaluate(witnesses, variables, prime):
    width = len(variables)
    result = (0,) * (width + 1)
    one = (1,) + (0,) * width
    for witness in witnesses:
        term = one
        for bit in range(width):
            if (witness >> bit) & 1:
                term = multiply(term, variables[bit], prime)
        result = add(result, term, prime)
    return result


def decode(vector):
    if vector[0] != 1 or any(value not in (0, 1) for value in vector[1:]):
        return None
    return sum(value << bit for bit, value in enumerate(vector[1:]))


def direct_tag_sum(witnesses, width, prime):
    result = (0,) * (width + 1)
    for witness in witnesses:
        result = add(result, tag(witness, width, prime), prime)
    return result


def cylinder_value(witnesses, width, zeros, ones, prime):
    """Inclusion-exclusion using only substitutions in the original C_S."""
    result = (0,) * (width + 1)
    one = (1,) + (0,) * width
    zero = (0,) * (width + 1)
    ones = tuple(ones)
    for selected in range(1 << len(ones)):
        forbidden = set(zeros)
        forbidden.update(ones[index] for index in range(len(ones))
                         if (selected >> index) & 1)
        variables = [
            zero if bit in forbidden else add(one, (0,) + tuple(
                int(index == bit) for index in range(width)), prime)
            for bit in range(width)
        ]
        value = evaluate(witnesses, variables, prime)
        result = add(result, scale(value, (-1) ** selected.bit_count(), prime), prime)
    return result


def cylinders(width, maximum):
    for size in range(min(width, maximum) + 1):
        for coordinates in combinations(range(width), size):
            for assignment in range(1 << size):
                ones = tuple(coordinates[index] for index in range(size)
                             if (assignment >> index) & 1)
                zeros = tuple(coordinates[index] for index in range(size)
                              if not ((assignment >> index) & 1))
                yield zeros, ones


def cylinder_witness(witnesses, width, bound_bits, prime, verify_all=False):
    checked = 0
    first = None
    witness_set = set(witnesses)
    for zeros, ones in cylinders(width, bound_bits):
        result = cylinder_value(witnesses, width, zeros, ones, prime)
        expected = [witness for witness in witnesses
                    if all(not ((witness >> bit) & 1) for bit in zeros)
                    and all((witness >> bit) & 1 for bit in ones)]
        assert result == direct_tag_sum(expected, width, prime)
        checked += 1
        candidate = decode(result)
        if candidate is not None and candidate in witness_set:
            first = candidate if first is None else min(first, candidate)
            if not verify_all:
                break
    return first, checked


def polynomial_remainder(value, divisor):
    while value and value.bit_length() >= divisor.bit_length():
        value ^= divisor << (value.bit_length() - divisor.bit_length())
    return value


class BinaryField:
    def __init__(self, degree):
        self.degree, self.size = degree, 1 << degree
        self.polynomial = next(
            candidate
            for candidate in range(self.size | 1, 2 * self.size, 2)
            if all(polynomial_remainder(candidate, divisor)
                   for part_degree in range(1, degree // 2 + 1)
                   for divisor in range(1 << part_degree, 1 << (part_degree + 1)))
        )

    def multiply(self, first, second):
        value = 0
        while second:
            if second & 1:
                value ^= first
            second >>= 1
            first <<= 1
            if first & self.size:
                first ^= self.polynomial
        return value

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            value = self.multiply(value, value)
            exponent >>= 1
        return result

    def hash_basis(self, point, width):
        powers = [1]
        for _ in range(ceil_div(width, self.degree)):
            powers.append(self.multiply(powers[-1], point))
        return tuple(self.multiply(1 << (bit % self.degree),
                                   powers[bit // self.degree])
                     for bit in range(width))


def field_for_bound(width, bound):
    degree = 1
    while (1 << degree) <= (bound - 1) * ceil_div(width, degree):
        degree += 1
    return BinaryField(degree)


def hash_horner(word, point, field, width):
    chunks = ceil_div(width, field.degree)
    value = 0
    for index in range(chunks - 1, -1, -1):
        value = field.multiply(value, point)
        value ^= (word >> (index * field.degree)) & (field.size - 1)
    return value


def hash_word(word, basis):
    result = 0
    for bit, value in enumerate(basis):
        if (word >> bit) & 1:
            result ^= value
    return result


def walsh_bucket(witnesses, basis, bucket, prime, field_size):
    width = len(basis)
    size = field_size
    result = (0,) * (width + 1)
    one = (1,) + (0,) * width
    for character in range(size):
        variables = []
        for bit, image in enumerate(basis):
            value = add(one, (0,) + tuple(
                int(index == bit) for index in range(width)), prime)
            sign = -1 if (character & image).bit_count() % 2 else 1
            variables.append(scale(value, sign, prime))
        value = evaluate(witnesses, variables, prime)
        sign = -1 if (character & bucket).bit_count() % 2 else 1
        result = add(result, scale(value, sign, prime), prime)
    return scale(result, pow(size, -1, prime), prime)


class PrimeExtension:
    """Polynomial-basis finite field, constructed by exhaustive factor testing."""

    def __init__(self, prime, degree):
        self.prime, self.degree = prime, degree
        self.size = prime ** degree

        def remainder(value, divisor):
            value = list(value)
            for position in range(len(value) - 1, len(divisor) - 2, -1):
                scalar = value[position] % prime
                offset = position - len(divisor) + 1
                for index, coefficient in enumerate(divisor):
                    value[offset + index] = (value[offset + index] - scalar * coefficient) % prime
            return any(value[:len(divisor) - 1])

        self.polynomial = next(
            self.coordinates(candidate, degree) + (1,)
            for candidate in range(self.size)
            if all(remainder(self.coordinates(candidate, degree) + (1,),
                             self.coordinates(divisor, part) + (1,))
                   for part in range(1, degree // 2 + 1)
                   for divisor in range(prime ** part)))

    def coordinates(self, value, length=None):
        length = self.degree if length is None else length
        return tuple(value // self.prime ** index % self.prime for index in range(length))

    def encode(self, coefficients):
        return sum((value % self.prime) * self.prime ** index
                   for index, value in enumerate(coefficients))

    def add(self, first, second):
        if self.prime == 2:
            return first ^ second
        return self.encode(a + b for a, b in zip(self.coordinates(first), self.coordinates(second)))

    def multiply(self, first, second):
        if self.degree == 1:
            return first * second % self.prime
        coefficients = [0] * (2 * self.degree - 1)
        for i, a in enumerate(self.coordinates(first)):
            for j, b in enumerate(self.coordinates(second)):
                coefficients[i + j] += a * b
        for position in range(len(coefficients) - 1, self.degree - 1, -1):
            scalar = coefficients[position] % self.prime
            for index, value in enumerate(self.polynomial):
                coefficients[position - self.degree + index] -= scalar * value
        return self.encode(coefficients[:self.degree])

    def power(self, value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            value = self.multiply(value, value)
            exponent >>= 1
        return result

    def dot(self, first, second):
        return sum(a * b for a, b in zip(self.coordinates(first), self.coordinates(second))) % self.prime

    def bit_basis(self, point, width):
        result, value = [], 1
        for _ in range(width):
            result.append(value)
            value = self.multiply(value, point)
        return tuple(result)

    def hash_word(self, witness, basis):
        value = 0
        for bit, image in enumerate(basis):
            if witness >> bit & 1:
                value = self.add(value, image)
        return value


def extension_field_evaluator(witnesses, variables, field, add_table, multiply_table):
    """Evaluate the original multilinear polynomial in square-zero coordinates."""
    width = len(variables)
    result = [0] * (width + 1)
    for witness in witnesses:
        term = [1] + [0] * width
        for bit, variable in enumerate(variables):
            if witness >> bit & 1:
                scalar = multiply_table[term[0]][variable[0]]
                tags = [add_table[multiply_table[term[0]][variable[index]]]
                                [multiply_table[variable[0]][term[index]]]
                        for index in range(1, width + 1)]
                term = [scalar] + tags
        result = [add_table[a][b] for a, b in zip(result, term)]
    return tuple(result)


def additive_bucket(witnesses, basis, bucket, control, coefficient, root,
                    add_table, multiply_table):
    width = len(basis)
    result = [0] * (width + 1)
    normalizer = pow(control.size % coefficient.prime, -1, coefficient.prime)
    phases = [coefficient.power(root, exponent) for exponent in range(control.prime)]
    for character in range(control.size):
        variables = []
        for bit, image in enumerate(basis):
            phase = phases[control.dot(character, image)]
            variables.append((phase,) + tuple(phase if index == bit else 0
                                              for index in range(width)))
        value = extension_field_evaluator(
            witnesses, variables, coefficient, add_table, multiply_table)
        factor = multiply_table[normalizer][phases[-control.dot(character, bucket) % control.prime]]
        result = [add_table[old][multiply_table[factor][new]] for old, new in zip(result, value)]
    return tuple(result)


def odd_primes(count):
    result, candidate = [], 3
    while len(result) < count:
        if all(candidate % divisor
               for divisor in range(3, isqrt(candidate) + 1, 2)):
            result.append(candidate)
        candidate += 2
    return result


def cyclic_recover(witnesses, width, prime):
    if not witnesses:
        return set(), 0
    count = width * (len(witnesses) - 1) + 1
    accepted, buckets_seen = set(), 0
    witness_set = set(witnesses)
    for modulus in odd_primes(count):
        buckets = {}
        for witness in witnesses:
            residue = witness % modulus
            buckets[residue] = add(
                buckets.get(residue, (0,) * (width + 1)),
                tag(witness, width, prime), prime)
        for value in buckets.values():
            buckets_seen += 1
            candidate = decode(value)
            if candidate in witness_set:
                accepted.add(candidate)
    return accepted, buckets_seen


def factor_with_cylinders(number, prime, parameter=None):
    if number == 1:
        return [], 0
    width = number.bit_length()
    logarithm = ceil_log2(width + 2)
    bound_bits = (max(1, isqrt(width // logarithm)) if parameter is None
                  else min(width, parameter))
    threshold_bits = ceil_div(width, bound_bits)
    factors, tested = [], 0
    candidate = 2
    while candidate <= min(1 << threshold_bits, number):
        while number % candidate == 0:
            factors.append(candidate)
            number //= candidate
        candidate += 1
    if number == 1:
        return sorted(factors), tested

    pending = [number]
    while pending:
        current = pending.pop()
        proper = [value for value in range(2, current)
                  if current % value == 0]
        assert len(proper) <= 1 << bound_bits
        found, calls = cylinder_witness(proper, width, bound_bits, prime)
        tested += calls
        if found is None:
            assert all(current % divisor
                       for divisor in range(2, isqrt(current) + 1))
            factors.append(current)
        else:
            assert 1 < found < current and current % found == 0
            pending.extend((found, current // found))
    return sorted(factors), tested


def verify():
    algebra_cases = 0
    for prime in (2, 3, 5):
        for width in range(1, 7):
            one = (1,) + (0,) * width
            variables = [add(one, (0,) + tuple(
                int(index == bit) for index in range(width)), prime)
                for bit in range(width)]
            for witness in range(1 << width):
                assert evaluate([witness], variables, prime) == tag(witness, width, prime)
                algebra_cases += 1

    cylinder_cases = cylinder_entries = 0
    for mask in range(1, 1 << 8):
        witnesses = [value for value in range(8) if (mask >> value) & 1]
        bound_bits = ceil_log2(len(witnesses))
        for prime in (2, 3, 5):
            found, checked = cylinder_witness(
                witnesses, 3, bound_bits, prime, verify_all=True)
            assert found in witnesses
            cylinder_cases += 1
            cylinder_entries += checked
    for prime in (2, 3, 5):
        found, _ = cylinder_witness([], 4, 2, prime, verify_all=True)
        assert found is None

    # A small cylinder need not isolate each prescribed point.
    star = [0, 1, 2, 4, 8]
    isolated_zero = False
    for zeros, ones in cylinders(4, 3):
        subset = [value for value in star
                  if all(not ((value >> bit) & 1) for bit in zeros)
                  and all((value >> bit) & 1 for bit in ones)]
        isolated_zero |= subset == [0]
    assert not isolated_zero

    cyclic_cases = cyclic_buckets = 0
    for number in range(2, 256):
        witnesses = [value for value in range(1, number + 1) if number % value == 0]
        for prime in (2, 3, 5):
            accepted, buckets = cyclic_recover(witnesses, number.bit_length(), prime)
            assert accepted == set(witnesses)
            cyclic_cases += 1
            cyclic_buckets += buckets

    hash_cases = projection_cases = 0
    fields = {}
    examples = [
        (7, [1, 18, 67]),
        (7, [3, 35, 99]),
        (7, [0, 64, 96]),
        (8, [1, 17, 65, 129]),
        (8, [5, 37, 101, 245]),
    ]
    for width, witnesses in examples:
        field = field_for_bound(width, len(witnesses))
        fields[field.degree] = field.polynomial
        isolated, recovered = set(), {prime: set() for prime in (3, 5, 7)}
        collision_counts = {pair: 0 for pair in combinations(witnesses, 2)}
        for point in range(field.size):
            basis = field.hash_basis(point, width)
            assert max(basis).bit_length() == field.degree
            for witness in witnesses:
                assert hash_word(witness, basis) == hash_horner(witness, point, field, width)
            for first, second in collision_counts:
                collision_counts[first, second] += (
                    hash_word(first, basis) == hash_word(second, basis))
            buckets = {}
            for witness in witnesses:
                buckets.setdefault(hash_word(witness, basis), []).append(witness)
            for values in buckets.values():
                if len(values) == 1:
                    isolated.add(values[0])
            for bucket in range(field.size):
                for prime in (3, 5, 7):
                    value = walsh_bucket(witnesses, basis, bucket, prime, field.size)
                    assert value == direct_tag_sum(buckets.get(bucket, []), width, prime)
                    candidate = decode(value)
                    if candidate in witnesses:
                        recovered[prime].add(candidate)
                    projection_cases += 1
            hash_cases += 1
        assert isolated == set(witnesses)
        assert all(found == set(witnesses) for found in recovered.values())
        assert all(count <= ceil_div(width, field.degree) - 1
                   for count in collision_counts.values())

    additive_parameters = additive_buckets = additive_characteristic_two = 0
    additive_field_records = []
    additive_examples = [(3, [1, 3, 6]), (4, [1, 5, 9, 14]), (5, [0, 17])]
    for prime, control_prime, coefficient_degree in ((2, 3, 2), (3, 2, 1), (5, 3, 2)):
        coefficient = PrimeExtension(prime, coefficient_degree)
        root = next(value for value in range(2, coefficient.size)
                    if coefficient.power(value, control_prime) == 1)
        assert root != 1
        # These tables have fixed size, independent of n. They implement only the
        # constant-degree coefficient extension, never the large control field.
        add_table = [[coefficient.add(a, b) for b in range(coefficient.size)]
                     for a in range(coefficient.size)]
        multiply_table = [[coefficient.multiply(a, b) for b in range(coefficient.size)]
                          for a in range(coefficient.size)]
        for width, witnesses in additive_examples:
            degree = 1
            while control_prime ** degree <= width * (len(witnesses) - 1):
                degree += 1
            control = PrimeExtension(control_prime, degree)
            additive_field_records.append({
                "coefficient_prime": prime,
                "coefficient_polynomial": coefficient.polynomial,
                "root": root,
                "control_prime": control_prime,
                "control_degree": degree,
                "control_polynomial": control.polynomial,
            })
            isolated, recovered = set(), set()
            collisions = {pair: 0 for pair in combinations(witnesses, 2)}
            for point in range(control.size):
                basis = control.bit_basis(point, width)
                buckets = {}
                for witness in witnesses:
                    # Compare the linear bit-basis implementation with Horner.
                    horner = 0
                    for bit in range(width - 1, -1, -1):
                        horner = control.add(control.multiply(horner, point), witness >> bit & 1)
                    image = control.hash_word(witness, basis)
                    assert image == horner
                    buckets.setdefault(image, []).append(witness)
                for pair in collisions:
                    collisions[pair] += control.hash_word(pair[0], basis) == (
                        control.hash_word(pair[1], basis))
                for members in buckets.values():
                    if len(members) == 1:
                        isolated.add(members[0])
                for bucket in range(control.size):
                    value = additive_bucket(witnesses, basis, bucket, control,
                                            coefficient, root, add_table, multiply_table)
                    assert value == direct_tag_sum(buckets.get(bucket, []), width, prime)
                    candidate = decode(value)
                    if candidate in witnesses:
                        recovered.add(candidate)
                    additive_buckets += 1
                    additive_characteristic_two += prime == 2
                additive_parameters += 1
            assert isolated == recovered == set(witnesses)
            assert all(count <= width - 1 for count in collisions.values())

    balancing_cases = 0
    for width in range(2, 100001):
        logarithm = ceil_log2(width + 2)
        bound_bits = max(1, isqrt(width // logarithm))
        threshold = ceil_div(width, bound_bits)
        root = isqrt(width * logarithm)
        root += root * root < width * logarithm
        assert threshold + bound_bits * logarithm <= 3 * root + 1
        balancing_cases += 1

    recovery_cases = recovery_cylinders = 0
    for number in range(2, 256):
        for prime in (2, 3, 5):
            factors, count = factor_with_cylinders(number, prime, parameter=3)
            value = 1
            for factor in factors:
                assert all(factor % divisor for divisor in range(2, isqrt(factor) + 1))
                value *= factor
            assert value == number
            recovery_cases += 1
            recovery_cylinders += count
    # Force nontrivial rough residuals at somewhat larger finite inputs.
    for number in (17 * 19, 17 ** 2, 19 * 23, 257):
        factors, count = factor_with_cylinders(number, 3, parameter=3)
        value = 1
        for factor in factors:
            value *= factor
        assert value == number
        recovery_cases += 1
        recovery_cylinders += count

    # Integer exponentiation of one tag does not encode digit bits.
    epsilon_zero = (1, 1, 0, 0)
    wrong_digit_tag = multiply(epsilon_zero, epsilon_zero, 2)
    assert wrong_digit_tag != tag(2, 3, 2)
    # Four distinct witnesses can imitate a one-witness tag in F_3.
    collision = direct_tag_sum([1, 2, 4, 8], 4, 3)
    assert decode(collision) == 15 and 15 not in (1, 2, 4, 8)
    # Binary Walsh signs coincide in characteristic two and their sum is zero.
    assert sum(((-1) ** ((value & 1).bit_count() % 2)) % 2
               for value in range(8)) % 2 == 0

    return {
        "scope": "exact recovery transforms with explicit finite witness-polynomial evaluator",
        "square_zero_tag_cases": algebra_cases,
        "cylinder_witness_sets": cylinder_cases,
        "cylinder_inclusion_exclusion_entries": cylinder_entries,
        "cyclic_all_divisor_cases": cyclic_cases,
        "cyclic_buckets_examined": cyclic_buckets,
        "binary_polynomial_hash_parameters": hash_cases,
        "odd_characteristic_walsh_buckets": projection_cases,
        "constructed_binary_fields": fields,
        "additive_hash_parameters": additive_parameters,
        "different_characteristic_additive_buckets": additive_buckets,
        "characteristic_two_additive_buckets": additive_characteristic_two,
        "constructed_additive_fields": additive_field_records,
        "complete_factor_recovery_cases": recovery_cases,
        "factor_recovery_cylinders": recovery_cylinders,
        "polynomial_space_time_balancing_cases": balancing_cases,
        "negative_controls": {
            "digit_power_tagging_is_wrong": True,
            "collision_candidate_needs_verification": 15,
            "small_cylinders_need_not_isolate_every_witness": True,
            "characteristic_two_walsh_projection_collapses": True,
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
