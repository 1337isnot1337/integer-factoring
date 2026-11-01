#!/usr/bin/env python3
"""Exact finite checks for the growing-number-field transfer.

This verifies the coefficient-ring construction and short multinomial
algorithms, not an asymptotic factoring improvement or a complete factorer.
Only Python's standard library is used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from math import factorial, gcd, isqrt, prod
from pathlib import Path
import random


class TensorQuotient:
    """Dense, exact tensor quotient by separate monic polynomials.

    A None modulus gives rational coordinates. Polynomials are low-degree
    first, with the final monic coefficient included. All basis products
    are compiled independently by repeated monic polynomial reduction.
    """

    def __init__(self, polynomials, modulus=None):
        self.polynomials = tuple(tuple(f) for f in polynomials)
        assert self.polynomials and all(len(f) >= 2 and f[-1] == 1
                                        for f in self.polynomials)
        self.modulus = modulus
        assert modulus is None or modulus >= 2
        self.degrees = tuple(len(f) - 1 for f in self.polynomials)
        self.basis = tuple(product(*(range(e) for e in self.degrees)))
        self.position = {powers: i for i, powers in enumerate(self.basis)}
        self.dimension = prod(self.degrees)
        self.zero = (0,) * self.dimension
        self.one = (1,) + (0,) * (self.dimension - 1)
        self.basis_products = {}
        reductions = []
        for f in self.polynomials:
            e = len(f) - 1
            per_power = []
            for power in range(2 * e - 1):
                coefficients = [0] * (max(power, e - 1) + 1)
                coefficients[power] = 1
                for leading in range(len(coefficients) - 1, e - 1, -1):
                    value = coefficients[leading]
                    for k in range(e + 1):
                        coefficients[leading - e + k] -= value * f[k]
                per_power.append(tuple((i, a) for i, a in enumerate(coefficients[:e]) if a))
            reductions.append(per_power)
        for i, first in enumerate(self.basis):
            for j, second in enumerate(self.basis):
                options = [reductions[t][a + b]
                           for t, (a, b) in enumerate(zip(first, second))]
                entries = []
                for terms in product(*options):
                    powers = tuple(term[0] for term in terms)
                    value = prod(term[1] for term in terms)
                    if value:
                        entries.append((self.position[powers], value))
                self.basis_products[i, j] = tuple(entries)

    def reduce_scalar(self, value):
        value = Fraction(value)
        if self.modulus is None:
            return value
        assert gcd(value.denominator, self.modulus) == 1
        return value.numerator * pow(value.denominator, -1, self.modulus) % self.modulus

    def scalar(self, value):
        return (self.reduce_scalar(value),) + (0,) * (self.dimension - 1)

    def add(self, a, b):
        return tuple(self.reduce_scalar(x + y) for x, y in zip(a, b))

    def neg(self, a):
        return tuple(self.reduce_scalar(-x) for x in a)

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        output = [0] * self.dimension
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if not y:
                    continue
                for k, z in self.basis_products[i, j]:
                    output[k] += x * y * z
        return tuple(self.reduce_scalar(x) for x in output)

    def power(self, a, exponent):
        result = self.one
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            exponent //= 2
        return result

    def variable(self, t):
        e = self.degrees[t]
        if e == 1:
            return self.scalar(-self.polynomials[t][0])
        powers = [0] * len(self.degrees)
        powers[t] = 1
        result = [0] * self.dimension
        result[self.position[tuple(powers)]] = 1
        return tuple(result)

    def embed(self, local_element, t):
        assert len(local_element) == self.degrees[t]
        result = [0] * self.dimension
        for j, coefficient in enumerate(local_element):
            powers = [0] * len(self.degrees)
            powers[t] = j
            result[self.position[tuple(powers)]] = self.reduce_scalar(coefficient)
        return tuple(result)


def inverse_in_local_field(field, value):
    """Independent rational Gaussian elimination for one small field inverse."""
    assert len(field.polynomials) == 1 and field.modulus is None
    e = field.dimension
    columns = []
    for j in range(e):
        vector = tuple(int(i == j) for i in range(e))
        columns.append(field.mul(value, vector))
    matrix = [[Fraction(columns[j][i]) for j in range(e)] + [Fraction(i == 0)]
              for i in range(e)]
    for column in range(e):
        pivot = next(row for row in range(column, e) if matrix[row][column])
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        divisor = matrix[column][column]
        matrix[column] = [x / divisor for x in matrix[column]]
        for row in range(e):
            if row != column:
                multiple = matrix[row][column]
                matrix[row] = [x - multiple * y
                               for x, y in zip(matrix[row], matrix[column])]
    result = tuple(row[-1] for row in matrix)
    assert field.mul(value, result) == field.one
    return result


def make_local_pieces(polynomial, label_count=2):
    """A scale-normalized crossing with algebraic, genuinely fractional factors."""
    field = TensorQuotient([polynomial])
    x = field.variable(0)
    scale = field.add(field.one, x)
    inverse = inverse_in_local_field(field, scale)
    vectors = []
    for a, b in ((2, 1), (3, -1), (1, 2), (-2, 3)):
        values = [
            field.add(field.scalar(Fraction(a, 3)), field.mul(field.scalar(Fraction(b, 5)), x)),
            field.add(field.scalar(Fraction(b, 5)), field.mul(field.scalar(Fraction(a, 3)), x)),
        ]
        for j in range(2, label_count):
            values.append(field.add(field.scalar(Fraction(a + j, 3)),
                field.mul(field.scalar(Fraction(b + 2 * j, 5)), x)))
        vectors.append(values)
    a, b, c, d = vectors
    tensor = {}
    row_left, row_right, col_left, col_right = {}, {}, {}, {}
    for first, second in product(range(label_count), repeat=2):
        row_left[first, second, 0] = field.mul(inverse, field.mul(a[first], b[second]))
        row_right[first, second, 0] = field.mul(c[first], d[second])
        col_left[first, second, 0] = field.mul(inverse, field.mul(a[first], c[second]))
        col_right[first, second, 0] = field.mul(b[first], d[second])
    for x0, y0, u0, v0 in product(range(label_count), repeat=4):
        tensor[x0, y0, u0, v0] = field.mul(row_left[x0, y0, 0], row_right[u0, v0, 0])
    first = dict(tensor=tensor, row_channels=1, col_channels=1,
                 row_left=row_left, row_right=row_right,
                 col_left=col_left, col_right=col_right)
    residual = {}
    for labels, value in tensor.items():
        x0, y0, u0, v0 = labels
        residual[labels] = field.sub(field.scalar(x0 == v0 and y0 == u0), value)
    row_left, row_right, col_left, col_right = {}, {}, {}, {}
    for a0, b0, channel in product(range(label_count), range(label_count), range(label_count ** 2)):
        row_left[a0, b0, channel] = field.scalar(channel == label_count * a0 + b0)
        col_left[a0, b0, channel] = field.scalar(channel == label_count * a0 + b0)
        row_right[a0, b0, channel] = residual[channel // label_count, channel % label_count, a0, b0]
        col_right[a0, b0, channel] = residual[channel // label_count, a0, channel % label_count, b0]
    second = dict(tensor=residual, row_channels=label_count ** 2, col_channels=label_count ** 2,
                  row_left=row_left, row_right=row_right,
                  col_left=col_left, col_right=col_right)
    return field, (first, second)


def piece_checks(ring, pieces, label_count=2):
    factorization_entries = crossing_entries = 0
    for labels in product(range(label_count), repeat=4):
        x, y, u, v = labels
        total = ring.zero
        for piece in pieces:
            row = ring.zero
            for channel in range(piece['row_channels']):
                row = ring.add(row, ring.mul(piece['row_left'][x, y, channel],
                                             piece['row_right'][u, v, channel]))
            column = ring.zero
            for channel in range(piece['col_channels']):
                column = ring.add(column, ring.mul(piece['col_left'][x, u, channel],
                                                   piece['col_right'][y, v, channel]))
            assert row == column == piece['tensor'][labels]
            total = ring.add(total, row)
            factorization_entries += 2
        assert total == ring.scalar(x == v and y == u)
        crossing_entries += 1
    return factorization_entries, crossing_entries


def embed_pieces(ring, pieces, t):
    result = []
    for piece in pieces:
        new = {key: value for key, value in piece.items() if key.endswith('channels')}
        for name in ('tensor', 'row_left', 'row_right', 'col_left', 'col_right'):
            new[name] = {key: ring.embed(value, t) for key, value in piece[name].items()}
        result.append(new)
    return tuple(result)


def algebra_checks():
    rng = random.Random(20261030)
    totals = dict(local_factorization_entries=0, local_crossing_entries=0,
                  reduced_factorization_entries=0, reduced_crossing_entries=0,
                  homomorphism_products=0, mixed_size_crossing_entries=0,
                  repeated_field_zero_divisors=0, defining_polynomial_checks=0,
                  forbidden_denominator_controls=0)
    records = []
    # X^e-2 is Eisenstein at two. The second case repeats one field;
    # the tensor quotient still has two independent variables.
    for degrees in ((2, 3, 4), (2, 3, 2)):
        polynomials = [(-2,) + (0,) * (e - 1) + (1,) for e in degrees]
        local = [make_local_pieces(f) for f in polynomials]
        denominator = 1
        for field, pieces in local:
            checks = piece_checks(field, pieces)
            totals['local_factorization_entries'] += checks[0]
            totals['local_crossing_entries'] += checks[1]
            for piece in pieces:
                for name in ('row_left', 'row_right', 'col_left', 'col_right'):
                    for value in piece[name].values():
                        for coefficient in value:
                            denominator *= Fraction(coefficient).denominator
        assert denominator % 3 == 0 and denominator % 5 == 0
        moduli = (49, 77, 121, 2 * denominator + 1)
        records.append(dict(degrees=degrees, tensor_dimension=prod(degrees),
                            denominator_bits=denominator.bit_length(),
                            constructed_modulus_bits=(2 * denominator + 1).bit_length(),
                            small_composite_moduli=moduli[:-1]))
        for modulus in moduli:
            assert gcd(denominator, modulus) == 1
            ring = TensorQuotient(polynomials, modulus)
            assert ring.one != ring.zero
            assert ring.mul(ring.scalar(2), ring.scalar(Fraction(1, 2))) == ring.one
            embedded = []
            for t, (field, pieces) in enumerate(local):
                reduced = embed_pieces(ring, pieces, t)
                embedded.append(reduced)
                checks = piece_checks(ring, reduced)
                totals['reduced_factorization_entries'] += checks[0]
                totals['reduced_crossing_entries'] += checks[1]
                theta = ring.variable(t)
                assert ring.power(theta, degrees[t]) == ring.scalar(2)
                totals['defining_polynomial_checks'] += 1
                for _ in range(24):
                    first = tuple(Fraction(rng.randrange(-4, 5), 15) for _ in range(degrees[t]))
                    second = tuple(Fraction(rng.randrange(-4, 5), 15) for _ in range(degrees[t]))
                    assert ring.embed(field.mul(first, second), t) == ring.mul(
                        ring.embed(first, t), ring.embed(second, t))
                    totals['homomorphism_products'] += 1
            # Two blocks from different source fields share one ring.
            for x, y, u, v in product(range(4), repeat=4):
                total = ring.zero
                first_labels = tuple(z // 2 for z in (x, y, u, v))
                second_labels = tuple(z % 2 for z in (x, y, u, v))
                for first, second in product(embedded[0], embedded[1]):
                    total = ring.add(total, ring.mul(first['tensor'][first_labels],
                                                     second['tensor'][second_labels]))
                assert total == ring.scalar(x == v and y == u)
                totals['mixed_size_crossing_entries'] += 1
            if degrees[0] == degrees[2]:
                difference = ring.sub(ring.variable(0), ring.variable(2))
                sum_value = ring.add(ring.variable(0), ring.variable(2))
                assert difference != ring.zero and sum_value != ring.zero
                assert ring.mul(difference, sum_value) == ring.zero
                totals['repeated_field_zero_divisors'] += 1
        # Reducing a denominator divisible by three is deliberately rejected.
        bad_ring = TensorQuotient(polynomials, 9)
        try:
            bad_ring.scalar(Fraction(1, 3))
        except AssertionError:
            totals['forbidden_denominator_controls'] += 1
        else:
            raise AssertionError('Nonunit denominator was accepted')
    return dict(totals=totals, cases=records)


def tail_checks():
    raw_factorizations = capped_factorizations = crossing_entries = 0
    wrong_padding_controls = 0
    for degree in (2, 3, 4):
        polynomial = (-2,) + (0,) * (degree - 1) + (1,)
        field, raw_pieces = make_local_pieces(polynomial, label_count=4)
        checked, _ = piece_checks(field, raw_pieces, label_count=4)
        raw_factorizations += checked
        for modulus in (49, 77):
            ring = TensorQuotient([polynomial], modulus)
            raw = embed_pieces(ring, raw_pieces, 0)
            restricted = []
            for original in raw:
                # One trailing zero is appended to every external label.
                tensor = {labels: original['tensor'][tuple(2 * x for x in labels)]
                          for labels in product(range(2), repeat=4)}
                row_left, row_right, col_left, col_right = {}, {}, {}, {}
                for a, b, channel in product(range(2), range(2), range(4)):
                    row_left[a, b, channel] = ring.scalar(channel == 2 * a + b)
                    col_left[a, b, channel] = ring.scalar(channel == 2 * a + b)
                    row_right[a, b, channel] = tensor[channel // 2, channel % 2, a, b]
                    col_right[a, b, channel] = tensor[channel // 2, a, channel % 2, b]
                restricted.append(dict(tensor=tensor, row_channels=4, col_channels=4,
                    row_left=row_left, row_right=row_right,
                    col_left=col_left, col_right=col_right))
            checked, crossings = piece_checks(ring, restricted)
            capped_factorizations += checked
            crossing_entries += crossings
            mismatched = ring.zero
            for piece in raw:
                mismatched = ring.add(mismatched, piece['tensor'][2, 0, 0, 1])
            assert mismatched == ring.zero != ring.one
            wrong_padding_controls += 1
    return dict(raw_factorization_entries=raw_factorizations,
                restricted_and_capped_factorization_entries=capped_factorizations,
                restricted_crossing_entries=crossing_entries,
                inconsistent_padding_negative_controls=wrong_padding_controls,
                raw_label_bits=2, tail_label_bits=1, same_raw_field_reused=True)


def short_multinomial(counts):
    """Return exact count and largest arithmetic intermediate, without factorials."""
    total = sum(counts)
    remaining, result, largest = total, 1, max(1, total)
    for count in counts:
        choose = min(count, remaining - count)
        coefficient = 1
        for j in range(1, choose + 1):
            numerator = coefficient * (remaining - j + 1)
            largest = max(largest, numerator)
            assert numerator % j == 0
            coefficient = numerator // j
        result *= coefficient
        largest = max(largest, result)
        remaining -= count
    assert remaining == 0
    assert largest <= result * max(1, total)
    return result, largest


def unrank_type_word(counts, ordinal):
    counts = list(counts)
    remaining = sum(counts)
    number, largest = short_multinomial(counts)
    assert 0 <= ordinal < number
    result = []
    while remaining:
        for letter, count in enumerate(counts):
            if not count:
                continue
            numerator = number * count
            largest = max(largest, numerator)
            assert numerator % remaining == 0
            first_count = numerator // remaining
            if ordinal >= first_count:
                ordinal -= first_count
            else:
                result.append(letter)
                counts[letter] -= 1
                remaining -= 1
                number = first_count
                break
        else:
            raise AssertionError('No first letter was selected')
    assert ordinal == 0
    return tuple(result), largest


def multinomial_checks():
    types = words = 0
    largest_bits = 0
    for alphabet in range(1, 5):
        for length in range(8):
            reference = {}
            for word in product(range(alphabet), repeat=length):
                counts = tuple(word.count(i) for i in range(alphabet))
                reference.setdefault(counts, []).append(word)
            for counts, expected_words in reference.items():
                number, largest = short_multinomial(counts)
                expected = factorial(length) // prod(factorial(x) for x in counts)
                assert number == expected == len(expected_words)
                assert largest <= number * max(1, length)
                largest_bits = max(largest_bits, largest.bit_length())
                types += 1
                for ordinal, expected_word in enumerate(expected_words):
                    actual, largest = unrank_type_word(counts, ordinal)
                    assert actual == expected_word
                    assert largest <= number * max(1, length)
                    words += 1
    # A long exact control distinguishes the short count from a factorial.
    counts = (2048, 2048)
    number, largest = short_multinomial(counts)
    factorial_bits = factorial(4096).bit_length()
    assert factorial_bits > 10 * number.bit_length()
    assert largest.bit_length() <= number.bit_length() + (4096).bit_length()
    first, _ = unrank_type_word(counts, 0)
    last, _ = unrank_type_word(counts, number - 1)
    assert first == (0,) * 2048 + (1,) * 2048
    assert last == (1,) * 2048 + (0,) * 2048
    return dict(type_counts=types, complete_unranking_checks=words,
                largest_small_intermediate_bits=largest_bits,
                long_case=dict(counts=counts, count_bits=number.bit_length(),
                               intermediate_bits=largest.bit_length(),
                               avoided_factorial_bits=factorial_bits))


def ceil_log2(value):
    assert value >= 1
    return (value - 1).bit_length()


def block_length(label_length, block_constant):
    return max(1, min(label_length, isqrt(ceil_log2(label_length + 2)) // block_constant))


def parameter_checks():
    cutoffs = envelopes = degree_products = 0
    for coefficient in range(1, 17):
        block_constant = isqrt(1024 * coefficient - 1) + 1
        assert block_constant * block_constant >= 1024 * coefficient
        for exponent in (4, 16, 64, 256, 1024, 4096, 16384):
            n = 1 << exponent
            envelope = (1 << 16) * (n + 2) ** 3
            log_bound = ceil_log2(envelope + 2)
            maximum = max(1, isqrt(log_bound) // block_constant)
            assert (maximum + 1) ** 2 * block_constant ** 2 <= (
                8 * block_constant ** 2 + 2 * log_bound)
            # Exact sufficient form of G(m) <= const_C * n^(1/64).
            assert 64 * coefficient * (maximum + 1) ** 2 <= (
                512 * coefficient + 67 + exponent)
            envelopes += 1
            for label in (1, 2, isqrt(isqrt(n)), n, n * n, envelope):
                assert block_length(label, block_constant) <= maximum
                cutoffs += 1
    for degree_rate in range(1, 9):
        for maximum in range(1, 17):
            degrees = [1 << (degree_rate * t) for t in range(1, maximum + 1)]
            dimension = prod(degrees)
            assert dimension == 1 << (degree_rate * maximum * (maximum + 1) // 2)
            coefficient_bits = 1 << (degree_rate * maximum)
            tensor_bits = dimension ** 3 * coefficient_bits
            # Deliberately loose model including polynomial arithmetic
            # overhead; the degree product is never replaced by max(e_t).
            operation_budget = dimension ** 6 * coefficient_bits ** 8 * maximum ** 4
            family_constant = 8 * degree_rate + 16
            global_budget = 1 << (family_constant * (maximum + 1) ** 2)
            assert tensor_bits <= operation_budget <= global_budget
            if maximum >= 2:
                assert dimension > max(degrees)
            degree_products += 1
    return dict(global_envelopes=envelopes, block_and_tail_cutoffs=cutoffs,
                complete_degree_product_checks=degree_products,
                degree_growth_control=dict(size_count=16, degree_rate=2,
                    maximum_field_degree_log2=32,
                    tensor_dimension_log2=272,
                    multiplication_tensor_bits_log2=848),
                bit_budget='G(m(n)) <= 2^(8*C+67/64) * n^(1/64) on checked dyadic n')


def integer_root(value, degree):
    assert value >= 0 and degree >= 1
    lower, upper = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while lower < upper:
        middle = (lower + upper + 1) // 2
        if middle ** degree <= value:
            lower = middle
        else:
            upper = middle - 1
    assert lower ** degree <= value < (lower + 1) ** degree
    return lower


def global_block_checks():
    preprocessing = active_nodes = roots = degree_growth = depth_products = 0
    powers = (1, 2, 3, 4, 5, 8, 16)
    for power in powers:
        for root in range(65):
            for value in (max(0, root ** power - 1), root ** power,
                          (root + 1) ** power - 1):
                computed = integer_root(value, power)
                assert computed ** power <= value < (computed + 1) ** power
                roots += 1
    for coefficient in range(1, 17):
        block_constant = 128 * coefficient
        for power in powers:
            for exponent in (64, 256, 1024, 4096, 16384, 65536):
                n = 1 << exponent
                raw_length = max(1, integer_root(ceil_log2(n + 2) // block_constant,
                                                power))
                # Exact sufficient form for G(m)=2^(C*(m^h+1)):
                # G(m)^128 <= 2^(256*C+1) * n on these dyadic n.
                assert 128 * coefficient * (raw_length ** power + 1) <= (
                    256 * coefficient + 1 + exponent)
                preprocessing += 1
                minimum_label = (1 << (exponent // 4)) // 8
                maximum_label = (1 << 16) * (n + 2) ** 3
                for label_length in (minimum_label, n, maximum_label):
                    assert raw_length <= label_length
                    blocks, tail = divmod(label_length, raw_length)
                    assert blocks * raw_length + tail == label_length
                    assert 0 <= tail < raw_length
                    assert label_length.bit_length() <= 3 * exponent + 20
                    if exponent >= 1024 * coefficient:
                        assert coefficient * (raw_length ** power + 1) <= exponent // 4 - 3
                    active_nodes += 1
            # Symbolically take log2(n)=2^u. Avoid materializing n.
            # The reciprocal-log raw error has a bounded depth product.
            threshold = 2 * ceil_log2(block_constant) + 4 * power
            for depth_parameter in (threshold, threshold + 1, threshold + 16,
                                    threshold + 64, threshold + 256):
                log_length = 1 << depth_parameter
                raw_length = max(1, integer_root((log_length + 1) // block_constant,
                                                power))
                log_raw_lower = (raw_length + 2).bit_length() - 1
                depth = depth_parameter + 1
                assert depth <= 4 * power * log_raw_lower
                assert Fraction(log_raw_lower + 1, log_raw_lower) ** depth <= 3 ** (4 * power)
                depth_products += 1
    for power in powers:
        for degree_rate in range(1, 9):
            for raw_length in (1, 2, 3, 4, 8, 16, 32):
                exponent = raw_length ** power
                degree_log2 = degree_rate * exponent
                coordinate_bits_log2 = degree_rate * exponent
                tensor_bits_log2 = 3 * degree_log2 + coordinate_bits_log2
                # Include a polynomial piece count (m+1)^8 and all four
                # factor tables, plus dense exact algebra overhead.
                raw_tables_log2 = (2 + 4 * raw_length + degree_log2 +
                                  coordinate_bits_log2 + 8 * ceil_log2(raw_length + 1))
                arithmetic_log2 = (6 * degree_log2 + 8 * coordinate_bits_log2 +
                                   4 * ceil_log2(raw_length + 1))
                global_log2 = (16 * degree_rate + 32) * (exponent + 1)
                assert max(tensor_bits_log2, raw_tables_log2, arithmetic_log2) <= global_log2
                degree_growth += 1
    return dict(global_preprocessing_budgets=preprocessing,
                active_node_divisions=active_nodes,
                integer_root_boundary_checks=roots,
                fixed_cost_powers=list(powers),
                single_field_degree_growth_budgets=degree_growth,
                bounded_reciprocal_log_depth_products=depth_products,
                bit_budget='G(m(n)) <= 2^(2*C+1/128) * n^(1/128) on checked dyadic n',
                raw_cost='2^(O(t^h)) for a fixed integer h>=1',
                raw_width_error='O(t/log(t+2)); bounded product across O(loglog n) levels',
                no_per_tail_field=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='PASS',
                  scope='Growing-field coefficient-ring and indexing interfaces; no new factoring score.',
                  algebra=algebra_checks(),
                  global_tail=tail_checks(),
                  multinomial=multinomial_checks(),
                  parameters=parameter_checks(),
                  global_blocks=global_block_checks())
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
