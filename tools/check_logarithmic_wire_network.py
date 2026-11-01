#!/usr/bin/env python3
"""Exact checks of logarithmic-wire multiplication and bilinear block queries."""

from itertools import combinations, product
import json
from math import isqrt
from random import Random

from check_polynomial_top_queries import Network, interpolate


def is_prime(number):
    return number >= 2 and all(number % divisor for divisor in range(2, isqrt(number) + 1))


def parameters(width, coefficient_prime=2):
    digit_bits = (width + 1).bit_length()
    size = 2
    while size * size // 2 * digit_bits < width:
        size *= 2
    base = 1 << digit_bits
    bound = size * size // 2 * (base - 1) ** 2 + 2 * size + coefficient_prime
    return digit_bits, size, base, bound


def interpolation_matrix(size, modulus):
    result = [[0] * size for _ in range(size)]
    for selected in range(size):
        polynomial, denominator = [1], 1
        for other in range(size):
            if other == selected:
                continue
            following = [0] * (len(polynomial) + 1)
            for index, value in enumerate(polynomial):
                following[index] -= other * value
                following[index + 1] += value
            polynomial = [value % modulus for value in following]
            denominator = denominator * (selected - other) % modulus
        inverse = pow(denominator, -1, modulus)
        for power, value in enumerate(polynomial):
            result[power][selected] = value * inverse % modulus
    return result


class BivariateNetwork(Network):
    def __init__(self, width, coefficient_prime=2, toy_base=None):
        self.d, self.a, source_base, bound = parameters(width, coefficient_prime)
        self.h = self.a // 2
        self.modulus = bound + 1
        while not is_prime(self.modulus):
            self.modulus += 1
        self.base = source_base if toy_base is None else toy_base
        self.digit_bits = self.base.bit_length() - 1
        self.input_bits = self.a * self.h * self.digit_bits
        self.f = self.modulus.bit_length()
        self.carry_limit = self.modulus
        self.inverse_interpolation = interpolation_matrix(self.a, self.modulus)
        self.gates = [(kind, index) for kind in ("p", "q", "mid")
                      for index in range(self.a)]
        self.gates += [("ib", index) for index in range(2 * self.a)]
        self.gates += [(kind, index) for kind in ("add", "carry")
                      for index in range(self.a ** 2)]
        self.order = {gate: index for index, gate in enumerate(self.gates)}
        self.edges = [((kind, row), ("mid", column)) for kind in ("p", "q")
                      for row in range(self.a) for column in range(self.a)]
        self.edges += [(("mid", row), ("ib", column))
                       for row in range(self.a) for column in range(2 * self.a)]
        self.edges += [(("ib", family * self.a + index % self.a), ("add", index))
                       for family in range(2) for index in range(self.a ** 2)]
        self.edges += [(("add", index), ("carry", index)) for index in range(self.a ** 2)]
        self.edges += [(("carry", index), ("carry", index + 1))
                       for index in range(self.a ** 2 - 1)]
        self.edges.sort(key=lambda edge: (self.order[edge[0]], self.order[edge[1]]))
        assert len(self.gates) == 2 * self.a ** 2 + 5 * self.a
        assert len(self.edges) == 8 * self.a ** 2 - 1
        assert self.modulus > bound and self.modulus <= 4 * (bound + 2) ** 2

    def fa_outputs(self, row, digits):
        return tuple(sum(value * pow(point, position, self.modulus)
                         for position, value in enumerate(digits)) % self.modulus
                     for point in range(self.a))

    def mid_outputs(self, column, first, second):
        return tuple(sum(first[row] * second[power - row] for row in range(self.a)
                         if 0 <= power - row < self.a) % self.modulus
                     for power in range(2 * self.a))

    def ib_outputs(self, incoming):
        return tuple(sum(value * coefficient for value, coefficient in zip(incoming, row))
                     % self.modulus for row in self.inverse_interpolation)

    def outputs(self, gate, setting):
        if gate[0] == "add":
            return (sum(setting) % self.modulus,)
        return super().outputs(gate, setting)

    def message(self, source, destination, settings, outgoing):
        kind, index = source
        other_kind, other_index = destination
        if kind == "ib" and other_kind == "add":
            quotient, residue = divmod(other_index, self.a)
            assert index % self.a == residue
            if index < self.a:
                return outgoing[source][quotient]
            return outgoing[source][quotient - 1] if quotient else 0
        if kind == "add":
            if other_kind == "ib":
                return settings[source][other_index >= self.a]
            assert other_kind == "carry" and index == other_index
            return outgoing[source][0]
        if kind == "carry" and other_kind == "add":
            assert index == other_index
            return settings[source][0]
        return super().message(source, destination, settings, outgoing)

    def grouping(self, count):
        assert 2 <= count <= self.a and self.a % count == 0
        block = self.a // count
        result = {}
        for gate in self.gates:
            kind, index = gate
            coordinate = index // self.a if kind in ("add", "carry") else index % self.a
            result[gate] = coordinate // block
        return result

    def global_setting(self, first, second):
        settings, outgoing = {}, {}
        for kind, digits in (("p", first), ("q", second)):
            for row in range(self.a):
                gate = (kind, row)
                settings[gate] = tuple(digits[row + self.a * i] for i in range(self.h))
                outgoing[gate] = self.outputs(gate, settings[gate])
        for column in range(self.a):
            gate = ("mid", column)
            settings[gate] = tuple(tuple(outgoing[(kind, row)][column]
                                        for row in range(self.a)) for kind in ("p", "q"))
            outgoing[gate] = self.outputs(gate, settings[gate])
        # Independent two-variable coefficient convolution.
        bivariate = [[0] * self.a for _ in range(2 * self.a)]
        ordinary = [0] * self.a ** 2
        for first_index, left in enumerate(first):
            for second_index, right in enumerate(second):
                x_power = first_index % self.a + second_index % self.a
                y_power = first_index // self.a + second_index // self.a
                bivariate[x_power][y_power] += left * right
                ordinary[first_index + second_index] += left * right
        for row in range(2 * self.a):
            gate = ("ib", row)
            settings[gate] = tuple(outgoing[("mid", column)][row] for column in range(self.a))
            outgoing[gate] = self.outputs(gate, settings[gate])
            assert outgoing[gate] == tuple(bivariate[row])
        carry = 0
        for index in range(self.a ** 2):
            quotient, residue = divmod(index, self.a)
            gate = ("add", index)
            settings[gate] = (outgoing[("ib", residue)][quotient],
                              outgoing[("ib", self.a + residue)][quotient - 1]
                              if quotient else 0)
            outgoing[gate] = self.outputs(gate, settings[gate])
            assert outgoing[gate][0] == ordinary[index]
        for index in range(self.a ** 2):
            gate = ("carry", index)
            settings[gate] = (outgoing[("add", index)][0], carry)
            outgoing[gate] = self.outputs(gate, settings[gate])
            carry = outgoing[gate][0]
            assert carry < self.modulus
        assert carry == 0
        first_value = sum(value * self.base ** i for i, value in enumerate(first))
        second_value = sum(value * self.base ** i for i, value in enumerate(second))
        number = first_value * second_value
        for index in range(self.a ** 2):
            coefficient, incoming = settings[("carry", index)]
            assert (coefficient + incoming) % self.base == number // self.base ** index % self.base
        return settings, outgoing, number, first_value

    def propagate(self, cluster, boundary, first_settings, number, cutoff, prefix):
        expected_keys = {(first, second) if first in cluster else (second, first)
                         for first, second in self.edges
                         if (first in cluster) != (second in cluster)}
        if boundary.keys() != expected_keys:
            return None
        if any(not 0 <= value < 1 << self.f for value in boundary.values()):
            return None
        settings, outgoing = {}, {}

        def read(destination, source):
            return (self.message(source, destination, settings, outgoing) if source in cluster
                    else boundary[(destination, source)])

        for gate in self.gates:
            if gate not in cluster:
                continue
            kind, index = gate
            if kind in ("p", "q"):
                setting = first_settings[gate]
                if len(setting) != self.h or any(not 0 <= value < self.base for value in setting):
                    return None
                if kind == "p" and not self.prefix_ok(index, setting, cutoff, prefix):
                    return None
            elif kind == "mid":
                setting = tuple(tuple(read(gate, (input_kind, row)) for row in range(self.a))
                                for input_kind in ("p", "q"))
                if any(not 0 <= value < self.modulus for vector in setting for value in vector):
                    return None
            elif kind == "ib":
                setting = tuple(read(gate, ("mid", row)) for row in range(self.a))
                if any(not 0 <= value < self.modulus for value in setting):
                    return None
            elif kind == "add":
                setting = tuple(read(gate, ("ib", family * self.a + index % self.a))
                                for family in range(2))
                if any(not 0 <= value < self.modulus for value in setting):
                    return None
            else:
                coefficient = read(gate, ("add", index))
                carry = 0 if index == 0 else read(gate, ("carry", index - 1))
                if not 0 <= coefficient < self.modulus or not 0 <= carry < self.modulus:
                    return None
                total = coefficient + carry
                if total % self.base != number // self.base ** index % self.base:
                    return None
                if index == self.a ** 2 - 1 and total // self.base:
                    return None
                setting = (coefficient, carry)
            settings[gate] = setting
            outgoing[gate] = self.outputs(gate, setting)
            for (inside, outside), value in boundary.items():
                if inside == gate and self.message(inside, outside, settings, outgoing) != value:
                    return None
        for first, second in self.edges:
            if first in cluster and second in cluster:
                assert self.message(first, second, settings, outgoing) == (
                    self.message(second, first, settings, outgoing))
        return settings

    def query(self, cluster, boundary, number, cutoff, prefix):
        assert sum(kind == "mid" for kind, _ in cluster) <= self.h
        first_settings = {}
        for gate in self.gates:
            if gate not in cluster or gate[0] not in ("p", "q"):
                continue
            selected = [column for column in range(self.a) if ("mid", column) not in cluster][:self.h]
            if any((gate, ("mid", column)) not in boundary for column in selected):
                return None
            values = [boundary[(gate, ("mid", column))] for column in selected]
            if any(not 0 <= value < self.modulus for value in values):
                return None
            first_settings[gate] = interpolate(selected, values, self.modulus)
        return self.propagate(cluster, boundary, first_settings, number, cutoff, prefix)


class CharacterAlgebra:
    """F_p[X]/(1+X+...+X^(M-1)), with no irreducibility assumption."""

    def __init__(self, modulus, prime):
        self.modulus, self.prime = modulus, prime
        self.zero = (0,) * (modulus - 1)

    def add(self, first, second):
        return tuple((a + b) % self.prime for a, b in zip(first, second))

    def scale(self, value, scalar):
        return tuple(scalar * entry % self.prime for entry in value)

    def phase(self, exponent):
        exponent %= self.modulus
        if exponent == self.modulus - 1:
            return (self.prime - 1,) * (self.modulus - 1)
        return tuple(int(index == exponent) for index in range(self.modulus - 1))

    def multiply(self, first, second):
        coefficients = [0] * (2 * self.modulus - 3)
        for i, a in enumerate(first):
            for j, b in enumerate(second):
                coefficients[i + j] += a * b
        for power in range(len(coefficients) - 1, self.modulus - 2, -1):
            scalar = coefficients[power] % self.prime
            for index in range(self.modulus):
                coefficients[power - self.modulus + 1 + index] -= scalar
        return tuple(value % self.prime for value in coefficients[:self.modulus - 1])


def linear_solve(matrix, target, modulus):
    rows, columns = len(matrix), len(matrix[0])
    augmented = [[value % modulus for value in row] + [rhs % modulus]
                 for row, rhs in zip(matrix, target)]
    pivot_columns, active = [], 0
    for column in range(columns):
        pivot = next((row for row in range(active, rows) if augmented[row][column]), None)
        if pivot is None:
            continue
        augmented[active], augmented[pivot] = augmented[pivot], augmented[active]
        inverse = pow(augmented[active][column], -1, modulus)
        augmented[active] = [value * inverse % modulus for value in augmented[active]]
        for row in range(rows):
            if row != active:
                scalar = augmented[row][column]
                augmented[row] = [(value - scalar * chosen) % modulus
                                  for value, chosen in zip(augmented[row], augmented[active])]
        pivot_columns.append(column)
        active += 1
        if active == rows:
            break
    if any(not any(row[:columns]) and row[columns] for row in augmented):
        return len(pivot_columns), None
    solution = [0] * columns
    for row, column in enumerate(pivot_columns):
        solution[column] = augmented[row][columns]
    return len(pivot_columns), tuple(solution)


def bilinear_closed(matrix, alpha, beta, algebra):
    modulus = algebra.modulus
    rank, solution = linear_solve(matrix, [-value for value in alpha], modulus)
    if solution is None:
        return algebra.zero
    transpose = [list(row) for row in zip(*matrix)]
    _, dual = linear_solve(transpose, beta, modulus)
    if dual is None:
        return algebra.zero
    phase = algebra.phase(sum(value * chosen for value, chosen in zip(beta, solution)))
    exponent = len(alpha) + len(beta) - rank
    return algebra.scale(phase, pow(modulus, exponent, algebra.prime))


def bilinear_direct(matrix, alpha, beta, algebra):
    value = algebra.zero
    for first in product(range(algebra.modulus), repeat=len(alpha)):
        for second in product(range(algebra.modulus), repeat=len(beta)):
            exponent = sum(a * b for a, b in zip(alpha, first))
            exponent += sum(a * b for a, b in zip(beta, second))
            exponent += sum(first[i] * coefficient * second[j]
                            for i, row in enumerate(matrix) for j, coefficient in enumerate(row))
            value = algebra.add(value, algebra.phase(exponent))
    return value


def carry_prefix_query(digits, frequencies, incoming_frequency, outgoing_frequency,
                       base, initial, final, algebra):
    modulus = algebra.modulus
    states = [algebra.phase(incoming_frequency * carry) if not initial or carry == 0
              else algebra.zero for carry in range(modulus)]
    for digit, frequency in zip(digits, frequencies):
        prefix = [algebra.zero]
        for carry in range(modulus):
            weighted = algebra.multiply(states[carry], algebra.phase(-frequency * carry))
            prefix.append(algebra.add(prefix[-1], weighted))
        following = []
        for outgoing in range(modulus):
            lower = max(0, digit + base * outgoing - modulus + 1)
            upper = min(modulus - 1, digit + base * outgoing)
            if lower > upper:
                following.append(algebra.zero)
            else:
                interval = algebra.add(prefix[upper + 1], algebra.scale(prefix[lower], -1))
                following.append(algebra.multiply(
                    algebra.phase(frequency * (digit + base * outgoing)), interval))
        states = following
    if final:
        return states[0]
    result = algebra.zero
    for carry in range(modulus):
        result = algebra.add(result, algebra.multiply(
            states[carry], algebra.phase(outgoing_frequency * carry)))
    return result


def carry_direct_query(digits, frequencies, incoming_frequency, outgoing_frequency,
                       base, initial, final, algebra):
    result = algebra.zero
    for incoming in range(1 if initial else algebra.modulus):
        for coefficients in product(range(algebra.modulus), repeat=len(digits)):
            carry = incoming
            for digit, coefficient in zip(digits, coefficients):
                total = coefficient + carry
                if total % base != digit:
                    break
                carry = total // base
                assert carry < algebra.modulus
            else:
                if final and carry:
                    continue
                phase = incoming_frequency * incoming + outgoing_frequency * carry
                phase += sum(frequency * coefficient
                             for frequency, coefficient in zip(frequencies, coefficients))
                result = algebra.add(result, algebra.phase(phase))
    return result


def verify():
    rng = Random(20261030)
    parameter_cases = 0
    for width in range(2, 100001):
        digit_bits, size, base, bound = parameters(width)
        assert digit_bits <= width
        assert size * size * digit_bits < 8 * width
        assert bound + 2 < 64 * width ** 3
        upper_modulus = 4 * (bound + 2) ** 2
        assert upper_modulus.bit_length() <= 16 * digit_bits
        assert 13 * size * size * upper_modulus.bit_length() < 1664 * width
        parameter_cases += 1

    pipeline_cases = edge_cases = group_cases = 0
    arbitrary_half_evaluations = 0
    prime_records = []
    for width in (2, 3, 9, 32, 64, 257):
        network = BivariateNetwork(width)
        prime_records.append({"n": width, "A": network.a, "d": network.d,
                              "M": network.modulus, "wire_bits": network.f})
        for _ in range(12):
            first = tuple(rng.randrange(network.base) for _ in range(network.a * network.h))
            second = tuple(rng.randrange(network.base) for _ in range(network.a * network.h))
            network.global_setting(first, second)
            pipeline_cases += 1
        digits = tuple(rng.randrange(network.base) for _ in range(network.h))
        outputs = network.fa_outputs(0, digits)
        selected_sets = (list(combinations(range(network.a), network.h)) if network.a <= 8
                         else [tuple(sorted(rng.sample(range(network.a), network.h))) for _ in range(24)])
        for selected in selected_sets:
            assert interpolate(selected, [outputs[column] for column in selected],
                               network.modulus) == digits
            arbitrary_half_evaluations += 1
        for power in range(1, network.a.bit_length()):
            count = 1 << power
            groups = network.grouping(count)
            actual = {(a, b): 0 for a in range(count) for b in range(a + 1, count)}
            for first, second in network.edges:
                a, b = sorted((groups[first], groups[second]))
                if a != b:
                    actual[a, b] += 1
            block = network.a // count
            for (a, b), edges in actual.items():
                assert edges == 12 * block ** 2 + (b == a + 1)
                edge_cases += 1
            label = (12 * block ** 2 + 1) * network.f
            assert count ** 2 * label < 1664 * width
            group_cases += 1

    cluster_cases = malformed = prefixes = weighted = 0
    network = BivariateNetwork(2)
    grouping = network.grouping(2)
    clusters = [{gate for gate in network.gates if grouping[gate] == index} for index in range(2)]

    def compare(model, cluster, boundary, number, cutoff, prefix):
        nonlocal cluster_cases, weighted
        direct = model.enumerate_query(cluster, boundary, number, cutoff, prefix)
        fast = model.query(cluster, boundary, number, cutoff, prefix)
        assert direct == ([] if fast is None else [fast])
        if fast is not None:
            assert len(direct) == 1
            expected_weight = sum(value * model.base ** (row + model.a * i)
                                  for (kind, row), digits in direct[0].items() if kind == "p"
                                  for i, value in enumerate(digits))
            actual_weight = sum(value * model.base ** (row + model.a * i)
                                for (kind, row), digits in fast.items() if kind == "p"
                                for i, value in enumerate(digits))
            assert expected_weight == actual_weight
            weighted += 1
        cluster_cases += 1

    product_witnesses = {}
    for first_value, second_value in product(range(16), repeat=2):
        first = (first_value % 4, first_value // 4)
        second = (second_value % 4, second_value // 4)
        settings, outgoing, number, value = network.global_setting(first, second)
        product_witnesses.setdefault(number, []).append(value)
        for cluster in clusters:
            boundary = network.boundary(cluster, settings, outgoing)
            for cutoff in (0, 2, 4):
                compare(network, cluster, boundary, number, cutoff, value >> cutoff)
                prefixes += 1
            if (first_value + 3 * second_value) % 16 == 0:
                for key in boundary:
                    changed = dict(boundary)
                    changed[key] = (changed[key] + 1) % network.modulus
                    compare(network, cluster, changed, number, 4, 0)
                    malformed += 1
    witness_cases = 0
    for number in range(2, 16):
        assert sorted(product_witnesses[number]) == [
            value for value in range(1, number + 1) if number % value == 0]
        witness_cases += 1

    # A=4 with radix four permits independent enumeration of every local digit.
    toy = BivariateNetwork(9, toy_base=4)
    groups = toy.grouping(4)
    for _ in range(6):
        first = tuple(rng.randrange(4) for _ in range(toy.a * toy.h))
        second = tuple(rng.randrange(4) for _ in range(toy.a * toy.h))
        settings, outgoing, number, value = toy.global_setting(first, second)
        for group in range(4):
            cluster = {gate for gate in toy.gates if groups[gate] == group}
            boundary = toy.boundary(cluster, settings, outgoing)
            compare(toy, cluster, boundary, number, toy.input_bits, 0)
            changed = dict(boundary)
            key = list(changed)[rng.randrange(len(changed))]
            changed[key] = (changed[key] + 1) % toy.modulus
            compare(toy, cluster, changed, number, toy.input_bits, 0)
            malformed += 1

    character_cases = gauss_cases = singular_cases = rectangular_cases = 0
    for modulus, coefficient_prime in ((3, 2), (5, 2), (5, 3), (7, 2)):
        algebra = CharacterAlgebra(modulus, coefficient_prime)
        for frequency in range(modulus):
            total = algebra.zero
            for variable in range(modulus):
                total = algebra.add(total, algebra.phase(frequency * variable))
            expected = algebra.scale(algebra.phase(0), modulus) if frequency == 0 else algebra.zero
            assert total == expected
            character_cases += 1
        for a, b in product(range(modulus), repeat=2):
            assert algebra.multiply(algebra.phase(a), algebra.phase(b)) == algebra.phase(a + b)
            character_cases += 1
        for rows, columns in ((1, 1), (2, 2), (2, 3)):
            for trial in range(12):
                if trial == 0:
                    matrix = [[0] * columns for _ in range(rows)]
                elif trial == 1:
                    matrix = [[1] * columns for _ in range(rows)]
                elif trial % 2:
                    gamma = [rng.randrange(modulus) for _ in range(rows + columns - 1)]
                    matrix = [[gamma[i + j] for j in range(columns)] for i in range(rows)]
                else:
                    matrix = [[rng.randrange(modulus) for _ in range(columns)] for _ in range(rows)]
                for choice in range(4):
                    alpha = ([0] * rows if choice == 0 else [rng.randrange(modulus) for _ in range(rows)])
                    beta = ([0] * columns if choice == 0 else
                            [rng.randrange(modulus) for _ in range(columns)])
                    assert bilinear_closed(matrix, alpha, beta, algebra) == (
                        bilinear_direct(matrix, alpha, beta, algebra))
                    rank, _ = linear_solve(matrix, [0] * rows, modulus)
                    singular_cases += rank < min(rows, columns)
                    rectangular_cases += rows != columns
                    gauss_cases += 1

    # A genuine zero-divisor coefficient algebra still supports every identity.
    algebra = CharacterAlgebra(7, 2)
    first = tuple(int(i in (0, 1, 3)) for i in range(6))
    second = tuple(int(i in (0, 2, 3)) for i in range(6))
    assert first != algebra.zero and second != algebra.zero
    assert algebra.multiply(first, second) == algebra.zero
    # The all-ones quotient for a composite arithmetic modulus is not enough.
    composite = CharacterAlgebra(9, 2)
    nonprincipal_sum = composite.zero
    for value in range(9):
        nonprincipal_sum = composite.add(nonprincipal_sum, composite.phase(3 * value))
    assert nonprincipal_sum != composite.zero

    carry_cases = 0
    for modulus, coefficient_prime, base in ((3, 2, 2), (5, 2, 4), (5, 3, 2), (7, 2, 4)):
        algebra = CharacterAlgebra(modulus, coefficient_prime)
        for length in (1, 2, 3):
            for initial, final in product((False, True), repeat=2):
                for trial in range(5):
                    digits = [rng.randrange(base) for _ in range(length)]
                    frequencies = [0] * length if trial == 0 else [
                        rng.randrange(modulus) for _ in range(length)]
                    incoming_frequency = rng.randrange(modulus)
                    outgoing_frequency = rng.randrange(modulus)
                    arguments = (digits, frequencies, incoming_frequency, outgoing_frequency,
                                 base, initial, final, algebra)
                    assert carry_prefix_query(*arguments) == carry_direct_query(*arguments)
                    carry_cases += 1

    # Contract a complete two-middle arithmetic core independently, including
    # interpolation and the shifted bivariate-to-univariate coefficient addition.
    # Incoming boundary values remain arbitrary; no FA consistency is assumed.
    core_cases = 0
    modulus, size = 3, 2
    algebra = CharacterAlgebra(modulus, 2)
    inverse = interpolation_matrix(size, modulus)
    for trial in range(8):
        alpha = [[rng.randrange(modulus) for _ in range(size)] for _ in range(size)]
        beta = [[rng.randrange(modulus) for _ in range(size)] for _ in range(size)]
        output_frequency = [rng.randrange(modulus) for _ in range(size ** 2)]
        if trial == 0:
            alpha = beta = [[0] * size for _ in range(size)]
            output_frequency = [0] * (size ** 2)
        coefficient_frequency = [[0] * size for _ in range(2 * size)]
        for index, frequency in enumerate(output_frequency):
            quotient, residue = divmod(index, size)
            coefficient_frequency[residue][quotient] += frequency
            if quotient:
                coefficient_frequency[size + residue][quotient - 1] += frequency
        middle_frequency = [[sum(coefficient_frequency[u][degree] * inverse[degree][point]
                                 for degree in range(size)) % modulus
                             for u in range(2 * size)] for point in range(size)]
        closed = algebra.phase(0)
        for point in range(size):
            matrix = [[-middle_frequency[point][row + column] % modulus
                       for column in range(size)] for row in range(size)]
            closed = algebra.multiply(closed, bilinear_closed(
                matrix, alpha[point], beta[point], algebra))
        direct = algebra.zero
        for inputs in product(range(modulus), repeat=2 * size ** 2):
            first = [inputs[2 * size * point:2 * size * point + size] for point in range(size)]
            second = [inputs[2 * size * point + size:2 * size * (point + 1)]
                      for point in range(size)]
            middle = [[sum(first[point][row] * second[point][power - row]
                           for row in range(size) if 0 <= power - row < size) % modulus
                       for power in range(2 * size)] for point in range(size)]
            coefficients = [[sum(middle[point][power] * inverse[degree][point]
                                 for point in range(size)) % modulus
                             for degree in range(size)] for power in range(2 * size)]
            univariate = []
            for index in range(size ** 2):
                quotient, residue = divmod(index, size)
                value = coefficients[residue][quotient]
                if quotient:
                    value += coefficients[size + residue][quotient - 1]
                univariate.append(value % modulus)
            phase = sum(alpha[point][row] * first[point][row]
                        + beta[point][row] * second[point][row]
                        for point in range(size) for row in range(size))
            phase -= sum(frequency * value for frequency, value in zip(output_frequency, univariate))
            direct = algebra.add(direct, algebra.phase(phase))
        assert direct == closed
        core_cases += 1

    return {
        "scope": "exact bivariate network and restricted bilinear contraction checks; no improved logarithmic saving exponent",
        "parameter_information_bounds": parameter_cases,
        "ordinary_vs_bivariate_pipelines": pipeline_cases,
        "arbitrary_half_evaluation_recoveries": arbitrary_half_evaluations,
        "exact_intergroup_edge_counts": edge_cases,
        "grouping_volume_cases": group_cases,
        "independent_cluster_enumerations": cluster_cases,
        "malformed_boundary_cases": malformed,
        "prefix_cases": prefixes,
        "weighted_unique_extensions": weighted,
        "complete_small_divisor_bijections": witness_cases,
        "character_algebra_identities": character_cases,
        "bilinear_gauss_cases": gauss_cases,
        "singular_gauss_cases": singular_cases,
        "rectangular_gauss_cases": rectangular_cases,
        "complete_arithmetic_core_contractions": core_cases,
        "streamed_carry_block_contractions": carry_cases,
        "constructed_arithmetic_primes": prime_records,
        "negative_controls": {
            "coefficient_zero_divisors_are_permitted": True,
            "composite_all_ones_quotient_fails_orthogonality": True,
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
