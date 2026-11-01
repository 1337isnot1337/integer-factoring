#!/usr/bin/env python3
"""Complete normalized Bezout pairs on circuit-isometric separated nodes.

The index maps are i+2*g(i), where g is a supplied integer arithmetic
circuit of unrestricted formal degree. Separation and isometry certify
the global denominator exponent. No node list or factorial is expanded.
"""

from __future__ import annotations

from math import comb

from circuit_box_products import (
    _weighted_packet_mod, circuit_box_moments, circuit_box_product,
    evaluate_product_circuit, factorial_valuation, validate_circuit,
)


def valuation_two(value):
    if not value:
        raise ValueError("The valuation of zero is not a finite integer")
    value = abs(value)
    return (value & -value).bit_length() - 1


def _multiply(first, second, modulus, size):
    result = [0] * size
    for i, a in enumerate(first[:size]):
        if a:
            for j, b in enumerate(second[:size - i]):
                result[i + j] = (result[i + j] + a * b) % modulus
    return result


def _inverse(series, modulus, size):
    inverse_constant = pow(series[0], -1, modulus)
    result = [inverse_constant]
    for index in range(1, size):
        result.append(-inverse_constant * sum(
            series[j] * result[index - j]
            for j in range(1, min(index + 1, len(series)))) % modulus)
    return result


def product_jet(circuit, interval, precision, degree):
    """Return prod_(i in interval)(1-q(i)T) through T^degree modulo 2^p.

    The root circuit need only have the correct values modulo 2^p.
    Reversing the global parity packets avoids a factorial guard for
    the elementary functions. No polynomial of numerical degree M is
    expanded, even when degree is much smaller than the precision.
    """
    if precision < 1 or degree < 0:
        raise ValueError("Require positive precision and nonnegative jet degree")
    count = interval[1] - interval[0]
    retained = min(count, degree, precision - 1)
    working = max(1, precision - 1)
    data = circuit_box_moments(circuit, (interval,), retained, working)
    modulus = 1 << precision
    result = [1] + [0] * degree
    for parity in (0, 1):
        packet_count = data["counts"][parity]
        packet_degree = min(packet_count, degree, precision - 1)
        packet = _weighted_packet_mod(data["moments"][parity], packet_degree,
                                      working, precision)
        if parity == 0:
            reverse = list(packet) + [0] * (degree + 1 - len(packet))
        else:
            reverse = [sum(packet[j] * (-1) ** (index - j)
                           * comb(packet_count - j, index - j)
                           for j in range(min(index, packet_degree) + 1)
                           if index - j <= packet_count - j) % modulus
                       for index in range(degree + 1)]
        result = _multiply(result, reverse, modulus, degree + 1)
    return result


class _ScalarRing:
    def __init__(self, precision):
        self.modulus = 1 << precision

    def constant(self, value):
        return value % self.modulus

    def add(self, first, second):
        return (first + second) % self.modulus

    def multiply(self, first, second):
        return first * second % self.modulus


def _root_product_value(circuit, count, precision):
    encoded = circuit_box_product(circuit, ((0, count),), precision)
    return evaluate_product_circuit(encoded, _ScalarRing(precision), 0)


def odd_factorial_unit(number, precision):
    """Odd part of number! mod 2^p from logarithmically many odd buckets."""
    if number < 0 or precision < 1:
        raise ValueError("Invalid factorial-unit parameters")
    odd_circuit = {"nodes": [["var", 0], ["const", 2], ["mul", 0, 1],
                              ["const", 1], ["add", 2, 3]], "output": 4}
    modulus = 1 << precision
    result = 1
    remaining = number
    while remaining:
        odd_count = (remaining + 1) // 2
        value = _root_product_value(odd_circuit, odd_count, precision)
        result = result * ((-1) ** odd_count * value) % modulus
        remaining //= 2
    return result


class Builder:
    def __init__(self):
        self.nodes = [["var", 0]]

    def constant(self, value):
        self.nodes.append(["const", value])
        return len(self.nodes) - 1

    def operation(self, kind, first, second):
        self.nodes.append([kind, first, second])
        return len(self.nodes) - 1

    def scale(self, index, value):
        return self.operation("mul", index, self.constant(value))

    def affine(self, index, start, step):
        return self.operation("add", self.constant(start), self.scale(index, step))

    def import_circuit(self, circuit):
        mapped = []
        for instruction in circuit["nodes"]:
            operation = instruction[0]
            if operation == "var":
                assert instruction[1] == 0
                mapped.append(0)
            elif operation == "const":
                mapped.append(self.constant(instruction[1]))
            else:
                mapped.append(self.operation(operation, mapped[instruction[1]], mapped[instruction[2]]))
        return mapped[circuit["output"]]

    def power(self, index, exponent):
        result = self.constant(1)
        while exponent:
            if exponent & 1:
                result = self.operation("mul", result, index)
            exponent >>= 1
            if exponent:
                index = self.operation("mul", index, index)
        return result

    def finish(self, output):
        return {"nodes": self.nodes, "output": output}


def circuit_value(circuit, index, modulus):
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = index
        elif operation == "const":
            value = instruction[1]
        elif operation == "add":
            value = values[instruction[1]] + values[instruction[2]]
        else:
            value = values[instruction[1]] * values[instruction[2]]
        values.append(value % modulus)
    return values[circuit["output"]]


class IsometricContext:
    def __init__(self, number, first_a, step_a, first_b, step_b,
                 count, height, circuit_a, circuit_b):
        values = number, first_a, step_a, first_b, step_b, count, height
        if any(not isinstance(value, int) for value in values) or min(values) < 1:
            raise ValueError("Require positive integer parameters")
        if height >= number.bit_length():
            raise ValueError("Separating height exceeds the encoded number")
        residue_modulus = 1 << (height + 1)
        if first_a % residue_modulus != 1 or first_b % residue_modulus != 1:
            raise ValueError("Both initial nodes must be one modulo 2^(h+1)")
        if valuation_two(step_a) != height + 1 or valuation_two(step_b) != height + 1:
            raise ValueError("Both step valuations must equal h+1")
        if number % residue_modulus != 1 + (1 << height):
            raise ValueError("The number does not separate the two root classes")
        self.degree_a = validate_circuit(circuit_a, 1)
        self.degree_b = validate_circuit(circuit_b, 1)
        self.number = number
        self.first_a = first_a
        self.step_a = step_a
        self.first_b = first_b
        self.step_b = step_b
        self.count = count
        self.height = height
        self.circuit_a = circuit_a
        self.circuit_b = circuit_b
        self._z_cache = {}
        self.a_circuit = self._domain_circuit(first_a, step_a, circuit_a)
        self.b_circuit = self._domain_circuit(first_b, step_b, circuit_b)

    @staticmethod
    def _domain_circuit(start, step, supplied):
        builder = Builder()
        perturbation = builder.import_circuit(supplied)
        index_map = builder.operation("add", 0, builder.scale(perturbation, 2))
        return builder.finish(builder.affine(index_map, start, step))

    def z_circuits(self, precision):
        if precision in self._z_cache:
            return self._z_cache[precision]
        even = self._domain_circuit((self.first_a - 1) >> self.height,
                                    self.step_a >> self.height, self.circuit_a)
        builder = Builder()
        perturbation = builder.import_circuit(self.circuit_b)
        index_map = builder.operation("add", 0, builder.scale(perturbation, 2))
        denominator = builder.affine(index_map, self.first_b, self.step_b)
        numerator = builder.affine(index_map, (self.number - self.first_b) >> self.height,
                                   -(self.step_b >> self.height))
        # The denominator is structurally odd. Euler's theorem modulo
        # 2^precision replaces its inverse by an ordinary binary power.
        inverse = builder.power(denominator, (1 << (precision - 1)) - 1)
        odd = builder.finish(builder.operation("mul", numerator, inverse))
        self._z_cache[precision] = even, odd
        return even, odd

    def node_jet(self, node_count, precision, degree, root_shifts=False):
        even, odd = self.z_circuits(precision)
        first = product_jet(even, (0, (node_count + 1) // 2), precision, degree)
        second = product_jet(odd, (0, node_count // 2), precision, degree)
        result = _multiply(first, second, 1 << precision, degree + 1)
        if root_shifts:
            result = [(value << (self.height * j)) % (1 << precision)
                      if self.height * j < precision else 0 for j, value in enumerate(result)]
        return result

    def root_shift(self, index, modulus):
        if index % 2 == 0:
            return (circuit_value(self.a_circuit, index // 2, modulus) - 1) % modulus
        denominator = circuit_value(self.b_circuit, index // 2, modulus)
        return (self.number * pow(denominator, -1, modulus) - 1) % modulus

    def domain_reversals(self, precision, degree):
        even, odd = self.z_circuits(precision)
        first = product_jet(even, (0, self.count), precision, degree)
        second = product_jet(odd, (0, self.count), precision, degree)
        modulus = 1 << precision
        leading_q = _root_product_value(self.b_circuit, self.count, precision)
        first = [(value << (self.height * j)) % modulus
                 if self.height * j < precision else 0 for j, value in enumerate(first)]
        second = [(leading_q * value << (self.height * j)) % modulus
                  if self.height * j < precision else 0 for j, value in enumerate(second)]
        return first, second


def normalized_parity_difference(context, order, precision):
    """The normalized huge-order parity divided difference from short jets."""
    # Each term already has a factor 2^j. Its net division loss is
    # v2((k+1)...(k+j))-j = popcount(k)-popcount(k+j).
    guard = max(order.bit_count() - (order + j).bit_count()
                for j in range(precision))
    working = precision + guard
    modulus = 1 << working
    inverse_nodes = _inverse(context.node_jet(order + 1, working, precision - 1),
                             modulus, precision)
    identity = {"nodes": [["var", 0]], "output": 0}
    falling = product_jet(identity, (0, order), working, precision - 1)
    denominator, answer = 1, 0
    output_modulus = 1 << precision
    for j in range(precision):
        if j:
            last = order + j - 1
            falling = [(falling[i] - (last * falling[i - 1] if i else 0)) % modulus
                       for i in range(precision)]
            denominator *= order + j
        coefficient = sum(falling[i] * inverse_nodes[j - i] for i in range(j + 1)) % modulus
        loss = valuation_two(denominator)
        if coefficient % (1 << loss):
            raise ArithmeticError("A divided-difference numerator lost its proved divisibility")
        # Only precision-j quotient bits are required; the following
        # shift restores the output precision. No unknown high bits of
        # the quotient are being asserted.
        quotient = (coefficient >> loss) * pow(denominator >> loss, -1, output_modulus) % output_modulus
        answer = (answer + (-1) ** j * (quotient << j)) % output_modulus
    if not answer & 1:
        raise ArithmeticError("The normalized parity divided difference must be odd")
    return {"unit": answer, "precision": precision, "order": order,
            "working_precision": working, "factorial_ratio_guard_bits": guard}


def circuit_isometric_bezout(number, first_a, step_a, first_b, step_b,
                            count, height, precision, circuit_a, circuit_b):
    """Complete U,V mod 2^p with UP+VQ=2^r and exactly optimal scale r."""
    if precision < 1:
        raise ValueError("Require positive precision")
    context = IsometricContext(number, first_a, step_a, first_b, step_b,
                                count, height, circuit_a, circuit_b)
    return _compile_context(context, precision)


def _compile_context(context, precision, difference_evaluator=None):
    if difference_evaluator is None:
        difference_evaluator = normalized_parity_difference
    count, height = context.count, context.height
    top = 2 * count - 1
    exponent = height * top - (count - 1).bit_count()
    retained = []
    for offset in range(min(top, 2 * precision)):
        order = top - offset
        gap = (height - 1) * offset + factorial_valuation(top) - factorial_valuation(order)
        if gap >= precision:
            break
        retained.append((order, gap))
    assert retained and retained[0] == (top, 0)
    modulus = 1 << precision
    factorial_unit = odd_factorial_unit(top, precision)
    root_degree = min(top, (precision - 1) // height)
    prefix = context.node_jet(top, precision, root_degree, root_shifts=True)
    f_degree = min(top, top - retained[-1][0] + root_degree)
    f_reverse = [0] * (f_degree + 1)
    certificates = []
    for order, gap in retained:
        difference = difference_evaluator(context, order, precision - gap)
        coefficient = ((-1) ** (order + 1) * difference["unit"]
                       * pow(factorial_unit, -1, modulus) << gap) % modulus
        offset = top - order
        for j, value in enumerate(prefix):
            if offset + j <= f_degree:
                f_reverse[offset + j] = (f_reverse[offset + j] + coefficient * value) % modulus
        certificates.append({**difference, "valuation_gap": gap,
                             "scaled_newton_coefficient": coefficient})
        if order > retained[-1][0]:
            shift = context.root_shift(order - 1, modulus)
            smaller = [prefix[0]]
            for j in range(1, len(prefix)):
                smaller.append((prefix[j] + shift * smaller[-1]) % modulus)
            prefix = smaller
            factorial_unit = factorial_unit * pow(order >> valuation_two(order), -1, modulus) % modulus
    degree = count - 1
    u_depth = min(degree, f_degree + (precision - 1) // (height + 1))
    v_depth = min(degree, f_degree + (precision - 1) // height)
    p_reverse, q_reverse = context.domain_reversals(precision, max(u_depth, v_depth))
    u_reverse = _multiply(f_reverse, _inverse(p_reverse, modulus, u_depth + 1), modulus, u_depth + 1)
    v_reverse = [(-value) % modulus for value in _multiply(
        f_reverse, _inverse(q_reverse, modulus, v_depth + 1), modulus, v_depth + 1)]
    assert u_reverse[0] & 1 and v_reverse[0] & 1
    assert len(u_reverse) < 4 * precision and len(v_reverse) < 4 * precision
    return {"parameters": {"number": context.number, "first_a": context.first_a,
                            "step_a": context.step_a, "first_b": context.first_b,
                            "step_b": context.step_b, "count": count,
                            "height": height, "precision": precision},
            "input_circuit_nodes": [len(context.circuit_a["nodes"]), len(context.circuit_b["nodes"])],
            "perturbation_degree_bounds": [context.degree_a, context.degree_b],
            "bezout_exponent": exponent, "polynomial_degree": degree, "basis_center": 1,
            "u_reversed_coefficients": u_reverse, "v_reversed_coefficients": v_reverse,
            "retained_newton_coefficients": certificates,
            "maximum_working_precision": max(item["working_precision"] for item in certificates),
            "direct_boundary_queries": max(0, len(retained) - 1),
            "scope": "Complete globally normalized Bezout pairs on certified circuit-isometric separated domains; no arbitrary factor-root search."}
