#!/usr/bin/env python3
"""Compressed polynomial-circuit moments of an implicit pole Cauchy inverse.

The circuit degree is unrestricted. The algorithm uses existing unit-
progression arithmetic and reduces every circuit node in two monic
quotients of degree at most precision-1. It does not evaluate rational
weights with nonunit denominators.
"""

from interval_unit_jets import QuadraticJets
from modular_interval_jets import unit_progression_jet


def _validate(number, left, right, displacement, precision):
    if number < 1 or number % 2 == 0 or displacement < 1 or precision < 1:
        raise ValueError("Require positive odd N and positive displacement/precision")
    sizes = []
    for domain in (left, right):
        if not domain:
            raise ValueError("Require nonempty node domains")
        previous, size = 0, 0
        for first, length in domain:
            if first <= previous or first % 2 == 0 or length < 1:
                raise ValueError("Require ordered disjoint positive odd progressions")
            previous = first + 2 * (length - 1)
            size += length
        sizes.append(size)
    if sizes[0] != sizes[1]:
        raise ValueError("Require equal domain cardinalities")
    return sizes[0]


def compressed_inverse_context(number, left, right, displacement, precision):
    size = _validate(number, left, right, displacement, precision)
    degree = min(size, precision - 1)
    ring = QuadraticJets(precision, degree)
    variable = list(ring.constant(0))
    if degree:
        variable[1] = (1, 0)
    variable = tuple(variable)
    shift = (1 << displacement) if displacement < precision else 0
    z = ring.constant(number, shift)
    variable_over_z = ring.multiply(variable, ring.inverse(z))
    normalized_column = ring.constant(1)
    normalized_row = ring.constant(1)
    for first, length in right:
        part = unit_progression_jet(
            ring, (0, length),
            ring.add(ring.constant(first), ring.scale(variable, -1)),
            ring.constant(2))
        part = ring.multiply(part, ring.inverse(ring.constant(*part[0])))
        normalized_column = ring.multiply(normalized_column, part)
    for first, length in left:
        part = unit_progression_jet(
            ring, (0, length),
            ring.add(ring.constant(1), ring.scale(variable_over_z, -first)),
            ring.scale(variable_over_z, -2))
        normalized_row = ring.multiply(normalized_row, part)
    common_inverse = ring.inverse(ring.power(
        ring.add(ring.constant(1), ring.scale(variable, -1)), size - degree))
    column = ring.multiply(normalized_column, common_inverse)
    row = ring.multiply(normalized_row, common_inverse)
    column = ring.multiply(column, ring.inverse(ring.constant(*column[-1])))
    row = ring.multiply(row, ring.inverse(ring.constant(*row[-1])))
    assert column[-1] == row[-1] == (1, 0)
    origin = ring.multiply(column, ring.inverse(row))
    infinity = ring.multiply(tuple(reversed(row)),
                             ring.inverse(tuple(reversed(column))))
    scalar = QuadraticJets(precision, 0)
    pair = scalar.pair_product
    z = (number % scalar.modulus, shift)
    inverse_row_constant = scalar.inverse(scalar.constant(*row[0]))[0]
    row_quotient = [(0, 0)] * (degree + 1)
    z_power = (1, 0)
    for j in range(degree + 1):
        row_quotient[degree - j] = pair(pair(row[j], z_power), inverse_row_constant)
        z_power = pair(z_power, z)
    assert row_quotient[-1] == (1, 0)
    moments = []
    for e in range(degree):
        line = []
        factor = pair(z, (-1 % scalar.modulus, -1 % scalar.modulus))
        for d in range(degree):
            value = [int(e == d), 0]
            for j in range(min(d, e) + 1):
                term = pair(origin[d - j], infinity[e - j])
                value[0] -= term[0]
                value[1] -= term[1]
            line.append(pair(tuple(value), factor))
            factor = pair(factor, z)
        moments.append(tuple(line))
    return {"number": number, "left": left, "right": right,
            "displacement": displacement, "precision": precision,
            "matrix_dimension": size, "quotient_degree": degree,
            "column_packet": column, "reciprocal_row_packet": row,
            "row_quotient": tuple(row_quotient), "column_quotient": column,
            "normalized_moments": tuple(moments)}


class CircuitMomentQuotient:
    def __init__(self, context):
        self.context = context
        self.degree = context["quotient_degree"]
        self.ring = QuadraticJets(context["precision"], 0)
        self.modulus = self.ring.modulus
        self.row = context["row_quotient"]
        self.column = context["column_quotient"]

    def pair_add(self, x, y):
        return ((x[0] + y[0]) % self.modulus, (x[1] + y[1]) % self.modulus)

    def constant(self, value):
        value = value % self.modulus
        return {(0, 0): (value, 0)} if self.degree and value else {}

    def reduce(self, polynomial):
        if not self.degree:
            return {}
        result = {key: (value[0] % self.modulus, value[1] % self.modulus)
                  for key, value in polynomial.items() if value != (0, 0)}
        maximum_row = max((d for d, e in result), default=-1)
        maximum_column = max((e for d, e in result), default=-1)
        for d in range(maximum_row, self.degree - 1, -1):
            for e in range(maximum_column + 1):
                value = result.pop((d, e), (0, 0))
                if value == (0, 0):
                    continue
                for j in range(self.degree):
                    term = self.ring.pair_product(value, self.row[j])
                    key = (d - self.degree + j, e)
                    result[key] = self.pair_add(result.get(key, (0, 0)), (-term[0], -term[1]))
        for e in range(maximum_column, self.degree - 1, -1):
            for d in range(self.degree):
                value = result.pop((d, e), (0, 0))
                if value == (0, 0):
                    continue
                for j in range(self.degree):
                    term = self.ring.pair_product(value, self.column[j])
                    key = (d, e - self.degree + j)
                    result[key] = self.pair_add(result.get(key, (0, 0)), (-term[0], -term[1]))
        return {key: value for key, value in result.items() if value != (0, 0)}

    def add(self, left, right):
        result = left.copy()
        for key, value in right.items():
            result[key] = self.pair_add(result.get(key, (0, 0)), value)
        return {key: value for key, value in result.items() if value != (0, 0)}

    def multiply(self, left, right):
        result = {}
        for (d, e), a in left.items():
            for (f, g), b in right.items():
                key = (d + f, e + g)
                term = self.ring.pair_product(a, b)
                result[key] = self.pair_add(result.get(key, (0, 0)), term)
        return self.reduce(result)

    def evaluate(self, circuit):
        values = []
        for index, node in enumerate(circuit["nodes"]):
            operation = node[0]
            if operation == "row":
                value = self.reduce({(1, 0): (1, 0)})
            elif operation == "column":
                value = self.reduce({(0, 1): (1, 0)})
            elif operation == "const":
                if len(node) != 2 or not isinstance(node[1], int):
                    raise ValueError("Constants must be explicit integers")
                value = self.constant(node[1])
            elif operation in ("add", "mul"):
                if len(node) != 3 or any(not isinstance(i, int) or not 0 <= i < index
                                         for i in node[1:]):
                    raise ValueError("Circuit edges must point to earlier nodes")
                function = self.add if operation == "add" else self.multiply
                value = function(values[node[1]], values[node[2]])
            else:
                raise ValueError("Only row, column, const, add and mul nodes are allowed")
            values.append(value)
        output = circuit["output"]
        if not isinstance(output, int) or not 0 <= output < len(values):
            raise ValueError("Invalid output node")
        return values[output]

    def contract(self, polynomial):
        result = (0, 0)
        for (d, e), value in polynomial.items():
            term = self.ring.pair_product(value, self.context["normalized_moments"][e][d])
            result = self.pair_add(result, term)
        return result


def inverse_circuit_moment(context, circuit):
    quotient = CircuitMomentQuotient(context)
    reduced = quotient.evaluate(circuit)
    residue = quotient.contract(reduced)
    return {"scale_exponent": -context["displacement"], "residue": list(residue),
            "precision": context["precision"],
            "quotient_degree": context["quotient_degree"],
            "circuit_nodes": len(circuit["nodes"]),
            "reduced_terms": len(reduced)}


def monomial_circuit(row_exponent, column_exponent):
    if row_exponent < 0 or column_exponent < 0:
        raise ValueError("Require nonnegative exponents")
    nodes = [["row"], ["column"], ["const", 1]]

    def binary(operation, left, right):
        nodes.append([operation, left, right])
        return len(nodes) - 1

    def power(node, exponent):
        result = 2
        while exponent:
            if exponent & 1:
                result = binary("mul", result, node)
            exponent >>= 1
            if exponent:
                node = binary("mul", node, node)
        return result

    row = power(0, row_exponent)
    column = power(1, column_exponent)
    return {"nodes": nodes, "output": binary("mul", row, column)}
