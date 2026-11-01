#!/usr/bin/env python3
"""Global modular products of circuit values on binary-described boxes.

For a fixed number of variables, product_(a in box)(Z-q(a)) modulo 2^p
has a polynomial-time constructor even when q's formal degree is huge.
The input circuit uses var, const, add and mul nodes. The output circuit
uses x, const, add and mul nodes; binary powering is fully expanded.

The algorithm computes ordinary polynomial coefficients, not only values
at integer points. It does not divide away singular product content.
"""

from __future__ import annotations

from itertools import product
from math import comb, prod


def factorial_valuation(number: int) -> int:
    if number < 0:
        raise ValueError("Require a nonnegative factorial argument")
    return number - number.bit_count()


def prefix_power_sums(stop: int, degree: int) -> tuple[int, ...]:
    """Exact sums over 0 <= x < stop; no numerical interval enumeration."""
    if stop < 0 or degree < 0:
        raise ValueError("Require a nonnegative endpoint and degree")
    result = [stop]
    for exponent in range(1, degree + 1):
        numerator = stop ** (exponent + 1) - sum(
            comb(exponent + 1, lower) * result[lower]
            for lower in range(exponent)
        )
        assert numerator % (exponent + 1) == 0
        result.append(numerator // (exponent + 1))
    return tuple(result)


def interval_power_sums(interval, degree):
    lower, upper = interval
    if lower < 0 or upper < lower:
        raise ValueError("Require a nonnegative half-open interval")
    low = prefix_power_sums(lower, degree)
    high = prefix_power_sums(upper, degree)
    return tuple(b - a for a, b in zip(low, high))


def validate_circuit(circuit, variables):
    """Check the written straight-line program and return a degree bound."""
    if variables < 1:
        raise ValueError("Require at least one variable")
    nodes = circuit["nodes"]
    degrees = []
    for index, instruction in enumerate(nodes):
        operation = instruction[0]
        if operation == "var" and len(instruction) == 2:
            if not isinstance(instruction[1], int) or not 0 <= instruction[1] < variables:
                raise ValueError("Variable index outside the supplied box")
            degrees.append(1)
        elif operation == "const" and len(instruction) == 2:
            if not isinstance(instruction[1], int):
                raise ValueError("Circuit constants must be integers")
            degrees.append(0)
        elif operation in ("add", "mul") and len(instruction) == 3:
            left, right = instruction[1:]
            if not all(isinstance(i, int) and 0 <= i < index for i in (left, right)):
                raise ValueError("Circuit inputs must refer to preceding nodes")
            degrees.append(max(degrees[left], degrees[right]) if operation == "add"
                           else degrees[left] + degrees[right])
        else:
            raise ValueError("Unsupported circuit instruction")
    output = circuit["output"]
    if not isinstance(output, int) or not 0 <= output < len(nodes):
        raise ValueError("Circuit output index is invalid")
    return degrees[output]


class IntegerFunctionQuotient:
    """Pointwise-valid integer functions modulo 2^precision.

    Q(X)=X(X-1)...(X-d+1), with v2(d!) >= precision, vanishes at
    every integer. Reduction by Q in every variable therefore preserves
    every circuit value to the requested precision.
    """

    def __init__(self, variables: int, precision: int):
        if variables < 1 or precision < 1:
            raise ValueError("Require positive dimension and precision")
        self.variables = variables
        self.precision = precision
        self.modulus = 1 << precision
        self.zero_key = (0,) * variables
        degree = 0
        while factorial_valuation(degree) < precision:
            degree += 1
        self.degree = degree
        polynomial = [1]
        for root in range(degree):
            following = [0] * (len(polynomial) + 1)
            for exponent, coefficient in enumerate(polynomial):
                following[exponent] -= root * coefficient
                following[exponent + 1] += coefficient
            polynomial = [value % self.modulus for value in following]
        assert polynomial[-1] == 1 and degree <= 2 * precision
        self.null_polynomial = tuple(polynomial)
        self.nonzero_reduction = tuple((i, value) for i, value in
                                       enumerate(polynomial[:-1]) if value)

    def constant(self, value):
        value %= self.modulus
        return {self.zero_key: value} if value else {}

    def variable(self, axis):
        key = tuple(int(i == axis) for i in range(self.variables))
        return {key: 1}

    def add(self, left, right):
        result = dict(left)
        for key, value in right.items():
            new = (result.get(key, 0) + value) % self.modulus
            if new:
                result[key] = new
            else:
                result.pop(key, None)
        return result

    def scale(self, value, coefficient):
        coefficient %= self.modulus
        if not coefficient:
            return {}
        result = {key: number * coefficient % self.modulus
                  for key, number in value.items()}
        return {key: number for key, number in result.items() if number}

    def multiply(self, left, right):
        if not left or not right:
            return {}
        if set(left) == {self.zero_key}:
            return self.scale(right, left[self.zero_key])
        if set(right) == {self.zero_key}:
            return self.scale(left, right[self.zero_key])
        result = {}
        for first, a in left.items():
            for second, b in right.items():
                key = tuple(x + y for x, y in zip(first, second))
                result[key] = (result.get(key, 0) + a * b) % self.modulus
        result = {key: value for key, value in result.items() if value}
        for axis in range(self.variables):
            for exponent in range(2 * self.degree - 2, self.degree - 1, -1):
                terms = [(key, value) for key, value in result.items()
                         if key[axis] == exponent]
                for key, value in terms:
                    del result[key]
                    for target, coefficient in self.nonzero_reduction:
                        reduced_key = list(key)
                        reduced_key[axis] = exponent - self.degree + target
                        reduced_key = tuple(reduced_key)
                        new = (result.get(reduced_key, 0) - value * coefficient) % self.modulus
                        if new:
                            result[reduced_key] = new
                        else:
                            result.pop(reduced_key, None)
        return result

    def parity_normalized_circuit(self, circuit, parities):
        """Return c and (q(e+2X)-c)/2 without a nonunit ring division."""
        values = []
        for instruction in circuit["nodes"]:
            operation = instruction[0]
            if operation == "var":
                value = parities[instruction[1]], self.variable(instruction[1])
            elif operation == "const":
                parity = instruction[1] % 2
                value = parity, self.constant((instruction[1] - parity) // 2)
            else:
                first, left = values[instruction[1]]
                second, right = values[instruction[2]]
                if operation == "add":
                    parity = (first + second) % 2
                    reduced = self.add(self.add(left, right),
                                       self.constant((first + second - parity) // 2))
                else:
                    assert operation == "mul"
                    parity = first * second
                    reduced = self.add(self.scale(left, second), self.scale(right, first))
                    reduced = self.add(reduced, self.scale(self.multiply(left, right), 2))
                value = parity, reduced
            values.append(value)
        return values[circuit["output"]]

    def box_sum(self, polynomial, sums):
        return sum(value * prod(sums[axis][exponent]
                                for axis, exponent in enumerate(key))
                   for key, value in polynomial.items()) % self.modulus


def circuit_box_moments(circuit, intervals, degree: int, precision: int):
    """Count root parities exactly; compute normalized moments modulo 2^P."""
    intervals = tuple(tuple(interval) for interval in intervals)
    variables = len(intervals)
    formal_degree = validate_circuit(circuit, variables)
    if degree < 0:
        raise ValueError("Require nonnegative moment degree")
    for lower, upper in intervals:
        if lower < 0 or upper < lower:
            raise ValueError("Require nonnegative half-open intervals")
    ring = IntegerFunctionQuotient(variables, precision)
    sections = []
    section_sums = []
    for lower, upper in intervals:
        axis_sections = [((lower - e + 1) // 2, (upper - e + 1) // 2)
                         for e in (0, 1)]
        sections.append(axis_sections)
        section_sums.append([tuple(value % ring.modulus for value in
                                  interval_power_sums(section, ring.degree - 1))
                             for section in axis_sections])
    moments = [[0] * (degree + 1), [0] * (degree + 1)]
    counts = [0, 0]
    nonempty_boxes = 0
    for parities in product((0, 1), repeat=variables):
        count = prod(sections[axis][e][1] - sections[axis][e][0]
                     for axis, e in enumerate(parities))
        if not count:
            continue
        nonempty_boxes += 1
        parity, reduced = ring.parity_normalized_circuit(circuit, parities)
        counts[parity] += count
        sums = [section_sums[axis][e] for axis, e in enumerate(parities)]
        power = ring.constant(1)
        for exponent in range(degree + 1):
            moments[parity][exponent] = (moments[parity][exponent]
                                        + ring.box_sum(power, sums)) % ring.modulus
            if exponent < degree:
                power = ring.multiply(power, reduced)
    assert sum(counts) == prod(upper - lower for lower, upper in intervals)
    assert all(moments[c][0] == counts[c] % ring.modulus for c in (0, 1))
    return {"counts": counts, "moments": moments, "precision": precision,
            "quotient_degree": ring.degree, "quotient_dimension": ring.degree ** variables,
            "nonempty_parity_boxes": nonempty_boxes, "formal_degree_bound": formal_degree}


def _newton_elementaries(moments, degree, working_precision):
    """Return e_j with its varying precision P-v2(j!), never inverting two."""
    if working_precision - factorial_valuation(degree) < 1:
        raise ValueError("Insufficient precision for the Newton recurrence")
    result = [1]
    for index in range(1, degree + 1):
        known_precision = working_precision - factorial_valuation(index - 1)
        numerator = sum((-1) ** (j - 1) * result[index - j] * moments[j]
                        for j in range(1, index + 1)) % (1 << known_precision)
        valuation = (index & -index).bit_length() - 1
        assert numerator % (1 << valuation) == 0
        retained_precision = known_precision - valuation
        modulus = 1 << retained_precision
        result.append((numerator >> valuation) * pow(index >> valuation, -1, modulus) % modulus)
    return tuple(result)


def _elementaries_mod(moments, degree, working_precision, output_precision):
    """Generic elementary residues, with a guard for every Newton division."""
    if working_precision - factorial_valuation(degree) < output_precision:
        raise ValueError("Insufficient guard precision for Newton identities")
    return tuple(value % (1 << output_precision) for value in
                 _newton_elementaries(moments, degree, working_precision))


def _weighted_packet_mod(moments, degree, working_precision, output_precision):
    """Compute (-2)^j e_j using only P=p-1 moment bits.

    The coefficient's error has valuation P+j-v2(j!)=P+popcount(j).
    Thus every positive-degree packet coefficient has p correct bits.
    The constant coefficient is exactly one.
    """
    if degree >= output_precision:
        raise ValueError("Only packet degrees below the output precision are needed")
    if degree and working_precision < output_precision - 1:
        raise ValueError("Insufficient precision for weighted packet coefficients")
    elementary = _newton_elementaries(moments, degree, working_precision)
    return tuple(((-2) ** j * value) % (1 << output_precision)
                 for j, value in enumerate(elementary))


def circuit_box_product(circuit, intervals, precision: int):
    """Construct the entire root product polynomial modulo 2^precision."""
    if precision < 1:
        raise ValueError("Require positive output precision")
    cardinality = prod(upper - lower for lower, upper in intervals)
    degree = min(cardinality, precision - 1)
    if cardinality < 0:
        raise ValueError("Invalid interval")
    working_precision = max(1, precision - 1)
    data = circuit_box_moments(circuit, intervals, degree, working_precision)
    modulus = 1 << precision
    nodes = [["x"]]

    def constant(value):
        nodes.append(["const", value % modulus])
        return len(nodes) - 1

    def binary(operation, left, right):
        nodes.append([operation, left, right])
        return len(nodes) - 1

    def power(base, exponent):
        result = constant(1)
        while exponent:
            if exponent & 1:
                result = binary("mul", result, base)
            exponent >>= 1
            if exponent:
                base = binary("mul", base, base)
        return result

    output = constant(1)
    packets = []
    for parity in (0, 1):
        count = data["counts"][parity]
        retained_degree = min(count, precision - 1)
        coefficients = list(_weighted_packet_mod(data["moments"][parity], retained_degree,
                                                  working_precision, precision))
        packets.append({"parity": parity, "count": count,
                        "retained_degree": retained_degree, "coefficients": coefficients})
        if not count:
            continue
        base = 0 if parity == 0 else binary("add", 0, constant(-1))
        series = constant(1)
        for coefficient in coefficients[1:]:
            series = binary("add", binary("mul", series, base), constant(coefficient))
        factor = binary("mul", power(base, count - retained_degree), series)
        output = binary("mul", output, factor)
    return {"intervals": [list(interval) for interval in intervals], "precision": precision,
            "input_degree": cardinality, "input_circuit_nodes": len(circuit["nodes"]),
            "formal_degree_bound": data["formal_degree_bound"],
            "working_precision": working_precision, "quotient_degree": data["quotient_degree"],
            "quotient_dimension": data["quotient_dimension"], "packets": packets,
            "nodes": nodes, "output": output}


def evaluate_product_circuit(circuit, ring, argument):
    """Evaluate the emitted global circuit in any commutative target algebra."""
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "x":
            value = argument
        elif operation == "const":
            value = ring.constant(instruction[1])
        elif operation == "add":
            value = ring.add(values[instruction[1]], values[instruction[2]])
        else:
            assert operation == "mul"
            value = ring.multiply(values[instruction[1]], values[instruction[2]])
        values.append(value)
    return values[circuit["output"]]


def circuit_box_unit_trace(circuit, intervals, ring, argument):
    """Complete resolvent sum when argument minus each root parity is a unit."""
    data = circuit_box_moments(circuit, intervals, ring.precision - 1,
                              max(1, ring.precision - 1))
    result = ring.constant(0)
    for parity in (0, 1):
        if not data["counts"][parity]:
            continue
        inverse = ring.inverse(ring.add(argument, ring.constant(-parity)))
        power = inverse
        for exponent, moment in enumerate(data["moments"][parity]):
            if exponent == 0:
                moment = data["counts"][parity]
            result = ring.add(result, ring.scale(power, (1 << exponent) * moment))
            if exponent + 1 < ring.precision:
                power = ring.multiply(power, inverse)
    return result


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                        help="JSON with circuit, intervals, and precision fields")
    parser.add_argument("--output", required=True)
    arguments = parser.parse_args()
    supplied = json.loads(Path(arguments.input).read_text())
    result = circuit_box_product(supplied["circuit"], supplied["intervals"], supplied["precision"])
    Path(arguments.output).write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {len(result['nodes'])} nodes for degree {result['input_degree']}.")
