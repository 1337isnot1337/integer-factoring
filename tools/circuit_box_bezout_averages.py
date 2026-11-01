#!/usr/bin/env python3
"""Count-independent precision for circuit-isometric Bezout normalization.

Known-length averages absorb the unique high-valuation factor in each
short factorial ratio. The remaining division loss is less than the
requested precision, even when the implicit node count has enormous
2-adic valuation. These are the same certified separated domains as in
circuit_box_bezout, not a factor-root search oracle.
"""

from __future__ import annotations

from math import factorial

from circuit_box_averages import AverageQuotient, validate_regular_circuit
from circuit_box_bezout import (
    Builder, IsometricContext, _compile_context, _multiply, valuation_two,
)
from circuit_box_products import factorial_valuation, prefix_power_sums


def regular_z_circuits(context):
    """Exact globally regular circuits for the two interleaved node classes."""
    if hasattr(context, "_regular_z_cache"):
        return context._regular_z_cache
    even = context._domain_circuit((context.first_a - 1) >> context.height,
                                   context.step_a >> context.height, context.circuit_a)
    builder = Builder()
    perturbation = builder.import_circuit(context.circuit_b)
    index_map = builder.operation("add", 0, builder.scale(perturbation, 2))
    denominator = builder.affine(index_map, context.first_b, context.step_b)
    numerator = builder.affine(index_map, (context.number - context.first_b) >> context.height,
                               -(context.step_b >> context.height))
    builder.nodes.append(["inv", denominator])
    inverse = len(builder.nodes) - 1
    odd = builder.finish(builder.operation("mul", numerator, inverse))
    validate_regular_circuit(even, 1)
    validate_regular_circuit(odd, 1)
    context._regular_z_cache = even, odd
    return even, odd


def regular_value(circuit, point, modulus):
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = point
        elif operation == "const":
            value = instruction[1]
        elif operation == "add":
            value = values[instruction[1]] + values[instruction[2]]
        elif operation == "mul":
            value = values[instruction[1]] * values[instruction[2]]
        else:
            assert operation == "inv"
            value = pow(values[instruction[1]], -1, modulus)
        values.append(value % modulus)
    return values[circuit["output"]]


def prefix_average_moments(circuit, length, degree, precision):
    """All 2/length sum q(i)^l using one regular-circuit quotient evaluation."""
    if length < 1 or degree < 0 or precision < 1:
        raise ValueError("Invalid normalized moment parameters")
    validate_regular_circuit(circuit, 1)
    ring = AverageQuotient(1, precision)
    sums = prefix_power_sums(length, ring.degree - 1)
    loss = valuation_two(length)
    inverse_odd = pow(length >> loss, -1, ring.modulus)
    functional = []
    for value in sums:
        assert (2 * value) % (1 << loss) == 0
        functional.append(((2 * value) >> loss) * inverse_odd % ring.modulus)
    reduced = ring.evaluate(circuit, (0,))
    power = ring.constant(1)
    result = []
    for order in range(degree + 1):
        result.append(sum(coefficient * functional[key[0]]
                          for key, coefficient in power.items()) % ring.modulus)
        if order < degree:
            power = ring.multiply(power, reduced)
    return result


def normalized_moment_differences(context, length, degree, precision):
    """Integral (1/length) sum (z_i^l-i^l), with no count-valuation guard."""
    if length < 1 or degree < 0 or precision < 1:
        raise ValueError("Invalid normalized difference parameters")
    if not hasattr(context, "_average_difference_cache"):
        context._average_difference_cache = {}
    key = length, degree, precision
    if key in context._average_difference_cache:
        return context._average_difference_cache[key]
    # Two extra bits cover the division by four when length is even.
    average_precision = precision + 2
    average_modulus = 1 << average_precision
    output_modulus = 1 << precision
    even, odd = regular_z_circuits(context)
    even_count, odd_count = (length + 1) // 2, length // 2
    first = prefix_average_moments(even, even_count, degree, average_precision)
    second = (prefix_average_moments(odd, odd_count, degree, average_precision)
              if odd_count else [0] * (degree + 1))
    identity = {"nodes": [["var", 0]], "output": 0}
    reference = prefix_average_moments(identity, length, degree, average_precision)
    result = []
    for order in range(degree + 1):
        if length % 2 == 0:
            numerator = (first[order] + second[order] - 2 * reference[order]) % average_modulus
            if numerator % 4:
                raise ArithmeticError("Isometry failed to make a normalized moment integral")
            value = (numerator >> 2) % output_modulus
        else:
            numerator = (even_count * first[order] + odd_count * second[order]
                         - length * reference[order]) % average_modulus
            if numerator % 2:
                raise ArithmeticError("Normalized odd-length moment lost divisibility")
            value = (numerator >> 1) * pow(length, -1, output_modulus) % output_modulus
        result.append(value)
    context._average_difference_cache[key] = result
    return result


def normalized_parity_difference_by_average(context, order, precision):
    """Normalized parity divided difference using at most about 2p bits.

    For each j, select the most divisible integer M' in M,...,M+j-1,
    with M=k+1. Known-length moments absorb M'. All other factorial
    factors have total valuation v2(s!)+v2((j-1-s)!), independent of k.
    """
    if order < 1 or precision < 1:
        raise ValueError("Require positive divided-difference order and precision")
    output_modulus = 1 << precision
    start = order + 1
    even, odd = regular_z_circuits(context)
    answer = 1
    certificates = []
    boundary_queries = 0
    for degree in range(1, precision):
        selected = max(range(degree), key=lambda offset: valuation_two(start + offset))
        length = start + selected
        denominator = 1
        for offset in range(degree):
            if offset != selected:
                denominator *= start + offset
        loss = valuation_two(denominator)
        assert loss == factorial_valuation(selected) + factorial_valuation(degree - 1 - selected)
        assert loss <= factorial_valuation(degree - 1)
        working = precision + loss
        modulus = 1 << working
        # Every exponential summand has at least one explicit factor
        # of two, so moments need one fewer bit than the numerator.
        moments = normalized_moment_differences(context, length, degree, working - 1)
        generator = [0] * (degree + 1)
        for index in range(1, degree + 1):
            index_loss = valuation_two(index)
            generator[index] = (moments[index] * pow(index >> index_loss, -1, modulus)
                                << (index - 1 - index_loss)) % modulus
        exponential = [0] * (degree + 1)
        power = [1] + [0] * degree
        for index in range(1, degree + 1):
            power = _multiply(power, generator, modulus, degree + 1)
            odd_factorial = factorial(index) >> factorial_valuation(index)
            weight = (pow(length, index - 1, modulus)
                      * pow(odd_factorial, -1, modulus) << index.bit_count()) % modulus
            exponential = [(a + weight * b) % modulus for a, b in zip(exponential, power)]
        correction = [1] + [0] * degree
        for point in range(length, start + degree - 1):
            correction = _multiply(correction, [1, -2 * point], modulus, degree + 1)
        for node_index in range(start, length):
            circuit = even if node_index % 2 == 0 else odd
            value = regular_value(circuit, node_index // 2, modulus)
            correction = _multiply(correction, [1, -2 * value], modulus, degree + 1)
            boundary_queries += 1
        numerator = sum(correction[index] * exponential[degree - index]
                        for index in range(degree + 1)) % modulus
        if numerator % (1 << loss):
            raise ArithmeticError("Weighted divided-difference numerator lost divisibility")
        contribution = (numerator >> loss) * pow(denominator >> loss, -1, output_modulus) % output_modulus
        if contribution % (1 << degree):
            raise ArithmeticError("The weighted contribution must vanish to order j")
        answer = (answer + (-1) ** degree * contribution) % output_modulus
        certificates.append({"degree": degree, "selected_offset": selected,
                             "absorbed_count_valuation": valuation_two(length),
                             "remaining_denominator_valuation": loss,
                             "numerator_precision": working,
                             "moment_precision": working - 1,
                             "average_precision": working + 1,
                             "weighted_contribution": contribution})
    assert answer & 1
    bound = precision + factorial_valuation(max(0, precision - 2)) + 1
    maximum = max((item["average_precision"] for item in certificates), default=precision)
    assert maximum <= bound
    return {"unit": answer, "precision": precision, "order": order,
            "working_precision": maximum,
            "count_independent_precision_bound": bound,
            "short_ratio_certificates": certificates,
            "average_boundary_queries": boundary_queries,
            "method": "known-count normalized averages"}


def circuit_isometric_bezout_by_average(number, first_a, step_a, first_b, step_b,
                                       count, height, precision, circuit_a, circuit_b):
    """The complete normalized Bezout compiler with count-independent guards."""
    if precision < 1:
        raise ValueError("Require positive precision")
    context = IsometricContext(number, first_a, step_a, first_b, step_b,
                                count, height, circuit_a, circuit_b)
    result = _compile_context(context, precision, normalized_parity_difference_by_average)
    result["average_boundary_queries"] = sum(item["average_boundary_queries"]
                                               for item in result["retained_newton_coefficients"])
    result["scope"] = ("Complete globally normalized Bezout pairs on certified circuit-isometric "
                       "separated domains, with count-independent coefficient precision; "
                       "no arbitrary factor-root search.")
    return result
