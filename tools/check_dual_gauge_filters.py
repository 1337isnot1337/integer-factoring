#!/usr/bin/env python3
"""Exact prime-field dual gauges and divided-translation mode filters.

All source evaluations retain their integer denominator inside the
modulus. No denominator divisible by the characteristic is inverted.
This checks identities, not an asymptotic crossing-rank improvement.
"""

from collections import defaultdict
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import argparse
import json

from check_partition_channels import partitions, partition_denominator, rank_mod
from check_modular_crossing import modular_operator
from check_integral_schur_crossing import scaled_duals, piece_matrices
from check_schur_fusion_gram import characters


def binomial(upper, lower):
    if upper >= 0:
        return comb(upper, lower) if lower <= upper else 0
    return (-1) ** lower * comb(lower - upper - 1, lower)


@lru_cache(None)
def scaled_state(weight, index):
    denominator = factorial(weight)
    return tuple((mu, denominator * value // partition_denominator(mu))
                 for mu, value in zip(partitions(weight), characters(weight)[index])
                 if value)


@lru_cache(None)
def raise_monomial(parts, degree):
    """Coefficient of z^degree in exp(zU)p_parts, in integer coordinates."""
    if not parts:
        return (((), 1),) if degree == 0 else ()
    result = defaultdict(int)
    first = parts[0]
    for shift in range(degree + 1):
        coefficient = comb(first + shift - 1, shift)
        for tail, value in raise_monomial(parts[1:], degree - shift):
            result[tuple(sorted((first + shift,) + tail))] += coefficient * value
    return tuple((mu, value) for mu, value in result.items() if value)


def raise_polynomial(polynomial, degree):
    result = defaultdict(int)
    for mu, coefficient in polynomial:
        for nu, value in raise_monomial(mu, degree):
            result[nu] += coefficient * value
    return tuple((mu, value) for mu, value in result.items() if value)


def schur_coordinates(polynomial, denominator, weight):
    values = dict(polynomial)
    result = []
    for row in characters(weight):
        value = sum(values.get(mu, 0) * coefficient
                    for mu, coefficient in zip(partitions(weight), row))
        assert value % denominator == 0
        result.append(value // denominator)
    return result


def source_entry(polynomial, denominator, input_weight, input_index,
                 output_weight, output_index, prime):
    input_denominator = factorial(input_weight)
    scale = denominator * input_denominator
    modulus = prime * scale
    output_character = characters(output_weight)[output_index]
    value = sum(c * d * character * modular_operator(modulus, a, rho, b)
                for a, c in polynomial
                for b, d in scaled_state(input_weight, input_index)
                for rho, character in zip(partitions(output_weight), output_character)
                if character) % modulus
    assert value % scale == 0
    return value // scale


def check_projection(prime, weight):
    parameters = [(w, i) for w in range(weight + 1)
                  for i in range(len(partitions(w)))]
    dimension = len(partitions(weight))
    projections = []
    vacuum = []
    translated = []
    integrality = 0
    for w, index in parameters:
        polynomial = scaled_state(w, index)
        raised = raise_polynomial(polynomial, weight - w) if w else ()
        coordinates = schur_coordinates(raised, factorial(w), weight)
        projections.append([value % prime for value in coordinates])
        integrality += len(coordinates)
        translated.append(raised)
        vacuum.append([source_entry(polynomial, factorial(w), weight, v, 0, 0, prime)
                       for v in range(dimension)])
    gram = vacuum[-dimension:]
    assert rank_mod(gram, prime) == dimension
    for projected, actual in zip(projections, vacuum):
        assert [sum(projected[j] * gram[j][v] for j in range(dimension)) % prime
                for v in range(dimension)] == actual
    assert projections[-dimension:] == [
        [int(i == j) for j in range(dimension)] for i in range(dimension)]

    # The displayed kernel basis has one unique lower-degree coordinate.
    kernel = []
    for index in range(len(parameters) - dimension):
        row = [0] * len(parameters)
        row[index] = 1
        for j, value in enumerate(projections[index]):
            row[len(parameters) - dimension + j] = -value % prime
        assert all(sum(row[k] * vacuum[k][v] for k in range(len(parameters))) % prime == 0
                   for v in range(dimension))
        kernel.append(row)
    assert rank_mod(kernel, prime) == len(parameters) - dimension

    mode_entries = 0
    for (w, index), raised in zip(parameters, translated):
        original = scaled_state(w, index)
        for input_weight in (0, 1):
            for input_index in range(len(partitions(input_weight))):
                for output_index in range(len(partitions(input_weight + weight))):
                    arguments = (factorial(w), input_weight, input_index,
                                 input_weight + weight, output_index, prime)
                    assert source_entry(original, *arguments) == source_entry(raised, *arguments)
                    mode_entries += 1
        for input_index in range(len(partitions(weight + 1))):
            arguments = (factorial(w), weight + 1, input_index, 1, 0, prime)
            assert source_entry(original, *arguments) == source_entry(raised, *arguments)
            mode_entries += 1
    return {"prime": prime, "weight": weight, "dimension": dimension,
            "parameter_dimension": len(parameters), "kernel_dimension": len(kernel),
            "integral_projection_coordinates": integrality,
            "direct_vacuum_entries": len(parameters) * dimension,
            "direct_exceptional_mode_entries": mode_entries}


def filter_coefficients(prime, period, locality):
    residue = -locality % period
    return [((-1) ** (k - residue) * comb(k, residue)) % prime
            if k >= residue else 0 for k in range(period)]


def check_newton_filters():
    mode_cases = laurent_cases = 0
    for prime, periods in ((2, (1, 2, 4, 8)), (3, (1, 3, 9)), (5, (1, 5, 25))):
        for period in periods:
            for locality in (0, 1, 2, 3, 5, 8, 13):
                coefficients = filter_coefficients(prime, period, locality)
                for index in range(-3 * period, 3 * period + 1):
                    value = sum(c * binomial(index, k)
                                for k, c in enumerate(coefficients)) % prime
                    assert value == int((index + locality) % period == 0)
                    mode_cases += 1
                # Direct divided Leibniz rule on z^s, with no Taylor-series shortcut.
                leibniz = [sum(coefficients[k] * binomial(-locality, k - j)
                               for k in range(j, period)) % prime
                           for j in range(period)]
                for exponent in range(-2 * period, 2 * period + 1):
                    result = defaultdict(int)
                    for j, coefficient in enumerate(leibniz):
                        coefficient *= binomial(exponent, j)
                        for k in range(j + 1):
                            result[exponent - j + k] += (
                                coefficient * comb(j, k) * (-1) ** (j - k))
                    result = {s: value % prime for s, value in result.items() if value % prime}
                    assert result == {period * (exponent // period): 1}
                    laurent_cases += 1
    return {"periodic_mode_filter_cases": mode_cases,
            "independent_divided_Leibniz_cases": laurent_cases}


def filter_polynomial(polynomial, prime, period, state_weight):
    result = defaultdict(int)
    for mu, c in polynomial:
        coefficients = filter_coefficients(prime, period, sum(mu) + state_weight)
        for k, coefficient in enumerate(coefficients):
            for nu, value in raise_monomial(mu, k):
                result[nu] += c * coefficient * value
    return tuple((mu, value) for mu, value in result.items() if value)


def source_pieces(weight, duals, denominator, prime=2):
    """Generalized source angle for duals with unrestricted finite weights."""
    basis = partitions(weight)
    dimension = len(basis)
    maximum = max(sum(mu) for dual in duals for mu, _ in dual)
    scale = denominator * factorial(weight) ** 2
    modulus = prime * scale
    states = [scaled_state(weight, i) for i in range(dimension)]
    chars = characters(weight)
    shifts = range(-maximum - weight, 2 * weight + 1)
    pieces = {s: [[0] * dimension ** 2 for _ in range(dimension ** 2)] for s in shifts}
    for energy in range(maximum + 3 * weight + 1):
        indices = [(i, energy - i) for i in range(maximum + weight + 1)
                   if energy - i in pieces]
        for intermediate in partitions(energy):
            right = [sum(c * chars[u][j] * modular_operator(modulus, mu, rho, intermediate)
                         for mu, c in states[v] for j, rho in enumerate(basis)) % modulus
                     for v in range(dimension) for u in range(dimension)]
            right = [(j, value) for j, value in enumerate(right) if value]
            if not right:
                continue
            left = []
            for dual in duals:
                for state in states:
                    values = defaultdict(int)
                    for mu, c in dual:
                        values[sum(mu)] += c * sum(
                            d * modular_operator(modulus, mu, intermediate, nu)
                            for nu, d in state)
                    left.append(values)
            for i, shift in indices:
                matrix = pieces[shift]
                for row, values in enumerate(left):
                    coefficient = sum(binomial(weight + w, i) * value
                                      for w, value in values.items()) % modulus
                    if coefficient:
                        for column, value in right:
                            matrix[row][column] = (matrix[row][column]
                                                  + (-1) ** i * coefficient * value) % modulus
    integral_entries = 0
    for matrix in pieces.values():
        for row in matrix:
            for j, value in enumerate(row):
                assert value % scale == 0
                row[j] = value // scale
                integral_entries += 1
    total = [[sum(matrix[i][j] for matrix in pieces.values()) % prime
              for j in range(dimension ** 2)] for i in range(dimension ** 2)]
    assert total == [[int(i == j) for j in range(dimension ** 2)]
                     for i in range(dimension ** 2)]
    return pieces, integral_entries


def check_filtered_crossings():
    records = []
    for weight, period in ((1, 4), (2, 4), (2, 8)):
        duals = scaled_duals(weight, homogeneous=True)
        filtered = tuple(filter_polynomial(dual, 2, period, weight) for dual in duals)
        actual, integral_entries = source_pieces(weight, filtered, factorial(weight))
        original = piece_matrices(weight, homogeneous=True)
        dimension = len(partitions(weight)) ** 2
        expected = defaultdict(lambda: [[0] * dimension for _ in range(dimension)])
        for shift, matrix in original.items():
            target = expected[period * (shift // period)]
            for i, row in enumerate(matrix):
                for j, value in enumerate(row):
                    target[i][j] ^= value
        zero = [[0] * dimension for _ in range(dimension)]
        for shift in set(actual) | set(expected):
            assert actual.get(shift, zero) == expected.get(shift, zero)
        nonzero = [s for s, matrix in actual.items() if any(any(row) for row in matrix)]
        records.append({"weight": weight, "period": period,
                        "maximum_parameter_weight": max(sum(mu) for dual in filtered for mu, _ in dual),
                        "nonzero_shifts": nonzero,
                        "complete_integrality_entries": integral_entries,
                        "interval_identity_entries": len(set(actual) | set(expected)) * dimension ** 2})
        modular_operator.cache_clear()
    return records


def verify():
    projections = [check_projection(prime, weight) for prime, weight in
                   ((2, 1), (2, 2), (2, 4), (3, 1), (3, 3), (5, 1), (5, 5))]
    # The prime-power hypothesis is necessary for this projection formula.
    original = scaled_state(1, 0)
    raised = raise_polynomial(original, 1)
    counterexample = [source_entry(poly, 1, 2, 0, 0, 0, 3)
                      for poly in (original, raised)]
    assert counterexample == [2, 0]
    result = {"scope": "exact dual parametrization and coefficient coarsening; no improved asymptotic rank",
              "bounded_weight_projection_cases": projections,
              "non_prime_power_negative_control": {"prime": 3, "weight": 2,
                                                     "original_and_projected_vacuum": counterexample}}
    result.update(check_newton_filters())
    result["direct_filtered_schur_crossings"] = check_filtered_crossings()
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = verify()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(rendered, end="")
