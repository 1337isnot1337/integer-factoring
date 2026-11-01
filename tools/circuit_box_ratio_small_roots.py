"""Exact classical small-root endpoint for dyadic factor residue classes.

The implementation uses rational LLL and Sturm isolation, with explicit
integer certificates for the parameter bound and the short polynomial.
It is a reference algorithm, not an optimized factoring implementation.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd, isqrt, lcm


def _gram_schmidt(basis):
    count, width = len(basis), len(basis[0])
    coefficients = [[Fraction(0) for _ in range(count)] for _ in range(count)]
    orthogonal, norms = [], []
    for i, row in enumerate(basis):
        vector = list(map(Fraction, row))
        for j in range(i):
            coefficients[i][j] = sum(Fraction(row[k]) * orthogonal[j][k] for k in range(width)) / norms[j]
            vector = [vector[k] - coefficients[i][j] * orthogonal[j][k] for k in range(width)]
        orthogonal.append(vector)
        norms.append(sum(value * value for value in vector))
        if not norms[-1]:
            raise ValueError("independent lattice rows required")
    return coefficients, norms


def _short_lattice_vector(basis, threshold):
    original = [row[:] for row in basis]
    basis = [row[:] for row in basis]
    dimension = len(basis)
    transform = [[int(i == j) for j in range(dimension)] for i in range(dimension)]
    coefficients, norms = _gram_schmidt(basis)
    index, swaps, reductions = 1, 0, 0
    while index < dimension:
        for j in range(index - 1, -1, -1):
            value = coefficients[index][j]
            quotient = (2 * value.numerator + value.denominator) // (2 * value.denominator)
            if quotient:
                basis[index] = [x - quotient * y for x, y in zip(basis[index], basis[j])]
                transform[index] = [x - quotient * y for x, y in zip(transform[index], transform[j])]
                coefficients, norms = _gram_schmidt(basis)
                reductions += 1
        if norms[index] >= (Fraction(3, 4) - coefficients[index][index - 1] ** 2) * norms[index - 1]:
            index += 1
        else:
            basis[index], basis[index - 1] = basis[index - 1], basis[index]
            transform[index], transform[index - 1] = transform[index - 1], transform[index]
            coefficients, norms = _gram_schmidt(basis)
            index = max(index - 1, 1)
            swaps += 1
        if sum(map(abs, basis[0])) < threshold:
            break
    vector = basis[0]
    if not any(vector) or sum(map(abs, vector)) >= threshold:
        raise AssertionError("the proved LLL parameter bound did not produce a short polynomial")
    if any(sum(transform[0][i] * original[i][j] for i in range(dimension)) != vector[j] for j in range(dimension)):
        raise AssertionError("incorrect short-vector lattice certificate")
    return vector, transform[0], swaps, reductions


def small_root_parameters(bits, factor_exponent, precision):
    """Choose the first dimension whose worst-case LLL bound is certified."""
    if not 1 <= factor_exponent < bits or precision > factor_exponent:
        raise ValueError("a nontrivial dyadic residue interval is required")
    length_exponent = factor_exponent + 1 - precision
    for dimension in range(2, 2 * bits + 1):
        overhead = (dimension ** (2 * dimension)).bit_length()
        for multiplicity in range(1, dimension):
            budget = (
                4 * factor_exponent * multiplicity * dimension
                - dimension * (dimension - 1)
                - 2 * bits * multiplicity * (multiplicity + 1)
                - 2 * length_exponent * dimension * (dimension - 1)
            )
            if budget >= overhead:
                return dimension, multiplicity
    raise ValueError("ratio precision is below the proved small-root threshold")


def residue_short_polynomial(number, factor_exponent, precision, residue):
    bits = number.bit_length()
    dimension, multiplicity = small_root_parameters(bits, factor_exponent, precision)
    modulus = 1 << precision
    bound = 1 << (factor_exponent + 1 - precision)
    constant = residue * pow(modulus, -1, number) % number
    if constant > number // 2:
        constant -= number
    rows = []
    for degree in range(dimension):
        power = min(degree, multiplicity)
        shift = degree - power
        scale = number ** max(multiplicity - degree, 0)
        row = [0] * dimension
        for j in range(power + 1):
            row[j + shift] = scale * comb(power, j) * constant ** (power - j) * bound ** (j + shift)
        rows.append(row)
    threshold = 1 << (factor_exponent * multiplicity)
    vector, certificate, swaps, reductions = _short_lattice_vector(rows, threshold)
    polynomial = []
    for degree, value in enumerate(vector):
        scale = bound**degree
        if value % scale:
            raise AssertionError("a scaled polynomial coefficient lost integrality")
        polynomial.append(value // scale)
    while polynomial and not polynomial[-1]:
        polynomial.pop()
    return {
        "polynomial": polynomial,
        "lattice_combination": certificate,
        "dimension": dimension,
        "multiplicity": multiplicity,
        "root_bound": bound,
        "linear_constant": constant,
        "weighted_l1_norm": sum(map(abs, vector)),
        "vanishing_threshold": threshold,
        "lattice_swaps": swaps,
        "lattice_reductions": reductions,
    }


def _trim(polynomial):
    while polynomial and not polynomial[-1]:
        polynomial.pop()
    return polynomial


def _primitive(polynomial):
    denominator = 1
    for value in polynomial:
        denominator = lcm(denominator, value.denominator)
    integers = [int(value * denominator) for value in polynomial]
    content = 0
    for value in integers:
        content = gcd(content, value)
    return _trim([value // content for value in integers]) if content else []


def _remainder(first, second):
    result = list(map(Fraction, first))
    while len(result) >= len(second):
        shift = len(result) - len(second)
        multiplier = result[-1] / second[-1]
        for j, value in enumerate(second):
            result[j + shift] -= multiplier * value
        _trim(result)
    return result


def _evaluate(polynomial, argument):
    value = 0
    for coefficient in reversed(polynomial):
        value = value * argument + coefficient
    return value


def bounded_integer_roots(polynomial, lower, upper):
    """All integer roots in [lower,upper), without interval enumeration."""
    polynomial = _trim(list(polynomial))
    if not polynomial:
        raise ValueError("a nonzero polynomial is required")
    if upper < lower:
        raise ValueError("ordered interval required")
    if len(polynomial) == 1 or lower == upper:
        return {"roots": [], "sturm_queries": 0, "visited_intervals": 0}
    derivative = [i * polynomial[i] for i in range(1, len(polynomial))]
    sequence = [_primitive(list(map(Fraction, polynomial))), _primitive(list(map(Fraction, derivative)))]
    while sequence[-1]:
        remainder = _remainder(sequence[-2], sequence[-1])
        following = _primitive([-value for value in remainder])
        if not following:
            break
        sequence.append(following)
    # Every endpoint t-1/D has reduced denominator D greater than the
    # leading coefficient, so it is not a root of the input polynomial.
    denominator = abs(polynomial[-1]) + 1
    queries = 0

    def variation(endpoint):
        nonlocal queries
        queries += 1
        argument = Fraction(endpoint * denominator - 1, denominator)
        signs = []
        for entry in sequence:
            value = _evaluate(entry, argument)
            if value:
                signs.append(1 if value > 0 else -1)
        return sum(left != right for left, right in zip(signs, signs[1:]))

    stack = [(lower, upper, variation(lower), variation(upper))]
    roots, visited = [], 0
    while stack:
        left, right, first_variation, last_variation = stack.pop()
        visited += 1
        if first_variation < last_variation:
            raise AssertionError("negative Sturm count")
        if first_variation == last_variation:
            continue
        if right - left == 1:
            if _evaluate(polynomial, left) == 0:
                roots.append(left)
            continue
        middle = (left + right) // 2
        middle_variation = variation(middle)
        stack.append((middle, right, middle_variation, last_variation))
        stack.append((left, middle, first_variation, middle_variation))
    return {"roots": sorted(roots), "sturm_queries": queries, "visited_intervals": visited}


def dyadic_residue_factors(number, factor_exponent, residue, guard_bits=2):
    if number < 3 or not number & 1 or guard_bits < 2:
        raise ValueError("odd input >=3 and at least two guard bits required")
    bits = number.bit_length()
    if not 0 <= factor_exponent < bits:
        raise ValueError("dyadic factor exponent outside input range")
    precision = (bits + 3) // 4 + guard_bits
    modulus = 1 << precision
    residue %= modulus
    if not residue & 1:
        raise ValueError("odd factor residue required")
    left, right = 1 << factor_exponent, min(1 << (factor_exponent + 1), number + 1)
    root_lower = -((residue - left) // modulus)
    root_upper = -((residue - right) // modulus)
    result = {"factors": [], "precision": precision, "dimension": 0, "lattice_swaps": 0, "lattice_reductions": 0, "sturm_queries": 0, "visited_intervals": 0}
    if root_upper <= root_lower:
        return result
    if root_upper - root_lower == 1:
        candidate = residue + modulus * root_lower
        if number % candidate == 0:
            result["factors"] = [candidate]
        return result
    witness = residue_short_polynomial(number, factor_exponent, precision, residue)
    root_data = bounded_integer_roots(witness["polynomial"], root_lower, root_upper)
    for root in root_data["roots"]:
        candidate = residue + modulus * root
        if left <= candidate < right and number % candidate == 0:
            result["factors"].append(candidate)
    for key in ("dimension", "lattice_swaps", "lattice_reductions"):
        result[key] = witness[key]
    for key in ("sturm_queries", "visited_intervals"):
        result[key] = root_data[key]
    result["root_polynomial_degree"] = len(witness["polynomial"]) - 1
    result["weighted_l1_norm_bits"] = witness["weighted_l1_norm"].bit_length()
    return result


def smallest_factor_by_small_roots(number):
    if number < 2:
        raise ValueError("input must be at least two")
    totals = {"residue_queries": 0, "lattice_queries": 0, "maximum_lattice_dimension": 0, "sturm_queries": 0}
    if number % 2 == 0:
        return {"factor": 2, **totals}
    bits = number.bit_length()
    precision = (bits + 3) // 4 + 2
    modulus = 1 << precision
    stop = isqrt(number)
    for exponent in range(1, bits):
        if 1 << exponent > stop:
            break
        least = None
        for residue in range(1, modulus, 2):
            result = dyadic_residue_factors(number, exponent, residue)
            totals["residue_queries"] += 1
            totals["lattice_queries"] += bool(result["dimension"])
            totals["maximum_lattice_dimension"] = max(totals["maximum_lattice_dimension"], result["dimension"])
            totals["sturm_queries"] += result["sturm_queries"]
            for factor in result["factors"]:
                if 2 <= factor <= stop:
                    least = factor if least is None else min(least, factor)
        if least is not None:
            return {"factor": least, **totals}
    return {"factor": number, **totals}
