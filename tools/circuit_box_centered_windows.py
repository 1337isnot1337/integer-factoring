"""Centered Coppersmith queries for dyadic divisor windows and low residues.

The exact lattice and Sturm engines are shared with the quarter-bit endpoint.
This wrapper trades ordinary interval width against required low-bit precision.
"""

from math import isqrt

from circuit_box_ratio_small_roots import bounded_integer_roots, residue_short_polynomial


def window_parameters(number, left, right):
    if number < 3 or not number & 1:
        raise ValueError("an odd integer at least three is required")
    if not 1 <= left <= right <= number + 1:
        raise ValueError("the positive divisor window must lie below N+1")
    factor_exponent = left.bit_length() - 1
    if right > 1 << (factor_exponent + 1):
        raise ValueError("the window must lie within one dyadic factor band")
    width = right - left
    width_exponent = (width - 1).bit_length() if width else 0
    precision = max(0, width_exponent - factor_exponent * factor_exponent // number.bit_length())
    return {
        "factor_exponent": factor_exponent,
        "width_exponent": width_exponent,
        "required_precision": precision,
    }


def centered_window_factors(number, left, right, residue=0, *, precision=None):
    """Return every divisor in the window and the selected residue class."""
    parameters = window_parameters(number, left, right)
    required = parameters["required_precision"]
    precision = required if precision is None else precision
    if not required <= precision <= number.bit_length() + 1:
        raise ValueError("precision is below the proved threshold or above the input bound")
    modulus = 1 << precision
    residue %= modulus
    result = {
        **parameters, "precision": precision, "residue": residue,
        "factors": [], "candidate_count": 0, "direct_divisibility_tests": 0,
        "root_divisibility_tests": 0, "lattice_queries": 0,
        "dimension": 0, "lattice_swaps": 0, "lattice_reductions": 0,
        "sturm_queries": 0, "visited_intervals": 0,
    }
    if left == right:
        return result
    first = left + (residue - left) % modulus
    count = max(0, (right - 1 - first) // modulus + 1)
    result["candidate_count"] = count
    if precision and not residue & 1:
        return result
    if count <= 3:
        for index in range(count):
            candidate = first + modulus * index
            result["direct_divisibility_tests"] += 1
            if number % candidate == 0:
                result["factors"].append(candidate)
        return result

    center_index = count // 2
    center = first + modulus * center_index
    lower, upper = -center_index, count - center_index
    radius_exponent = (center_index - 1).bit_length()
    exponent = parameters["factor_exponent"]
    if radius_exponent > exponent * exponent // number.bit_length() - 1:
        raise AssertionError("centered radius exceeded the proved lattice budget")

    # The old engine uses radius 2**(k+1-effective_precision). Choosing
    # this scaled formal residue preserves the desired constant a0/R mod N.
    effective_precision = exponent + 1 - radius_exponent
    shift = effective_precision - precision
    if shift < 0:
        raise AssertionError("the centered lattice scaling has a negative shift")
    witness = residue_short_polynomial(number, exponent, effective_precision, center << shift)
    if witness["linear_constant"] % number != center * pow(modulus, -1, number) % number:
        raise AssertionError("the centered linear congruence changed")
    roots = bounded_integer_roots(witness["polynomial"], lower, upper)
    for root in roots["roots"]:
        candidate = center + modulus * root
        result["root_divisibility_tests"] += 1
        if left <= candidate < right and number % candidate == 0:
            result["factors"].append(candidate)
    result.update({
        "center": center, "root_interval": [lower, upper],
        "effective_precision": effective_precision,
        "radius_exponent": radius_exponent,
        "lattice_queries": 1, "witness": witness,
    })
    for key in ("dimension", "lattice_swaps", "lattice_reductions"):
        result[key] = witness[key]
    for key in ("sturm_queries", "visited_intervals"):
        result[key] = roots[key]
    return result


def _totals():
    return {
        "residue_queries": 0, "lattice_queries": 0,
        "direct_divisibility_tests": 0, "root_divisibility_tests": 0,
        "maximum_lattice_dimension": 0, "sturm_queries": 0,
        "largest_precision": 0,
    }


def _add(totals, result):
    totals["residue_queries"] += 1
    for key in ("lattice_queries", "direct_divisibility_tests", "root_divisibility_tests", "sturm_queries"):
        totals[key] += result[key]
    totals["maximum_lattice_dimension"] = max(totals["maximum_lattice_dimension"], result["dimension"])
    totals["largest_precision"] = max(totals["largest_precision"], result["precision"])


def centered_factor_rectangle(number, row_interval, column_interval, *, collect_pairs=False):
    """Exact positive-rectangle factor pairs, with all class queries charged."""
    if number < 3 or not number & 1:
        raise ValueError("an odd integer at least three is required")
    row_left, row_right = row_interval
    col_left, col_right = column_interval
    if row_left < 1 or col_left < 1 or row_right < row_left or col_right < col_left:
        raise ValueError("ordered positive intervals are required")
    totals, factors, count, least = _totals(), [], 0, None
    row_right = min(row_right, number + 1)
    left = row_left
    while left < row_right:
        right = min(row_right, 1 << left.bit_length())
        precision = window_parameters(number, left, right)["required_precision"]
        residues = range(1, 1 << precision, 2) if precision else (0,)
        for residue in residues:
            result = centered_window_factors(number, left, right, residue)
            _add(totals, result)
            for factor in result["factors"]:
                partner = number // factor
                if col_left <= partner < col_right:
                    count += 1
                    least = factor if least is None else min(least, factor)
                    if collect_pairs:
                        factors.append((factor, partner))
        left = right
    result = {"count": count, "smallest_row": least, **totals}
    if collect_pairs:
        result["factor_pairs"] = sorted(factors)
    return result


def smallest_factor_by_centered_windows(number):
    if number < 2:
        raise ValueError("input must be at least two")
    totals = _totals()
    if number % 2 == 0:
        return {"factor": 2, **totals}
    left, stop = 2, isqrt(number) + 1
    while left < stop:
        right = min(stop, 1 << left.bit_length())
        precision = window_parameters(number, left, right)["required_precision"]
        residues = range(1, 1 << precision, 2) if precision else (0,)
        least = None
        for residue in residues:
            result = centered_window_factors(number, left, right, residue)
            _add(totals, result)
            if result["factors"]:
                current = min(result["factors"])
                least = current if least is None else min(least, current)
        if least is not None:
            return {"factor": least, **totals}
        left = right
    return {"factor": number, **totals}
