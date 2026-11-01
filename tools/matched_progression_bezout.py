#!/usr/bin/env python3
"""Primitive modular Bezout pairs for deeply matched reflected progressions."""

from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt, prod

from balanced_bezout_modular import _inverse, _linear, _multiply, _scalar_coefficients, _unit_product
from interval_unit_jets import QuadraticJets
from pole_cauchy_conditioning import affine_product, valuation_two


def _factorial_order(n):
    return n - n.bit_count()


@lru_cache(maxsize=128)
def _factorial_unit(n, precision):
    ring = QuadraticJets(precision, 0)
    order, value = affine_product(ring, 1, 1, n, None)
    if order != _factorial_order(n):
        raise ArithmeticError("The factorial valuation disagrees with Legendre")
    return value[0][0]


def _choose(n, k):
    return comb(n, k) if n >= 0 else (-1) ** k * comb(k - n - 1, k)


def _stirling_rows(stop):
    rows = [[1]]
    for _ in range(stop):
        old = rows[-1]
        rows.append([0] + [old[j - 1] + (j * old[j] if j < len(old) else 0)
                           for j in range(1, len(old) + 1)])
    return rows


def _ratio_order(m, r, v, t):
    top = m + r - v - 2
    numerator = _factorial_order(top) - _factorial_order(r - 1) - _factorial_order(m - v - 1)
    denominator = t + _factorial_order(m - 1) + _factorial_order(r - v - 1)
    return numerator - denominator


def _ratio_jet(m, r, v, t, a, eta, precision, degree):
    """Unit jet of the Chu ratio at c+2^t Z; its scale is _ratio_order."""
    ring = QuadraticJets(precision, degree)
    lower, upper = 1 - m, r - v - 1
    length = upper - lower + 1
    value = ring.constant(1)
    covered = 0
    for order in range(t):
        stride, residue = 1 << (order + 1), 1 << order
        first = residue + stride * (-((residue - lower) // stride))
        count = max(0, (upper - first) // stride + 1)
        if not count:
            continue
        base = (m * eta + a * first) >> order
        unit = _linear(ring, base, a << (t - order))
        value = ring.multiply(value, _unit_product(
            ring, (0, count), unit, ring.constant(2 * a)))
        covered += count
    if lower <= 0 <= upper:
        value = ring.multiply(value, _linear(ring, eta, a))
        covered += 1
    if covered != length:
        raise ArithmeticError("The short-ratio valuation buckets did not cover its interval")
    value = ring.scale(value, pow(a, -length, ring.modulus))
    top = m + r - v - 2
    numerator_unit = (_factorial_unit(top, precision)
                      * pow(_factorial_unit(r - 1, precision), -1, ring.modulus)
                      * pow(_factorial_unit(m - v - 1, precision), -1, ring.modulus)) % ring.modulus
    return _scalar_coefficients(ring.scale(ring.inverse(value), (-1) ** (r - 1) * numerator_unit))


def _literal_coefficient(m, r, a, eta, height, exponent, precision):
    """Polynomial-size fallback when m is at most the declared series budget."""
    step = 1 << height
    c, lam = Fraction(m * eta, a), Fraction(step, a)
    nodes = [(j + c) / (1 - lam * j) for j in range(r)]
    coefficient = sum((1 / (prod(node - i for i in range(m))
                                  * prod(node - other for j, other in enumerate(nodes) if j != index))
                       for index, node in enumerate(nodes)), Fraction(0))
    if not coefficient:
        return 0
    scale = exponent + valuation_two(coefficient.numerator) - valuation_two(coefficient.denominator)
    if scale < 0:
        raise ArithmeticError("A normalized matched Newton coefficient is nonintegral")
    if scale >= precision:
        return 0
    numerator = coefficient.numerator >> valuation_two(coefficient.numerator)
    denominator = coefficient.denominator >> valuation_two(coefficient.denominator)
    return ((numerator * pow(denominator, -1, 1 << precision)) << scale) % (1 << precision)


def matched_newton_coefficient(m, t, height, a, eta, r, precision):
    """Return the normalized Newton coefficient at k=m+r-1."""
    top, order = 2 * m - 1, m + r - 1
    exponent = height * top + 2 * _factorial_order(m - 1)
    local_scale = exponent - height * order
    terms = (precision + t + height - t - 1) // (height - t)
    if m <= 2 * (terms - 1) or r <= 2 * (terms - 1):
        if m > 4 * precision:
            raise ValueError("A small-prefix literal fallback is outside the terminal coefficient contract")
        return {"coefficient": _literal_coefficient(
                    m, r, a, eta, height, local_scale, precision),
                "perturbation_terms": terms, "working_precision": precision,
                "fallback_nodes": m + r}
    stirling = _stirling_rows(2 * (terms - 1))
    orders = [_ratio_order(m, r, v, t) for v in range(2 * (terms - 1) + 1)]
    schedule = []
    for s in range(terms):
        for q in range(s + 1):
            first = _choose(r - 2 - q, s - q)
            if not first:
                continue
            for ell in range(s + 1):
                for derivative in range(q + 1):
                    scalar = ((-1) ** (derivative + ell) * first * comb(s, ell)
                              * comb(ell + q, q - derivative))
                    for v, number in enumerate(stirling[ell + derivative]):
                        if not number:
                            continue
                        scale = (local_scale + orders[v] - t * derivative
                                 + height * s + t * (s - ell))
                        if scale < precision:
                            schedule.append((v, derivative, scale, scalar * number, s - ell, 2 * s - ell))
    guard = max(0, -min((item[2] for item in schedule), default=0))
    working, modulus = precision + guard, 1 << (precision + guard)
    jets = {v: _ratio_jet(m, r, v, t, a, eta, working, terms - 1)
            for v in {item[0] for item in schedule}}
    answer = 0
    for v, derivative, scale, scalar, eta_power, a_power in schedule:
        unit = (scalar * jets[v][derivative] * pow(eta, eta_power, modulus)
                * pow(a, -a_power, modulus)) % modulus
        answer = (answer + (unit << (scale + guard))) % modulus
    ring = QuadraticJets(working, 0)
    lam = (1 << height) * pow(a, -1, modulus) % modulus
    c = m * eta * pow(a, -1, modulus) % modulus
    prefactor = _unit_product(ring, (0, r), ring.constant(1), ring.constant(-lam))[0][0]
    prefactor = prefactor * pow((1 + lam * c) % modulus, -(r - 1), modulus) % modulus
    answer = answer * prefactor % modulus
    if answer % (1 << guard):
        raise ArithmeticError("The complete normalized coefficient lost its dyadic integrality")
    return {"coefficient": answer >> guard, "perturbation_terms": terms,
            "working_precision": working, "factorial_ratio_jets": len(jets),
            "retained_series_summands": len(schedule), "fallback_nodes": 0}


def matched_progression_bezout(count, height, center, eta, precision):
    """Complete primitive U,V modulo 2^p for the declared matched family.

A_i=a+2^h i, B_i=a-2^h i, N=a^2+2^h*m*eta. The output basis is
(X-a)^(m-1-j). The common cleared denominator has its exact binary value.
"""
    if any(not isinstance(x, int) for x in (count, height, center, eta, precision)):
        raise ValueError("Require encoded integer parameters")
    if count < 2 or count & (count - 1) or precision < 1:
        raise ValueError("Require a power-of-two count at least two and positive precision")
    t = count.bit_length() - 1
    if height < 2 * t + 1 or center < 1 or center % 2 == 0 or eta % 2 == 0:
        raise ValueError("Require h>=2t+1, an odd positive center, and odd eta")
    if height >= center.bit_length():
        raise ValueError("The reflected progression cannot stay positive at this height")
    step = 1 << height
    if center <= step * (count - 1):
        raise ValueError("Every reflected B node must be positive")
    n = center * center + step * count * eta
    if n <= 0:
        raise ValueError("N must be positive")
    m, a, h, p = count, center, height, precision
    top = 2 * m - 1
    exponent = h * top + 2 * _factorial_order(m - 1)
    modulus = 1 << p
    node_degree = min(m - 1, (p - 1) // h)
    ring = QuadraticJets(p, node_degree)
    numerator = _unit_product(ring, (0, m - 1),
                              _linear(ring, a, -step * m * eta),
                              _linear(ring, -step, -step * a))
    denominator = _unit_product(ring, (0, m - 1), ring.constant(a), ring.constant(-step))
    prefix = _scalar_coefficients(ring.multiply(numerator, ring.inverse(denominator)))
    maximum_offset = min(m - 1, 2 * p - 1)
    u_depth = min(m - 1, maximum_offset + node_degree)
    u_reverse = [0] * (u_depth + 1)
    records = []
    for offset in range(maximum_offset + 1):
        order, r = top - offset, m - offset
        base_gap = h * offset + _factorial_order(offset) - valuation_two(offset + 1)
        perturbation_floor = h - 2 * t + h * offset + _factorial_order(top) - _factorial_order(order)
        if min(base_gap, perturbation_floor) < p:
            record = matched_newton_coefficient(m, t, h, a, eta, r, p)
            coefficient = record["coefficient"]
            records.append({**record, "newton_order": order, "translated_gap": base_gap,
                            "perturbation_floor": perturbation_floor})
            for j, value in enumerate(prefix):
                if offset + j <= u_depth:
                    u_reverse[offset + j] = (u_reverse[offset + j] + coefficient * value) % modulus
        if offset < maximum_offset:
            last = r - 2
            shift = step * (a * last + m * eta) * pow(a - step * last, -1, modulus) % modulus
            shorter = [prefix[0]]
            for j in range(1, len(prefix)):
                shorter.append((prefix[j] + shift * shorter[-1]) % modulus)
            prefix = shorter
    v_depth = min(m - 1, u_depth + node_degree)
    ring = QuadraticJets(p, v_depth)
    p_reverse = _scalar_coefficients(_unit_product(
        ring, (0, m), _linear(ring, 1), _linear(ring, 0, -step)))
    q_reverse = _scalar_coefficients(_unit_product(
        ring, (0, m), _linear(ring, -a, step * m * eta), _linear(ring, step, step * a)))
    ratio = _multiply(p_reverse, _inverse(q_reverse, modulus, v_depth + 1), modulus, v_depth + 1)
    v_reverse = [(-x) % modulus for x in _multiply(u_reverse, ratio, modulus, v_depth + 1)]
    if u_reverse[0] % 2 != 1 or v_reverse[0] % 2 != 1:
        raise ArithmeticError("The complete primitive pair lost its odd leading coefficients")
    while len(u_reverse) > 1 and u_reverse[-1] == 0:
        u_reverse.pop()
    while len(v_reverse) > 1 and v_reverse[-1] == 0:
        v_reverse.pop()
    return {"parameters": {"count": m, "height": h, "center": a, "eta": eta, "precision": p},
            "n": n, "matched_residue_depth": h + t,
            "bezout_exponent": exponent, "polynomial_degree": m - 1, "basis_center": a,
            "u_reversed_coefficients": u_reverse, "v_reversed_coefficients": v_reverse,
            "u_terms": [[m - 1 - j, value] for j, value in enumerate(u_reverse) if value],
            "v_terms": [[m - 1 - j, value] for j, value in enumerate(v_reverse) if value],
            "newton_records": records,
            "maximum_working_precision": max(record["working_precision"] for record in records),
            "enumerated_nodes": max(record["fallback_nodes"] for record in records),
            "scope": "Complete modular primitive pair for the stable deeply matched reflected family; no general inverse-graph aggregation."}


def reflected_root_incidence(count, height, center, eta):
    """Exact common-root incidence on the wider reflected family, with arbitrary eta."""
    if count < 1 or count & (count - 1) or height < count.bit_length() - 1 or center % 2 == 0:
        raise ValueError("Require a power-of-two count, D divisible by count, and odd center")
    if center < 1 or height >= center.bit_length():
        raise ValueError("The reflected progression cannot stay positive at this height")
    step = 1 << height
    if center <= step * (count - 1) or center * center + step * count * eta <= 0:
        raise ValueError("Require positive nodes and positive N")
    square_numerator = -count * eta
    index = None
    if square_numerator >= 0 and square_numerator % step == 0:
        square = square_numerator // step
        candidate = isqrt(square)
        if candidate < count and candidate * candidate == square:
            index = candidate
    return {"gcd_degree": int(index is not None), "matched_index": index,
            "common_root": None if index is None else center + step * index,
            "enumerated_nodes": 0}
