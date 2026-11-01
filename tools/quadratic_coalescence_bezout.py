#!/usr/bin/env python3
"""Normalized quadratic binomial traces and a two-depth Bezout testbed."""

from interval_unit_jets import QuadraticJets
from pole_cauchy_conditioning import affine_product, valuation_two


def _factorial_unit(n, precision):
    ring = QuadraticJets(precision, 0)
    order, value = affine_product(ring, 1, 1, n, None)
    if order != n - n.bit_count():
        raise ArithmeticError("The factorial valuation must agree with its binary formula")
    return value[0][0]


def _odd(value):
    return value >> valuation_two(value)


def quadratic_binomial_trace(t, w, precision):
    """Return 2^(2-t) sum_j binom(2^t-1,j)^2/(j^2+w) modulo2^p.

    The residue w=1 modulo4 makes the result an odd 2-adic integer.
    Only O(p) initial Mahler coefficients are needed; no domain is listed.
    """
    if (any(not isinstance(x, int) for x in (t, w, precision))
            or t < 1 or precision < 1 or w % 4 != 1):
        raise ValueError("Require t,p>=1 and w=1 modulo four")
    n = (1 << t) - 1
    last = min(n, 4 * precision + 7)
    guard = max(0, max(valuation_two(k + 1) - 1 for k in range(last + 1)))
    working = precision + guard
    modulus, output_modulus = 1 << working, 1 << precision
    values = []
    for j in range(last + 1):
        denominator = j * j + w
        loss = valuation_two(denominator)
        if loss not in (0, 1):
            raise ArithmeticError("The quadratic denominator has an unexpected depth")
        values.append((1 << (1 - loss)) * pow(denominator >> loss, -1, modulus) % modulus)
    mahler = []
    while values:
        mahler.append(values[0])
        values = [(values[j + 1] - values[j]) % modulus for j in range(len(values) - 1)]
    n_unit = _factorial_unit(n, working)
    moment_unit = _factorial_unit(2 * n, working) * pow(n_unit, -2, modulus) % modulus
    total, records = 0, []
    for k, coefficient in enumerate(mahler):
        shift = 1 - valuation_two(k + 1)
        if shift < 0:
            if coefficient % (1 << -shift):
                raise ArithmeticError("A normalized quadratic Mahler term lost integrality")
            term = (coefficient >> -shift) * moment_unit
        else:
            term = (coefficient * moment_unit) << shift
        total = (total + term) % output_modulus
        records.append({"index": k, "moment_valuation": t - valuation_two(k + 1),
                        "normalized_shift": shift, "normalized_term": term % output_modulus})
        if k < last:
            moment_unit = (moment_unit * _odd(n - k) ** 2
                           * pow(_odd(k + 1) * _odd(2 * n - k), -1, modulus)) % modulus
    if total % 2 != 1:
        raise ArithmeticError("The normalized quadratic trace must be odd")
    return {"t": t, "w": w, "precision": precision, "coefficient": total,
            "working_precision": working, "mahler_terms": last + 1,
            "factorial_unit": n_unit % output_modulus,
            "normalized_moment_records": records, "enumerated_domain_nodes": 0}


def quadratic_coalescence_bezout(t, h, a, w, precision):
    """Complete primitive pair in its rigorously charged precision window.

    A_i=a+2^h i, B_i=a-2^h i, N=a^2+2^(2h)w, with 2^t nodes each.
    Both complete modular polynomials have only their leading coefficient.
    """
    if any(not isinstance(x, int) for x in (t, h, a, w, precision)) or t < 1:
        raise ValueError("Use integer parameters and t>=1")
    if h < 2 * t - 1 or not 1 <= precision <= h - 2 * t + 2:
        raise ValueError("Require h>=2t-1 and 1<=p<=h-2t+2")
    m, step = 1 << t, 1 << h
    if a % 2 != 1 or a <= step * (m - 1) or w % 4 != 1:
        raise ValueError("Require positive reflected nodes, odd a, and w=1 modulo four")
    n = a * a + step * step * w
    if n <= 0:
        raise ValueError("The encoded integer N must be positive")
    trace = quadratic_binomial_trace(t, w, precision)
    modulus = 1 << precision
    scalar = trace["coefficient"] * pow(trace["factorial_unit"], -2, modulus) % modulus
    u = a * scalar % modulus
    v = -pow(pow(a, -1, modulus), m - 1, modulus) * scalar % modulus
    return {"parameters": {"t": t, "h": h, "a": a, "w": w, "precision": precision},
            "n": n, "count": m, "basis_center": a,
            "bezout_exponent": 2 * h * m + 2 * m - 3 * t,
            "u_degree": m - 1, "v_degree": m - 1,
            "u_terms": [[m - 1, u]], "v_terms": [[m - 1, v]],
            "u_reversed_coefficients": [u], "v_reversed_coefficients": [v],
            "trace": trace, "enumerated_domain_nodes": 0,
            "scope": "Complete primitive pairs in a two-depth reflected testbed; incidence already has a square test."}
