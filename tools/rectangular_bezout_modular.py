#!/usr/bin/env python3
"""Complete normalized Bezout polynomials on arbitrary separated rectangles."""

from functools import lru_cache

from balanced_bezout_modular import _inverse, _linear, _multiply, _scalar_coefficients
from interval_unit_jets import QuadraticJets, interval_elementaries
from pole_cauchy_conditioning import affine_product, valuation_two
from rectangular_bezout_scale import rectangular_newton_window


@lru_cache(maxsize=32)
def _elementaries(length, degree):
    return interval_elementaries(0, length, degree)


def unit_progression_general(ring, length, unit, step):
    """Unit progression with even constant step and arbitrary higher coefficients.

A coefficient of degree at most d in W^j uses at most d nonconstant
factors. If W(0) has valuation c>=1, j>d+(P-1)//c vanishes mod2^P.
A zero constant gives the exact nilpotence cutoff d.
"""
    if length < 0:
        raise ValueError("Require a nonnegative encoded length")
    if not length:
        return ring.constant(1)
    increment = ring.multiply(step, ring.inverse(unit))
    if all(x == 0 for pair in increment for x in pair):
        return ring.power(unit, length)
    constants = [x for x in increment[0] if x]
    degree = ring.size - 1
    if not constants:
        cutoff = degree
    else:
        order = min(valuation_two(x) for x in constants)
        if order < 1:
            raise ValueError("The progression increment must have an even constant coefficient")
        cutoff = degree + (ring.precision - 1) // order
    cutoff = min(length, cutoff)
    elementary = _elementaries(length, cutoff)
    value = power = ring.constant(1)
    for j in range(1, cutoff + 1):
        power = ring.multiply(power, increment)
        value = ring.add(value, ring.scale(power, elementary[j]))
    return ring.multiply(ring.power(unit, length), value)


def _child_product(ring, count, n, first, step, reciprocal):
    """Product (1-x_i T) on even-child or reciprocal odd-child indices."""
    if not reciprocal:
        return unit_progression_general(ring, count,
            _linear(ring, 1, -(first - 1) // 4), _linear(ring, 0, -step // 4))
    numerator = unit_progression_general(ring, count,
        _linear(ring, first, (3 * first - n) // 4), _linear(ring, step, 3 * step // 4))
    denominator = unit_progression_general(ring, count, ring.constant(first), ring.constant(step))
    return ring.multiply(numerator, ring.inverse(denominator))


def low_surplus_binomials(n, first, step, count, precision, *, reciprocal=False):
    """C_(a,a+i)=a! [x_0,...,x_a] binom(X,a+i), for 0<=i<p."""
    if count < 1 or precision < 1:
        raise ValueError("Require a nonempty isometric prefix and positive precision")
    a = count - 1
    denominator = 1
    for i in range(1, precision):
        denominator *= a + i
    guard = valuation_two(denominator)
    ring = QuadraticJets(precision + guard, precision - 1)
    inverse_nodes = _scalar_coefficients(ring.inverse(_child_product(
        ring, count, n, first, step, reciprocal)))
    elementary = interval_elementaries(0, a, precision - 1)
    falling = [(-1) ** i * value % ring.modulus for i, value in enumerate(elementary)]
    denominator, values = 1, []
    modulus = 1 << precision
    for i in range(precision):
        if i:
            last = a + i - 1
            falling = [(falling[j] - (last * falling[j - 1] if j else 0)) % ring.modulus
                       for j in range(precision)]
            denominator *= a + i
        coefficient = sum(falling[j] * inverse_nodes[i - j] for j in range(i + 1)) % ring.modulus
        loss = valuation_two(denominator)
        if coefficient % (1 << loss):
            raise ArithmeticError("An isometry binomial divided difference lost integrality")
        values.append((coefficient >> loss) * pow(denominator >> loss, -1, modulus) % modulus)
    if values[0] != 1:
        raise ArithmeticError("The leading low-surplus coefficient must equal one")
    return {"coefficients": values, "precision": precision,
            "working_precision": ring.precision, "factorial_ratio_guard": guard,
            "enumerated_nodes": 0}


def normalized_rectangle_unit(a, b, first_coefficients, second_coefficients, precision):
    """d / [(-1)^b binom(a+b,a)/prod_(v=-a)^b(1+2v)] modulo2^p."""
    modulus = 1 << precision
    left = right = rising = 1
    left_inverses, right_inverses, risings = [1], [1], [1]
    for i in range(1, precision):
        left = left * (1 - 2 * a - 2 * i) % modulus
        right = right * (1 + 2 * b + 2 * i) % modulus
        rising = rising * (a + b + i) % modulus
        left_inverses.append(pow(left, -1, modulus))
        right_inverses.append(pow(right, -1, modulus))
        risings.append(rising)
    answer = 0
    for i in range(precision):
        for j in range(precision - i):
            term = ((-1) ** j * risings[i + j] * first_coefficients[i]
                    * second_coefficients[j] * left_inverses[i] * right_inverses[j])
            answer = (answer + (term << (i + j))) % modulus
    if answer % 2 != 1:
        raise ArithmeticError("The normalized rectangular divided difference must be odd")
    return answer


def _factorial_unit(n, precision):
    ring = QuadraticJets(precision, 0)
    order, value = affine_product(ring, 1, 1, n, None)
    if order != n - n.bit_count():
        raise ArithmeticError("The factorial unit query failed its independent valuation check")
    return value[0][0]


def rectangular_bezout_modular(n, first_a, step_a, rows, first_b, step_b, columns, precision):
    """Complete optimal primitive U,V mod2^p; degrees columns-1 and rows-1."""
    parameters = (n, first_a, step_a, rows, first_b, step_b, columns, precision)
    if any(not isinstance(x, int) for x in parameters) or min(parameters) < 1:
        raise ValueError("Require positive encoded domain data and numerical precision")
    if n % 4 != 3 or first_a % 4 != 1 or first_b % 4 != 1:
        raise ValueError("Require separated residue-one domains and N=3 modulo four")
    if valuation_two(step_a) != 2 or valuation_two(step_b) != 2:
        raise ValueError("Both progression steps must have valuation exactly two")
    window = rectangular_newton_window(rows, columns, precision)
    gaps = {term["index"]: term["valuation_gap"] for term in window["retained_newton_terms"]}
    first_index, exponent = window["first_newton_index"], window["exponent"]
    p, modulus = precision, 1 << precision
    b, top_a = columns - 1, rows - 1
    scalar_ring = QuadraticJets(p, 0)
    column_product = unit_progression_general(
        scalar_ring, columns, scalar_ring.constant(first_b), scalar_ring.constant(step_b))[0][0]
    q_leading_inverse = pow((-1) ** columns * column_product % modulus, -1, modulus)
    odd_product = unit_progression_general(scalar_ring, rows + columns - 1,
        scalar_ring.constant(1 - 2 * top_a), scalar_ring.constant(2))[0][0]
    odd_product_inverse = pow(odd_product, -1, modulus)
    binomial_unit = (_factorial_unit(top_a + b, p)
                     * pow(_factorial_unit(top_a, p), -1, modulus)
                     * pow(_factorial_unit(b, p), -1, modulus)) % modulus
    second = low_surplus_binomials(n, first_b, step_b, columns, p, reciprocal=True)
    node_degree = min(rows - 1, (p - 1) // 2)
    ring = QuadraticJets(p, node_degree)
    prefix = _scalar_coefficients(unit_progression_general(
        ring, rows - 1, _linear(ring, 1, 1 - first_a), _linear(ring, 0, -step_a)))
    v_reverse = [0] * (window["reverse_v_depth"] + 1)
    records = []
    for offset in range(rows - first_index):
        k = rows - 1 - offset
        if k in gaps:
            gap = gaps[k]
            local_precision = p - gap
            first = low_surplus_binomials(n, first_a, step_a, k + 1, local_precision)
            unit = normalized_rectangle_unit(
                k, b, first["coefficients"], second["coefficients"], local_precision)
            coefficient = ((-(-1) ** b * binomial_unit * odd_product_inverse
                            * unit * q_leading_inverse) << gap) % modulus
            for j, value in enumerate(prefix):
                if offset + j < len(v_reverse):
                    v_reverse[offset + j] = (v_reverse[offset + j] + coefficient * value) % modulus
            records.append({"index": k, "valuation_gap": gap, "scaled_newton_coefficient": coefficient,
                            "local_precision": local_precision, "normalized_difference_unit": unit,
                            "row_working_precision": first["working_precision"]})
        if k > first_index:
            shift = first_a - 1 + step_a * (k - 1)
            shorter = [prefix[0]]
            for j in range(1, len(prefix)):
                shorter.append((prefix[j] + shift * shorter[-1]) % modulus)
            prefix = shorter
            numerator, denominator = k, k + b
            numerator >>= valuation_two(numerator)
            denominator >>= valuation_two(denominator)
            binomial_unit = binomial_unit * numerator * pow(denominator, -1, modulus) % modulus
            odd_product_inverse = odd_product_inverse * (1 - 2 * k) % modulus
    depth = window["reverse_u_depth"]
    # The complete ratio Q*/P* has degree-j valuation at least j.
    # Its tail from degree p vanishes, even when the output defect is large.
    input_depth = min(depth, p - 1)
    ring = QuadraticJets(p, input_depth)
    p_reverse = _scalar_coefficients(unit_progression_general(
        ring, rows, _linear(ring, 1, 1 - first_a), _linear(ring, 0, -step_a)))
    q_reverse = _scalar_coefficients(unit_progression_general(
        ring, columns, _linear(ring, -first_b, n - first_b), _linear(ring, -step_b, -step_b)))
    ratio = _multiply(q_reverse, _inverse(p_reverse, modulus, input_depth + 1),
                      modulus, input_depth + 1)
    u_reverse = [(-x) % modulus for x in _multiply(v_reverse, ratio, modulus, depth + 1)]
    if not any(x % 2 for x in u_reverse) or not any(x % 2 for x in v_reverse):
        raise ArithmeticError("The optimal modular Bezout pair must retain primitive coefficients")
    while len(u_reverse) > 1 and u_reverse[-1] == 0:
        u_reverse.pop()
    while len(v_reverse) > 1 and v_reverse[-1] == 0:
        v_reverse.pop()
    return {"parameters": {"n": n, "first_a": first_a, "step_a": step_a, "rows": rows,
                            "first_b": first_b, "step_b": step_b, "columns": columns,
                            "precision": p},
            "bezout_exponent": exponent, "scale_and_window": window, "basis_center": 1,
            "u_degree": columns - 1, "v_degree": rows - 1,
            "u_reversed_coefficients": u_reverse, "v_reversed_coefficients": v_reverse,
            "u_terms": [[columns - 1 - j, x] for j, x in enumerate(u_reverse) if x],
            "v_terms": [[rows - 1 - j, x] for j, x in enumerate(v_reverse) if x],
            "newton_records": records, "column_surplus_working_precision": second["working_precision"],
            "maximum_working_precision": max([second["working_precision"]]
                                              + [x["row_working_precision"] for x in records]),
            "enumerated_progression_nodes": 0,
            "scope": "Complete primitive modular Bezout pairs on independent separated arithmetic rectangles; no unknown-root search."}
