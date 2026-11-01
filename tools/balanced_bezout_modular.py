#!/usr/bin/env python3
"""Complete modular Bezout polynomials on separated arithmetic progressions.

The optimal dyadic exponent and polynomial degrees use binary encoding.
Only O(precision) highest powers of X-1 can survive in either polynomial.
No progression node, factorial, or unscaled Bezout polynomial is expanded.
"""

from functools import lru_cache

from interval_unit_jets import QuadraticJets, interval_elementaries
from pole_cauchy_conditioning import affine_product, valuation_two


@lru_cache(maxsize=32)
def _elementaries(start, stop, degree):
    return interval_elementaries(start, stop, degree)


def _unit_product(ring, interval, unit, step):
    """The unit-progression identity with its nilpotent-degree cutoff charged.

If the increment's constant has valuation c and all coefficients have
valuation at least v, a power k has coefficientwise valuation at least
c*k-(c-v)*degree. A zero constant gives the usual degree+1 nilpotence.
"""
    start, stop = interval
    length = stop - start
    if not length:
        return ring.constant(1)
    increment = ring.multiply(step, ring.inverse(unit))
    nonzero = [x for pair in increment for x in pair if x]
    if not nonzero:
        return ring.power(unit, length)
    minimum = min(valuation_two(x) for x in nonzero)
    if minimum < 1:
        raise ArithmeticError("A unit progression step is not even coefficientwise")
    constants = [x for x in increment[0] if x]
    if not constants:
        cutoff = ring.size - 1
    else:
        constant_order = min(valuation_two(x) for x in constants)
        cutoff = (ring.precision - 1 + (constant_order - minimum) * (ring.size - 1)) // constant_order
    cutoff = min(cutoff, (ring.precision - 1) // minimum, length)
    elementary = _elementaries(start, stop, cutoff)
    series = power = ring.constant(1)
    for j in range(1, cutoff + 1):
        power = ring.multiply(power, increment)
        series = ring.add(series, ring.scale(power, elementary[j]))
    return ring.multiply(ring.power(unit, length), series)


def _linear(ring, constant, slope=0):
    value = list(ring.constant(constant))
    if ring.size > 1:
        value[1] = (slope % ring.modulus, 0)
    return tuple(value)


def _scalar_coefficients(jet):
    if any(b for _, b in jet):
        raise ArithmeticError("A scalar construction acquired a quadratic component")
    return [a for a, _ in jet]


def _multiply(a, b, modulus, size):
    result = [0] * size
    for i, x in enumerate(a[:size]):
        if x:
            for j, y in enumerate(b[:size - i]):
                result[i + j] = (result[i + j] + x * y) % modulus
    return result


def _inverse(a, modulus, size):
    inverse_constant = pow(a[0], -1, modulus)
    result = [inverse_constant]
    for k in range(1, size):
        result.append(-inverse_constant * sum(
            a[j] * result[k - j] for j in range(1, min(k + 1, len(a)))) % modulus)
    return result


def _validate(n, first_a, step_a, first_b, step_b, count, height, precision):
    values = (n, first_a, step_a, first_b, step_b, count, height, precision)
    if any(not isinstance(value, int) for value in values):
        raise ValueError("All parameters must be encoded integers")
    if min(n, first_a, first_b, count, height, precision) < 1:
        raise ValueError("Require positive data, count, height, and numerical precision")
    if step_a <= 0 or step_b <= 0:
        raise ValueError("Require positive progression steps")
    if height >= n.bit_length():
        raise ValueError("The separating height exceeds the encoded size of N")
    modulus = 1 << (height + 1)
    if first_a % modulus != 1 or first_b % modulus != 1:
        raise ValueError("Both progressions must start at one modulo 2^(height+1)")
    if valuation_two(step_a) != height + 1 or valuation_two(step_b) != height + 1:
        raise ValueError("Both step valuations must equal height+1")
    if n % modulus != 1 + (1 << height):
        raise ValueError("N must separate the two root classes at the declared height")


def _node_product(ring, node_count, n, first_a, step_a, first_b, step_b, height,
                  *, root_shifts=False):
    """Product (1-z_i*T), or (1-(x_i-1)*T), on an interleaved prefix.

Here x_(2i)=A_i, x_(2i+1)=N/B_i, and z_i=(x_i-1)/2^height.
All constant terms are units. Every progression step is even as a jet.
"""
    a_count, b_count = (node_count + 1) // 2, node_count // 2
    divisor = 1 if root_shifts else 1 << height
    a_unit = _linear(ring, 1, -(first_a - 1) // divisor)
    a_step = _linear(ring, 0, -step_a // divisor)
    a_product = _unit_product(ring, (0, a_count), a_unit, a_step)
    # 1-z_(2i+1)*T = (B_i + (B_i-N)*T/2^height)/B_i.
    b_unit = _linear(ring, first_b, (first_b - n) // divisor)
    b_step = _linear(ring, step_b, step_b // divisor)
    numerator = _unit_product(ring, (0, b_count), b_unit, b_step)
    denominator = _unit_product(
        ring, (0, b_count), ring.constant(first_b), ring.constant(step_b))
    return ring.multiply(a_product, ring.multiply(numerator, ring.inverse(denominator)))


def normalized_parity_difference(n, first_a, step_a, first_b, step_b, height,
                                 order, precision):
    """Return (-1)^(k+1) k! f[z_0,...,z_k]/2^(k-1) modulo 2^p.

The nodes form the declared interleaved isometry sequence and f is parity.
This odd unit needs only p Mahler terms, even when k is enormous.
"""
    if order < 1 or precision < 1:
        raise ValueError("Require positive divided-difference order and precision")
    denominator = 1
    for j in range(1, precision):
        denominator *= order + j
    guard = valuation_two(denominator)
    working = precision + guard
    ring = QuadraticJets(working, precision - 1)
    inverse_nodes = _scalar_coefficients(ring.inverse(_node_product(
        ring, order + 1, n, first_a, step_a, first_b, step_b, height)))
    elementary = interval_elementaries(0, order, precision - 1)
    falling = [(-1) ** j * value % ring.modulus for j, value in enumerate(elementary)]
    denominator, answer = 1, 0
    modulus = 1 << precision
    for j in range(precision):
        if j:
            last = order + j - 1
            falling = [(falling[i] - (last * falling[i - 1] if i else 0)) % ring.modulus
                       for i in range(precision)]
            denominator *= order + j
        coefficient = sum(falling[i] * inverse_nodes[j - i] for i in range(j + 1)) % ring.modulus
        loss = valuation_two(denominator)
        if coefficient % (1 << loss):
            raise ArithmeticError("The isometry divided-difference integrality was lost")
        quotient = (coefficient >> loss) * pow(denominator >> loss, -1, modulus) % modulus
        answer = (answer + (-1) ** j * (quotient << j)) % modulus
    if answer % 2 != 1:
        raise ArithmeticError("The normalized parity divided difference is not a unit")
    return {"unit": answer, "precision": precision, "order": order,
            "working_precision": working, "mahler_terms": precision,
            "factorial_ratio_guard_bits": guard, "enumerated_nodes": 0}


def _root_shift(n, first_a, step_a, first_b, step_b, index, modulus):
    if index % 2 == 0:
        return (first_a + step_a * (index // 2) - 1) % modulus
    denominator = first_b + step_b * (index // 2)
    return (n * pow(denominator, -1, modulus) - 1) % modulus


def _domain_reversals(n, first_a, step_a, first_b, step_b, count, precision, degree):
    """Reverse P(1+Y), Q(1+Y), retaining their true leading coefficients."""
    ring = QuadraticJets(precision, degree)
    p = _unit_product(ring, (0, count),
                             _linear(ring, 1, 1 - first_a), _linear(ring, 0, -step_a))
    # Q(1+Y)=prod_i((N-B_i)-B_i*Y).
    q = _unit_product(ring, (0, count),
                             _linear(ring, -first_b, n - first_b),
                             _linear(ring, -step_b, -step_b))
    return _scalar_coefficients(p), _scalar_coefficients(q)


def balanced_bezout_modular(n, first_a, step_a, first_b, step_b, count, height, precision):
    """Complete optimal U,V mod 2^p, with UP+VQ=2^r over Z_(2)[X].

Output coefficient index j means the coefficient of (X-1)^(count-1-j).
All omitted coefficients are zero modulo 2^precision. The returned r is
binary and is never used as a literal working-precision or allocation size.
"""
    _validate(n, first_a, step_a, first_b, step_b, count, height, precision)
    top = 2 * count - 1
    factorial_valuation = top - top.bit_count()
    exponent = (height - 1) * top + 1 + factorial_valuation
    if exponent != height * top - (count - 1).bit_count():
        raise ArithmeticError("The two exact denominator formulas disagree")
    retained = []
    for offset in range(min(top, 2 * precision)):
        order = top - offset
        gap = exponent - ((height - 1) * order + 1 + order - order.bit_count())
        if gap >= precision:
            break
        retained.append((order, gap))
    if not retained or retained[0] != (top, 0):
        raise ArithmeticError("The leading normalized Newton coefficient must survive")
    if len(retained) == 2 * precision and top > 2 * precision:
        following = top - 2 * precision
        following_gap = exponent - ((height - 1) * following + 1
                                      + following - following.bit_count())
        if following_gap < precision:
            raise ArithmeticError("The terminal-window bound failed")
    modulus = 1 << precision
    scalar_ring = QuadraticJets(precision, 0)
    actual_valuation, factorial_unit_jet = affine_product(scalar_ring, 1, 1, top, None)
    if actual_valuation != factorial_valuation:
        raise ArithmeticError("The normalized factorial oracle disagrees with Legendre")
    factorial_unit = _scalar_coefficients(factorial_unit_jet)[0]
    root_degree = min(top, (precision - 1) // height)
    root_ring = QuadraticJets(precision, root_degree)
    prefix = _scalar_coefficients(_node_product(
        root_ring, top, n, first_a, step_a, first_b, step_b, height, root_shifts=True))
    f_degree = min(top, top - retained[-1][0] + root_degree)
    f_reverse = [0] * (f_degree + 1)
    certificates = []
    for order, gap in retained:
        difference = normalized_parity_difference(
            n, first_a, step_a, first_b, step_b, height, order, precision - gap)
        coefficient = ((-1) ** (order + 1) * difference["unit"]
                       * pow(factorial_unit, -1, modulus) << gap) % modulus
        offset = top - order
        for j, value in enumerate(prefix):
            if offset + j <= f_degree:
                f_reverse[offset + j] = (f_reverse[offset + j] + coefficient * value) % modulus
        certificates.append({**difference, "valuation_gap": gap,
                             "scaled_newton_coefficient": coefficient})
        if order > retained[-1][0]:
            # N_(k-1)=N_k/(Y-(x_(k-1)-1)); invert the reversed linear factor.
            shift = _root_shift(n, first_a, step_a, first_b, step_b, order - 1, modulus)
            smaller = [prefix[0]]
            for j in range(1, len(prefix)):
                smaller.append((prefix[j] + shift * smaller[-1]) % modulus)
            prefix = smaller
            factorial_unit = factorial_unit * pow(order >> valuation_two(order), -1, modulus) % modulus
    degree = count - 1
    u_depth = min(degree, f_degree + (precision - 1) // (height + 1))
    v_depth = min(degree, f_degree + (precision - 1) // height)
    p_reverse, q_reverse = _domain_reversals(
        n, first_a, step_a, first_b, step_b, count, precision, max(u_depth, v_depth))
    u_reverse = _multiply(f_reverse, _inverse(p_reverse, modulus, u_depth + 1), modulus, u_depth + 1)
    v_reverse = [(-value) % modulus for value in _multiply(
        f_reverse, _inverse(q_reverse, modulus, v_depth + 1), modulus, v_depth + 1)]
    if not u_reverse[0] % 2 or not v_reverse[0] % 2:
        raise ArithmeticError("The optimal normalized Bezout polynomials must be primitive")
    return {"parameters": {"n": n, "first_a": first_a, "step_a": step_a,
                            "first_b": first_b, "step_b": step_b, "count": count,
                            "height": height, "precision": precision},
            "bezout_exponent": exponent, "polynomial_degree": degree, "basis_center": 1,
            "u_reversed_coefficients": u_reverse, "v_reversed_coefficients": v_reverse,
            "u_terms": [[degree - j, value] for j, value in enumerate(u_reverse) if value],
            "v_terms": [[degree - j, value] for j, value in enumerate(v_reverse) if value],
            "normalized_interpolant_degree": top,
            "normalized_interpolant_reversed_coefficients": f_reverse,
            "retained_newton_coefficients": certificates,
            "maximum_working_precision": max(item["working_precision"] for item in certificates),
            "enumerated_progression_nodes": 0,
            "scope": "Complete modular optimal Bezout polynomials on the declared separated progressions; no unknown-root search."}
