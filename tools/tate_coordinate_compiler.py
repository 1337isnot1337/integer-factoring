"""Algebraic q-torsion coordinates in v, before the composition v=Delta(t).

The compiler gives polynomial-size rational diagonal data for coefficients
in v. It does NOT compute a high-index coefficient after Delta substitution.
"""

from functools import lru_cache
from math import isqrt

from cyclotomic_character_decoder import (
    modular_inverse, modular_product, polynomial_product, polynomial_remainder,
)
from lambert_theta_sections import AdditivePhaseAlgebra


def binary_add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for degree, value in polynomial.items():
            coefficient = result.get(degree, 0) ^ value
            if coefficient:
                result[degree] = coefficient
            else:
                result.pop(degree, None)
    return result


def binary_product(left, right):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            value = result.get(i+j, 0) ^ polynomial_product(a, b)
            if value:
                result[i+j] = value
            else:
                result.pop(i+j, None)
    return result


def binary_square(polynomial):
    return {2*j: sum(1 << (2*i) for i in range(value.bit_length()) if value >> i & 1)
            for j, value in polynomial.items()}


def division_polynomial(index):
    """psi_index on y^2+xy=x^3+v, represented by {v_degree: x_bitmask}."""
    if not isinstance(index, int) or index < 0:
        raise ValueError("a nonnegative division index is required")

    @lru_cache(None)
    def psi(n):
        if n == 0:
            return {}
        if n == 1:
            return {0: 1}
        if n == 2:
            return {0: 2}
        if n == 3:
            return {0: (1 << 4) ^ (1 << 3), 1: 1}
        if n == 4:
            return {0: 1 << 6, 1: 1 << 2}
        m = n//2
        if n & 1:
            return binary_add(
                binary_product(psi(m+2), binary_product(psi(m), binary_square(psi(m)))),
                binary_product(psi(m-1), binary_product(psi(m+1), binary_square(psi(m+1)))))
        numerator = binary_product(psi(m), binary_add(
            binary_product(psi(m+2), binary_square(psi(m-1))),
            binary_product(psi(m-2), binary_square(psi(m+1)))))
        if any(value & 1 for value in numerator.values()):
            raise ValueError("inexact division by psi_2=x")
        return {j: value >> 1 for j, value in numerator.items()}

    return psi(index)


def ring_polynomial_product(left, right, algebra):
    result = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            position = i+k, j+ell
            value = result.get(position, 0) ^ algebra.product(a, b)
            if value:
                result[position] = value
            else:
                result.pop(position, None)
    return result


def shifted_equation(psi, constant, algebra):
    """F(v,w)=psi(constant+w,v), using the binary Pascal coefficients."""
    degree = max(value.bit_length()-1 for value in psi.values())
    powers = [1]
    for _ in range(degree):
        powers.append(algebra.product(powers[-1], constant))
    result = {}
    for v_degree, mask in psi.items():
        for x_degree in range(mask.bit_length()):
            if not mask >> x_degree & 1:
                continue
            w_degree = x_degree
            while True:
                position = v_degree, w_degree
                value = result.get(position, 0) ^ powers[x_degree-w_degree]
                if value:
                    result[position] = value
                else:
                    result.pop(position, None)
                if not w_degree:
                    break
                w_degree = (w_degree-1) & x_degree
    return result


class PolynomialCoefficientAlgebra:
    def __init__(self, modulus):
        self.modulus = modulus

    def product(self, left, right):
        return modular_product(left, right, self.modulus)


def plan_algebra(plan):
    return PolynomialCoefficientAlgebra(plan['coefficient_modulus'])


def compile_coordinate(prime, intrinsic=False):
    phase_algebra = AdditivePhaseAlgebra(prime)
    psi = division_polynomial(prime)
    nodal_shift = prime*(prime-1)//2
    nodal_value = psi.get(0, 0)
    if nodal_value & ((1 << nodal_shift)-1):
        raise ValueError("unexpected nodal multiplicity")
    nodal_modulus = nodal_value >> nodal_shift
    if nodal_modulus.bit_length()-1 != (prime-1)//2:
        raise ValueError("unexpected nonsingular nodal degree")
    if intrinsic:
        algebra = PolynomialCoefficientAlgebra(nodal_modulus)
        constant = polynomial_remainder(2, nodal_modulus)
    else:
        algebra = phase_algebra
        inverse = algebra.inverse_theta_constant(1)
        constant = algebra.product(algebra.phase(1), algebra.product(inverse, inverse))
    equation = shifted_equation(psi, constant, algebra)
    if equation.get((0, 0), 0):
        raise ValueError("the selected nodal constant is not a root")
    derivative = equation.get((0, 1), 0)
    derivative_inverse = modular_inverse(derivative, algebra.modulus)
    numerator, denominator = {}, {}
    for (i, j), value in equation.items():
        if i+j < 1:
            raise ValueError("invalid shifted equation constant")
        denominator[i, i+j-1] = value
        if j & 1:
            numerator[i, i+j] = value
    assert denominator.get((0, 0), 0) == derivative
    degree_v = max(i for i, _ in numerator.keys() | denominator.keys())
    degree_w = max(j for _, j in numerator.keys() | denominator.keys())
    return {"prime": prime, "intrinsic": intrinsic, "coefficient_modulus": algebra.modulus,
            "nodal_modulus": nodal_modulus,
            "constant": constant, "psi": psi, "equation": equation,
            "derivative": derivative, "derivative_inverse": derivative_inverse,
            "numerator": numerator, "denominator": denominator,
            "degree_v": degree_v, "degree_w": degree_w,
            "linear_space_bound": (degree_v+1)*(degree_w+1)}


def embed_coefficient(value, plan):
    """Embed the intrinsic nodal algebra into the additive phase algebra."""
    if not plan['intrinsic']:
        return value
    algebra = AdditivePhaseAlgebra(plan['prime'])
    inverse = algebra.inverse_theta_constant(1)
    constant = algebra.product(algebra.phase(1), algebra.product(inverse, inverse))
    result = 0
    for i in range(value.bit_length()-1, -1, -1):
        result = algebra.product(result, constant) ^ (value >> i & 1)
    return result


def algebraic_coefficient(index, plan):
    """Coefficient of x_q(v); poly(prime)*log(index+1) binary diagonal steps."""
    if not isinstance(index, int) or index < 0:
        raise ValueError("a nonnegative coefficient index is required")
    if index == 0:
        return {"value": plan['constant'], "steps": [], "ring_products": 0}
    algebra = plan_algebra(plan)
    numerator, denominator = plan['numerator'], plan['denominator']
    target, records = index, []
    while target:
        bit = target & 1
        child = {}
        products = 0
        for (i, j), a in numerator.items():
            for (k, ell), b in denominator.items():
                if (i+k) & 1 != bit or (j+ell) & 1 != bit:
                    continue
                products += 1
                position = (i+k)//2, (j+ell)//2
                value = child.get(position, 0) ^ algebra.product(a, b)
                if value:
                    child[position] = value
                else:
                    child.pop(position, None)
        records.append({"target_bits": target.bit_length(), "bit": bit,
                        "numerator_terms": len(numerator), "child_terms": len(child),
                        "denominator_terms": len(denominator), "ring_products": products})
        numerator = child
        denominator = {position: algebra.product(value, value)
                       for position, value in denominator.items()}
        target //= 2
    inverse = modular_inverse(denominator[0, 0], algebra.modulus)
    value = algebra.product(numerator.get((0, 0), 0), inverse)
    return {"value": value, "steps": records,
            "ring_products": sum(record['ring_products'] for record in records)}


def series_product(left, right, limit, algebra):
    result = [0]*(limit+1)
    for i, a in enumerate(left):
        if not a or i > limit:
            continue
        for j, b in enumerate(right[:limit-i+1]):
            if b:
                result[i+j] ^= algebra.product(a, b)
    return result


def evaluate_psi(psi, series, limit, algebra, derivative=False):
    degree = max(value.bit_length()-1 for value in psi.values())
    powers = [[1]+[0]*limit]
    for _ in range(degree):
        powers.append(series_product(powers[-1], series, limit, algebra))
    result = [0]*(limit+1)
    for v_degree, mask in psi.items():
        if v_degree > limit:
            continue
        for x_degree in range(mask.bit_length()):
            if not mask >> x_degree & 1 or derivative and not x_degree & 1:
                continue
            power = powers[x_degree-int(derivative)]
            for k, value in enumerate(power[:limit-v_degree+1]):
                result[k+v_degree] ^= value
    return result


def hensel_coefficients(plan, limit):
    """Independent Newton/Hensel prefix, for controls of the diagonal compiler."""
    algebra = plan_algebra(plan)
    result = [plan['constant']]
    precision = 1
    while precision <= limit:
        next_precision = min(2*precision, limit+1)
        degree = next_precision-1
        residual = evaluate_psi(plan['psi'], result, degree, algebra)
        derivative = evaluate_psi(plan['psi'], result, degree, algebra, derivative=True)
        inverse = [modular_inverse(derivative[0], algebra.modulus)]+[0]*degree
        for n in range(1, next_precision):
            value = 0
            for k in range(1, n+1):
                value ^= algebra.product(derivative[k], inverse[n-k])
            inverse[n] = algebra.product(inverse[0], value)
        correction = series_product(residual, inverse, degree, algebra)
        assert not any(correction[:precision])
        result = [result[i] if i < len(result) else 0 for i in range(next_precision)]
        result = [a ^ b for a, b in zip(result, correction)]
        precision = next_precision
    return result[:limit+1]


def delta_composition(coefficients, limit):
    """Dense finite control of x_q(Delta(t)); this retains the original gap."""
    delta = sum(1 << (r*r) for r in range(1, isqrt(limit)+1, 2))
    cutoff = (1 << (limit+1))-1
    power, result = 1, [0]*(limit+1)
    for coefficient in coefficients[:limit+1]:
        if coefficient:
            mask = power
            while mask:
                bit = mask & -mask
                result[bit.bit_length()-1] ^= coefficient
                mask ^= bit
        power = polynomial_product(power, delta) & cutoff
    return result


def q3_algebraic_coefficient(index):
    """Independent Fuss-Catalan parity in v, not the original nome t."""
    if not isinstance(index, int) or index < 0:
        raise ValueError("a nonnegative coefficient index is required")
    if index == 0:
        return 1
    binomial_valuation = ((index-1).bit_count()+(3*index-1).bit_count()
                          -(4*index-2).bit_count())
    divisor_valuation = (index & -index).bit_length()-1
    assert binomial_valuation >= divisor_valuation
    return int(binomial_valuation == divisor_valuation)


def nodal_trace(value, plan):
    """Algebra trace over F_2, without factoring the nodal modulus."""
    if not plan['intrinsic']:
        raise ValueError("the intrinsic coefficient algebra is required")
    algebra = plan_algebra(plan)
    degree = algebra.modulus.bit_length()-1
    result = 0
    for i in range(degree):
        result ^= algebra.product(value, 1 << i) >> i & 1
    return result


def trace_one_element(plan):
    """A fixed coefficient eta with Tr(eta)=1, including even algebra degrees."""
    degree = plan['coefficient_modulus'].bit_length()-1
    for i in range(degree):
        if nodal_trace(1 << i, plan):
            return 1 << i
    raise ValueError("the supposedly etale trace pairing is degenerate")


def trace_tate_coefficient(integer, prime):
    """Cheap original-nome trace coefficient; the full phase value is still open."""
    if not isinstance(integer, int) or integer < 1:
        raise ValueError("a positive original index is required")
    AdditivePhaseAlgebra(prime)
    while integer % 2 == 0:
        integer //= 2
    while integer % prime == 0:
        integer //= prime
    return int(isqrt(integer)**2 == integer)
