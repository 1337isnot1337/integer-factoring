"""Binary trace characters through classical Artin--Hasse splitting.

Reuse the existing exact Artin--Hasse acquisition and canonical ring carry /
scalar inverse compiler. Fourier projection spends r precision bits, not a
divided-power coefficient algebra. No native integer digit circuit is handled.
"""
from functools import lru_cache

from precision_rational_carries import rational_sum, factor_degree
from teichmuller_cyclic_carries import UnramifiedRing
from witt_artin_hasse import artin_hasse_coefficients


def residue_bits(ring, value):
    return sum((c & 1) << i for i, c in enumerate(ring.coefficients(value)))


@lru_cache(None)
def splitting_coefficients(precision):
    """Acquire E_2(pi X) modulo 2^P; pi is the even root E_2(pi)=-1.

    Terms of degree >=P vanish because v_2(pi)=1. The imported recurrence
    computes rational coefficients with odd denominators, not a field table.
    """
    if precision < 1:
        raise ValueError('Positive splitting precision required')
    modulus = 1 << precision
    rational = artin_hasse_coefficients(precision-1)
    coefficients = tuple(c.numerator*pow(c.denominator, -1, modulus) % modulus
                         for c in rational)

    def evaluate(value):
        answer = 0
        for c in reversed(coefficients):
            answer = (answer*value+c) % modulus
        return answer

    pi = 0
    for k in range(1, precision):
        if (evaluate(pi)+1) % (1 << (k+1)):
            pi += 1 << k
    assert (evaluate(pi)+1) % modulus == 0
    theta = tuple(c*pow(pi, n, modulus) % modulus
                  for n, c in enumerate(coefficients))
    return pi, theta


def acquire_character_factors(ring, width, row_modulus, variables, polynomial):
    """Split a polynomial's residue trace character MONOMIAL BY MONOMIAL.

    A polynomial sum need not be Teichmuller. Applying theta to that sum is
    invalid. Each residue coefficient is freshly Teichmuller-lifted, so every
    monomial value to which the splitting lemma is applied is Teichmuller.
    """
    t = row_modulus
    if width < 1 or ring.degree % width or t < 1 or width % t or variables < 1:
        raise ValueError('Incompatible character dimensions')
    columns = width//t
    pi, theta = splitting_coefficients(ring.precision)
    factors = []
    degree_total = 0
    monomials = 0
    for coefficient, exponents in polynomial:
        if len(exponents) != variables or any(e < 0 for e in exponents):
            raise ValueError('Ordinary polynomial exponents required')
        if ring.mod_two(ring.frobenius(coefficient, width)) != ring.mod_two(coefficient):
            raise ValueError('Coefficient residue outside the summation subfield')
        residue = residue_bits(ring, coefficient)
        if not residue:
            continue
        coefficient = ring.teichmuller(residue)
        monomials += 1
        degree_total += sum(exponents)
        for shift in range(width):
            row = shift % t
            column = shift//t
            base = ring.frobenius(coefficient, shift)
            power = 1
            terms = []
            for n, scalar in enumerate(theta):
                if scalar:
                    counts = [0]*(variables*columns)
                    for a, e in enumerate(exponents):
                        counts[a*columns+column] = n*e
                    terms.append((ring.mul(ring.integer(scalar), power), tuple(counts)))
                power = ring.mul(power, base)
            factors.append((row, terms))
    return factors, dict(monomials=monomials, monomial_degree_total=degree_total,
                         factors=len(factors), splitting_degree=ring.precision-1,
                         pi=pi)


def trace_character_sum(ring, width, row_modulus, variables, predicates,
                        numerators=(), denominators=(), nonzero=False):
    """Weighted trace-constrained sum modulo 2^(P-r), given precision P.

    Predicates are (ordinary polynomial, target bit). Their written length
    r is charged, without an inferred rank reduction for nonlinear functions.
    Weights are scalar synchronized polynomials / total unit inverses over
    the supplied precision-P ring. Return a value packed in the precision
    P-r ring with the SAME defining polynomial, and a resource invoice.
    """
    rank = len(predicates)
    precision = ring.precision
    output_precision = precision-rank
    if output_precision < 1:
        raise ValueError('Fourier projection requires r extra precision bits')
    groups = []
    acquisition = []
    targets = []
    for polynomial, bit in predicates:
        if bit not in (0, 1):
            raise ValueError('Binary trace target required')
        factors, metadata = acquire_character_factors(
            ring, width, row_modulus, variables, polynomial)
        groups.append(factors)
        acquisition.append(metadata)
        targets.append(bit)
    numerator = 0
    backend_queries = 0
    actual_degree = 0
    for mask in range(1 << rank):
        factors = list(numerators)
        sign = 0
        for i, group in enumerate(groups):
            if (mask >> i) & 1:
                factors.extend(group)
                sign ^= targets[i]
        value, backend = rational_sum(ring, width, row_modulus, variables,
                                      factors, denominators, nonzero)
        numerator = ring.sub(numerator, value) if sign else ring.add(numerator, value)
        backend_queries += backend['queries']
        actual_degree = max(actual_degree, backend['max_degree_budget'])
    divisor = 1 << rank
    coefficients = ring.coefficients(numerator)
    assert all(c % divisor == 0 for c in coefficients)
    output_ring = UnramifiedRing(output_precision, ring.polynomial)
    answer = output_ring.from_polynomial([c//divisor for c in coefficients])
    degree_total = sum(a['monomial_degree_total'] for a in acquisition)
    numerator_degree = sum(map(factor_degree, numerators))
    denominator_degree = sum(map(factor_degree, denominators))
    budget = (numerator_degree+((2*precision+1)*ring.degree-1)*denominator_degree
              +(precision-1)*width*degree_total)
    assert actual_degree <= budget
    invoice = dict(input_precision=precision, output_precision=output_precision,
                   predicates=rank, fourier_queries=1 << rank,
                   backend_queries=backend_queries, precision_loss=rank,
                   division_coefficientwise_exact=True,
                   coefficient_degree=ring.degree, summation_width=width,
                   relative_columns=width//row_modulus,
                   max_degree_budget=budget, actual_max_degree=actual_degree,
                   acquisition=acquisition)
    if row_modulus == width:
        # A complete Frobenius orbit spreads each monomial over all rows.
        # Both the row degree and rotated initial carry obey this sharper cap.
        invoice['full_row_degree_and_carry_bound'] = (
            numerator_degree+(2*precision+1)*(ring.degree//width)*denominator_degree
            +(precision-1)*degree_total)
    return answer, invoice
