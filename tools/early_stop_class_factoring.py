"""Factor from complete quadratic-class generators and a supplied annihilator.

For discriminant -4*N with N odd, at most one nonidentity ambiguous class
fails to give a proper coprime factorization.  The two-group insertion
algorithm can therefore stop before retaining a second cyclic generator.
This module does not construct class-group generators or an annihilator.
"""

from math import gcd

from arithmetic_two_group import order_log_two, power
from quadratic_form_arithmetic import (
    compose_primitive, form_discriminant, inverse_primitive, reduce_primitive,
)


def ambiguous_factorization(form, odd_n):
    """Read the coprime factors of odd_n from a reduced ambiguous class."""
    if odd_n <= 1 or odd_n % 2 == 0 or form_discriminant(form) != -4*odd_n:
        raise ValueError("an odd integer greater than one and discriminant -4*N are required")
    reduced = reduce_primitive(form)
    a, b, c = reduced
    if b == 0:
        factors, shape = (a, c), "zero_middle"
    elif b == a:
        assert a % 2 == 0
        factors, shape = (a//2, 2*c-a//2), "middle_equals_leading"
    elif a == c:
        assert b % 2 == 0
        factors, shape = (a-b//2, a+b//2), "equal_outer_coefficients"
    else:
        raise ValueError("the supplied class is not ambiguous")
    factors = tuple(sorted(factors))
    assert factors[0] > 0 and factors[0]*factors[1] == odd_n
    assert gcd(*factors) == 1
    proper = factors[0] if factors[0] > 1 else None
    if proper is None:
        identity = (1, 0, odd_n)
        exceptional = (2, 2, (odd_n+1)//2)
        assert reduced == identity or (odd_n % 4 == 1 and reduced == exceptional)
    return {"form": reduced, "shape": shape, "coprime_factors": factors,
            "proper_factor": proper}


def factor_from_annihilated_generators(odd_n, generators, annihilator, with_trace=False):
    """Return a factor, or report that no supplied two-part direction gave one.

    Completeness of the supplied class-group generators is a promise.  If
    they generate the entire group and N has two distinct prime factors,
    the result is guaranteed to contain a proper factor.  Annihilation is
    checked for every generator actually processed.  No factorization of
    N or of the odd part of the annihilator is used.
    """
    if odd_n <= 1 or odd_n % 2 == 0 or not annihilator:
        raise ValueError("an odd integer greater than one and a nonzero annihilator are required")
    annihilator = abs(annihilator)
    odd_part = annihilator
    while odd_part % 2 == 0:
        odd_part //= 2
    generators = tuple(generators)
    identity = (1, 0, odd_n)
    order_bound = (4*odd_n).bit_length()
    basis = None
    stats = {"generator_count": len(generators), "processed_generators": 0,
             "group_multiplications": 0, "group_inversions": 0,
             "top_checks": 0, "dependent_top_cancellations": 0,
             "order_swaps": 0, "maximum_retained_basis_rank": 0,
             "maximum_cyclic_order_log": 0, "order_bit_bound": order_bound}
    trace = []

    def multiply(left, right):
        stats["group_multiplications"] += 1
        return compose_primitive(left, right)

    def inverse(value):
        stats["group_inversions"] += 1
        return inverse_primitive(value)

    def result(factor_data=None):
        answer = {"odd_n": odd_n, "annihilator": annihilator,
                  "odd_annihilator_part": odd_part,
                  "status": "factored" if factor_data else "no_factor_in_supplied_subgroup",
                  "factor": factor_data["proper_factor"] if factor_data else None,
                  "factor_data": factor_data, "statistics": stats}
        if with_trace:
            answer["trace"] = trace
        return answer

    for index, generator in enumerate(generators):
        if form_discriminant(generator) != -4*odd_n:
            raise ValueError("all generators must have discriminant -4*N")
        generator = reduce_primitive(generator)
        stats["processed_generators"] += 1
        if power(generator, annihilator, identity, multiply) != identity:
            raise ValueError("the supplied integer does not annihilate a processed generator")
        candidate = power(generator, odd_part, identity, multiply)
        candidate_order = order_log_two(candidate, identity, multiply, order_bound)
        while candidate_order:
            stats["maximum_cyclic_order_log"] = max(
                stats["maximum_cyclic_order_log"], candidate_order)
            top = power(candidate, 1 << (candidate_order-1), identity, multiply)
            assert top != identity
            factor_data = ambiguous_factorization(top, odd_n)
            stats["top_checks"] += 1
            if with_trace:
                trace.append({"generator_index": index, "candidate": candidate,
                              "candidate_order_log": candidate_order, "top": top,
                              "proper_factor": factor_data["proper_factor"]})
            if factor_data["proper_factor"] is not None:
                return result(factor_data)
            # Every nonidentity top that reaches this point is the unique
            # exceptional form (2,2,(N+1)/2).  There is no subset-product search.
            assert odd_n % 4 == 1 and top == (2, 2, (odd_n+1)//2)
            if basis is None:
                basis = (candidate, candidate_order)
                stats["maximum_retained_basis_rank"] = 1
                break
            if basis[1] < candidate_order:
                basis, (candidate, candidate_order) = (candidate, candidate_order), basis
                stats["order_swaps"] += 1
            previous_order = candidate_order
            correction = power(basis[0], 1 << (basis[1]-candidate_order), identity, multiply)
            candidate = multiply(candidate, inverse(correction))
            candidate_order = order_log_two(candidate, identity, multiply, order_bound)
            stats["dependent_top_cancellations"] += 1
            assert candidate_order < previous_order
    return result()
