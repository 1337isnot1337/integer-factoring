"""Polynomial early extraction when at most one nonidentity top is useless.

The group algorithm uses only the supplied operations and annihilator.
The quadratic-form wrapper applies it at discriminant -4N. See
research/polynomial-ambiguous-extraction.md for the proof and scope.
"""

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from quadratic_ideal_hnf import compose_forms_hnf


def first_useful_top(generators, annihilator, identity, multiply, inverse,
                     order_bit_bound, inspect_top):
    """Find a useful order-two element, retaining at most one cyclic direction.

    ``inspect_top(top)`` returns a non-None certificate on useful elements.
    Among nonidentity order-two elements, it may reject at most one element.
    The annihilator must kill every generator. The supplied bound must bound
    log2 of the order of each element in their Sylow two-subgroup.

    Return ``(top, certificate, stats)`` or ``(None, None, stats)``. The latter
    means that the generated subgroup has no useful order-two element. Full
    generation is an application hypothesis, not tested by this routine.
    """
    annihilator = abs(annihilator)
    if annihilator == 0:
        raise ValueError("the annihilator must be nonzero")
    if order_bit_bound < 0:
        raise ValueError("the order bound must be nonnegative")
    two_valuation = (annihilator & -annihilator).bit_length()-1
    odd_part = annihilator >> two_valuation
    stats = {
        "generators_processed": 0,
        "group_multiplications": 0,
        "group_inversions": 0,
        "tested_tops": 0,
        "cancellations": 0,
        "order_swaps": 0,
        "peak_retained_cyclic_directions": 0,
    }

    def mul(left, right):
        stats["group_multiplications"] += 1
        return multiply(left, right)

    def inv(value):
        stats["group_inversions"] += 1
        return inverse(value)

    def order_and_top(value):
        if value == identity:
            return 0, identity
        current = value
        for order in range(1, order_bit_bound+1):
            top = current
            current = mul(current, current)
            if current == identity:
                return order, top
        raise ValueError("element order exceeds the supplied two-power bound")

    retained = None
    retained_order = 0
    retained_top = None
    for generator in generators:
        stats["generators_processed"] += 1
        candidate = power(generator, odd_part, identity, mul)
        candidate_order, top = order_and_top(candidate)
        while candidate_order:
            stats["tested_tops"] += 1
            certificate = inspect_top(top)
            if certificate is not None:
                return top, certificate, stats
            if retained is None:
                retained, retained_order, retained_top = candidate, candidate_order, top
                stats["peak_retained_cyclic_directions"] = 1
                break
            if top != retained_top:
                raise ValueError("the top predicate rejected two distinct nonidentity tops")
            if candidate_order > retained_order:
                candidate, retained = retained, candidate
                candidate_order, retained_order = retained_order, candidate_order
                stats["order_swaps"] += 1
            previous_order = candidate_order
            correction = retained
            for _ in range(retained_order-candidate_order):
                correction = mul(correction, correction)
            candidate = mul(candidate, inv(correction))
            stats["cancellations"] += 1
            candidate_order, top = order_and_top(candidate)
            assert candidate_order < previous_order
    return None, None, stats


def factor_from_annihilator(integer, generators, annihilator):
    """Split odd N from generators of C_{-4N} and any supplied annihilator.

    A proper factor is guaranteed when N has two distinct prime factors and
    the generators span the full class group. Inputs are checked for the
    stated discriminant; generation and annihilation are supplied promises.
    No factorization of N, the discriminant, or the annihilator is used.
    """
    if integer <= 1 or integer % 2 == 0:
        raise ValueError("the target must be odd and greater than one")
    discriminant = -4*integer
    canonical = []
    for form in generators:
        a, middle, c = form
        if middle*middle-4*a*c != discriminant:
            raise ValueError("a generator has the wrong discriminant")
        canonical.append(reduce_form(form))
    identity = (1, 0, integer)
    # There are at most (floor(sqrt(D))+1)^2 reduced coefficient pairs.
    # Thus |C_Delta| < 4D, and bits(D)+2 bounds all two-power order logs.
    order_bit_bound = (-discriminant).bit_length()+2
    top, factor, stats = first_useful_top(
        canonical, annihilator, identity, compose_forms_hnf, inverse_form,
        order_bit_bound, lambda form: ambiguous_factor(form, integer),
    )
    return {"factor": factor, "ambiguous_form": top, "stats": stats}
