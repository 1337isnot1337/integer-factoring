"""Recover a fixed smooth-square witness from variable-cap decision queries.

The decision callback is an explicit unresolved oracle.  Every query keeps
a square norm and a squared target class.  No efficient decider is claimed.
"""

from math import gcd

from quadratic_form_arithmetic import (
    compose_primitive, form_discriminant, power_primitive, reduce_primitive,
)
from universal_smooth_square import representation_from_square_word


def recover_uniform_square_word(form, primes, roots, exponent_cap, decide, with_trace=False):
    """Recover Q(x,y)=M^2 from a promised positive fixed-square instance.

    Fixed primes and roots define the canonical prime forms g_p.  The callback
    decide(target, caps) must exactly decide whether target represents
    (product p**caps[p])**2, including nonprimitive representations.  It must
    handle unequal caps, zero caps, and negative instances.  The initial
    instance is promised positive; it is not tested here.
    """
    discriminant = form_discriminant(form)
    primes = tuple(primes)
    if (exponent_cap < 1 or len(set(primes)) != len(primes)
            or any(prime < 2 or prime not in roots or gcd(prime, discriminant) != 1
                   for prime in primes)):
        raise ValueError("distinct unramified primes, their roots, and a positive cap are required")
    generators = []
    for prime in primes:
        root = roots[prime]
        if (root*root-discriminant) % (4*prime):
            raise ValueError("an invalid splitting root was supplied")
        generators.append(reduce_primitive((prime, root,
                                            (root*root-discriminant)//(4*prime))))
    squares = [power_primitive(generator, 2) for generator in generators]
    parity = discriminant % 2
    identity = (1, parity, (parity-discriminant)//4)
    target = reduce_primitive(form)
    caps = [exponent_cap]*len(primes)
    offsets = [0]*len(primes)
    trace, query_count = [], 0

    def query(candidate, candidate_caps, index, kind):
        nonlocal query_count
        answer = decide(candidate, tuple(candidate_caps))
        if not isinstance(answer, bool):
            raise ValueError("the representation decider must return an exact boolean")
        query_count += 1
        if with_trace:
            trace.append({"coordinate": index, "kind": kind,
                          "target": candidate, "caps": tuple(candidate_caps),
                          "answer": answer})
        return answer

    for index, square in enumerate(squares):
        while caps[index] > 1:
            previous_radius = caps[index]
            shift, radius = previous_radius//2, (previous_radius+1)//2
            candidate_caps = caps[:]
            candidate_caps[index] = radius
            positive_shift = power_primitive(square, shift)
            left_target = compose_primitive(target, positive_shift)
            if query(left_target, candidate_caps, index, "left_interval"):
                target = left_target
                offsets[index] -= shift
            else:
                target = compose_primitive(target, power_primitive(square, -shift))
                offsets[index] += shift
            caps = candidate_caps
            assert radius < previous_radius
            assert abs(offsets[index])+radius <= exponent_cap
        if caps[index] == 1:
            candidate_caps = caps[:]
            candidate_caps[index] = 0
            left_target = compose_primitive(target, square)
            if query(left_target, candidate_caps, index, "left_singleton"):
                target = left_target
                offsets[index] -= 1
            elif query(target, candidate_caps, index, "middle_singleton"):
                pass
            else:
                target = compose_primitive(target, power_primitive(square, -1))
                offsets[index] += 1
            caps = candidate_caps
        assert abs(offsets[index]) <= exponent_cap
    if target != identity:
        raise ValueError("decision responses are inconsistent with the promised representation")
    value = identity
    for square, exponent in zip(squares, offsets):
        value = compose_primitive(value, power_primitive(square, exponent))
    assert value == reduce_primitive(form)
    witness = representation_from_square_word(form, primes, roots, offsets, exponent_cap)
    query_bound = len(primes)*((exponent_cap-1).bit_length()+2)
    assert query_count <= query_bound
    result = {"half_word": offsets, "query_count": query_count,
              "query_bound": query_bound, "representation": witness}
    if with_trace:
        result["queries"] = trace
    return result
