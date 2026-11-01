"""Primitive representation of a supplied factored norm by orientation matching.

Prime roots are supplied factor-base data. The norm is fixed by the input;
this routine does not choose a smooth norm or enumerate a class group.
"""

from bisect import bisect_left
from math import gcd, prod

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form, reduce_form
from quadratic_ideal_hnf import (
    NormRepresentation, compose_forms_hnf, reduce_with_transform,
    represented_norm_from_word,
)


def represent_factored_norm(form, prime_factors):
    """Return a primitive representation or certify absence for the fixed norm.

    Each input factor is (p, e, root), with distinct certified primes p,
    e>0, gcd(p, Delta)=1, 0<=root<2p, and root^2=Delta mod 4p.
    Thus all factors split, and the root has O(log p) bits.
    The primality certificates/prime sieve are external factor-base data.
    Prime powers are retained; the norm need not be squarefree.

    The deterministic meet-in-the-middle table uses sorting and binary
    search, not an expected-time hash-table bound. Return a dictionary
    with ``representation`` (a NormRepresentation or None), ``prime_terms``
    (a signed word or None), and exact finite work counters.
    """
    canonical, transform = reduce_with_transform(form)
    discriminant = form[1]*form[1]-4*form[0]*form[2]
    parity = discriminant % 2
    identity = reduce_form((1, parity, (parity*parity-discriminant)//4))
    factors = tuple(prime_factors)
    if len({prime for prime, _, _ in factors}) != len(factors):
        raise ValueError("combine repeated copies of each prime into one exponent")
    for prime, exponent, root in factors:
        if prime < 2 or exponent <= 0 or gcd(prime, discriminant) != 1:
            raise ValueError("expected a positive power of a prime coprime to Delta")
        if not 0 <= root < 2*prime:
            raise ValueError("supply the root as a representative in [0, 2p)")
        if (root*root-discriminant) % (4*prime):
            raise ValueError("invalid split-prime root")
    norm = prod(prime**exponent for prime, exponent, _ in factors)
    stats = {
        "distinct_primes": len(factors),
        "multiplicity": sum(exponent for _, exponent, _ in factors),
        "group_multiplications": 0,
        "left_products": 0,
        "right_products_examined": 0,
    }

    def multiply(left, right):
        stats["group_multiplications"] += 1
        return compose_forms_hnf(left, right)

    choices = []
    for prime, exponent, root in factors:
        prime_form = reduce_form((prime, root, (root*root-discriminant)//(4*prime)))
        plus = power(prime_form, exponent, identity, multiply)
        choices.append((plus, inverse_form(plus)))
    split = len(choices)//2

    def products(start, stop):
        # The explicit stack also works beyond Python's recursion limit.
        stack = [(start, identity, 0)]
        while stack:
            index, value, mask = stack.pop()
            if index == stop:
                yield value, mask
                continue
            plus, minus = choices[index]
            stack.append((index+1, multiply(value, minus), mask | (1 << index)))
            stack.append((index+1, multiply(value, plus), mask))

    left = sorted(products(0, split))
    left_keys = [value for value, _ in left]
    stats["left_products"] = len(left)
    for right, right_mask in products(split, len(choices)):
        stats["right_products_examined"] += 1
        needed = multiply(canonical, inverse_form(right))
        position = bisect_left(left_keys, needed)
        if position == len(left) or left[position][0] != needed:
            continue
        mask = left[position][1] | right_mask
        terms = tuple((prime, -exponent if mask >> index & 1 else exponent, root)
                      for index, (prime, exponent, root) in enumerate(factors))
        represented = represented_norm_from_word(discriminant, terms)
        assert represented.form == canonical
        assert represented.norm == norm and represented.removed_integer_content == 1
        x, y = represented.vector
        original_vector = (transform[0][0]*x+transform[0][1]*y,
                           transform[1][0]*x+transform[1][1]*y)
        x, y = original_vector
        assert gcd(x, y) == 1
        assert form[0]*x*x+form[1]*x*y+form[2]*y*y == norm
        result = NormRepresentation(tuple(form), original_vector, norm, norm, 1)
        return {"representation": result, "prime_terms": terms, "stats": stats}
    return {"representation": None, "prime_terms": None, "stats": stats}
