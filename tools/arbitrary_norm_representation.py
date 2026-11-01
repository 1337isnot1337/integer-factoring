"""Complete fixed-norm representation, allowing imprimitive coordinates.

Split prime powers use all allocations between conjugate ideals. Odd
inert primes have forced integer content. The norm is supplied and
factored; this module performs no norm selection or class-group listing.
"""

from bisect import bisect_left
from math import gcd, prod

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form, reduce_form
from quadratic_ideal_hnf import (
    NormRepresentation, compose_forms_hnf, reduce_with_transform,
    represented_norm_from_word,
)


def balanced_radix_partition(radices):
    """Partition radices with product ratio at most their largest value."""
    if any(not isinstance(radix, int) or radix < 2 for radix in radices):
        raise ValueError('Every radix must be an integer at least two.')
    buckets, sizes = [[], []], [1, 1]
    for index in sorted(range(len(radices)), key=lambda i: (-radices[i], i)):
        side = 0 if sizes[0] <= sizes[1] else 1
        buckets[side].append(index)
        sizes[side] *= radices[index]
    if sizes[0] > sizes[1]:
        buckets.reverse()
        sizes.reverse()
    assert sizes[1] <= max(radices, default=1)*sizes[0]
    return tuple(map(tuple, buckets)), tuple(sizes)


def represent_arbitrary_factored_norm(form, prime_factors, *, maximum_primitive_norm=None):
    """Find any representation of a supplied coprime factored norm.

    Each factor is (p,e,root), with distinct certified primes p, e>0,
    and gcd(p,Delta)=1. A split prime has a bounded root 0<=root<2p
    satisfying root^2=Delta mod4p. For an odd inert prime, root=None;
    the routine verifies the nonresidue condition. Inert primes of odd
    exponent certify nonrepresentation. The prime2 inert case is outside
    this interface; a split prime2 with a valid root is allowed.

    The result contains both the full represented norm and a primitive
    represented divisor. The latter supplies the signed prime word.
    The two norms differ by integer_content**2. Prime certification is
    supplied by the caller's sieve or other certified factor-base data.
    An optional maximum_primitive_norm requires the recovered primitive
    divisor to satisfy that bound. This search is complete across all
    possible contents, not only the first unbounded representation.
    """
    if (maximum_primitive_norm is not None
            and (not isinstance(maximum_primitive_norm, int) or maximum_primitive_norm < 1)):
        raise ValueError('The optional primitive norm bound must be a positive integer.')
    if (len(form) != 3 or any(not isinstance(entry, int) for entry in form)
            or form[0] <= 0 or gcd(gcd(*form[:2]), form[2]) != 1):
        raise ValueError('Expected a primitive positive integral form.')
    discriminant = form[1]*form[1]-4*form[0]*form[2]
    if discriminant >= 0:
        raise ValueError('Expected a negative form discriminant.')
    canonical, transform = reduce_with_transform(form)
    parity = discriminant % 2
    identity = reduce_form((1, parity, (parity*parity-discriminant)//4))
    factors = tuple(prime_factors)
    if len({prime for prime, _, _ in factors}) != len(factors):
        raise ValueError('Combine repeated copies of each prime into one exponent.')
    for prime, exponent, root in factors:
        if (not isinstance(prime, int) or not isinstance(exponent, int)
                or prime < 2 or exponent <= 0 or gcd(prime, discriminant) != 1):
            raise ValueError('Expected a positive power of a certified coprime prime.')
        if root is None:
            if prime == 2 or pow(discriminant % prime, (prime-1)//2, prime) != prime-1:
                raise ValueError('A missing root must designate an odd inert prime.')
        elif (not isinstance(root, int) or not 0 <= root < 2*prime
              or (root*root-discriminant) % (4*prime)):
            raise ValueError('Invalid bounded split-prime root.')
    norm = prod(prime**exponent for prime, exponent, _ in factors)
    split_factors = tuple(item for item in factors if item[2] is not None)
    inert_factors = tuple(item for item in factors if item[2] is None)
    radices = tuple(exponent+1 for _, exponent, _ in split_factors)
    partition, capacities = balanced_radix_partition(radices)
    stats = {
        'distinct_primes': len(factors), 'split_primes': len(split_factors),
        'inert_primes': len(inert_factors),
        'multiplicity': sum(exponent for _, exponent, _ in factors),
        'choice_product': prod(radices), 'maximum_radix': max(radices, default=1),
        'left_capacity': capacities[0], 'right_capacity': capacities[1],
        'group_multiplications': 0, 'left_products': 0,
        'right_products_examined': 0,
        'maximum_primitive_norm': maximum_primitive_norm,
        'over_bound_matches': 0,
        'odd_inert_exponent': any(exponent % 2 for _, exponent, _ in inert_factors),
    }

    def absent():
        return {'representation': None, 'primitive_representation': None,
                'integer_content': None, 'prime_terms': None, 'stats': stats}

    if stats['odd_inert_exponent']:
        return absent()
    inert_content = prod(prime**(exponent//2) for prime, exponent, _ in inert_factors)

    def multiply(left, right):
        stats['group_multiplications'] += 1
        return compose_forms_hnf(left, right)

    choices = []
    for prime, exponent, root in split_factors:
        prime_form = reduce_form((prime, root, (root*root-discriminant)//(4*prime)))
        first = inverse_form(power(prime_form, exponent, identity, multiply))
        step = multiply(prime_form, prime_form)
        values = [first]
        for _ in range(exponent):
            values.append(multiply(values[-1], step))
        choices.append(tuple(values))

    def products(indices):
        strides, stride = [], 1
        for index in indices:
            strides.append(stride)
            stride *= radices[index]
        stack = [(0, identity, 0)]
        while stack:
            position, value, code = stack.pop()
            if position == len(indices):
                yield value, code
                continue
            index = indices[position]
            for allocation in reversed(range(radices[index])):
                stack.append((position+1, multiply(value, choices[index][allocation]),
                              code+allocation*strides[position]))

    left_indices, right_indices = partition
    if maximum_primitive_norm is None:
        left = sorted(products(left_indices))
        left_code_index = 1
    else:
        norm_choices = tuple(tuple(prime**abs(2*allocation-exponent)
                                   for allocation in range(exponent+1))
                             for prime, exponent, _ in split_factors)

        def partial_norm(indices, code):
            value = 1
            for index in indices:
                value *= norm_choices[index][code % radices[index]]
                code //= radices[index]
            assert code == 0
            return value

        # For each class, binary search returns its least primitive
        # partial norm. Larger left norms cannot make a bound feasible
        # for a right allocation when this least one already fails.
        left = sorted((value, partial_norm(left_indices, code), code)
                      for value, code in products(left_indices))
        left_code_index = 2
    left_keys = [entry[0] for entry in left]
    stats['left_products'] = len(left)
    assert len(left) == capacities[0]
    for right, right_code in products(right_indices):
        stats['right_products_examined'] += 1
        needed = multiply(canonical, inverse_form(right))
        position = bisect_left(left_keys, needed)
        if position == len(left) or left[position][0] != needed:
            continue
        if maximum_primitive_norm is not None:
            smallest_norm = left[position][1]*partial_norm(right_indices, right_code)
            if smallest_norm > maximum_primitive_norm:
                stats['over_bound_matches'] += 1
                continue
        allocations = [0]*len(split_factors)
        for indices, code in ((left_indices, left[position][left_code_index]), (right_indices, right_code)):
            for index in indices:
                allocations[index] = code % radices[index]
                code //= radices[index]
            assert code == 0
        terms, content = [], inert_content
        for allocation, (prime, exponent, root) in zip(allocations, split_factors):
            signed_exponent = 2*allocation-exponent
            content *= prime**min(allocation, exponent-allocation)
            if signed_exponent:
                terms.append((prime, signed_exponent, root))
        terms = tuple(terms)
        represented = represented_norm_from_word(discriminant, terms)
        assert represented.form == canonical and represented.removed_integer_content == 1
        assert represented.norm*content*content == norm
        assert maximum_primitive_norm is None or represented.norm <= maximum_primitive_norm
        x, y = represented.vector
        primitive_vector = (transform[0][0]*x+transform[0][1]*y,
                            transform[1][0]*x+transform[1][1]*y)
        full_vector = tuple(content*entry for entry in primitive_vector)
        assert gcd(*primitive_vector) == 1 and gcd(*full_vector) == content
        a, middle, c = form
        for vector, expected in ((primitive_vector, represented.norm), (full_vector, norm)):
            x, y = vector
            assert a*x*x+middle*x*y+c*y*y == expected
        primitive = NormRepresentation(tuple(form), primitive_vector,
                                       represented.norm, represented.norm, 1)
        full = NormRepresentation(tuple(form), full_vector, norm, norm, 1)
        return {'representation': full, 'primitive_representation': primitive,
                'integer_content': content, 'prime_terms': terms, 'stats': stats}
    return absent()
