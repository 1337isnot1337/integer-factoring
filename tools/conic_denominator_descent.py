"""Exact denominator cancellation in a known-square quadratic conic.

This parameterizes rational points and certifies the cancellation, not a
smooth-denominator search.  Outside a directly exposed gcd factor of N,
the residual denominator is a value of the original root norm form or
its product with the unique unhelpful ambiguous class.
"""

from math import gcd, isqrt

from early_stop_class_factoring import ambiguous_factorization
from quadratic_form_arithmetic import (
    compose_primitive, composition_data, extended_gcd, form_discriminant,
    inverse_primitive, reduce_primitive, reduce_with_transform,
)


def square_target_from_root(root_form, odd_n):
    """Prepare H=(a^2,2b,C), or expose a factor while making a odd and coprime."""
    if odd_n <= 1 or odd_n % 2 == 0 or form_discriminant(root_form) != -4*odd_n:
        raise ValueError("an odd N>1 and a root form of discriminant -4*N are required")
    reduced_root = reduce_primitive(root_form)
    a, middle, c = reduced_root
    if a % 2 == 0:
        a, middle, c = c, -middle, a
    divisor = gcd(a, odd_n)
    if divisor != 1:
        assert 1 < divisor < odd_n
        return {"status": "factored", "factor": divisor}
    assert a % 2 == 1 and gcd(a, 2*odd_n) == 1
    root = (a, middle, c)
    data = composition_data(root, root)
    assert data["content"] == 1
    target = data["primitive_form"]
    assert target[0] == a*a and target[1] % 2 == 0
    reduced_target, transform = reduce_with_transform(target)
    root_norm_form = (a, -target[1], a*target[2])
    assert reduce_primitive(root_norm_form) == inverse_primitive(reduced_root)
    return {"status": "prepared", "root": reduced_root,
            "root_with_odd_coprime_leading": root, "target": target,
            "reduced_target": reduced_target, "target_reduction_transform": transform,
            "initial_denominator": a, "root_norm_form": root_norm_form}


def recover_square_root(target, vector, denominator):
    """Recover a root class from a primitive unramified square representation."""
    discriminant = form_discriminant(target)
    x, y = vector
    a, b, c = target
    if denominator <= 0 or gcd(denominator, discriminant) != 1:
        raise ValueError("a positive denominator coprime to the discriminant is required")
    if a*x*x+b*x*y+c*y*y != denominator*denominator or gcd(x, y) != 1:
        raise ValueError("a primitive exact square representation is required")
    common, bezout_x, bezout_y = extended_gcd(x, y)
    assert common == 1
    complement_x, complement_y = -bezout_y, bezout_x
    middle = (2*a*x*complement_x+b*(x*complement_y+y*complement_x)
              +2*c*y*complement_y)
    tail = (middle*middle-discriminant)//(4*denominator*denominator)
    assert middle*middle-4*denominator*denominator*tail == discriminant
    root = reduce_primitive((denominator, middle, denominator*tail))
    assert compose_primitive(root, root) == reduce_primitive(target)
    return root


def parameter_point(target, r, s, known_root=None):
    """Evaluate a primitive conic parameter and its exact content certificate."""
    discriminant = form_discriminant(target)
    if discriminant % 4 or (-discriminant//4) % 2 == 0:
        raise ValueError("the discriminant must be -4*N with N odd")
    odd_n = -discriminant//4
    squared_a, twice_b, c = target
    a = isqrt(squared_a)
    if (a*a != squared_a or a % 2 == 0 or gcd(a, odd_n) != 1
            or twice_b % 2 or gcd(r, s) != 1):
        raise ValueError("a square leading coefficient, coprimality, and primitive parameters are required")
    b = twice_b//2
    w = a*r-b*s
    x = c*s*s-r*r
    y = 2*s*w
    z = a*r*r-2*b*r*s+a*c*s*s
    assert z > 0 and squared_a*x*x+twice_b*x*y+c*y*y == z*z
    content = gcd(gcd(x, y), z)
    linear_gcd = gcd(w, odd_n)
    parity_factor = 2 if z % 2 == 0 else 1
    assert content == parity_factor*linear_gcd
    vector, denominator = (x//content, y//content), z//content
    assert gcd(*vector) == 1
    answer = {"parameter": (r, s), "linear_value": w,
              "raw_point": (x, y, z), "content": content,
              "linear_gcd_with_N": linear_gcd, "parity_factor": parity_factor,
              "vector": vector, "denominator": denominator,
              "proper_factor": linear_gcd if 1 < linear_gcd < odd_n else None}
    if answer["proper_factor"] is not None:
        return answer

    norm_form = (a, -twice_b, a*c)
    norm_parameter = (r, s)
    original_norm_class = reduce_primitive(norm_form)
    if linear_gcd == odd_n:
        residue = b*pow(a, -1, odd_n) % odd_n
        assert (r-residue*s) % odd_n == 0
        norm_parameter = ((r-residue*s)//odd_n, s)
        norm_form = (a*odd_n, 2*a*residue-twice_b,
                     (a*residue*residue-twice_b*residue+a*c)//odd_n)
        assert reduce_primitive(norm_form) == original_norm_class
    if parity_factor == 2:
        leading, middle, tail = norm_form
        residue = tail % 2
        old_r, old_s = norm_parameter
        assert (old_r-residue*old_s) % 2 == 0
        norm_parameter = ((old_r-residue*old_s)//2, old_s)
        norm_form = (2*leading, 2*leading*residue+middle,
                     (leading*residue*residue+middle*residue+tail)//2)
        if odd_n % 4 == 3:
            assert all(coefficient % 2 == 0 for coefficient in norm_form)
            assert denominator % 2 == 0
        else:
            exceptional = (2, 2, (odd_n+1)//2)
            assert reduce_primitive(norm_form) == compose_primitive(original_norm_class, exceptional)
    nr, ns = norm_parameter
    leading, middle, tail = norm_form
    assert leading*nr*nr+middle*nr*ns+tail*ns*ns == denominator
    answer["residual_norm_form"] = norm_form
    answer["residual_norm_parameter"] = norm_parameter
    answer["residual_norm_form_primitive"] = gcd(gcd(leading, middle), tail) == 1
    if gcd(denominator, discriminant) == 1:
        recovered_root = recover_square_root(target, vector, denominator)
        answer["recovered_root"] = recovered_root
        if known_root is not None:
            assert compose_primitive(known_root, known_root) == reduce_primitive(target)
            difference = compose_primitive(recovered_root, inverse_primitive(known_root))
            classification = ambiguous_factorization(difference, odd_n)
            assert classification["proper_factor"] is None
            answer["root_difference"] = difference
    return answer
