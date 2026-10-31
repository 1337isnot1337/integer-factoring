"""Polynomial integer arithmetic for primitive positive binary forms.

Composition multiplies the two ideal lattices and uses extended gcd; it does
not search for a coprime represented leading coefficient.  The original
binary_quadratic_forms.py remains an independent small-input comparison.
"""

from math import gcd


def extended_gcd(left, right):
    """Return nonnegative g and x,y with x*left+y*right=g."""
    a, b = left, right
    x, y, u, v = 1, 0, 0, 1
    while b:
        quotient = a // b
        a, b, x, u, y, v = b, a-quotient*b, u, x-quotient*u, v, y-quotient*v
    if a < 0:
        return -a, -x, -y
    return a, x, y


def form_discriminant(form):
    a, b, c = form
    discriminant = b*b-4*a*c
    if a <= 0 or discriminant >= 0 or gcd(gcd(a, b), c) != 1:
        raise ValueError("a primitive positive definite form is required")
    return discriminant


def reduce_with_transform(form):
    """Return (reduced form, U), where reduced(z)=form(U*z), det(U)=1."""
    discriminant = form_discriminant(form)
    a, b, c = form
    u00, u01, u10, u11 = 1, 0, 0, 1
    while True:
        shift = (a-b)//(2*a)
        c, b = c+b*shift+a*shift*shift, b+2*a*shift
        u01, u11 = u00*shift+u01, u10*shift+u11
        if a > c:
            a, b, c = c, -b, a
            u00, u01, u10, u11 = u01, -u00, u11, -u10
            continue
        if a == c and b < 0:
            b = -b
            u00, u01, u10, u11 = u01, -u00, u11, -u10
        assert b*b-4*a*c == discriminant and -a < b <= a <= c
        assert u00*u11-u01*u10 == 1
        return (a, b, c), ((u00, u01), (u10, u11))


def reduce_primitive(form):
    return reduce_with_transform(form)[0]


def inverse_primitive(form):
    a, b, c = form
    return reduce_primitive((a, -b, c))


def composition_data(left, right):
    """Return the primitive product form and its exact ideal-lattice data."""
    discriminant = form_discriminant(left)
    if form_discriminant(right) != discriminant:
        raise ValueError("forms must have the same discriminant")
    a, b, c = left
    other_a, other_b, other_c = right
    parity = discriminant % 2
    constant = (parity*parity-discriminant)//4
    u, v = (b+parity)//2, (other_b+parity)//2
    vertical = parity-u-v
    first_gcd, first_x, first_y = extended_gcd(a, other_a)
    content, second_x, second_y = extended_gcd(first_gcd, vertical)
    coefficients = (second_x*first_x, second_x*first_y, second_y)
    horizontal = (coefficients[0]*(-a*v)
                  + coefficients[1]*(-other_a*u)
                  + coefficients[2]*(u*v-constant))
    assert horizontal % content == 0
    assert a*other_a % (content*content) == 0
    leading = a*other_a//(content*content)
    middle = -2*(horizontal//content)-parity
    middle += 2*leading*((leading-middle)//(2*leading))
    numerator = middle*middle-discriminant
    assert numerator % (4*leading) == 0
    primitive = (leading, middle, numerator//(4*leading))
    assert form_discriminant(primitive) == discriminant
    return {
        "primitive_form": primitive, "content": content,
        "ideal_generators": ((a*other_a, 0), (-a*v, a),
                             (-other_a*u, other_a), (u*v-constant, vertical)),
        "vertical_bezout": coefficients,
        "vertical_one_after_content": (horizontal//content, 1),
        "primitive_ideal_index": leading,
    }


def compose_primitive(left, right):
    """Gauss composition by integer ideal multiplication and extended gcd."""
    return reduce_primitive(composition_data(left, right)["primitive_form"])


def power_primitive(form, exponent):
    discriminant = form_discriminant(form)
    parity = discriminant % 2
    result = (1, parity, (parity-discriminant)//4)
    if exponent < 0:
        form, exponent = inverse_primitive(form), -exponent
    while exponent:
        if exponent & 1:
            result = compose_primitive(result, form)
        exponent >>= 1
        if exponent:
            form = compose_primitive(form, form)
    return result
