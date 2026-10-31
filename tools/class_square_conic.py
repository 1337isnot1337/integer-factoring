"""Parametrize Q(x,y)=z^2 from a supplied class square root.

The quadratic map is explicit ideal squaring.  It does not find a smooth
denominator or solve the fixed-square representation search problem.
"""

from math import gcd

from quadratic_form_arithmetic import (
    composition_data, form_discriminant, reduce_with_transform,
)


def evaluate(form, x, y):
    a, b, c = form
    return a*x*x+b*x*y+c*y*y


def square_parametrization(root):
    """Return the reduced square and three integral quadratic polynomials.

    If their values are X,Y,Z at u,v, then Q(X,Y)=Z^2 and Z=F(u,v).
    Polynomials use the coefficient order (u^2, uv, v^2).
    """
    discriminant = form_discriminant(root)
    a, b, _ = root
    parity = discriminant % 2
    constant = (parity*parity-discriminant)//4
    offset = (b+parity)//2
    data = composition_data(root, root)
    content = data["content"]
    raw = data["primitive_form"]
    leading = raw[0]
    horizontal_basis = (-raw[1]-parity)//2
    target, transform = reduce_with_transform(raw)
    (u00, u01), (u10, u11) = transform

    def point(u, v):
        horizontal = (a*u+offset*v)**2-constant*v*v
        vertical = -2*a*u*v-b*v*v
        assert horizontal % content == vertical % content == 0
        other = vertical//content
        numerator = horizontal//content-horizontal_basis*other
        assert numerator % leading == 0
        first, second = numerator//leading, -other
        return (u11*first-u01*second, -u10*first+u00*second,
                evaluate(root, u, v))

    at_zero, at_infinity, at_one = point(1, 0), point(0, 1), point(1, 1)
    polynomials = tuple((left, middle-left-right, right)
                        for left, right, middle in zip(at_zero, at_infinity, at_one))
    assert polynomials[2] == tuple(root)
    assert sum(polynomials[0][i]*cross(polynomials[1], polynomials[2])[i]
               for i in range(3)) == -discriminant
    return {"root": tuple(root), "target": target, "polynomials": polynomials}


def parametrized_point(data, u, v):
    """Evaluate the integral map; nonzero primitive parameters are required."""
    if gcd(u, v) != 1:
        raise ValueError("primitive nonzero parameters are required")
    point = tuple(evaluate(poly, u, v) for poly in data["polynomials"])
    content = gcd(gcd(point[0], point[1]), point[2])
    assert content > 0
    assert (-form_discriminant(data["root"])) % content == 0
    assert evaluate(data["target"], point[0], point[1]) == point[2]**2
    return {"parameters": (u, v), "point": point, "content": content,
            "primitive_point": tuple(value//content for value in point)}


def cross(left, right):
    a, b, c = left
    x, y, z = right
    return b*z-c*y, c*x-a*z, a*y-b*x


def parameters_from_point(data, point):
    """Invert the projective conic map using its integral adjugate."""
    x, y, z = point
    if z <= 0 or evaluate(data["target"], x, y) != z*z:
        raise ValueError("a conic point with positive last coordinate is required")
    point_content = gcd(gcd(x, y), z)
    primitive = tuple(value//point_content for value in point)
    columns = tuple(zip(*data["polynomials"]))
    adjugate_rows = (cross(columns[1], columns[2]), cross(columns[2], columns[0]),
                     cross(columns[0], columns[1]))
    rank_one = tuple(sum(a*b for a, b in zip(row, primitive)) for row in adjugate_rows)
    assert rank_one[0]*rank_one[2] == rank_one[1]**2
    if rank_one[0] == 0:
        parameters = (0, 1)
    else:
        u, v = rank_one[:2]
        common = gcd(u, v)
        assert common
        u, v = u//common, v//common
        parameters = (-u, -v) if u < 0 or (u == 0 and v < 0) else (u, v)
    result = parametrized_point(data, *parameters)
    assert result["primitive_point"] == primitive
    adjugate_content = gcd(gcd(rank_one[0], rank_one[1]), rank_one[2])
    assert adjugate_content*result["content"] == -form_discriminant(data["root"])
    result["input_content"] = point_content
    return result


def representation_from_denominator(data, u, v, square_root):
    """Scale a primitive point when its positive denominator divides M."""
    result = parametrized_point(data, u, v)
    x, y, denominator = result["primitive_point"]
    if square_root <= 0 or square_root % denominator:
        raise ValueError("the normalized denominator must divide the supplied M")
    scale = square_root//denominator
    result["representation"] = (scale*x, scale*y)
    assert evaluate(data["target"], scale*x, scale*y) == square_root*square_root
    return result


def content_factor(data, point, integer):
    """Return a proper gcd certificate, if inverse-map content provides one."""
    result = parameters_from_point(data, point)
    divisor = gcd(result["content"], integer)
    result["factor"] = divisor if 1 < divisor < integer else None
    return result
