"""Exact ratio classes via short lattice lines and integer quadratic roots.

At R**3 >= the maximum positive rectangle product, at most twenty
lattice lines are examined per ratio class. The complete reference
count has 2**(n/3)*poly(n) time and polynomial bit space.
"""

from __future__ import annotations

from math import isqrt

from circuit_box_kloosterman_blocks import odd_square_roots


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def reduced_affine_lattice(modulus: int, slope: int, width: int, height: int) -> dict:
    """A two-dimensional size reduction with a polynomial swap potential."""
    if min(modulus, width, height) <= 0:
        raise ValueError("positive lattice parameters required")
    first = (modulus, 0)
    second = (-slope, 1)
    scale_x, scale_y = height * height, width * width

    def inner(left, right):
        return scale_x * left[0] * right[0] + scale_y * left[1] * right[1]

    swaps = reductions = 0
    while True:
        norm = inner(first, first)
        quotient = (2 * inner(first, second) + norm) // (2 * norm)
        second = (second[0] - quotient * first[0], second[1] - quotient * first[1])
        reductions += 1
        if 4 * inner(second, second) >= 3 * norm:
            break
        first, second = second, first
        swaps += 1
    determinant = first[0] * second[1] - first[1] * second[0]
    if determinant < 0:
        second = (-second[0], -second[1])
        determinant = -determinant
    if determinant != modulus:
        raise AssertionError("lattice determinant changed")
    return {
        "first": first,
        "second": second,
        "first_norm": inner(first, first),
        "second_norm": inner(second, second),
        "inner_product": inner(first, second),
        "determinant": determinant,
        "swaps": swaps,
        "reductions": reductions,
    }


def _quadratic_integer_roots(quadratic: int, linear: int, constant: int) -> list[int]:
    if quadratic == 0:
        if linear == 0:
            if constant == 0:
                raise AssertionError("a positive product curve cannot contain this nonconstant line")
            return []
        return [-constant // linear] if constant % linear == 0 else []
    discriminant = linear * linear - 4 * quadratic * constant
    if discriminant < 0:
        return []
    root = isqrt(discriminant)
    if root * root != discriminant:
        return []
    denominator = 2 * quadratic
    return sorted({
        numerator // denominator
        for numerator in (-linear + root, -linear - root)
        if numerator % denominator == 0
    })


def _line_interval(origin, direction, first, second):
    lower = upper = None
    for base, step, (left, right) in zip(origin, direction, (first, second)):
        if step > 0:
            low = _ceil_div(left - base, step)
            high = (right - 1 - base) // step
        elif step < 0:
            low = _ceil_div(base - right + 1, -step)
            high = (base - left) // (-step)
        elif left <= base < right:
            continue
        else:
            return None
        lower = low if lower is None else max(lower, low)
        upper = high if upper is None else min(upper, high)
    if lower is None or upper is None:
        raise AssertionError("zero lattice direction")
    return (lower, upper) if lower <= upper else None


def ratio_lattice_factors(number, first, second, ratio_precision, ratio) -> dict:
    """All exact factor pairs in one ratio class; charge every lattice line."""
    if number < 1 or not (number & ratio & 1) or ratio_precision < 1:
        raise ValueError("positive odd target and ratio, positive precision required")
    if min(first[0], second[0]) < 1 or first[1] < first[0] or second[1] < second[0]:
        raise ValueError("positive ordered rectangle required")
    modulus = 1 << ratio_precision
    ratio %= modulus
    roots = odd_square_roots(ratio * number, ratio_precision)
    pairs = set()
    line_queries = maximum_lines = swaps = reductions = singleton_queries = 0
    for row in roots:
        column = number * pow(row, -1, modulus) % modulus
        rows = tuple(_ceil_div(endpoint - row, modulus) for endpoint in first)
        columns = tuple(_ceil_div(endpoint - column, modulus) for endpoint in second)
        width, height = rows[1] - rows[0], columns[1] - columns[0]
        if not width or not height:
            continue
        if width == 1:
            a = row + modulus * rows[0]
            b, remainder = divmod(number, a)
            if remainder == 0 and second[0] <= b < second[1] and (a - ratio * b) % modulus == 0:
                pairs.add((a, b))
            singleton_queries += 1
            continue
        if height == 1:
            b = column + modulus * columns[0]
            a, remainder = divmod(number, b)
            if remainder == 0 and first[0] <= a < first[1] and (a - ratio * b) % modulus == 0:
                pairs.add((a, b))
            singleton_queries += 1
            continue
        target = ((number - row * column) // modulus) * pow(column, -1, modulus) % modulus
        basis = reduced_affine_lattice(modulus, ratio, width, height)
        direction, transverse = basis["first"], basis["second"]
        swaps += basis["swaps"]
        reductions += basis["reductions"]
        corner_crosses = [
            direction[0] * y - direction[1] * (x - target)
            for x in (rows[0], rows[1] - 1)
            for y in (columns[0], columns[1] - 1)
        ]
        start = _ceil_div(min(corner_crosses), modulus)
        stop = max(corner_crosses) // modulus + 1
        maximum_lines = max(maximum_lines, stop - start)
        for line in range(start, stop):
            line_queries += 1
            origin = (target + line * transverse[0], line * transverse[1])
            interval = _line_interval(origin, direction, rows, columns)
            if interval is None:
                continue
            a0, b0 = row + modulus * origin[0], column + modulus * origin[1]
            da, db = modulus * direction[0], modulus * direction[1]
            candidates = _quadratic_integer_roots(da * db, da * b0 + db * a0, a0 * b0 - number)
            for candidate in candidates:
                if interval[0] <= candidate <= interval[1]:
                    pair = (a0 + da * candidate, b0 + db * candidate)
                    if pair[0] * pair[1] != number:
                        raise AssertionError("incorrect integer quadratic root")
                    pairs.add(pair)
    return {
        "pairs": sorted(pairs),
        "count": len(pairs),
        "ratio_precision": ratio_precision,
        "ratio": ratio,
        "root_count": len(roots),
        "line_queries": line_queries,
        "maximum_lines_per_root": maximum_lines,
        "singleton_queries": singleton_queries,
        "lattice_swaps": swaps,
        "lattice_reductions": reductions,
    }


def lattice_factor_rectangle(number, first, second) -> dict:
    """Complete exact count and least row, streaming the cube-root ratio grid."""
    if number < 1 or not number & 1:
        raise ValueError("positive odd target required")
    if min(first[0], second[0]) < 1 or first[1] < first[0] or second[1] < second[0]:
        raise ValueError("positive ordered rectangle required")
    result = {
        "count": 0,
        "least_row": None,
        "ratio_precision": 0,
        "ratio_queries": 0,
        "line_queries": 0,
        "maximum_lines_per_root": 0,
        "singleton_queries": 0,
        "lattice_swaps": 0,
        "lattice_reductions": 0,
    }
    if first[0] == first[1] or second[0] == second[1]:
        return result
    minimum_product = first[0] * second[0]
    maximum_product = (first[1] - 1) * (second[1] - 1)
    result["maximum_product"] = maximum_product
    if not minimum_product <= number <= maximum_product:
        return result
    precision = max(1, (maximum_product.bit_length() + 2) // 3)
    modulus = 1 << precision
    result["ratio_precision"] = precision
    ratios = range(number % 8, modulus, 8) if precision >= 3 else range(1, modulus, 2)
    for ratio in ratios:
        current = ratio_lattice_factors(number, first, second, precision, ratio)
        if current["maximum_lines_per_root"] > 5:
            raise AssertionError("cube-root lattice terminal exceeded its line bound")
        result["ratio_queries"] += 1
        result["count"] += current["count"]
        for key in ("line_queries", "singleton_queries", "lattice_swaps", "lattice_reductions"):
            result[key] += current[key]
        result["maximum_lines_per_root"] = max(result["maximum_lines_per_root"], current["maximum_lines_per_root"])
        if current["pairs"]:
            row = current["pairs"][0][0]
            result["least_row"] = row if result["least_row"] is None else min(result["least_row"], row)
    return result


def smallest_factor_by_ratio_lattice(number: int) -> dict:
    if number < 2:
        raise ValueError("input must be at least two")
    totals = {"rectangle_queries": 0, "ratio_queries": 0, "line_queries": 0}
    if number % 2 == 0:
        return {"factor": 2, **totals}
    stop, start = isqrt(number) + 1, 2
    while start < stop:
        first = (start, min(2 * start, stop))
        second = (number // (2 * start) + 1, number // start + 1)
        result = lattice_factor_rectangle(number, first, second)
        totals["rectangle_queries"] += 1
        totals["ratio_queries"] += result["ratio_queries"]
        totals["line_queries"] += result["line_queries"]
        if result["least_row"] is not None:
            return {"factor": result["least_row"], **totals}
        start *= 2
    return {"factor": number, **totals}
