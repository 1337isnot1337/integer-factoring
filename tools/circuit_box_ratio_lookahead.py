"""Positive ratio-class counts and an adaptive two-bit product search.

One ratio bit is added at each level, while the product congruence gains
two bits. The full search has a cube-root-scale worst-case bound; it is
not a replacement for the repository's stronger half-log factoring bound.
"""

from __future__ import annotations

from math import isqrt

from circuit_box_kloosterman_blocks import affine_interval_pairs, odd_square_roots


def _condition(interval: tuple[int, int], residue: int, modulus: int) -> tuple[int, int]:
    return tuple((endpoint - residue + modulus - 1) // modulus for endpoint in interval)


def ratio_lookahead_count(
    number: int,
    first: tuple[int, int],
    second: tuple[int, int],
    ratio_precision: int,
    ratio: int,
) -> dict:
    """Count ab=N mod 2**(2h), a=ratio*b mod 2**h on two intervals."""
    if not (number & ratio & 1) or ratio_precision < 1:
        raise ValueError("odd target and ratio, and positive ratio precision, required")
    if first[1] < first[0] or second[1] < second[0]:
        raise ValueError("reversed interval endpoints")
    modulus = 1 << ratio_precision
    ratio %= modulus
    roots = odd_square_roots(ratio * number, ratio_precision)
    answer = queries = 0
    for row in roots:
        column = number * pow(row, -1, modulus) % modulus
        rows = _condition(first, row, modulus)
        columns = _condition(second, column, modulus)
        if rows[0] == rows[1] or columns[0] == columns[1]:
            continue
        target = ((number - row * column) // modulus) * pow(column, -1, modulus)
        answer += affine_interval_pairs(rows, columns, ratio, target, modulus)
        queries += 1
    return {
        "count": answer,
        "ratio_precision": ratio_precision,
        "product_precision": 2 * ratio_precision,
        "ratio": ratio,
        "root_count": len(roots),
        "affine_pair_queries": queries,
    }


def adaptive_factor_rectangle(
    number: int,
    first: tuple[int, int],
    second: tuple[int, int],
) -> dict:
    """Exact factor-pair count by positive, nested ratio lookahead queries."""
    if number < 1 or not number & 1:
        raise ValueError("positive odd target required")
    if first[0] < 1 or second[0] < 1:
        raise ValueError("positive rectangle required")
    if first[1] < first[0] or second[1] < second[0]:
        raise ValueError("reversed interval endpoints")
    output = {
        "count": 0,
        "terminal_ratio_precision": 0,
        "visited_nodes": 0,
        "live_nodes": 0,
        "affine_pair_queries": 0,
        "maximum_stack_size": 0,
        "levels": [],
    }
    if first[0] == first[1] or second[0] == second[1]:
        output["early_exit"] = "empty interval"
        return output
    minimum_product = first[0] * second[0]
    maximum_product = (first[1] - 1) * (second[1] - 1)
    output["minimum_product"] = minimum_product
    output["maximum_product"] = maximum_product
    if not minimum_product <= number <= maximum_product:
        output["early_exit"] = "target outside positive product range"
        return output
    height = max(number - minimum_product, maximum_product - number)
    terminal = max(1, (height.bit_length() + 1) // 2)
    output["maximum_difference"] = height
    output["terminal_ratio_precision"] = terminal
    levels = [
        {"ratio_precision": h, "visited": 0, "live": 0, "congruence_pairs": 0}
        for h in range(1, terminal + 1)
    ]
    stack = [(1, 1)]
    while stack:
        output["maximum_stack_size"] = max(output["maximum_stack_size"], len(stack))
        precision, ratio = stack.pop()
        current = ratio_lookahead_count(number, first, second, precision, ratio)
        output["visited_nodes"] += 1
        output["affine_pair_queries"] += current["affine_pair_queries"]
        level = levels[precision - 1]
        level["visited"] += 1
        if not current["count"]:
            continue
        level["live"] += 1
        level["congruence_pairs"] += current["count"]
        output["live_nodes"] += 1
        if precision == terminal:
            output["count"] += current["count"]
        else:
            # The stack holds at most one pending sibling at each level.
            stack.append((precision + 1, ratio + (1 << precision)))
            stack.append((precision + 1, ratio))
    output["levels"] = levels
    return output


def smallest_factor_by_ratio_lookahead(number: int) -> dict:
    """Complete deterministic reference factoring, with every search counted."""
    if number < 2:
        raise ValueError("input must be at least two")
    if not number & 1:
        return {"factor": 2, "rectangle_queries": 0, "visited_nodes": 0, "affine_pair_queries": 0}
    total = {"rectangle_queries": 0, "visited_nodes": 0, "affine_pair_queries": 0}

    def query(first, second):
        result = adaptive_factor_rectangle(number, first, second)
        total["rectangle_queries"] += 1
        total["visited_nodes"] += result["visited_nodes"]
        total["affine_pair_queries"] += result["affine_pair_queries"]
        return result["count"]

    stop = isqrt(number) + 1
    start = 2
    while start < stop:
        first = (start, min(2 * start, stop))
        second = (number // (2 * start) + 1, number // start + 1)
        if query(first, second):
            left, right = first
            while right - left > 1:
                middle = (left + right) // 2
                if query((left, middle), second):
                    right = middle
                else:
                    left = middle
            if number % left:
                raise AssertionError("positive exact count did not recover a divisor")
            return {"factor": left, **total}
        start *= 2
    return {"factor": number, **total}
