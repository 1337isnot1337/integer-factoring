#!/usr/bin/env python3
"""Binary-size optimal scales and terminal windows for balanced rectangles."""


def binomial_valuation(left: int, right: int) -> int:
    """v2(binomial(left+right,left)), without expanding a factorial."""
    if not isinstance(left, int) or not isinstance(right, int) or min(left, right) < 0:
        raise ValueError("Require two nonnegative encoded integers")
    return left.bit_count() + right.bit_count() - (left + right).bit_count()


def rectangular_bezout_scale(rows: int, columns: int) -> dict:
    """Return max_{0<=k<rows} columns+k-v2(binomial(columns+k-1,k)).

    The maximizing index lies in a final binary-length window. This function
    does not evaluate any polynomial or decide whether input domains intersect.
    """
    if any(not isinstance(x, int) or x < 1 for x in (rows, columns)):
        raise ValueError("Require positive encoded row and column counts")
    top = rows + columns - 1
    bound = top.bit_length()
    first = max(0, rows - 1 - bound)
    best, index = -1, None
    for k in range(first, rows):
        value = columns + k - binomial_valuation(k, columns - 1)
        if value > best:
            best, index = value, k
    return {"exponent": best, "maximizing_index": index,
            "defect": top - best, "first_examined_index": first,
            "examined_indices": rows - first,
            "binary_window_bound": bound + 1}


def rectangular_newton_window(rows: int, columns: int, precision: int) -> dict:
    """Indices whose scaled inverse Newton coefficients can survive mod 2^p.

    Gaps are not monotone. The complete terminal interval is inspected before
    skipping any individually vanishing coefficients.
    """
    if not isinstance(precision, int) or precision < 1:
        raise ValueError("Require positive numerical precision")
    scale = rectangular_bezout_scale(rows, columns)
    exponent = scale["exponent"]
    first = max(0, exponent - columns - precision + 1)
    terms = []
    for k in range(first, rows):
        gap = exponent - columns - k + binomial_valuation(k, columns - 1)
        if not gap >= 0:
            raise ArithmeticError("The computed optimal scale failed")
        if gap < precision:
            terms.append({"index": k, "valuation_gap": gap})
    defect = scale["defect"]
    if rows - first > defect + precision:
        raise ArithmeticError("The terminal-window bound failed")
    # At reverse Newton offset ell the coefficient has valuation >=ell-defect.
    # Every root shift is even, so another j places in that monomial gains j.
    # The symmetric rectangle formula applies the same argument to both sides.
    reverse_v_depth = min(rows - 1, defect + precision - 1)
    reverse_u_depth = min(columns - 1, defect + precision - 1)
    return {**scale, "precision": precision, "first_newton_index": first,
            "examined_newton_indices": rows - first, "retained_newton_terms": terms,
            "reverse_v_depth": reverse_v_depth, "reverse_u_depth": reverse_u_depth}
