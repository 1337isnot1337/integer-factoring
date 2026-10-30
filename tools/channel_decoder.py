"""Modular rank-56 channel decoder for the proposed crossing pieces.

This implements a single factor evaluation, not the recursive evaluator or a
factoring algorithm. Labels use the encoding in research/integration.md.
"""

from tools.check_crossing import crossing_c, crossing_r, psi, residual_left, residual_right


def channel_factor(width: int, index: int, label: int, first: int, second: int,
                   modulus: int, *, left: bool) -> int:
    """Evaluate one side of a padded row channel modulo any modulus >= 2.

    For a column channel, transpose the inputs and use index t-i. Values
    outside the 16^t * 7^index encoding give zero.
    """
    if width < 0 or index < 0 or label < 0 or modulus < 2:
        raise ValueError("invalid width, index, label, or modulus")
    t = (width + 2) // 3
    if index > t or first < 0 or second < 0 or first.bit_length() > width or second.bit_length() > width:
        raise ValueError("index or input outside its range")
    mask = label % (1 << t)
    rest = label >> t
    classes = []
    for _ in range(t):
        classes.append(rest % 8)
        rest //= 8
    digits = []
    for _ in range(index):
        digits.append(rest % 7)
        rest //= 7
    selected_count = mask.bit_count()
    if rest or selected_count > index or any(digits[selected_count:]):
        return 0

    if left:
        # Pascal's recurrence avoids division, including for composite moduli.
        degree = index - selected_count
        row = [1] + [0] * degree
        for _ in range(t - selected_count):
            for j in range(degree, 0, -1):
                row[j] = (row[j] + row[j - 1]) % modulus
        value = row[degree] * pow(4, t - index, modulus) * pow(-3, degree, modulus) % modulus
    else:
        value = 1
    selected = 0
    for s, k in enumerate(classes):
        p, q = (first >> (3 * s)) & 7, (second >> (3 * s)) & 7
        if mask >> s & 1:
            a = digits[selected] + 1
            selected += 1
            factor = residual_left(k, a, p, q) if left else residual_right(k, a, p, q)
        else:
            factor = psi(k, p, q) if left else psi(k, q, p)
        value = value * factor % modulus
    return value


def piece_factor(width: int, index: int, label: int, first: int, second: int,
                 modulus: int, *, left: bool) -> int:
    """Evaluate the narrower row factor, including the trivial-width fallback.

    Column orientation uses transposed inputs and index t-index. The label
    has min(2*width, 4*t+ceil(73*index/26)) bits; padding evaluates to zero.
    """
    if width < 0 or index < 0 or label < 0 or modulus < 2:
        raise ValueError("invalid width, index, label, or modulus")
    t = (width + 2) // 3
    if index > t or first < 0 or second < 0 or first.bit_length() > width or second.bit_length() > width:
        raise ValueError("index or input outside its range")
    expanded_bits = 4 * t + (73 * index + 25) // 26
    if 2 * width > expanded_bits:
        return channel_factor(width, index, label, first, second, modulus, left=left)
    if label.bit_length() > 2 * width:
        return 0
    p, q = label >> width, label & ((1 << width) - 1)
    if left:
        return int(first == p and second == q)

    # Coefficient of prod(R + z*C), truncated to the required degree.
    coefficients = [1] + [0] * index
    for s in range(t):
        block = tuple((value >> (3 * s)) & 7 for value in (p, q, first, second))
        r, c = crossing_r(*block), crossing_c(*block)
        for j in range(index, 0, -1):
            coefficients[j] = (r * coefficients[j] + c * coefficients[j - 1]) % modulus
        coefficients[0] = r * coefficients[0] % modulus
    return coefficients[index] * pow(4, t, modulus) % modulus
