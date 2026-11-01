"""Binary operand prefixes with fixed words substituted before contraction."""
from itertools import product

from trace_constraint_contraction import constrained_trace


def trace_bit(field, width, value):
    answer = 0
    for _ in range(width):
        answer ^= value
        value = field.mul(value, value)
    if answer not in (0, 1):
        raise ValueError('Trace input outside subfield')
    return answer


def dual_basis(field, width, basis):
    if len(basis) != width or any(not 0 <= a < field.size or field.power(a, 1 << width) != a for a in basis):
        raise ValueError('Bad subfield basis')
    rows = []
    for i, a in enumerate(basis):
        row = sum(trace_bit(field, width, field.mul(a, b)) << j for j, b in enumerate(basis))
        rows.append(row | (1 << (width + i)))
    for j in range(width):
        pivot = next((i for i in range(j, width) if (rows[i] >> j) & 1), None)
        if pivot is None:
            raise ValueError('Dependent subfield basis')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(width):
            if i != j and (rows[i] >> j) & 1:
                rows[i] ^= rows[j]
    result = []
    for row in rows:
        value = 0
        for j, a in enumerate(basis):
            if (row >> (width + j)) & 1:
                value ^= a
        result.append(value)
    return result


def encode_word(value, basis):
    answer = 0
    for j, a in enumerate(basis):
        if (value >> j) & 1:
            answer ^= a
    return answer


def specialize(field, width, t, factors, fixed, active):
    """Fold fixed-variable powers into coefficients; preserve each row."""
    D = width // t
    out = []
    for row, terms in factors:
        combined = {}
        for coefficient, counts in terms:
            for i, value in fixed.items():
                exponent = sum(counts[i * D + j] << (row + t * j) for j in range(D))
                coefficient = field.mul(coefficient, field.power(value, exponent))
            kept = tuple(counts[i * D + j] for i in active for j in range(D))
            combined[kept] = combined.get(kept, 0) ^ coefficient
        out.append((row, [(c, h) for h, c in combined.items() if c]))
    return out


def prefix_trace(field, width, t, groups, basis, numerators, denominators, nonzero=False):
    """groups: (word_count, high_prefix_length, prefix_value), low words first.

    Groups occupy consecutive variable blocks. Each is an integer encoded
    by little-endian width-bit words in the supplied binary subfield basis.
    """
    if width < 1 or t < 1 or field.degree % width or width % t:
        raise ValueError('Incompatible dimensions')
    dual = dual_basis(field, width, basis)
    variables = sum(words for words, _, _ in groups)
    mask = [nonzero] * variables if isinstance(nonzero, bool) else list(nonzero)
    if len(mask) != variables:
        raise ValueError('Bad domain mask')
    fixed = {}
    partial = []
    offset = 0
    for words, length, value in groups:
        if words < 1 or not 0 <= length <= words * width or not 0 <= value < (1 << length):
            raise ValueError('Bad prefix')
        whole, remainder = divmod(length, width)
        for k in range(whole):
            digit = (value >> (length - (k + 1) * width)) & ((1 << width) - 1)
            fixed[offset + words - 1 - k] = encode_word(digit, basis)
        if remainder:
            i = offset + words - whole - 1
            for j in range(width - remainder, width):
                partial.append((i, j, (value >> (j - width + remainder)) & 1))
        offset += words
    if any(mask[i] and not value for i, value in fixed.items()):
        return 0, {'empty_fixed_domain': True, 'active_variables': variables - len(fixed)}
    active = [i for i in range(variables) if i not in fixed]
    nums = specialize(field, width, t, numerators, fixed, active)
    dens = specialize(field, width, t, denominators, fixed, active)
    if not active:
        answer = 1
        for _, terms in nums:
            answer = field.mul(answer, terms[0][0] if terms else 0)
        for _, terms in dens:
            value = terms[0][0] if terms else 0
            answer = field.mul(answer, field.inverse(value)) if value else 0
        return answer, {'active_variables': 0, 'rank': 0, 'summands': 1}
    positions = {old: new for new, old in enumerate(active)}
    constraints = []
    for old, bit_position, bit in partial:
        coefficients = [0] * len(active)
        coefficients[positions[old]] = dual[bit_position]
        constraints.append((tuple(coefficients), bit))
    answer, stats = constrained_trace(field, width, t, len(active), nums, dens,
                                      constraints, [mask[i] for i in active])
    stats['active_variables'] = len(active)
    stats['fixed_words'] = len(fixed)
    return answer, stats


def dyadic_prefixes(bits, lower, upper):
    """Disjoint binary prefixes for the half-open interval [lower, upper)."""
    if bits < 0 or not 0 <= lower <= upper <= (1 << bits):
        raise ValueError('Bad integer interval')
    result = []
    while lower < upper:
        available = (upper - lower).bit_length() - 1
        aligned = (lower & -lower).bit_length() - 1 if lower else bits
        free = min(available, aligned)
        result.append((bits - free, lower >> free))
        lower += 1 << free
    return result


def interval_trace(field, width, t, boxes, basis, numerators, denominators, nonzero=False):
    """boxes: (word_count, lower, upper), with half-open integer intervals."""
    covers = [[(words, length, value) for length, value in
               dyadic_prefixes(words * width, lower, upper)]
              for words, lower, upper in boxes]
    answer = 0
    cells = 0
    for groups in product(*covers):
        value, _ = prefix_trace(field, width, t, groups, basis, numerators, denominators, nonzero)
        answer ^= value
        cells += 1
    return answer, {'prefix_cells': cells}
