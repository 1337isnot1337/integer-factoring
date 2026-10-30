#!/usr/bin/env python3
"""Linear-bit-space regular F_2 factor queries with compressed partitions.

The field helper uses scalar arithmetic only. A partition is a sentinel
followed by Elias-gamma codes of its successive positive gaps. This module
never expands a partition or matrix into a list.
"""

def polynomial_remainder(left, right):
    while left and left.bit_length() >= right.bit_length():
        left ^= right << (left.bit_length() - right.bit_length())
    return left


def polynomial_gcd(left, right):
    while right:
        left, right = right, polynomial_remainder(left, right)
    return left


def field_multiply(left, right, modulus, degree):
    value = 0
    while right:
        if right & 1:
            value ^= left
        right >>= 1
        left <<= 1
        if left >> degree:
            left ^= modulus
    return value


def irreducible_binary(degree):
    if degree == 1:
        return 3
    for candidate in range((1 << degree) | 1, 1 << (degree + 1), 2):
        value = 2
        valid = True
        for index in range(1, degree + 1):
            value = field_multiply(value, value, candidate, degree)
            if index <= degree // 2 and polynomial_gcd(value ^ 2, candidate) != 1:
                valid = False
                break
        if valid and value == 2:
            return candidate
    raise AssertionError("An irreducible polynomial must exist")


class BinaryField:
    """A polynomial-basis field with no element/power lookup table."""

    def __init__(self, degree):
        self.degree = degree
        self.order = 1 << degree
        self.modulus = irreducible_binary(degree)

    def multiply(self, left, right):
        return field_multiply(left, right, self.modulus, self.degree)

    def power(self, element, exponent):
        value = 1
        while exponent:
            if exponent & 1:
                value = self.multiply(value, element)
            element = self.multiply(element, element)
            exponent >>= 1
        return value


def coefficient_from_values(field, degree_bound, index, evaluate):
    """Extract one coefficient; evaluate is a streamed scalar callback."""
    if index < 0 or index > degree_bound:
        return 0
    assert field.order - 1 > degree_bound
    value = 0
    exponent = (-index) % (field.order - 1)
    for point in range(1, field.order):
        value ^= field.multiply(evaluate(point), field.power(point, exponent))
    assert value in (0, 1)
    return value


def encode_partition(parts):
    """Convenience input adapter; the core algorithms accept the code directly."""
    code = 1
    previous = 0
    for part in parts:
        gap = part - previous
        assert gap > 0
        code = (code << (2 * gap.bit_length() - 1)) | gap
        previous = part
    return code


def coded_gaps(code):
    remaining = code.bit_length() - 1
    while remaining:
        zeros = 0
        while ((code >> (remaining - 1)) & 1) == 0:
            remaining -= 1
            zeros += 1
        size = zeros + 1
        assert remaining >= size
        gap = (code >> (remaining - size)) & ((1 << size) - 1)
        remaining -= size
        yield gap


def coded_parts(code):
    current = 0
    for gap in coded_gaps(code):
        current += gap
        yield current


def coded_stats(code):
    count = weight = last = final_gap = 0
    for gap in coded_gaps(code):
        count += 1
        last += gap
        weight += last
        final_gap = gap
    return count, weight, last, final_gap


def partition_stream(weight):
    """Iterative DFS: one gap code and scalar state, with no frame stack."""
    if weight < 0:
        return
    if weight == 0:
        yield 1
        return
    code = 1
    used = last = 0
    candidate = 1
    while True:
        if candidate <= weight - used:
            gap = candidate - last
            code = (code << (2 * gap.bit_length() - 1)) | gap
            used += candidate
            last = candidate
            if used == weight:
                yield code
                code >>= 2 * gap.bit_length() - 1
                used -= candidate
                last = coded_stats(code)[2]
                candidate += 2
            else:
                candidate = last + 2
        else:
            if code == 1:
                return
            _, _, old_last, gap = coded_stats(code)
            code >>= 2 * gap.bit_length() - 1
            used -= old_last
            last = coded_stats(code)[2]
            candidate = old_last + 2


def partition_count(weight, cutoff=None):
    count = 0
    for _ in partition_stream(weight):
        count += 1
        if cutoff is not None and count >= cutoff:
            break
    return count


def partition_unrank(weight, index):
    if index < 0:
        raise IndexError(index)
    for code in partition_stream(weight):
        if not index:
            return code
        index -= 1
    raise IndexError("compressed partition rank outside its domain")


def label_weight(bits):
    if not bits:
        return 0
    for weight in range(1, 8 * bits * bits + 1):
        if partition_count(weight, 1 << bits) >= 1 << bits:
            return weight
    raise AssertionError("The explicit block construction guarantees success")


def coded_contains(code, part):
    for candidate in coded_parts(code):
        if candidate >= part:
            return candidate == part
    return False


def coded_is_unaffected(output, mask, part):
    for position, candidate in enumerate(coded_parts(output)):
        if candidate >= part:
            return candidate == part and (mask >> position) & 1
    return False


def unaffected_masks(label, output, input_):
    required = coded_stats(output)[0] + coded_stats(input_)[0] - coded_stats(label)[0]
    if required < 0 or required & 1:
        return
    required //= 2
    for mask in range(1 << coded_stats(output)[0]):
        if mask.bit_count() != required:
            continue
        valid = True
        for position, part in enumerate(coded_parts(output)):
            if (mask >> position) & 1 and not coded_contains(input_, part):
                valid = False
                break
        if valid:
            yield mask


def signed_columns(output, input_, mask):
    for position, part in enumerate(coded_parts(output)):
        if not ((mask >> position) & 1):
            yield True, part
    for part in coded_parts(input_):
        if not coded_is_unaffected(output, mask, part):
            yield False, part


def top_value(label, output, input_, point, field):
    length, h, _, _ = coded_stats(label)
    value = 0
    for mask in unaffected_masks(label, output, input_):
        for subset in range(1 << length):
            product = 1
            for n in coded_parts(label):
                row_sum = 0
                for column, (creation, k) in enumerate(signed_columns(output, input_, mask)):
                    if (subset >> column) & 1:
                        exponent = (n & k) - 1 if creation else n & ~k
                        row_sum ^= field.power(point, exponent)
                product = field.multiply(product, row_sum)
                if not product:
                    break
            value ^= product
    return field.multiply(value, field.power(1 ^ point, h + length))


def top_factor(label, index, output, input_, field=None):
    _, h, _, _ = coded_stats(label)
    if not h:
        return int(index == 0 and output == input_)
    if index < 0 or index > 2 * h:
        return 0
    if field is None:
        field = BinaryField((2 * h + 1).bit_length())
    return coefficient_from_values(
        field, 2 * h, index,
        lambda point: top_value(label, output, input_, point, field),
    )


def bottom_factor(label, output, input_):
    length = coded_stats(label)[0]
    value = 0
    for mask in unaffected_masks(label, output, input_):
        for subset in range(1 << length):
            product = 1
            for n in coded_parts(label):
                row_sum = 0
                for column, (creation, k) in enumerate(signed_columns(output, input_, mask)):
                    if (subset >> column) & 1:
                        row_sum ^= int(((n - 1) & ~(k - 1) == 0) if creation
                                       else ((n - 1) & k == 0))
                product &= row_sum
                if not product:
                    break
            value ^= product
    return value


def channel_stream(h, shift, exchanged):
    for angle in range(2 * h + 1):
        weight = shift + angle if exchanged else 2 * h - shift - angle
        for code in partition_stream(weight):
            yield angle, code


def channel_count(h, shift, exchanged):
    return sum(1 for _ in channel_stream(h, shift, exchanged))


def channel_unrank(h, shift, exchanged, index):
    if index < 0:
        raise IndexError(index)
    for channel in channel_stream(h, shift, exchanged):
        if not index:
            return channel
        index -= 1
    raise IndexError("compressed channel rank outside its domain")


def piece_entry(x, y, u, v, shift, exchanged=False):
    h = coded_stats(x)[1]
    assert coded_stats(y)[1] == coded_stats(u)[1] == coded_stats(v)[1] == h
    field = BinaryField((2 * h + 1).bit_length())
    value = 0
    for angle, intermediate in channel_stream(h, shift, exchanged):
        if exchanged:
            lower = bottom_factor(v, u, intermediate)
            if lower:
                value ^= top_factor(x, angle, intermediate, y, field)
        else:
            lower = bottom_factor(v, intermediate, y)
            if lower:
                value ^= top_factor(x, angle, u, intermediate, field)
    return value
