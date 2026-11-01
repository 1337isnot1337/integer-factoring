#!/usr/bin/env python3
"""Near-linear-space queries for the original regular F_2 factors.

Production routines use scalar extension-field arithmetic and streamed
subset/partition loops. Polynomial arrays and independent source imports
occur only in verify(). The accompanying note proves the space bound.
"""

from itertools import combinations
import argparse
import json


from linear_modular_factors import BinaryField, coefficient_from_values


def unaffected_subsets(label, output, input_):
    count = len(output) + len(input_) - len(label)
    if count < 0 or count & 1:
        return
    count //= 2
    common = tuple(k for k in output if k in input_)
    if count <= len(common):
        yield from combinations(common, count)


def signed_columns(output, input_, unchanged):
    return (tuple((True, k) for k in output if k not in unchanged)
            + tuple((False, k) for k in input_ if k not in unchanged))


def top_value(label, output, input_, point, field):
    """Evaluate the folded angle profile without a polynomial/matrix array."""
    length = len(label)
    value = 0
    for unchanged in unaffected_subsets(label, output, input_):
        columns = signed_columns(output, input_, unchanged)
        assert len(columns) == length
        for subset in range(1 << length):
            product = 1
            for n in label:
                row_sum = 0
                for column, (creation, k) in enumerate(columns):
                    if (subset >> column) & 1:
                        exponent = (n & k) - 1 if creation else n & ~k
                        row_sum ^= field.power(point, exponent)
                product = field.multiply(product, row_sum)
                if not product:
                    break
            value ^= product
    return field.multiply(value, field.power(1 ^ point, sum(label) + length))


def top_factor(label, index, output, input_, field=None):
    h = sum(label)
    if not label:
        return int(index == 0 and tuple(output) == tuple(input_))
    if index < 0 or index > 2 * h:
        return 0
    if field is None:
        field = BinaryField((2 * h + 1).bit_length())
    return coefficient_from_values(
        field, 2 * h, index,
        lambda point: top_value(label, output, input_, point, field),
    )


def bottom_factor(label, output, input_):
    value = 0
    length = len(label)
    for unchanged in unaffected_subsets(label, output, input_):
        columns = signed_columns(output, input_, unchanged)
        for subset in range(1 << length):
            product = 1
            for n in label:
                row_sum = 0
                for column, (creation, k) in enumerate(columns):
                    if (subset >> column) & 1:
                        row_sum ^= int(((n - 1) & ~(k - 1) == 0) if creation
                                       else ((n - 1) & k == 0))
                product &= row_sum
                if not product:
                    break
            value ^= product
    return value


def regular_partitions_stream(weight, minimum=1):
    """Increasing distinct odd partitions; no memoization/count table."""
    if weight == 0:
        yield ()
        return
    for part in range(minimum | 1, weight + 1, 2):
        for tail in regular_partitions_stream(weight - part, part + 2):
            yield (part,) + tail


def partition_count(weight, cutoff=None):
    count = 0
    for _ in regular_partitions_stream(weight):
        count += 1
        if cutoff is not None and count >= cutoff:
            break
    return count


def partition_unrank(weight, index):
    if index < 0:
        raise IndexError(index)
    for partition in regular_partitions_stream(weight):
        if index == 0:
            return partition
        index -= 1
    raise IndexError("partition rank outside its domain")


def label_weight(bits):
    if bits == 0:
        return 0
    target = 1 << bits
    for weight in range(1, 8 * bits * bits + 1):
        if partition_count(weight, target) >= target:
            return weight
    raise AssertionError("The explicit equal-weight block bound guarantees success")


def channel_stream(h, shift, exchanged):
    for angle in range(2 * h + 1):
        weight = shift + angle if exchanged else 2 * h - shift - angle
        if weight >= 0:
            for partition in regular_partitions_stream(weight):
                yield angle, partition


def channel_count(h, shift, exchanged):
    return sum(1 for _ in channel_stream(h, shift, exchanged))


def channel_unrank(h, shift, exchanged, index):
    if index < 0:
        raise IndexError(index)
    for channel in channel_stream(h, shift, exchanged):
        if index == 0:
            return channel
        index -= 1
    raise IndexError("channel rank outside its domain")


def piece_entry(x, y, u, v, shift, exchanged=False):
    h = sum(x)
    assert sum(y) == sum(u) == sum(v) == h
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


def verify():
    # These reference tables/arrays are confined to finite verification.
    from check_modular_factor_space import streamed_top_profile, streamed_bottom_factor
    from check_modular_crossing import modular_angle, modular_operator, regular_partitions
    from itertools import product
    import random
    import linear_modular_factors as compressed

    field_cases = extraction_cases = 0
    for degree in range(1, 9):
        field = BinaryField(degree)
        for element in range(1, field.order):
            assert field.power(element, field.order - 1) == 1
            field_cases += 1
    rng = random.Random(931847)
    for degree_bound in range(1, 25):
        field = BinaryField((degree_bound + 1).bit_length())
        polynomials = list(range(1 << (degree_bound + 1))) if degree_bound <= 5 else [
            rng.getrandbits(degree_bound + 1) for _ in range(12)]
        for polynomial in polynomials:
            def evaluate(point):
                value = 0
                for k in range(degree_bound, -1, -1):
                    value = field.multiply(value, point) ^ ((polynomial >> k) & 1)
                return value
            for index in range(degree_bound + 1):
                assert coefficient_from_values(field, degree_bound, index, evaluate) == (
                    (polynomial >> index) & 1)
                extraction_cases += 1

    partitions_checked = unrankings = compressed_unrankings = 0
    maximum_code_bits = 0
    for weight in range(65):
        actual = tuple(regular_partitions_stream(weight))
        expected = {tuple(sorted(mu)) for mu in regular_partitions(2, weight)}
        assert len(actual) == len(set(actual))
        assert set(actual) == expected
        assert partition_count(weight) == len(expected)
        codes = tuple(compressed.partition_stream(weight))
        assert tuple(tuple(compressed.coded_parts(code)) for code in codes) == actual
        assert compressed.partition_count(weight) == len(expected)
        for index, partition in enumerate(actual):
            assert partition_unrank(weight, index) == partition
            unrankings += 1
            code = compressed.encode_partition(partition)
            assert compressed.partition_unrank(weight, index) == code == codes[index]
            assert compressed.coded_stats(code)[1] == weight
            bits = code.bit_length() - 1
            assert bits * bits <= 25 * weight
            maximum_code_bits = max(maximum_code_bits, bits)
            compressed_unrankings += 1
        partitions_checked += len(actual)

    weights = {}
    for bits in range(1, 7):
        weight = label_weight(bits)
        assert partition_count(weight) >= 1 << bits
        assert all(partition_count(k) < 1 << bits for k in range(weight))
        weights[str(bits)] = weight
        assert compressed.label_weight(bits) == weight

    states = tuple(mu for weight in range(10) for mu in regular_partitions_stream(weight))
    top_entries = bottom_entries = 0
    for h in range(1, 10):
        field = BinaryField((2 * h + 1).bit_length())
        for label in regular_partitions_stream(h):
            for output in states:
                for input_ in states:
                    profile = streamed_top_profile(label, output, input_)
                    for index in range(2 * h + 1):
                        actual = top_factor(label, index, output, input_, field)
                        assert actual == ((profile >> index) & 1)
                        assert actual == modular_angle(2, label, index, output, input_)
                        assert actual == compressed.top_factor(
                            compressed.encode_partition(label), index,
                            compressed.encode_partition(output), compressed.encode_partition(input_), field)
                        top_entries += 1
                    actual = bottom_factor(label, output, input_)
                    assert actual == streamed_bottom_factor(label, output, input_)
                    assert actual == modular_operator(2, label, output, input_)
                    assert actual == compressed.bottom_factor(compressed.encode_partition(label),
                        compressed.encode_partition(output), compressed.encode_partition(input_))
                    bottom_entries += 1

    # A four-particle label exercises more than the smallest finite cases.
    label = (1, 3, 5, 7)
    for output, input_ in [((1, 3, 5, 7), ()), ((), (1, 3, 5, 7)),
                           ((1, 3, 7, 9), ()), ((1, 3, 5), (7,))]:
        profile = streamed_top_profile(label, output, input_)
        for index in range(33):
            assert top_factor(label, index, output, input_) == ((profile >> index) & 1)
            assert compressed.top_factor(compressed.encode_partition(label), index,
                compressed.encode_partition(output), compressed.encode_partition(input_)) == ((profile >> index) & 1)
            top_entries += 1

    channels_checked = 0
    for h in range(1, 9):
        for shift in range(-2 * h, 2 * h + 1):
            for exchanged in (False, True):
                channels = tuple(channel_stream(h, shift, exchanged))
                expected_count = sum(len(regular_partitions(2,
                    shift + i if exchanged else 2 * h - shift - i))
                    for i in range(2 * h + 1)
                    if (shift + i if exchanged else 2 * h - shift - i) >= 0)
                assert channel_count(h, shift, exchanged) == len(channels) == expected_count
                assert compressed.channel_count(h, shift, exchanged) == expected_count
                for index, channel in enumerate(channels):
                    assert channel_unrank(h, shift, exchanged, index) == channel
                    angle, code = compressed.channel_unrank(h, shift, exchanged, index)
                    assert (angle, tuple(compressed.coded_parts(code))) == channel
                    channels_checked += 1

    pieces_checked = crossings_checked = 0
    for h in (1, 4, 8):
        labels = tuple(regular_partitions_stream(h))
        for x, y, u, v in product(labels, repeat=4):
            total = 0
            for shift in range(-2 * h, 2 * h + 1):
                first = piece_entry(x, y, u, v, shift)
                second = piece_entry(x, y, u, v, shift, True)
                assert first == second
                codes = tuple(compressed.encode_partition(mu) for mu in (x, y, u, v))
                assert compressed.piece_entry(*codes, shift) == first
                assert compressed.piece_entry(*codes, shift, True) == second
                total ^= first
                pieces_checked += 1
            assert total == int(x == v and y == u)
            crossings_checked += 1

    # The strict interpolation bound is necessary: S^(q-1) aliases with1.
    alias_field = BinaryField(2)
    aliased = 0
    for point in range(1, alias_field.order):
        aliased ^= alias_field.power(point, alias_field.order - 1)
    assert aliased == 1
    # Removing the common angle factor changes this source entry.
    assert [top_factor((1,), i, (1,), ()) for i in range(3)] == [1, 0, 1]

    return {
        "status": "pass",
        "scope": "exact local factors and streamed indexing; space bounds are proved in the accompanying note",
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "nonzero_field_elements_checked": field_cases,
        "independent_coefficient_extractions": extraction_cases,
        "regular_partitions_checked": partitions_checked,
        "partition_unrankings": unrankings,
        "compressed_partition_unrankings": compressed_unrankings,
        "maximum_partition_code_bits_through_weight_64": maximum_code_bits,
        "least_label_weights": weights,
        "top_angle_entries": top_entries,
        "bottom_factor_entries": bottom_entries,
        "channel_unrankings": channels_checked,
        "both_order_piece_entries": pieces_checked,
        "complete_crossing_entries": crossings_checked,
        "negative_controls": ["field size aliases the top degree", "omitted common angle factor"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    arguments = parser.parse_args()
    encoded = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    if arguments.output:
        from pathlib import Path
        Path(arguments.output).write_text(encoded)
    print(encoded, end="")
