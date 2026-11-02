#!/usr/bin/env python3
"""Finite controls for the trace-character ANF selection counterexample."""
import argparse
from hashlib import sha256
import json
from math import log2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/trace-anf-selection-counterexample-controls.json'

WIDTH = 3
SIZE = 1 << WIDTH
FIELDS = (
    ('x^3+x+1', 0b1011),
    ('x^3+x^2+1', 0b1101),
)


def gf_mul(left, right, width, modulus):
    result = 0
    for index in range(width):
        if (right >> index) & 1:
            result ^= left << index
    for index in range(2 * width - 2, width - 1, -1):
        if (result >> index) & 1:
            result ^= modulus << (index - width)
    return result


def gf_pow(value, exponent, width, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = gf_mul(result, value, width, modulus)
        value = gf_mul(value, value, width, modulus)
        exponent >>= 1
    return result


def is_primitive(width, modulus):
    value = 1
    for exponent in range(1, 1 << width):
        value = gf_mul(value, 2, width, modulus)
        if value == 1:
            return exponent == (1 << width) - 1
    return False


def field_trace(value, width, modulus):
    result = 0
    for _ in range(width):
        result ^= value
        value = gf_mul(value, value, width, modulus)
    return result


def anf_coefficients(values, modulus):
    coefficients = [value % modulus for value in values]
    size = len(values)
    width = size.bit_length() - 1
    for index in range(width):
        step = 1 << index
        for mask in range(size):
            if mask & step:
                coefficients[mask] = (
                    coefficients[mask] - coefficients[mask ^ step]
                ) % modulus
    return coefficients


def anf_degree(values, modulus):
    return max(
        (mask.bit_count()
         for mask, value in enumerate(anf_coefficients(values, modulus))
         if value),
        default=0,
    )


def linear_trace(coefficient, point, modulus):
    return field_trace(gf_mul(coefficient, point, WIDTH, modulus),
                       WIDTH, modulus)


def indicator(support, size=SIZE):
    selected = set(support)
    return [1 if point in selected else 0 for point in range(size)]


def fourier_sum(coefficients, targets, point, modulus):
    total = 0
    for mask in range(1 << len(coefficients)):
        sign = -1 if sum(
            ((mask >> index) & 1) * targets[index]
            for index in range(len(coefficients))
        ) & 1 else 1
        exponent = 0
        for index, coefficient in enumerate(coefficients):
            if (mask >> index) & 1:
                exponent ^= linear_trace(coefficient, point, modulus)
        total += sign * (-1 if exponent else 1)
    return total


def case_record(modulus, coefficients, targets):
    support = [
        point for point in range(SIZE)
        if all(linear_trace(coefficient, point, modulus) == target
               for coefficient, target in zip(coefficients, targets))
    ]
    values = indicator(support)
    degree_f2 = anf_degree(values, 2)
    degree_f3 = anf_degree(values, 3)
    entropy = WIDTH - log2(len(support)) if support else 0.0
    claimed_degree = log2(2)
    return dict(
        modulus=modulus,
        coefficients=list(coefficients),
        targets=list(targets),
        support=support,
        support_size=len(support),
        degree_f2=degree_f2,
        degree_f3=degree_f3,
        entropy=entropy,
        claimed_degree=claimed_degree,
        degree_bound_fails=degree_f3 > claimed_degree,
        entropy_bound_fails=entropy > claimed_degree,
    )


def counterexample_controls():
    posted = case_record(0b1011, (4, 7), (1, 0))
    assert posted['support'] == [3, 6]
    assert posted['support_size'] == 2
    assert posted['degree_f2'] == 2
    assert posted['degree_f3'] == 3
    assert posted['entropy'] == 2.0
    assert posted['claimed_degree'] == 1.0
    assert posted['degree_bound_fails']
    assert posted['entropy_bound_fails']

    single = case_record(0b1101, (1,), (1,))
    assert single['support'] == [1, 2, 4, 7]
    assert single['support_size'] == 4
    assert single['degree_f2'] == 1
    assert single['degree_f3'] == 3
    assert single['entropy'] == 1.0
    assert single['degree_bound_fails']
    assert not single['entropy_bound_fails']

    return dict(posted=posted, single_predicate=single)


def family_controls():
    records = []
    for name, modulus in FIELDS:
        assert is_primitive(WIDTH, modulus)
        for coefficient in range(1, SIZE):
            for target in (0, 1):
                records.append(dict(
                    field=name,
                    **case_record(modulus, (coefficient,), (target,)),
                ))
        for first in range(1, SIZE):
            for second in range(1, SIZE):
                for targets in ((0, 0), (0, 1), (1, 0), (1, 1)):
                    record = case_record(modulus, (first, second), targets)
                    if record['support']:
                        records.append(dict(field=name, **record))
    for record in records:
        if record['support']:
            assert record['degree_f3'] + 1e-9 >= record['entropy']
    degree_failures = [record for record in records
                       if record['degree_bound_fails']]
    entropy_failures = [record for record in records
                        if record['entropy_bound_fails']]
    assert degree_failures
    assert entropy_failures
    assert any(record['support'] == [3, 6] for record in records)
    return dict(
        field_records=len(records),
        degree_bound_failures=len(degree_failures),
        entropy_bound_failures=len(entropy_failures),
        maximum_degree_f3=max(record['degree_f3'] for record in records),
        minimum_entropy=min(record['entropy'] for record in records
                            if record['support']),
        posted_witness=posted_witness(records),
    )


def posted_witness(records):
    for record in records:
        if record['support'] == [3, 6]:
            return record
    raise AssertionError('posted witness missing')


def xor_degree_controls():
    records = []
    for variables in range(1, WIDTH + 1):
        values = [
            sum((point >> index) & 1 for index in range(variables)) & 1
            for point in range(SIZE)
        ]
        records.append(dict(
            variables=variables,
            degree_f2=anf_degree(values, 2),
            degree_f3=anf_degree(values, 3),
        ))
        assert anf_degree(values, 2) == 1
        assert anf_degree(values, 3) == variables
    trace_bits = [linear_trace(1, point, 0b1101) for point in range(SIZE)]
    assert trace_bits == [
        sum((point >> index) & 1 for index in range(WIDTH)) & 1
        for point in range(SIZE)
    ]
    return dict(records=records, trace_u_degree_f3=anf_degree(trace_bits, 3))


def fourier_loss_controls():
    modulus = 0b1101
    traces = [linear_trace(1, point, modulus) for point in range(SIZE)]
    characters = [1 - 2 * trace for trace in traces]
    selected = [1 - trace for trace in traces]
    integer_sum = [1 + character for character in characters]
    assert integer_sum == [2 * value for value in selected]
    assert all(value % 2 == 0 for value in integer_sum)
    assert any(selected)
    assert [value % 3 for value in integer_sum] == [
        (2 * value) % 3 for value in selected
    ]
    assert [value % 2 for value in integer_sum] == [0] * SIZE
    return dict(
        modulus=modulus,
        integer_sum=integer_sum,
        selected=selected,
        mod2=[value % 2 for value in integer_sum],
        mod3=[value % 3 for value in integer_sum],
        sign_degree_f3=anf_degree(characters, 3),
        trace_degree_f3=anf_degree(traces, 3),
    )


def f2_lemma_controls():
    records = []
    for name, modulus in FIELDS:
        for exponent in range(2 * SIZE):
            terms = [(1, (exponent,))]
            values = [
                field_trace(
                    gf_pow(point, exponent, WIDTH, modulus),
                    WIDTH, modulus,
                )
                for point in range(SIZE)
            ]
            degree = anf_degree(values, 2)
            bound = log2(exponent + 1) if exponent else 0.0
            assert degree + 1e-9 >= 0
            assert degree <= bound + 1e-9
            records.append(dict(field=name, exponent=exponent,
                                degree_f2=degree, bound=bound))
        for first in range(5):
            for second in range(5 - first):
                degree_total = first + second
                values = []
                for point in range(SIZE * SIZE):
                    value = gf_mul(
                        gf_pow(point & (SIZE - 1), first, WIDTH, modulus),
                        gf_pow(point >> WIDTH, second, WIDTH, modulus),
                        WIDTH, modulus,
                    )
                    values.append(field_trace(value, WIDTH, modulus))
                observed = anf_degree(values, 2)
                bound = 2 * log2(degree_total + 1)
                assert observed <= bound + 1e-9
                records.append(dict(field=name, variables=2,
                                    degree=degree_total,
                                    degree_f2=observed, bound=bound))
    return records


def report():
    sources = [
        'tools/check_trace_anf_counterexample.py',
        'research/trace-anf-selection-counterexample.md',
        'research/artin-hasse-trace-projection.md',
    ]
    result = dict(
        scope=(
            'Finite controls for the counterexample to the trace-character '
            'selection-entropy bound and for the corrected F_2-ANF degree '
            'lemma. They verify the stated finite identities and arithmetic, '
            'not an independent proof, an Artin-Hasse representation, or a '
            'factoring improvement.'
        ),
        base='dd84ac3',
        counterexamples=counterexample_controls(),
        family=family_controls(),
        xor_degrees=xor_degree_controls(),
        fourier_loss=fourier_loss_controls(),
        f2_lemma=f2_lemma_controls(),
        source_sha256={
            path: sha256((ROOT/path).read_bytes()).hexdigest()
            for path in sources
        },
    )
    return json.loads(json.dumps(result))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = report()
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    else:
        assert result == json.loads(RECEIPT.read_text())
    print(json.dumps(dict(
        scope=result['scope'],
        family_records=result['family']['field_records'],
        degree_bound_failures=result['family']['degree_bound_failures'],
        entropy_bound_failures=result['family']['entropy_bound_failures'],
        f2_lemma_records=len(result['f2_lemma']),
    )))


if __name__ == '__main__':
    main()
