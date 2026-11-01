"""Compare encoded-word prefix contractions with literal integer boxes."""
import itertools
import json
import random
from hashlib import sha256
from pathlib import Path

from check_carry_fourier_elimination import BinaryField
from check_frobenius_rational_trace import field_of_degree
from check_joint_frobenius_carries import evaluate
from check_trace_constraint_contraction import reindex
from frobenius_prefix_contraction import prefix_trace, dual_basis, trace_bit, dyadic_prefixes, interval_trace


def literal(field, width, groups, basis, num, den, mask, intervals=False):
    words = sum(g[0] for g in groups)
    answer = 0
    assignments = 0
    for digits in itertools.product(range(1 << width), repeat=words):
        assignments += 1
        if any(nz and not digit for nz, digit in zip(mask, digits)):
            continue
        offset = 0
        allowed = True
        for size, length, value in groups:
            integer = sum(digits[offset + j] << (j * width) for j in range(size))
            accepted = length <= integer < value if intervals else integer >> (size * width - length) == value
            if not accepted:
                allowed = False
                break
            offset += size
        if not allowed:
            continue
        point = []
        for digit in digits:
            element = 0
            for j, vector in enumerate(basis):
                if digit & (1 << j):
                    element ^= vector
            point.append(element)
        a = evaluate(field, words, width, num, point)
        b = evaluate(field, words, width, den, point)
        answer ^= field.mul(a, field.inverse(b)) if b else 0
    return answer, assignments


def main():
    rng = random.Random(318748)
    counts = dict(queries=0, direct_assignments=0, all_fixed=0,
                  empty_fixed_domains=0, nonzero_answers=0, dual_pairings=0,
                  dyadic_covers=0, interval_queries=0)
    for bits in range(7):
        for lower in range((1 << bits) + 1):
            for upper in range(lower, (1 << bits) + 1):
                prefixes = dyadic_prefixes(bits, lower, upper)
                assert len(prefixes) <= 2 * bits + 1
                expanded = [x for length, value in prefixes
                            for x in range(value << (bits - length), (value + 1) << (bits - length))]
                assert expanded == list(range(lower, upper))
                counts['dyadic_covers'] += 1
    for width, degree in [(1, 1), (2, 2), (2, 4), (3, 3), (4, 4)]:
        field = BinaryField(3) if degree == 1 else field_of_degree(degree)
        # Choose a basis by an independent exhaustive span check, only here.
        span = {0}
        basis = []
        for a in reversed(range(field.size)):
            if field.power(a, 1 << width) == a and a not in span:
                basis.append(a)
                span |= {x ^ a for x in span}
        assert len(basis) == width
        dual = dual_basis(field, width, basis)
        for i, a in enumerate(dual):
            for j, b in enumerate(basis):
                assert trace_bit(field, width, field.mul(a, b)) == (i == j)
                counts['dual_pairings'] += 1
        cases = []
        if width <= 2:
            cases.extend([(2, length, value)] for length in range(2 * width + 1)
                         for value in range(1 << length))
        else:
            cases.extend([(2, length, rng.randrange(1 << length))]
                         for length in range(2 * width + 1))
        for _ in range(12):
            groups = []
            for _ in range(2):
                length = rng.randrange(width + 1)
                groups.append((1, length, rng.randrange(1 << length)))
            cases.append(groups)
        for index, groups in enumerate(cases):
            num = (index % width, [(1, (1, 1)), (rng.randrange(field.size), (0, 0))])
            den = (0, [(1, (1, 0)), (rng.randrange(field.size), (0, 0))])
            if index % 13 == 0:
                den = (0, [])
            for mask in [(False, False), (True, False), (False, True), (True, True)]:
                expected, assignments = literal(field, width, groups, basis, num, den, mask)
                counts['direct_assignments'] += assignments
                counts['nonzero_answers'] += bool(expected)
                for t in range(1, width + 1):
                    if width % t:
                        continue
                    actual, stats = prefix_trace(field, width, t, groups, basis,
                                                 [reindex(num, width, t, 2)],
                                                 [reindex(den, width, t, 2)], mask)
                    assert actual == expected, (width, degree, groups, mask, t, actual, expected)
                    counts['queries'] += 1
                    counts['all_fixed'] += stats['active_variables'] == 0
                    counts['empty_fixed_domains'] += stats.get('empty_fixed_domain', False)
                    if 'rank' in stats:
                        assert stats['rank'] == sum(length % width for _, length, _ in groups)
        for trial in range(10):
            boxes = []
            for _ in range(2):
                lower = rng.randrange((1 << width) + 1)
                upper = rng.randrange(lower, (1 << width) + 1)
                boxes.append((1, lower, upper))
            num = (trial % width, [(1, (1, 1)), (1, (0, 0))])
            den = (0, [(1, (1, 0)), (1, (0, 0))])
            for mask in [(False, False), (True, False), (False, True), (True, True)]:
                expected, assignments = literal(field, width, boxes, basis, num, den, mask, intervals=True)
                actual, _ = interval_trace(field, width, width, boxes, basis, [num], [den], mask)
                assert actual == expected, (width, degree, boxes, mask, actual, expected)
                counts['direct_assignments'] += assignments
                counts['interval_queries'] += 1
    assert counts['nonzero_answers'] > 10 and counts['all_fixed'] > 0
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/frobenius_prefix_contraction.py', 'tools/check_frobenius_prefix_contraction.py',
               'tools/trace_constraint_contraction.py', 'tools/joint_frobenius_carries.py',
               'tools/square_zero_trace_contraction.py']
    receipt = {'scope': 'Exact encoded-word prefix controls; arithmetic phase representation remains open',
               'counts': counts,
               'sha256': {p: sha256((root / p).read_bytes()).hexdigest() for p in sources}}
    (root / 'evidence/frobenius-prefix-contraction.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
