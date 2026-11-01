"""Independent literal sums for the sparse binary-constraint interface."""
import itertools
import json
import random
from hashlib import sha256
from pathlib import Path

from check_carry_fourier_elimination import BinaryField
from check_frobenius_rational_trace import field_of_degree
from check_joint_frobenius_carries import evaluate
from joint_frobenius_carries import joint_trace
from trace_constraint_contraction import constrained_trace, trace_packets


def absolute_trace(field, width, value):
    answer = 0
    for _ in range(width):
        answer ^= value
        value = field.mul(value, value)
    assert answer in (0, 1)
    return answer


def direct(field, width, variables, numerator, denominator, constraints, mask):
    subfield = [u for u in range(field.size) if field.power(u, 1 << width) == u]
    answer = 0
    assignments = 0
    for point in itertools.product(*([u for u in subfield if u or not nz] for nz in mask)):
        assignments += 1
        accepted = True
        for coefficients, bit in constraints:
            value = 0
            for a, u in zip(coefficients, point):
                value ^= field.mul(a, u)
            if absolute_trace(field, width, value) != bit:
                accepted = False
                break
        if accepted:
            num = evaluate(field, variables, width, numerator, point)
            den = evaluate(field, variables, width, denominator, point)
            answer ^= field.mul(num, field.inverse(den)) if den else 0
    return answer, assignments


def reindex(factor, width, t, variables):
    row, terms = factor
    out = []
    for c, h in terms:
        counts = [0] * (variables * (width // t))
        for i in range(variables):
            counts[i * (width // t) + row // t] = h[i]
        out.append((c, tuple(counts)))
    return row % t, out


def main():
    rng = random.Random(197331)
    counts = dict(queries=0, direct_assignments=0, inconsistent=0,
                  redundant=0, unsplit_controls=0, nonzero_answers=0)
    for width, degree in [(1, 1), (2, 2), (3, 3), (2, 4), (4, 4)]:
        field = BinaryField(3) if degree == 1 else field_of_degree(degree)
        subfield = [u for u in range(field.size) if field.power(u, 1 << width) == u]
        variables = 2
        for trial in range(8):
            a = tuple(rng.choice(subfield) for _ in range(variables))
            b = tuple(rng.choice(subfield) for _ in range(variables))
            constraints = [(a, trial % 2)]
            if trial == 0:
                constraints = []
            elif trial == 1:
                constraints += [(a, trial % 2)]
            elif trial == 2:
                constraints += [(a, 1 - trial % 2)]
            elif trial >= 3:
                constraints += [(b, (trial // 2) % 2)]
                if trial == 4:
                    constraints += [(tuple(x ^ y for x, y in zip(a, b)),
                                     (trial % 2) ^ ((trial // 2) % 2))]
            numerator = (trial % width, [(1, (1, 1)), (rng.randrange(field.size), (0, 0))])
            denominator = (0, [(1, (1, 0)), (rng.randrange(field.size), (0, 0))])
            if trial == 7:
                denominator = (0, [(0, (0, 0))])
            for mask in [(False, False), (True, False), (False, True), (True, True)]:
                expected, assignments = direct(field, width, variables, numerator, denominator, constraints, mask)
                counts['direct_assignments'] += assignments
                counts['nonzero_answers'] += bool(expected)
                for t in range(1, width + 1):
                    if width % t:
                        continue
                    num = reindex(numerator, width, t, variables)
                    den = reindex(denominator, width, t, variables)
                    actual, stats = constrained_trace(field, width, t, variables, [num], [den], constraints, mask)
                    assert actual == expected, (width, degree, trial, mask, t, actual, expected)
                    counts['queries'] += 1
                    counts['inconsistent'] += stats['inconsistent']
                    if not stats['inconsistent']:
                        assert stats['summands'] <= t ** stats['rank']
                        counts['redundant'] += stats['rank'] < len(constraints)
                    if t == 1:
                        # No elimination: multiply the original full indicators.
                        full = [num]
                        for coefficients, bit in constraints:
                            packets = trace_packets(field, width, 1, coefficients, bit)
                            full.append(packets[0] if packets else (0, []))
                        assert joint_trace(field, degree, width, 1, variables, full, [den], mask) == expected
                        counts['unsplit_controls'] += 1
    assert counts['nonzero_answers'] > 10
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/trace_constraint_contraction.py', 'tools/check_trace_constraint_contraction.py',
               'tools/joint_frobenius_carries.py']
    receipt = {'scope': 'Exact finite controls; no stronger factoring bound', 'counts': counts,
               'sha256': {p: sha256((root / p).read_bytes()).hexdigest() for p in sources}}
    (root / 'evidence/trace-constraint-contraction.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
