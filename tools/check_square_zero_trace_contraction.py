"""Literal field sums and packet controls for square-zero trace aggregation."""
import json
import random
from hashlib import sha256
from pathlib import Path

from check_carry_fourier_elimination import BinaryField
from check_frobenius_rational_trace import field_of_degree
from check_trace_constraint_contraction import direct, reindex
from square_zero_trace_contraction import SquareZeroAlgebra, square_zero_constrained_trace
from trace_constraint_contraction import constrained_trace


def main():
    rng = random.Random(318749)
    counts = dict(queries=0, packet_comparisons=0, direct_assignments=0,
                  nonzero_answers=0, inconsistent=0, maximum_rank=0, nilpotent_controls=0)
    for width, degree in [(1, 1), (2, 2), (3, 3), (2, 4), (4, 4)]:
        field = BinaryField(3) if degree == 1 else field_of_degree(degree)
        algebra = SquareZeroAlgebra(field, 3)
        generators = [1 << (field.degree * (1 << k)) for k in range(3)]
        for z in generators:
            assert z != 0 and algebra.mul(z, z) == 0
            assert algebra.mul(1 ^ z, 1 ^ z) == 1
            counts['nilpotent_controls'] += 2
        assert algebra.mul(algebra.mul(generators[0], generators[1]), generators[2]) == 1 << (7 * field.degree)
        counts['nilpotent_controls'] += 1
        subfield = [u for u in range(field.size) if field.power(u, 1 << width) == u]
        for trial in range(8):
            constraints = [(tuple(rng.choice(subfield) for _ in range(2)), rng.randrange(2))
                           for _ in range(trial % 4)]
            if trial == 4:
                constraints = [((1, 1), 0), ((1, 1), 1)]
            if trial == 5:
                constraints = [((1, 1), 0), ((1, 1), 0)]
            if trial == 6 and width == degree == 2:
                constraints = [((1, 0), 0), ((2, 0), 1), ((0, 1), 1), ((0, 2), 0)]
            num = (trial % width, [(1, (1, 1)), (rng.randrange(field.size), (0, 0))])
            den = (0, [(1, (1, 0)), (rng.randrange(field.size), (0, 0))])
            if trial == 7:
                den = (0, [])
            for mask in [(False, False), (True, False), (False, True), (True, True)]:
                expected, assignments = direct(field, width, 2, num, den, constraints, mask)
                counts['direct_assignments'] += assignments
                counts['nonzero_answers'] += bool(expected)
                for t in sorted({1, width}):
                    nums = [reindex(num, width, t, 2)]
                    dens = [reindex(den, width, t, 2)]
                    actual, stats = square_zero_constrained_trace(field, width, t, 2, nums, dens, constraints, mask)
                    packet, _ = constrained_trace(field, width, t, 2, nums, dens, constraints, mask)
                    assert actual == packet == expected, (width, degree, trial, mask, t, actual, expected)
                    counts['queries'] += 1
                    counts['packet_comparisons'] += 1
                    counts['inconsistent'] += stats['inconsistent']
                    counts['maximum_rank'] = max(counts['maximum_rank'], stats['rank'] or 0)
    assert counts['maximum_rank'] == 4 and counts['nonzero_answers'] > 10
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/square_zero_trace_contraction.py', 'tools/check_square_zero_trace_contraction.py',
               'tools/joint_frobenius_carries.py', 'tools/trace_constraint_contraction.py']
    receipt = {'scope': 'Exact square-zero aggregation controls; no factoring score improvement',
               'counts': counts,
               'sha256': {p: sha256((root / p).read_bytes()).hexdigest() for p in sources}}
    (root / 'evidence/square-zero-trace-contraction.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
