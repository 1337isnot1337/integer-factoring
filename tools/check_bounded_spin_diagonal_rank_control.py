"""Control showing diagonal evaluation counts are not cut-rank bounds.

At length eight the full six-vertex coefficient span is small enough to
build the actual cut operators.  For several deterministic coefficient
combinations this prints the actual rank pair and the number of nonzero
diagonal evaluations.  The control asserts that both inequalities fail
in the sample, so the diagonal multiplicity diagnostic cannot be used as
a rank lower or upper bound without additional structure.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from bounded_spin_batched_blocks import coefficient_tensor  # noqa: E402
from check_bounded_spin_six_vertex import operators  # noqa: E402
from native_gf5_flint import rank  # noqa: E402


def main() -> int:
    n = 8
    prime = 5
    monomials, tensors = coefficient_tensor(n)
    q = tensors.shape[1]
    forms = [
        tuple(int(tensors[j, a, b, a, b]) % prime for j in range(len(monomials)))
        for a in range(q) for b in range(q)
    ]
    generator = random.Random(1337)
    samples = []
    for _ in range(12):
        weights = [generator.randrange(prime) for _ in monomials]
        tensor = np.tensordot(
            np.array(weights, dtype=np.int64), tensors, axes=(0, 0)
        ) % prime
        row, column = operators(tensor[None, :], n)
        actual = (rank(row[0]), rank(column[0]))
        diagonal = sum(
            int(sum(weights[j] * form[j] for j in range(len(monomials))) % prime != 0)
            for form in forms
        )
        samples.append((actual, diagonal))
    assert any(max(pair) > diagonal for pair, diagonal in samples)
    assert any(max(pair) < diagonal for pair, diagonal in samples)
    for actual, diagonal in samples:
        print(f"rank_pair={actual} nonzero_diagonal={diagonal}")
    print("control: diagonal count is neither a lower nor an upper rank bound")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
