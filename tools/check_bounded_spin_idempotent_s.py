"""Check the central-idempotent rank sum of the ten-site bounded-spin model.

The published ten-site model gives the coefficient algebra, its seven
orthogonal idempotents and the cut-reciprocity matrix.  This script uses
only the rank formula from ``check_bounded_spin_growth.check_strata`` and
the model JSON; it does not build a cut operator or a tensor.  It is a
finite diagnostic for the central decomposition, not a rank theorem.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "evidence" / "bounded-spin-ten-algebra-model.json"
PRIME = 5


def dot(a, b):
    return sum(int(x) * int(y) for x, y in zip(a, b)) % PRIME


def ranks(model, forms, vector):
    dimensions = model["dimensions"]
    nil = model["nil_component"]
    nil_rank = model["nil_rank"]
    vanishing = [dot(form, vector) == 0 for form in forms]
    result = []
    for cut in (0, 1):
        value = sum(d for j, d in enumerate(dimensions)
                    if not vanishing[cut * 8 + j])
        if vanishing[cut * 8 + nil] and not vanishing[cut * 8 + 7]:
            value += nil_rank
        result.append(value)
    return result


def main() -> int:
    model = json.loads(MODEL.read_text(encoding="utf-8"))
    characters = model["characters"]
    idempotents = model["idempotents"]
    reciprocity = model["optimization_reciprocity"]
    radical = model["radical"]
    pivot = next(i for i, x in enumerate(radical) if x)
    inverse = pow(radical[pivot], -1, PRIME)
    # The derivative form is  (e_p - sum_j e_j[p] chi_j) / radical[p].
    derivative = [
        (int(i == pivot)
         - sum(idempotents[j][pivot] * characters[j][i]
               for j in range(len(idempotents)))) * inverse % PRIME
        for i in range(8)
    ]
    row_forms = characters + [derivative]
    forms = row_forms + [[
        sum(form[i] * reciprocity[j][i] for i in range(8)) % PRIME
        for j in range(8)
    ] for form in row_forms]

    total = 0.0
    pairs = []
    for vector in idempotents:
        row, column = ranks(model, forms, vector)
        pairs.append((row, column))
        total += math.sqrt(row * column)

    expected_rows = model["dimensions"]
    expected_columns = [931, 1156, 356, 581, 771, 776, 715]
    assert [r for r, _ in pairs] == expected_rows
    assert [c for _, c in pairs] == expected_columns
    assert total > 34 ** math.sqrt(2)
    print("idempotent row/column rank pairs:", pairs)
    print(f"central square-root rank sum: {total:.12f}")
    print(f"ten-site target: {34 ** math.sqrt(2):.12f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
