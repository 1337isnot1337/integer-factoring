"""Count diagonal evaluation forms of the bounded-spin coefficient span.

For each cut product state ``(a,b)`` this evaluates the full six-vertex
coefficient monomials at ``(a,b,a,b)``.  Distinct projective forms and
their state multiplicities are compressed.  The computation uses only
the published monomial-entry constructor; it does not build a cut
operator and makes no rank claim.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from bounded_spin_crossing import paths  # noqa: E402
from bounded_spin_six_vertex import monomial_entry  # noqa: E402

PRIME = 5


def diagonal_forms(length: int):
    path_list = paths(length)
    half = length // 2
    monomials = [(a, b, half - a - b)
                 for a in range(half + 1) for b in range(half - a + 1)]
    multiplicities = {}
    for left in path_list:
        for right in path_list:
            coefficients = monomial_entry((left, right, left, right))
            vector = tuple(coefficients.get(m, 0) % PRIME for m in monomials)
            inverse = next((pow(x, -1, PRIME) for x in vector if x), None)
            if inverse is None:
                continue
            projective = tuple(x * inverse % PRIME for x in vector)
            multiplicities[projective] = multiplicities.get(projective, 0) + 1
    return {
        "length": length,
        "labels": len(path_list),
        "monomials": [list(m) for m in monomials],
        "distinct_forms": len(multiplicities),
        "multiplicities": sorted(multiplicities.values(), reverse=True),
        "square_root_sum": sum(math.sqrt(m) for m in multiplicities.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-evidence", action="store_true")
    args = parser.parse_args()
    records = [diagonal_forms(n) for n in (8, 10, 12)]
    for record in records:
        target = record["labels"] ** math.sqrt(2)
        print(f"n={record['length']}: labels={record['labels']} "
              f"forms={record['distinct_forms']} "
              f"sqrt-sum={record['square_root_sum']:.6f} "
              f"target={target:.6f}")
    if args.write_evidence:
        path = ROOT / "evidence" / "bounded-spin-diagonal-weights.json"
        path.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
