#!/usr/bin/env python3
"""Exact arithmetic checks for the generic Mat12/Mat16 dimension bounds.

These finite calculations check the dimension tables and the elementary
Jordan arithmetic.  The algebraic-geometric proposition is proved in the
companion note, not by this bounded execution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time


def partitions(n, maximum=None):
    if n == 0:
        yield ()
        return
    if maximum is None:
        maximum = n
    for first in range(min(n, maximum), 0, -1):
        for tail in partitions(n-first, first):
            yield (first,)+tail


def centralizer_by_blocks(blocks):
    return sum(min(a, b) for a in blocks for b in blocks)


def centralizer_by_columns(blocks):
    return sum(sum(a >= k for a in blocks)**2
               for k in range(1, max(blocks, default=0)+1))


def dimension_table(a, c):
    b = a*c
    strata = {}
    for u in range(1, a+1):
        for v in range(1, c+1):
            r = u*v
            dimension = u*(2*a-u)+v*(2*c-v)-2
            strata[r] = max(strata.get(r, -1), dimension)
    rows = []
    for r, dimension in sorted(strata.items()):
        if r == b:
            continue
        bound = (b-r)**2+r*r
        rows.append({"rank": r, "product_stratum_dimension": dimension,
                     "centralizer_bound": bound,
                     "incidence_bound": 2*dimension+bound})
    # A trace-nonzero four-dimensional factor in characteristic two has
    # at least two eigenvalues, hence every multiplicity is at most three.
    eigenvalue_cap = 3*a
    full_centralizer = max(sum(m*m for m in multiplicities)
                           for multiplicities in partitions(b)
                           if max(multiplicities) <= eigenvalue_cap)
    full = {"rank": b, "product_stratum_dimension": strata[b],
            "eigenvalue_multiplicity_cap": eigenvalue_cap,
            "centralizer_bound": full_centralizer,
            "incidence_bound": 2*strata[b]+full_centralizer}
    maximum = max(row["incidence_bound"] for row in rows+[full])
    assert maximum < b*b
    return {"factor_dimensions": [a, c], "auxiliary_dimension": b,
            "nonfull_rank_rows": rows, "full_rank_row": full,
            "maximum_incidence_dimension": maximum,
            "general_linear_group_dimension": b*b,
            "dimension_gap": b*b-maximum}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    tables = [dimension_table(3, 4), dimension_table(4, 4)]
    expected = {
        12: [(1,10,122,142),(2,15,104,134),(3,18,90,126),
             (4,19,80,118),(6,21,72,114),(8,22,80,124),
             (9,22,90,134)],
        16: [(1,12,226,250),(2,17,200,234),(3,20,178,218),
             (4,22,160,204),(6,25,136,186),(8,26,128,180),
             (9,28,130,186),(12,29,160,218)],
    }
    for table in tables:
        actual = [tuple(row[key] for key in
                        ("rank", "product_stratum_dimension",
                         "centralizer_bound", "incidence_bound"))
                  for row in table["nonfull_rank_rows"]]
        assert actual == expected[table["auxiliary_dimension"]]
    assert tables[0]["full_rank_row"]["incidence_bound"] == 136
    assert tables[1]["full_rank_row"]["incidence_bound"] == 220
    partition_checks = rank_bound_checks = 0
    for b in range(1, 17):
        for zero_size in range(b+1):
            for zero_blocks in partitions(zero_size):
                zc = centralizer_by_blocks(zero_blocks)
                assert zc == centralizer_by_columns(zero_blocks)
                partition_checks += 1
                rank = b-len(zero_blocks)
                # Any choice of nonzero-eigenvalue Jordan blocks has
                # centralizer at most the square of their total size.
                # Thus this verifies the worst remaining contribution.
                worst = zc+(b-zero_size)**2
                assert worst <= (b-rank)**2+rank**2
                rank_bound_checks += 1
    output = {"status": "passed", "dimension_tables": tables,
              "jordan_block_column_identity_checks": partition_checks,
              "rank_centralizer_bound_checks": rank_bound_checks,
              "maximum_matrix_order_in_jordan_checks": 16,
              "scope": "Finite arithmetic verification only; generic existence follows from the companion incidence-dimension proof.",
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "seconds": time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
