#!/usr/bin/env python3
"""Literal finite controls for the complete raw group-pair boundary minor.

The all-size argument is in research/joint-top-boundary-rank.md. These
checks do not establish a new factoring exponent or any rank bound after
vertex weights/internal constraints have been absorbed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
import random


class Source:
    def __init__(self, s: int):
        self.h = 1 << s
        self.a = 2 * self.h
        self.d = 2 * self.h * self.h
        self.base = 1 << self.d
        self.modulus = (1 << (4 * self.d)) + 1
        self.length = self.a * self.a
        self.root = 16
        self.inverse_root = pow(self.root, self.length - 1, self.modulus)
        self.inverse_length = pow(self.length, -1, self.modulus)
        self.gates = tuple((kind, i) for kind, length in
                           (("P", self.a), ("Q", self.a), ("M", self.a),
                            ("I", self.a), ("C", self.length))
                           for i in range(length))

    def zero(self, gate):
        kind, _ = gate
        if kind in ("P", "Q"):
            return (0,) * self.h
        if kind == "M":
            return ((0,) * self.a, (0,) * self.a)
        if kind == "I":
            return (0,) * self.a
        return (0, 0)

    def fa(self, r, digits, recipient):
        # The source's ext0 permits A summands; the last h summands are zero.
        return pow(16, r * recipient, self.modulus) * sum(
            digits[i] * pow(16, self.a * i * recipient, self.modulus)
            for i in range(self.h)) % self.modulus

    def middle(self, k, left, right, recipient):
        # Literal forward transforms, pointwise product, inverse stage.
        products = []
        for frequency in range(self.a):
            a = sum(left[r] * pow(16, self.a * r * frequency, self.modulus)
                    for r in range(self.a))
            b = sum(right[r] * pow(16, self.a * r * frequency, self.modulus)
                    for r in range(self.a))
            products.append(a * b % self.modulus)
        return (pow(self.inverse_root, k * recipient, self.modulus)
                * sum(products[f] * pow(self.inverse_root,
                                        self.a * f * recipient, self.modulus)
                      for f in range(self.a)) % self.modulus)

    def inverse(self, values, recipient):
        return self.inverse_length * sum(
            values[k] * pow(self.inverse_root, self.a * k * recipient,
                            self.modulus)
            for k in range(self.a)) % self.modulus

    def linked(self, left, right):
        a, i = left
        b, j = right
        if a in ("P", "Q") and b == "M":
            return True
        if a == "M" and b in ("P", "Q", "I"):
            return True
        if a == "I" and b == "M":
            return True
        if a == "I" and b == "C":
            return j % self.a == i
        if a == "C" and b == "I":
            return i % self.a == j
        return a == b == "C" and abs(i - j) == 1

    def message(self, left, right, setting):
        a, i = left
        b, j = right
        if a in ("P", "Q") and b == "M":
            return self.fa(i, setting, j)
        if a == "M" and b in ("P", "Q"):
            return setting[0 if b == "P" else 1][j]
        if a == "M" and b == "I":
            return self.middle(i, setting[0], setting[1], j)
        if a == "I" and b == "M":
            return setting[j]
        if a == "I" and b == "C":
            return self.inverse(setting, j // self.a)
        if a == "C" and b == "I":
            return setting[0]
        if a == b == "C":
            return ((setting[0] + setting[1]) // self.base
                    if j == i + 1 else setting[1])
        return 0

    def boundary(self, u, v, row, column):
        return int(all(self.message(a, b, row[a]) ==
                       self.message(b, a, column[b])
                       for a in u for b in v if self.linked(a, b)))


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, parts - 1):
                yield (first,) + rest


def minor_settings(source, u, v, digit_values):
    """Construct both sides independently from one common identity index."""
    qu = tuple(i for kind, i in u if kind == "Q")
    qv = tuple(i for kind, i in v if kind == "Q")
    mu = tuple(i for kind, i in u if kind == "M")
    mv = tuple(i for kind, i in v if kind == "M")
    tu, tv = min(len(mv), source.h), min(len(mu), source.h)
    width = len(qu) * tu + len(qv) * tv
    rows, columns = [], []
    for label in product(digit_values, repeat=width):
        cursor = 0
        source_digits = {}
        for indices, t in ((qu, tu), (qv, tv)):
            for a in indices:
                source_digits[a] = tuple(label[cursor:cursor + t]) + (0,) * (source.h - t)
                cursor += t
        row = {gate: source.zero(gate) for gate in u}
        col = {gate: source.zero(gate) for gate in v}
        for indices, settings in ((qu, row), (qv, col)):
            for a in indices:
                settings[("Q", a)] = source_digits[a]
        for indices, opposite_q, settings in ((mu, qv, row), (mv, qu, col)):
            for b in indices:
                values = [0] * source.a
                for a in opposite_q:
                    values[a] = source.fa(a, source_digits[a], b)
                settings[("M", b)] = ((0,) * source.a, tuple(values))
        rows.append(row)
        columns.append(col)
    return rows, columns, width


def check_minors(totals, records):
    rng = random.Random(0x51A7B0)
    for s, digit_values in ((0, tuple(range(4))), (1, (0, 1))):
        source = Source(s)
        core = tuple((kind, i) for kind in ("Q", "M") for i in range(source.a))
        others = tuple(g for g in source.gates if g not in core)
        # Every U/V/outside partition of the four s=0 core gates. The
        # s=1 suite adds deliberately bidirectional and randomized cases.
        if s == 0:
            assignments = list(product(range(3), repeat=len(core)))
        else:
            assignments = [tuple(rng.randrange(3) for _ in core) for _ in range(28)]
            assignments += [tuple((i % 2) for _, i in core),
                            tuple(0 if kind == "Q" else 1 for kind, _ in core)]
        for case, groups in enumerate(assignments):
            placement = dict(zip(core, groups))
            placement.update({g: rng.randrange(3) for g in others})
            u = tuple(g for g in source.gates if placement[g] == 0)
            v = tuple(g for g in source.gates if placement[g] == 1)
            if not u or not v:
                continue
            rows, columns, width = minor_settings(source, u, v, digit_values)
            # A full 256-by-256 minor is affordable; bound this suite's
            # largest declared rank explicitly rather than truncating rows.
            assert len(rows) <= 256
            links = [(a, b) for a in u for b in v if source.linked(a, b)]
            row_messages = [tuple(source.message(a, b, row[a]) for a, b in links)
                            for row in rows]
            col_messages = [tuple(source.message(b, a, col[b]) for a, b in links)
                            for col in columns]
            for i, row in enumerate(row_messages):
                for j, col in enumerate(col_messages):
                    actual = int(all(x == y for x, y in zip(row, col)))
                    assert actual == int(i == j), (s, case, i, j)
                    totals["literal_boundary_entries"] += 1
            # The cached tuples above use all original edge formulas.
            # Compare direct uncached entries at diagonal/corner controls.
            for i, j in {(0, 0), (len(rows) - 1, len(rows) - 1),
                         (0, len(rows) - 1)}:
                assert source.boundary(u, v, rows[i], columns[j]) == int(i == j)
                totals["uncached_boundary_entries"] += 1
            totals["complete_boundary_minors"] += 1
            if any(a[0] == "Q" and b[0] == "M" for a, b in links) and any(
                    a[0] == "M" and b[0] == "Q" for a, b in links):
                totals["bidirectional_minors"] += 1
            for a, b in links:
                totals["edge_type_" + "".join(sorted((a[0], b[0])))] += 1
            records.append({"s": s, "case": case, "groups": [list(u), list(v)],
                            "digit_alphabet": list(digit_values),
                            "free_digit_count": width, "identity_order": len(rows),
                            "all_edge_count": len(links)})
    assert totals["bidirectional_minors"] > 0
    assert all(totals["edge_type_" + kind] > 0 for kind in ("MP", "MQ", "IM", "CI", "CC"))


def check_tradeoff(totals):
    for h in (1, 2, 4):
        a, d = 2 * h, 2 * h * h
        capacity = a * h * d
        for p in (2, 3, 4):
            counts = list(compositions(a, p))
            for q in counts:
                for m in counts:
                    # Exact integer lower bound for log2 of the Cartesian
                    # setting domain; log2 M > 4d accounts for the 8Ad term.
                    local = max(h * d * qi + 8 * a * d * mi for qi, mi in zip(q, m))
                    channels = sum(d * (q[u] * min(m[v], h) + q[v] * min(m[u], h))
                                   for u in range(p) for v in range(u + 1, p))
                    assert 2 * local + channels >= 2 * capacity
                    if max(m) <= h:
                        assert channels == d * (a * a - sum(qi * mi for qi, mi in zip(q, m)))
                        assert 2 * h * d * max(q) + channels >= 2 * capacity
                    else:
                        assert local > 8 * capacity
                    totals["complete_count_compositions"] += 1


def check_units(totals):
    for s in range(5):
        source = Source(s)
        zeta = pow(16, source.a, source.modulus)
        assert pow(zeta, source.h, source.modulus) == source.modulus - 1
        roots = [pow(zeta, j, source.modulus) for j in range(source.a)]
        for i, j in combinations(range(source.a), 2):
            assert gcd(roots[i] - roots[j], source.modulus) == 1
            totals["source_unit_differences"] += 1
    source = Source(1)
    assert source.modulus == 641 * 6700417


def negative_controls(totals):
    source = Source(0)
    u, v = (("Q", 0), ("M", 1)), (("Q", 1), ("M", 0))
    rows, cols, _ = minor_settings(source, u, v, (0, 1, 2, 3))
    assert len(rows) == 16
    # Absorbing one local predicate kills rows. It is not part of the raw
    # boundary and cannot be dismissed by the identity-minor argument.
    kept = [i for i, row in enumerate(rows) if row[("Q", 0)][0] == 0]
    assert len(kept) == 4
    for i, row in enumerate(rows):
        for j, col in enumerate(cols):
            weighted = int(row[("Q", 0)][0] == 0) * source.boundary(u, v, row, col)
            assert weighted == int(i == j and i in kept)
            totals["local_predicate_negative_entries"] += 1
    # Zeroing one operand is essential: a nonzero first middle vector can
    # create a nonzero outgoing message and invalidate a zero IB state.
    left, right = (1, 0), (1, 0)
    assert source.middle(0, left, right, 0) == source.a
    assert source.middle(0, (0, 0), right, 0) == 0
    totals["nonzero_first_operand_negative_controls"] += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    totals, records = Counter(), []
    check_units(totals)
    check_minors(totals, records)
    check_tradeoff(totals)
    negative_controls(totals)
    result = {"status": "PASS", "source_revision": "c91cc3f4947caca8977e81a5d14244f73b5ae2b1",
              "scope": "raw group-pair edge constraints; local weights and internal constraints excluded",
              "factoring_score_improved": False, "counts": dict(sorted(totals.items())),
              "minor_records": records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "counts": result["counts"]}, indent=2))


if __name__ == "__main__":
    main()
