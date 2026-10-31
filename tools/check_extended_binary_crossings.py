#!/usr/bin/env python3
"""Replay complete rank-one-constant quadratics and parity-preserving cubics."""
import argparse
from collections import Counter
from itertools import product
import json

from check_quadratic_crossings import ROOT, certify, matrices, orders


def multiply(a, b):
    return sum((sum((a >> (2*i+k) & 1) * (b >> (2*k+j) & 1)
                    for k in range(2)) % 2) << (2*i+j)
               for i, j in product(range(2), repeat=2))


def constant_orbits():
    nonzero = range(1, 16)
    inverse = {a: b for a, b in product(nonzero, repeat=2)
               if multiply(a, b) == multiply(b, a) == 9}
    unseen, result = set(product(nonzero, repeat=2)), []
    while unseen:
        a, b = min(unseen)
        orbit = {(multiply(multiply(s, a), inverse[t]),
                  multiply(multiply(t, b), inverse[s])) for s, t in product(inverse, repeat=2)}
        assert orbit <= unseen
        unseen -= orbit
        result.append((a, b, len(orbit)))
    assert len(result) == 11 and sum(size for _, _, size in result) == 225
    return result


def parity_mask(index, degree=3):
    positions = [entry for entry in range(16) if entry.bit_count() % 2 == 0]
    return sum((index >> (8*j+bit) & 1) << (16*j+entry)
               for bit, entry in enumerate(positions) for j in range(degree-1))


def replay(data, full=False, schedule=None, sample_count=256):
    degree, constants = data["degree"], tuple(data.get("constant_factors", (9, 9)))
    rows = data["profiles"]
    expected = {(tuple(row["infinity"]), tuple(row["zero"])): row for row in rows}
    witnesses = {row["first_mask"]: key for key, row in expected.items()}
    assert len(expected) == len(rows)
    assert len(witnesses) == len(rows)
    assert sum(row["count"] for row in rows) == data["tensors"]
    masks = schedule if full else sorted({row["first_mask"] for row in rows} | set(schedule[:sample_count]))
    counts, first = Counter(), {}
    for mask in masks:
        key = tuple(orders(matrix) for matrix in matrices(mask, degree, constants))
        assert key in expected, (mask, key)
        counts[key] += 1
        first.setdefault(key, mask)
        if mask in witnesses:
            assert key == witnesses[mask]
    if full:
        assert all(counts[key] == row["count"] and first[key] == row["first_mask"]
                   for key, row in expected.items())
    certificates = certify(rows, degree)
    return {"degree": degree, "constant_factors": constants, "catalogued_tensors": data["tensors"],
            "independently_checked_tensors": len(masks), "paired_profiles": len(rows),
            "maximal_profiles": len(certificates), "pressure_floor": "256/27",
            "bernstein_leaves": sum(c["leaves"] for c in certificates),
            "max_bernstein_depth": max(c["max_depth"] for c in certificates)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="independently replay all 786,432 structured tensors")
    args = parser.parse_args()
    result = {"constant_orbits": constant_orbits(), "structured": []}
    data = json.loads((ROOT / "evidence/rank-one-constant-quadratic-catalogues.json").read_text())
    assert [tuple(row["constant_factors"]) for row in data] == [(a, b) for a, b, _ in constant_orbits()]
    for row in data:
        assert row["degree"] == 2 and row["tensors"] == 65536
        result["structured"].append(replay(row, args.full, range(65536)))
    data = json.loads((ROOT / "evidence/parity-cubic-crossing-catalogue.json").read_text())
    assert data["degree"] == 3 and data["tensors"] == 65536
    result["structured"].append(replay(data, args.full, [parity_mask(i) for i in range(65536)]))
    # The unstructured cubic run is a sample, separately marked and independently
    # checked only at every saved representative and the first 256 scheduled masks.
    data = json.loads((ROOT / "evidence/cubic-crossing-million-sample.json").read_text())
    assert data["degree"] == 3 and data["tensors"] == 1000000
    assert data["mask_schedule"] == "(i*0x9e3779b1+0x7f4a7c15) mod 2^32"
    schedule = [(i*0x9e3779b1+0x7f4a7c15) % 2**32 for i in range(256)]
    result["unstructured_cubic_sample"] = replay(data, False, schedule)
    data = json.loads((ROOT / "evidence/parity-quartic-crossing-catalogue.json").read_text())
    assert data["degree"] == 4 and data["tensors"] == 2**24
    schedule = [parity_mask((i*0x9e3779+0x7f4a7c) % 2**24, 4) for i in range(4096)]
    result["quartic_catalogue_sampled_replay"] = replay(data, False, schedule, 4096)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
