#!/usr/bin/env python3
"""Exact character-lift guards on previously certified Path C data.

This postprocessing performs no class compositions and enumerates no new
arithmetic cycle. Its prerequisite checks the original group certificates.
"""
import argparse
from hashlib import sha256
import json
import lzma
from math import gcd, isqrt
from pathlib import Path

from check_cofactor_orbit_moments import maximum_zero_run

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "evidence/clustered-power-relations"
EXPECTED = ROOT / "evidence/power-character-lifts.json"


def smooth_support(primes, height):
    result = []

    def visit(i, norm, exponents):
        if i == len(primes):
            result.append((norm, tuple(exponents)))
            return
        exponent = 0
        while norm <= height:
            visit(i + 1, norm, exponents + [exponent])
            norm *= primes[i]
            exponent += 1

    visit(0, 1, [])
    return result


def compute():
    rows = json.loads((DATA / "row-audit.json").read_text())
    pool = json.loads((DATA / "pool-audit.json").read_text())
    orbit = json.loads((DATA / "receipt.json").read_text())
    primes = rows["primes"]
    h = pool["period"]
    indices = [pool["derived_generator_indices"][str(p)] for p in primes]
    assert primes == [5, 7, 19] and indices[0] == 1
    height = isqrt(pool["N"] // 3)
    support = smooth_support(primes, height)
    certificate = next(c for c in orbit["certificates"] if c["N"] == pool["N"])
    base_flags = lzma.decompress((ROOT / certificate["flags_file"]).read_bytes())
    assert len(base_flags) == h
    results = []
    for case in rows["results"]:
        p = case["prime"]
        first, second = case["selected_rows"][:2]
        z = [first[1] * second[2] - first[2] * second[1],
             first[2] * second[0] - first[0] * second[2],
             first[0] * second[1] - first[1] * second[0]]
        content = gcd(*z)
        assert content
        z = [a // content for a in z]
        assert all((a - z[0] * j) % h == 0 for a, j in zip(z, indices))
        assert gcd(z[0], h) == 1
        budget, witness = max((sum(e * abs(a) for e, a in zip(es, z)), (n, es))
                              for n, es in support)
        coordinate = z[primes.index(p)]
        threshold = (h - budget + abs(coordinate) - 1) // abs(coordinate)
        flags = bytes(base_flags[(k * indices[primes.index(p)]) % h] for k in range(h))
        allowed = bytes(int(min((k * coordinate) % h, (-k * coordinate) % h) <= budget)
                        for k in range(h))
        assert all(not good or admissible for good, admissible in zip(flags, allowed))
        outside = []
        forced_rows = 0
        for row in case["rows"]:
            value = sum(a * b for a, b in zip(row["row"], z))
            assert value % h == 0
            if row["k"] * abs(coordinate) + budget < h:
                assert value == 0
                forced_rows += 1
            if value:
                outside.append((row["k"], value, row["row"]))
        assert outside and outside[0][0] >= threshold
        assert outside[0][0] == case["first_full_rank_exponent"]
        # Every base phase lies in the open +/-1/6 arc, so no distribution
        # on these three generators or their inverses has bias <=1/2 here.
        assert max(6 * abs(a) for a in z) < h
        results.append(dict(prime=p, normal=z, character_multiplier=z[0],
                            exact_word_budget=budget,
                            budget_witness=dict(norm=witness[0], exponents=list(witness[1])),
                            first_possible_escape_exponent=threshold,
                            first_actual_escape_exponent=outside[0][0],
                            escape_dot=outside[0][1], escape_row=outside[0][2],
                            forced_hyperplane_rows=forced_rows,
                            forced_rough_phases=h-sum(allowed),
                            longest_forced_rough_run=maximum_zero_run(allowed),
                            base_admits_bias_at_most_half=False))
    inputs = {"receipt.json":orbit, "row-audit.json":rows, "pool-audit.json":pool}
    semantic_hashes = {
        name:sha256(json.dumps({k:v for k,v in data.items() if k != "seconds"},
                               sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        for name, data in inputs.items()
    }
    return dict(scope="Integer postprocessing of certified saved cycles; no new class compositions, orbit or cutoff search.",
                N=pool["N"], primes=primes, period=h, height=height,
                generator_indices=indices, unsigned_norm_support_size=len(support),
                input_semantic_sha256=semantic_hashes,
                input_hash_scope="Canonical JSON with only the top-level runtime seconds omitted.",
                results=results, complete_factoring_exponent="1/2; unchanged")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    actual = compute()
    if args.output:
        args.output.write_text(json.dumps(actual, indent=2) + "\n")
    else:
        assert actual == json.loads(EXPECTED.read_text())
    print(json.dumps({"unsigned_norm_support_size":actual["unsigned_norm_support_size"],
                      "guards":[{k:r[k] for k in ("prime", "first_possible_escape_exponent",
                                                  "first_actual_escape_exponent", "longest_forced_rough_run")}
                                for r in actual["results"]],
                      "complete_factoring_exponent":"1/2; unchanged"}))


if __name__ == "__main__":
    main()
