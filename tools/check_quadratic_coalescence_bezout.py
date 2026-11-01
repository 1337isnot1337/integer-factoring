#!/usr/bin/env python3
"""Exact controls for the two-depth reflected primitive Bezout testbed."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from check_balanced_bezout_modular import literal_bezout, residue
from matched_progression_bezout import reflected_root_incidence
from quadratic_coalescence_bezout import quadratic_binomial_trace, quadratic_coalescence_bezout

ROOT = Path(__file__).resolve().parents[1]


def controls():
    traces = pairs = coefficients = moments = incidence = 0
    for t in range(1, 8):
        m, n = 1 << t, (1 << t) - 1
        for k in range(min(n, 20) + 1):
            direct = sum(comb(n, j) ** 2 * comb(j, k) for j in range(k, m))
            assert direct == comb(n, k) * comb(2 * n - k, n - k)
            order = (direct & -direct).bit_length() - 1
            assert order == t - ((k + 1) & -(k + 1)).bit_length() + 1
            moments += 1
        for w in (-7, -3, 1, 5, 9, 13):
            exact = Fraction(4, m) * sum(
                (Fraction(comb(n, j) ** 2, j * j + w) for j in range(m)), Fraction())
            assert exact.numerator % 2 and exact.denominator % 2
            for p in (1, 4, 9):
                result = quadratic_binomial_trace(t, w, p)
                assert result["coefficient"] == residue(exact, p)
                traces += 1
    for t in range(1, 5):
        m = 1 << t
        for h in (2 * t - 1, 2 * t + 2):
            step, a = 1 << h, (1 << h) * m + 1
            aa, bb = [a + step * i for i in range(m)], [a - step * j for j in range(m)]
            for w in (-3, 1, 5):
                for p in sorted(set((1, h - 2 * t + 2))):
                    result = quadratic_coalescence_bezout(t, h, a, w, p)
                    u, v = literal_bezout(result["n"], aa, bb, result["bezout_exponent"])
                    for key, poly in (("u_terms", u), ("v_terms", v)):
                        moved = [sum(poly[j] * comb(j, i) * (a - 1) ** (j - i)
                                     for j in range(i, m)) for i in range(m)]
                        expected = [0] * (m - 1) + [result[key][0][1]]
                        assert [residue(value, p) for value in moved] == expected
                        assert moved[-1].numerator % 2 and moved[-1].denominator % 2
                        coefficients += m
                    pairs += 1
        h, step = 2 * t + 1, 1 << (2 * t + 1)
        a = step * m + 1
        for w in sorted(set([1, 5] + [-i * i for i in range(m)] + [1 - i * i for i in range(m)])):
            n = a * a + step * step * w
            if n <= 0:
                continue
            result = reflected_root_incidence(m, h, a, (step // m) * w)
            matches = [(i, j) for i in range(m) for j in range(m)
                       if (a + step * i) * (a - step * j) == n]
            assert len(matches) == result["gcd_degree"]
            assert not matches or matches == [(result["matched_index"], result["matched_index"])]
            incidence += 1
    invalid = 0
    for args in ((1, 1, 5, 1, 2), (3, 4, 257, 1, 1), (2, 5, 129, 3, 1),
                 (2, 5, 95, 1, 1), (2, 5, 129, -31, 1)):
        try:
            quadratic_coalescence_bezout(*args)
        except ValueError:
            invalid += 1
        else:
            raise AssertionError(("An invalid parameter contract was accepted", args))
    return {"exact_trace_residues": traces, "exact_weighted_moment_identities": moments,
            "complete_rational_pairs": pairs, "complete_coefficients": coefficients,
            "factor_and_no_incidence_controls": incidence, "invalid_contracts": invalid}


def implicit_controls():
    results = []
    for t, p, w in ((200, 12, 1), (500, 8, (1 << 601) + 1), (1000, 5, -3)):
        m, h = 1 << t, 2 * t + p + 3
        a = (1 << h) * m + 1
        started = perf_counter()
        lower = quadratic_coalescence_bezout(t, h, a, w, p)
        seconds = perf_counter() - started
        higher = quadratic_coalescence_bezout(t, h, a, w, p + 3)
        assert lower["bezout_exponent"] == higher["bezout_exponent"]
        for key in ("u_terms", "v_terms"):
            assert lower[key][0][0] == higher[key][0][0]
            assert lower[key][0][1] == higher[key][0][1] % (1 << p)
        assert lower["trace"]["coefficient"] == higher["trace"]["coefficient"] % (1 << p)
        assert lower["enumerated_domain_nodes"] == higher["enumerated_domain_nodes"] == 0
        results.append({"lower_seconds": seconds, "lower": lower, "higher": higher})
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "A complete primitive pair in a two-depth reflected testbed",
              "controls": controls(), "factoring_exponent_improved": False}
    if not args.small_only:
        report["implicit_controls"] = implicit_controls()
    report["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "quadratic_coalescence_bezout.py", "interval_unit_jets.py",
                     "pole_cauchy_conditioning.py", "matched_progression_bezout.py")}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["controls"], indent=2))
    if not args.small_only:
        print(json.dumps([{"dimension_bits": x["lower"]["count"].bit_length(),
                           "N_bits": x["lower"]["n"].bit_length(),
                           "working_precision": x["lower"]["trace"]["working_precision"],
                           "mahler_terms": x["lower"]["trace"]["mahler_terms"],
                           "lower_seconds": x["lower_seconds"]}
                          for x in report["implicit_controls"]], indent=2))


if __name__ == "__main__":
    main()
