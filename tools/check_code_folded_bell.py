#!/usr/bin/env python3
"""Exact controls for the code-quotient duality and Frobenius no-gain example."""

import argparse
from fractions import Fraction
import json

from code_folded_bell import code_path, multiply, nullspace, span
from folded_bell import polynomials
from math import comb


def rational_projection(n, generators):
    c = span(generators)
    even_c = tuple(x for x in c if x.bit_count() % 2 == 0)
    v = span(nullspace(generators, n))
    h = max(x.bit_count() for x in v) // 2
    k = (n - 1) // 2
    f = polynomials(k)
    spectra = []
    for y in v:
        j = y.bit_count() // 2
        spectra.append(multiply(tuple(comb(h - j, i) for i in range(h - j + 1)),
                                tuple((-1) ** i * comb(j, i) for i in range(j + 1))))
    common = tuple(Fraction(comb(k - h, i), 2 ** (k - h)) for i in range(k - h + 1))
    for x in v:
        averaged = [Fraction(0)] * (k + 1)
        for y in even_c:
            w = (x ^ y).bit_count()
            ell = min(w, n - w)
            for i, value in enumerate(f[k - ell]):
                averaged[i + ell] += Fraction(value, len(even_c))
        normalized = tuple(Fraction(sum((-1) ** ((x & y).bit_count() % 2) * row[i]
                                        for y, row in zip(v, spectra)), 2 ** h)
                           for i in range(h + 1))
        assert multiply(common, normalized) == tuple(averaged)
        order = next(i for i, a in enumerate(normalized) if a)
        assert order == min((x ^ y).bit_count() for y in c)
    return len(v) * (k + 1)


def check() -> dict:
    cases = (
        ("antipodal7", 7, (127,)),
        ("three_odd_blocks", 9, (7, 56, 448)),
        ("one_tripled_coordinate", 9, (511, 3, 6)),
        ("mixed_lcd", 9, (511, 15, 24)),
        ("full_tripling", 9, (511, 3, 6, 24, 48, 192, 384)),
    )
    results = []
    rational_checks = 0
    for name, n, generators in cases:
        rational_checks += rational_projection(n, generators)
        for prime in (3, 5, 7):
            result = code_path(n, generators, prime)
            result["case"] = name
            # Preserve compact summary receipts; code+parameters reproduce
            # the full coefficients and exact dual identity internally.
            result.pop("coefficients_by_degree")
            results.append(result)
    triple = code_path(9, cases[-1][2], 3)
    base = code_path(3, (7,), 3)
    assert triple["degree_after_maximal_common_division"] == 3 * base["degree_after_maximal_common_division"]
    assert tuple(triple["valuation_histogram"][::3]) == tuple(base["valuation_histogram"])
    assert all(x == 0 for i, x in enumerate(triple["valuation_histogram"]) if i % 3)
    # Each code has 2^m distinct cosets and no field-dependent physical
    # dimension is silently discarded.
    for result in results:
        assert sum(result["valuation_histogram"]) == 1 << result["dual_dimension"]
        assert sum(result["distance_histogram"]) == 1 << result["dual_dimension"]
        m, h = result["dual_dimension"], result["degree_after_maximal_common_division"]
        for x in (Fraction(1), Fraction(3, 2), Fraction(3), Fraction(9)):
            a = sum(Fraction(count) / x ** i for i, count in enumerate(result["distance_histogram"]))
            assert a >= (1 + 1 / x) ** m
            assert x ** h * a ** 2 >= x ** (m // 2) * (1 + 1 / x) ** (2 * m)
    return {"status": "passed", "exact_quotient_cases": len(results),
            "rational_projection_coefficients": rational_checks,
            "full_tripling": "orders and degree scale by 3; pressure optimum unchanged",
            "cases": results,
            "factoring_score": "1/2; unchanged"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check()
    if args.output:
        with open(args.output, "w") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps({"status": result["status"], "exact_quotient_cases": result["exact_quotient_cases"],
                      "cases": [{"case": x["case"], "prime": x["prime"],
                                 "rho_pressure": x["rho_pressure_approx"],
                                 "raised_cosets": x["raised_valuation_cosets"]}
                                for x in result["cases"]]}, indent=2))


if __name__ == "__main__":
    main()
