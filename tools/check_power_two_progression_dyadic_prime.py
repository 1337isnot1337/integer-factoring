#!/usr/bin/env python3
"""Exact constants and integer-schedule controls for the BMOR prime lemma.

The prime existence claim uses the cited analytic theorem. This checker
does not substitute sampled primes for that theorem or search its enormous
explicit lower cutoff.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path


def schedule(n):
    u = (n + 1).bit_length()
    p = isqrt(n) // (16 * u * u)
    a = 1 << (p * u - 1).bit_length()
    total = a * a
    denominator = total - a - 1
    w = 2 * ((n + 1 + 2 * denominator - 1) // (2 * denominator))
    d = total // 2
    t = d.bit_length() - 1
    return u, p, a, total, w, 2 * w + t


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    r1 = Fraction(9645908801, 10**9)
    assert 4 * r1 < 40
    assert 32 * 14579 < 500000
    assert 32 * 1001 < 100000000
    assert 2**39 > (128 * 41**2 + 1)**2
    assert Fraction(42, 41)**4 < 2

    # Every nonprincipal real character modulo a power of two is one of
    # these three sign patterns. Check their complete multiplicativity
    # on representative powers of two.
    patterns = ([1, -1, 1, -1], [1, -1, -1, 1], [1, 1, -1, -1])
    character_products = 0
    for modulus in (8, 16, 32, 64, 128):
        units = list(range(1, modulus, 2))
        for pattern in patterns:
            char = lambda value: pattern[(value % 8) // 2]
            for first in units:
                for second in units:
                    assert char(first * second) == char(first) * char(second)
                    character_products += 1

    # These rational checks supplement the uniform monotonicity and
    # strict-convexity proof for all real positive exponents.
    positive_blocks = 0
    for sigma in range(1, 7):
        for block in range(20):
            a = 8 * block
            values = [Fraction(1, (a + j)**sigma) for j in (1, 3, 5, 7)]
            for pattern in patterns:
                assert sum(sign * value for sign, value in zip(pattern, values)) > 0
                positive_blocks += 1

    cases = sorted({2**exponent + offset
                    for exponent in range(40, 513)
                    for offset in (-2, -1, 0, 1, 2)
                    if 2**exponent + offset >= 2**40})
    examples = []
    for n in cases:
        u, p, a, total, w, s = schedule(n)
        assert u >= 41 and p >= 8
        assert total >= 100000
        assert a >= p * u and a < 2 * p * u
        assert 64 * u * u * total < n
        assert w > 64 * u * u and s > 128 * u * u
        assert total < n < 2**u
        assert total // 2 == 1 << ((total // 2).bit_length() - 1)
        if n in (2**40, 2**64, 2**128, 2**256, 2**512):
            examples.append({"n": str(n), "u": u, "P": str(p), "A": str(a),
                             "T": str(total), "w": str(w), "s": str(s)})
    report = {
        "status": "all exact controls passed",
        "scope": "rational analytic-error margins, real-character patterns, and the integer schedule; prime existence uses BMOR Lemmas 6.10--6.12",
        "primary_source": "https://personal.math.ubc.ca/~bennett/BeMaObRe-IJM-2018.pdf",
        "rational_constant_controls": 5,
        "complete_character_products": character_products,
        "positive_rational_character_blocks": positive_blocks,
        "integer_schedule_cases": len(cases),
        "schedule_examples": examples,
        "source_sha256": {Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()},
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ("source_sha256", "schedule_examples")}))


if __name__ == "__main__":
    main()
