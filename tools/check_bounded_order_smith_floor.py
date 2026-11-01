#!/usr/bin/env python3
"""Exact entropy floor when both bare Smith profiles are bounded by degree."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def digest_integer(value):
    return sha256(value.to_bytes((value.bit_length()+7)//8, "big")).hexdigest()


def entropy_certificate(q):
    n, denominator = q*q, 20*q*(q-1)
    counts = [(7*q+6)*(q-1)]+[6*(q-1)]*(q-1)+[7]*(n-q)
    assert len(counts) == n and sum(counts) == denominator
    assert all(value > 0 for value in counts)
    prefix_weight = Fraction(3, 10*q)
    tail_weight = Fraction(7, 20*q*(q-1))
    assert all(Fraction(value, denominator) == prefix_weight for value in counts[1:q])
    assert all(Fraction(value, denominator) == tail_weight for value in counts[q:])
    assert prefix_weight*Fraction(q, 2)+tail_weight*(n-q) == Fraction(1, 2)
    common = gcd(denominator, 625)
    multiplier = 625//common
    left = denominator**(denominator*multiplier)
    right = q**(884*denominator//common)
    for value, multiplicity in Counter(counts).items():
        right *= value**(value*multiplicity*multiplier)
    assert left > right
    return {
        "alphabet": q, "probability_denominator": denominator,
        "probability_counts": sorted(Counter(counts).items(), reverse=True),
        "prefix_weight": str(prefix_weight), "tail_weight": str(tail_weight),
        "weighted_order_budget": "1/2",
        "integer_comparison_multiplier": multiplier,
        "left_integer_bits": left.bit_length(), "right_integer_bits": right.bit_length(),
        "left_sha256": digest_integer(left), "right_sha256": digest_integer(right),
        "strict_entropy_inequality": True, "rho_lower_bound": "442/625",
    }


def eleven_letter_formal_witness():
    q, degree = 11, 20
    profile = (0,)+(11,)*10+(20,)*110
    assert len(profile) == q*q and max(profile) <= degree and profile[0] == 0
    assert sum(profile) <= q*q*degree
    for selected in range(2, q+1):
        assert 2*sum(profile[:selected]) <= selected*degree
    polynomial = 5**20+10*4**11*5**9+110*4**20
    left, right = polynomial**500, 11**707*5**5000*4**5000
    assert left < right
    assert 2*707**2 < 1000**2
    return {
        "alphabet": q, "degree": degree,
        "both_profiles": sorted(Counter(profile).items()),
        "scale": "5/4", "strict_formal_rho_upper_bound": "707/1000",
        "left_integer_bits": left.bit_length(), "right_integer_bits": right.bit_length(),
        "left_sha256": digest_integer(left), "right_sha256": digest_integer(right),
        "crossing_realization_claimed": False,
    }


def run():
    assert 2*442**2 > 625**2
    return {
        "scope": "all fields and degrees; both bare Smith profiles must have orders at most actual degree",
        "entropy_certificates": [entropy_certificate(q) for q in range(2, 11)],
        "strictly_above_one_over_sqrt_two": 2*442**2-625**2,
        "eleven_letter_formal_relaxation_witness": eleven_letter_formal_witness(),
        "bell_consequence": "eight-letter polynomial order/degree pressure only, every characteristic different from two",
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT/"evidence"/"bounded-order-smith-floor.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
