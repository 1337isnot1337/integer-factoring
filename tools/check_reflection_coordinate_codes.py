#!/usr/bin/env python3
"""Exact all-scale pressure floor for coordinate codes of two reflection sites.

--full covers every common coordinate subset of the sixteen physical words,
modulo exchanging the two sites.  No floating optimization enters the proof.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import time

import check_block_smith_crossings as smith
from check_cancellation_module_automata import pmul
import check_triangular_auxiliary_row_fiber as reflection


ROOT = Path(__file__).resolve().parent.parent
HEIGHT = tuple(reflection.HEIGHT[w & 3]+reflection.HEIGHT[w >> 2]
               for w in range(16))
PARITY = tuple(((w & 3).bit_count() & 1)+2*((w >> 2).bit_count() & 1)
               for w in range(16))
SWAP = tuple((w >> 2)+4*(w & 3) for w in range(16))


def ambient_tensor():
    original = dict(zip(reflection.LABELS, reflection.reflection_seed()))
    zero, infinity = {}, {}
    for (x, y, u, v), value in original.items():
        shift = reflection.HEIGHT[x]+reflection.HEIGHT[y]-4
        if shift < 0:
            assert value & ((1 << -shift)-1) == 0
            zero[x, y, u, v] = value >> -shift
        else:
            zero[x, y, u, v] = value << shift
        degree = reflection.HEIGHT[u]+reflection.HEIGHT[v]+4
        infinity[x, y, u, v] = smith.reverse_polynomial(value, degree)
    for x, y, u, v in original:
        assert infinity[x, y, u, v] == zero[x, u, y, v]

    tensor = [0]*(16**4)
    for first, second in product(reflection.LABELS, repeat=2):
        labels = tuple(a+4*b for a, b in zip(first, second))
        x, y, u, v = labels
        value = pmul(zero[first], zero[second])
        tensor[x+16*y+256*u+4096*v] = value
        assert value.bit_count() % 2 == int(x == v and y == u)
        if value:
            assert PARITY[x] ^ PARITY[y] ^ PARITY[u] ^ PARITY[v] == 0
    return tensor


def code_matrix(tensor, code):
    return [[tensor[x+16*y+256*u+4096*v] for u, v in product(code, repeat=2)]
            for x, y in product(code, repeat=2)]


def code_profile(tensor, code):
    groups = {}
    for x, y in product(code, repeat=2):
        groups.setdefault(PARITY[x] ^ PARITY[y], []).append((x, y))
    result = []
    for pairs in groups.values():
        matrix = [[tensor[x+16*y+256*u+4096*v] for u, v in pairs]
                  for x, y in pairs]
        degree = max(value.bit_length()-1 for row in matrix for value in row)
        orders, _ = smith.smith_profile(matrix, degree)
        result.extend(orders)
    assert len(result) == len(code)**2
    return tuple(sorted(result))


def cumulative_certificate(code, orders):
    """Prove g(x^5)^7 - (sum x^(7h))^10 >= 0 for 0<=x<=1."""
    q = len(code)
    counts, heights = Counter(orders), Counter(HEIGHT[w] for w in code)
    minimum_order, minimum_height = min(counts), min(heights)
    # A selected crossing diagonal supplies this simple entry valuation.
    assert minimum_order <= 2*minimum_height
    common = min(35*minimum_order, 70*minimum_height)
    # Nonnegative convolution coefficients have total at most q^14.
    # This slot width rules out every carry in Kronecker substitution.
    bits = (q**14).bit_length()+1
    mask = (1 << bits)-1
    left = sum(value << (bits*5*(order-minimum_order))
               for order, value in counts.items())**7
    right = sum(value << (bits*7*(height-minimum_height))
                for height, value in heights.items())**10
    left <<= bits*(35*minimum_order-common)
    right <<= bits*(70*minimum_height-common)
    degree = max(35*max(orders), 70*max(heights))-common
    first = second = 0
    first_min = second_min = 0
    for _ in range(degree+1):
        coefficient = (left & mask)-(right & mask)
        left >>= bits
        right >>= bits
        first += coefficient
        second += first
        first_min = min(first_min, first)
        second_min = min(second_min, second)
    assert not left and not right
    assert first == q**14-q**10 > 0
    assert second >= 0
    assert second_min >= 0, (code, orders, first_min, second_min)
    return {
        "cumulative_depth": 1 if first_min >= 0 else 2,
        "polynomial_degree_after_monomial_removal": degree,
        "removed_monomial": common,
        "first_terminal_sum": first,
        "second_terminal_sum": second,
    }


def swapped_mask(mask):
    return sum(1 << SWAP[w] for w in range(16) if mask >> w & 1)


def control_cases(tensor):
    codes = ((0, 1), (0, 1, 2), (0, 1, 2, 3),
             (0, 1, 2, 3, 12), tuple(range(16)))
    answer = []
    for code in codes:
        orders = code_profile(tensor, code)
        certificate = cumulative_certificate(code, orders)
        matrix = code_matrix(tensor, code)
        for precision in range(1, 5):
            assert smith.toeplitz_rank(matrix, precision) == sum(
                max(0, precision-order) for order in orders)
        answer.append({"code": code, "profile": sorted(Counter(orders).items()),
                       "heights": sorted(Counter(HEIGHT[w] for w in code).items()),
                       **certificate})
    full = Counter({0: 1})
    for value in (1, 3, 5, 7)*2:
        nxt = full.copy()
        for order, multiplicity in full.items():
            nxt[order+value] += multiplicity
        full = nxt
    assert answer[-1]["profile"] == sorted(full.items())
    return answer


def run(full, progress):
    start = time.monotonic()
    tensor = ambient_tensor()
    controls = control_cases(tensor)
    result = {
        "scope": "two-site reflection common coordinate codes with inherited natural gauges",
        "field": "F2", "full_enumeration": full,
        "pressure_lower_bound": "5/7",
        "strictly_above_one_over_sqrt_two": 2*5**2-7**2,
        "literal_local_duality_entries": 4**4,
        "literal_two_site_endpoint_entries": 16**4,
        "independent_toeplitz_ranks": 4*len(controls),
        "controls": controls,
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }
    if full:
        orbit_counts, covered_counts, depth_counts = Counter(), Counter(), Counter()
        digest = sha256()
        maximum_degree = 0
        for mask in range(1 << 16):
            q = mask.bit_count()
            if q < 2:
                continue
            swapped = swapped_mask(mask)
            if mask > swapped:
                continue
            code = tuple(w for w in range(16) if mask >> w & 1)
            orders = code_profile(tensor, code)
            certificate = cumulative_certificate(code, orders)
            multiplicity = 1 if mask == swapped else 2
            orbit_counts[q] += 1
            covered_counts[q] += multiplicity
            depth_counts[certificate["cumulative_depth"]] += 1
            maximum_degree = max(maximum_degree,
                                 certificate["polynomial_degree_after_monomial_removal"])
            digest.update(json.dumps([mask, sorted(Counter(orders).items()), certificate],
                                     sort_keys=True, separators=(",", ":")).encode()+b"\n")
            if progress and sum(orbit_counts.values()) % 1000 == 0:
                print(json.dumps({"processed_orbits": sum(orbit_counts.values()),
                                  "seconds": round(time.monotonic()-start, 2)}), flush=True)
        assert covered_counts == Counter({q: comb(16, q) for q in range(2, 17)})
        assert sum(orbit_counts.values()) == 33269
        assert sum(covered_counts.values()) == 65519
        result.update({
            "orbit_representatives_by_alphabet": dict(sorted(orbit_counts.items())),
            "coordinate_codes_by_alphabet": dict(sorted(covered_counts.items())),
            "total_orbit_representatives": sum(orbit_counts.values()),
            "total_coordinate_codes": sum(covered_counts.values()),
            "certificate_depth_counts": dict(sorted(depth_counts.items())),
            "maximum_certificate_polynomial_degree": maximum_degree,
            "ordered_profile_and_certificate_sha256": digest.hexdigest(),
        })
    result["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    result["seed_checker_sha256"] = sha256(Path(reflection.__file__).read_bytes()).hexdigest()
    result["seconds"] = round(time.monotonic()-start, 3)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="enumerate all 33,269 site-exchange orbits")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--output", type=Path,
                        default=ROOT/"evidence"/"reflection-coordinate-codes.json")
    args = parser.parse_args()
    data = run(args.full, args.progress)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2)+"\n")
    print(json.dumps(data, indent=2))
