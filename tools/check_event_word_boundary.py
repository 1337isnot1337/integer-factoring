#!/usr/bin/env python3
"""Read-only saved-word controls for the exact p-boundary event formula."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import prod
from pathlib import Path

from binary_quadratic_forms import compose_forms
from check_relation_change_events import event_flags

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "evidence/event-word-boundary.json"


def norm(primes, word):
    return prod(p ** abs(e) for p, e in zip(primes, word))


def contributed_euler_control():
    """Discussion #51 comments 19356441/19356480; retain the supplied panel.

    Adapted from /tmp/factoring-boundary-euler-control.py, source SHA256
    37fef094a7110f3c78f2b60f4f8baf72c415c821e1e23b8637037f0cece256e9.
    """
    primes, p = (3, 5, 7), 3
    words = [(v, norm(primes, v)) for v in product(range(-6, 7), repeat=3)
             if norm(primes, v) <= 200]

    def counters(radius):
        raw = Counter(v for v, n in words if n <= radius)
        boundary = Counter()
        for u, n in words:
            if u[0] or n > radius:
                continue
            a = 0
            while n * p ** (a + 1) <= radius:
                a += 1
            boundary[(-a, u[1], u[2])] += 1
        return raw, boundary

    def image(v):
        return ((v[0] + v[1] + 2 * v[2]) % 3, (v[1] + 2 * v[2]) % 5)

    for radius in range(1, 201):
        raw, boundary = counters(radius)
        small_raw, small_boundary = counters(radius // p)
        left = boundary + small_boundary
        right = raw.copy()
        for v, count in small_raw.items():
            right[(v[0] + 1, v[1], v[2])] -= count
        assert all(left[v] == right[v] for v in set(left) | set(right))
        projected_left, projected_right = Counter(), Counter()
        for v, count in left.items():
            projected_left[image(v)] += count
        for v, count in right.items():
            projected_right[image(v)] += count
        assert projected_left == projected_right
        expected = sum(1 for u, n in words if u[0] == 0 and n <= radius and image(u)[1] == 0)
        assert sum(count for v, count in boundary.items() if image(v)[1] == 0) == expected
    return dict(primes=list(primes), entry_prime=p, heights=200, exact_checks=600,
                words=len(words), group_orders=[3, 5],
                original_source_sha256="37fef094a7110f3c78f2b60f4f8baf72c415c821e1e23b8637037f0cece256e9",
                source_discussion_comments=[19356441, 19356480])


def compute():
    data = json.loads((ROOT / "evidence/clustered-power-relations/row-audit.json").read_text())
    fixed = json.loads((ROOT / "evidence/small-norm-event-separation.json").read_text())
    primes = data["primes"]
    results = []
    for case in data["results"]:
        p, h, k0 = case["prime"], case["period"], case["K"]
        i = primes.index(p)
        radius = next(c["radius"] for c in fixed["controls"] if c["prime"] == p)
        words = {k: tuple(k if j == i else 0 for j in range(len(primes))) for k in range(k0)}
        for row in case["rows"]:
            k = row["k"]
            words[k] = tuple((k if j == i else 0) - r for j, r in enumerate(row["row"]))
        low = {k: v for k, v in words.items() if norm(primes, v) <= radius}
        hits = bytes(int(k in low) for k in range(h))
        events = event_flags(hits, low, i)
        boundary = bytes(int(k in low and low[k][i] <= 0 and p * norm(primes, low[k]) > radius)
                         for k in range(h))
        assert boundary == events
        pfree = {k: v for k, v in low.items() if v[i] == 0}
        assert len(pfree) == sum(events)
        covered, starts, fibers = set(), set(), []
        for phase, u in pfree.items():
            n = norm(primes, u)
            a = 0
            while n * p ** (a + 1) <= radius:
                a += 1
            event_phase = (phase - a) % h
            assert events[event_phase] and event_phase not in starts
            starts.add(event_phase)
            for z in range(-a, a + 1):
                k = (phase + z) % h
                expected = tuple(e + (z if j == i else 0) for j, e in enumerate(u))
                assert low[k] == expected and k not in covered
                covered.add(k)
            fibers.append((phase, n, a))
        assert covered == set(low)
        assert starts == {k for k, flag in enumerate(events) if flag}
        rows = {tuple((k if j == i else 0) - e for j, e in enumerate(low[k])) for k in starts}
        assert len(rows) == len(starts)
        fourier = []
        for modulus in (1, 2, 5, 10):
            assert h % modulus == 0
            direct = [sum(events[k] for k in range(r, h, modulus)) for r in range(modulus)]
            annular = [0] * modulus
            largest_a = max(a for _, _, a in fibers)
            for a in range(largest_a + 1):
                upper, lower = radius // p ** a, radius // p ** (a + 1)
                for phase, n, _ in fibers:
                    if lower < n <= upper:
                        annular[(phase - a) % modulus] += 1
            assert direct == annular
            fourier.append(dict(character_order=modulus, residue_counts=direct))
        results.append(dict(prime=p, period=h, radius=radius, hits=len(low),
                            p_free_words=len(pfree), events=sum(events),
                            disjoint_fibers=len(fibers), distinct_integer_event_rows=len(rows),
                            event_flags_sha256=sha256(events).hexdigest(),
                            exact_annular_character_controls=fourier))

    # With an excessive radius the word map ceases to be injective.
    g = (2, 1, 3)
    assert compose_forms(compose_forms(g, g), g) == (1, 1, 6)
    canonical = {0: [0], 1: [1], 2: [-1]}
    old = event_flags(bytes([1, 1, 1]), canonical, 0)
    excessive = bytes(int(v[0] <= 0 and 2 * 2 ** abs(v[0]) > 4) for v in canonical.values())
    assert old == bytes([0, 0, 1]) and excessive == bytes(3)
    return dict(scope="Saved-word boundary identities and a fixed height counterexample; no new orbit search.",
                controls=results,
                contributed_euler_control=contributed_euler_control(),
                excessive_height_negative_control=dict(discriminant=-23, prime=2, radius=4,
                                                        actual_events=1, incorrect_boundary_events=0),
                long_range_character_estimate_proved=False,
                complete_factoring_exponent="1/2; unchanged")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compute()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads(EXPECTED.read_text())
    print(json.dumps(result))


if __name__ == "__main__":
    main()
