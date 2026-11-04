#!/usr/bin/env python3
"""Exact saved-row controls for the principal-word gap and event separation."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd, isqrt, prod
from pathlib import Path

from binary_quadratic_forms import compose_forms
from check_event_cumulant_windows import joint_cumulant
from check_relation_change_events import event_flags
from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "evidence/small-norm-event-separation.json"


def capped_norm(primes, word, cap):
    result = 1
    for p, e in zip(primes, word):
        e = abs(e)
        if e >= cap.bit_length():
            return cap
        result *= p ** e
        if result >= cap:
            return cap
    return result


def compute():
    source = ROOT / "evidence/clustered-power-relations/row-audit.json"
    data = json.loads(source.read_text())
    primes = data["primes"]
    results = []
    for case in data["results"]:
        p, h, n, initial_k = case["prime"], case["period"], case["N"], case["K"]
        delta = -n if n % 4 == 3 else -4 * n
        d = abs(delta)
        assert gcd(prod(primes), d) == 1
        coordinate = primes.index(p)
        words = {k: [k if j == coordinate else 0 for j in range(len(primes))]
                 for k in range(initial_k)}
        for record in case["rows"]:
            k = record["k"]
            word = [(k if j == coordinate else 0) - r for j, r in enumerate(record["row"])]
            assert prod(q ** abs(e) for q, e in zip(primes, word)) == record["leading"]
            words[k] = word
        phases = sorted(words)
        checked, equal = 0, 0
        cap = (d + 3) // 4
        for i, k in enumerate(phases):
            for ell in phases[i + 1:]:
                difference = [(ell - k if j == coordinate else 0) - words[ell][j] + words[k][j]
                              for j in range(len(primes))]
                if not any(difference):
                    equal += 1
                    continue
                assert 4 * capped_norm(primes, difference, cap) >= d
                checked += 1

        # One prescribed separation length, no cutoff or orbit search.
        separation = 2
        radius = min(isqrt(d) // 4, isqrt((d - 1) // (4 * p ** separation)))
        assert radius >= 1 and 4 * p ** separation * radius ** 2 < d
        low_words = {k: v for k, v in words.items()
                     if prod(q ** abs(e) for q, e in zip(primes, v)) <= radius}
        hits = bytes(int(k in low_words) for k in range(h))
        events = event_flags(hits, low_words, coordinate)
        k_radius = least_power_above(p, radius)
        assert sum(hits) <= (2 * k_radius - 1) * sum(events)
        positions = [k for k, event in enumerate(events) if event]
        assert positions
        gaps = [(positions[(j + 1) % len(positions)] - k) % h or h
                for j, k in enumerate(positions)]
        assert min(gaps) > separation
        for shift in range(1, separation + 1):
            assert not any(events[k] and events[(k + shift) % h] for k in positions)
        density = F(sum(events), h)
        samples = [list(events[j:] + events[:j]) for j in range(separation + 1)]
        assert joint_cumulant(samples[:2]) == -density ** 2
        assert joint_cumulant(samples) == 2 * density ** 3
        results.append(dict(prime=p, period=h, unramified_word_pairs=checked,
                            equal_row_pairs=equal, separation=separation, radius=radius,
                            low_norm_hits=sum(hits), low_norm_events=sum(events),
                            row_multiplicity_cap=2 * k_radius - 1,
                            minimum_cyclic_event_distance=min(gaps),
                            covariance=str(-density ** 2), third_cumulant=str(2 * density ** 3),
                            event_flags_sha256=sha256(events).hexdigest()))

    # A ramified square is a rational principal ideal and must be excluded.
    ramified = (2, 2, 3)
    assert compose_forms(ramified, ramified) == (1, 0, 5)
    assert 4 * 2 ** 2 < 20
    semantic = {k: v for k, v in data.items() if k != "seconds"}
    return dict(scope="Principal-word norm gap and fixed T=2 low-height event separation; no new orbit.",
                source_row_semantic_sha256=sha256(json.dumps(semantic, sort_keys=True,
                                                            separators=(",", ":")).encode()).hexdigest(),
                controls=results,
                ramified_negative_control=dict(discriminant=-20, prime=2, word_exponent=2,
                                                word_norm=4, violates_unqualified_gap=True),
                long_range_cumulant_estimate_proved=False,
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
