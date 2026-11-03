#!/usr/bin/env python3
"""Exact controls for the conditional connected-event window criterion.

No group compositions: the arithmetic controls only read existing flags whose
equality to relation-change events was checked by the upstream event checker.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import argparse
import itertools
import json
import lzma
from math import comb, factorial, prod
from pathlib import Path

from check_entry_moment_cycles import cyclic_window_counts, entry_flags

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "evidence/event-cumulant-windows.json"


def cumulants(moments):
    """moments[0] = 1; return ordinary cumulants, including the mean."""
    result = [F(0)]
    for r in range(1, len(moments)):
        result.append(moments[r] - sum(comb(r - 1, j - 1) * result[j]
                                      * moments[r - j] for j in range(1, r)))
    return result


def moments_from_cumulants(values):
    result = [F(1)]
    for r in range(1, len(values)):
        result.append(sum(comb(r - 1, j - 1) * values[j] * result[r - j]
                          for j in range(1, r + 1)))
    return result


def centered_moments(moments):
    mean = moments[1]
    return [sum(comb(r, j) * (-mean) ** (r - j) * moments[j]
                for j in range(r + 1)) for r in range(len(moments))]


def minimum_integer_scale(values):
    mean = values[1]
    assert mean > 0

    def works(a):
        return all(abs(values[r]) <= mean * factorial(r) * a ** (r - 1)
                   for r in range(2, len(values)))

    low, high = 0, 1
    while not works(high):
        high *= 2
    while high - low > 1:
        mid = (low + high) // 2
        if works(mid):
            high = mid
        else:
            low = mid
    assert works(high) and (high == 1 or not works(high - 1))
    return high


def distribution_record(counts, order):
    hist = Counter(counts)
    h = len(counts)
    moments = [sum(F(number * value ** r, h) for value, number in hist.items())
               for r in range(2 * order + 1)]
    values = cumulants(moments)
    assert moments_from_cumulants(values) == moments
    centered = centered_moments(moments)
    shifted = values.copy()
    shifted[1] = 0
    assert moments_from_cumulants(shifted) == centered
    a = minimum_integer_scale(values)
    # The finite-order theorem must never certify a flag with an empty window.
    if moments[1] >= 36 * a * order:
        assert F(hist[0], h) <= F(1, 2 ** order)
        if h < 2 ** order:
            assert hist[0] == 0
    exact = [[str(x.numerator), str(x.denominator)] for x in values]
    return dict(mean=str(moments[1]), minimum_count=min(counts), maximum_count=max(counts),
                empty_windows=hist[0], checked_cumulant_orders=2 * order,
                minimum_integer_scale=a,
                mean_threshold_met=moments[1] >= 36 * a * order,
                exact_cumulants_sha256=sha256(json.dumps(exact).encode()).hexdigest())


def partitions(items):
    if not items:
        yield []
        return
    first, *rest = items
    for blocks in partitions(rest):
        yield [(first,)] + blocks
        for j in range(len(blocks)):
            yield blocks[:j] + [(first,) + blocks[j]] + blocks[j + 1:]


def joint_cumulant(samples):
    h = len(samples[0])
    value = F(0)
    for blocks in partitions(list(range(len(samples)))):
        term = F((-1) ** (len(blocks) - 1) * factorial(len(blocks) - 1))
        for block in blocks:
            term *= F(sum(prod(samples[i][x] for i in block) for x in range(h)), h)
        value += term
    return value


def compute():
    # Independent Bernoulli controls, with no finite-cycle interpretation.
    independent = []
    for p in (F(1, 8), F(1, 2), F(7, 8)):
        q, n = 12, 4096
        one = cumulants([F(1)] + [p] * (2 * q))
        values = [F(0)] + [n * x for x in one[1:]]
        a = minimum_integer_scale(values)
        assert a == 1 and values[1] >= 36 * a * q
        centered_values = values.copy()
        centered_values[1] = 0
        moment = moments_from_cumulants(centered_values)[2 * q]
        # This checks the actual even-moment bound, not only the void probability.
        assert moment * 2 ** q <= values[1] ** (2 * q)
        assert (1 - p) ** n <= F(1, 2 ** q)
        independent.append(dict(trials=n, probability=str(p), order=q, scale=a,
                                even_moment_void_bound_checked=True))

    # Summing mixed cumulants must include coincident time indices.
    flags = bytes([1, 0, 1, 1, 0, 0, 0])
    w = 3
    samples = [list(flags[j:] + flags[:j]) for j in range(w)]
    counts = cyclic_window_counts(flags, w)
    moments = [F(sum(x ** r for x in counts), len(flags)) for r in range(6)]
    values = cumulants(moments)
    for r in range(2, 6):
        total = sum(joint_cumulant([samples[i] for i in indices])
                    for indices in itertools.product(range(w), repeat=r))
        assert total == values[r]
    assert values[2] != sum(joint_cumulant([samples[i], samples[j]])
                            for i in range(w) for j in range(w) if i != j)

    # Pairwise independence does not give the higher connected bound.
    pairwise = [[(x >> 0) & 1 for x in range(4)],
                [(x >> 1) & 1 for x in range(4)],
                [((x >> 0) ^ (x >> 1)) & 1 for x in range(4)]]
    assert all(joint_cumulant([pairwise[i], pairwise[j]]) == 0
               for i in range(3) for j in range(i))
    assert joint_cumulant(pairwise) == -F(1, 8)

    regular = distribution_record(cyclic_window_counts(bytes([0, 1]) * 4096, 2048), 14)
    assert regular["mean_threshold_met"] and regular["empty_windows"] == 0
    clustered = distribution_record(cyclic_window_counts(bytes([1]) * 128 + bytes(128), 64), 9)
    assert clustered["empty_windows"] > 0 and clustered["minimum_integer_scale"] > 1

    # Check the numerical inequalities used in the proof by rational arithmetic.
    for q in range(1, 129):
        # At mu=36 A q, t=4q/mu, the positive-series exponent is q/2.
        assert F(16 * q * q, 36 * q) / (1 - F(4 * q, 36 * q)) == F(q, 2)
        assert factorial(2 * q) <= (2 * q) ** (2 * q)
        if q >= 2:
            assert q ** (q - 2) <= factorial(q) * 3 ** (q - 1)
    assert sum(F(1, 2 ** j * factorial(j)) for j in range(20)) < 2
    # A geometric tail bounds the omitted exp(1/2) terms strictly below 2.
    assert sum(F(1, 2 ** j * factorial(j)) for j in range(20)) + F(2, 2 ** 20 * factorial(20)) < 2

    saved = []
    receipt = json.loads((ROOT / "evidence/clustered-power-relations/receipt.json").read_text())
    pool = json.loads((ROOT / "evidence/clustered-power-relations/pool-audit.json").read_text())
    source = next(x for x in receipt["certificates"] if x["N"] == pool["N"])
    base = lzma.decompress((ROOT / source["flags_file"]).read_bytes())
    expected_events = json.loads((ROOT / "evidence/relation-change-events.json").read_text())
    for p in (5, 7):
        stride = pool["derived_generator_indices"][str(p)]
        raw = bytes(base[(j * stride) % len(base)] for j in range(len(base)))
        events = entry_flags(raw)
        old = next(x for x in expected_events["results"] if x["prime"] == p)
        assert sha256(events).hexdigest() == old["event_flags_sha256"]
        h = len(events)
        full_q = h.bit_length()
        assert sum(events) < 36 * full_q
        records = []
        for window in (3260, h // 4):
            records.append(dict(window=window,
                                **distribution_record(cyclic_window_counts(events, window), 6)))
        saved.append(dict(prime=p, period=h, events=sum(events),
                          event_flags_sha256=sha256(events).hexdigest(),
                          all_phase_order=full_q,
                          threshold_impossible_for_nonwrapping_windows=True,
                          controls=records))
    return dict(scope="Exact cumulant algebra and saved-event postprocessing; no new arithmetic walks.",
                independent_controls=independent, repeated_index_orders=[2, 3, 4, 5],
                pairwise_independent_nonzero_third_cumulant="-1/8",
                regular_control=regular, clustered_negative_control=clustered,
                saved_event_controls=saved,
                arithmetic_tree_bound_proved=False, target_indices_bounded=False,
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
