#!/usr/bin/env python3
"""Integer-only replay of periodic events when a smooth relation changes."""
import argparse
from hashlib import sha256
import json
import lzma
from math import isqrt, prod
from pathlib import Path

from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "evidence/clustered-power-relations"
EXPECTED = ROOT / "evidence/relation-change-events.json"


def event_flags(smooth, words, coordinate, factors=None):
    h = len(smooth)
    factors = factors or bytes(h)
    result = []
    for k in range(h):
        previous = (k - 1) % h
        same = (smooth[k] and smooth[previous] and
                all(a - b == int(i == coordinate)
                    for i, (a, b) in enumerate(zip(words[k], words[previous]))))
        result.append(int(bool(factors[k] or (smooth[k] and not same))))
    return bytes(result)


def segments(hits, events):
    h = len(hits)
    result = []
    for k, event in enumerate(events):
        if not event:
            continue
        length = 0
        for offset in range(h):
            j = (k + offset) % h
            if not hits[j] or (offset and events[j]):
                break
            length += 1
        result.append(length)
    assert result and sum(result) == sum(hits)
    return result


def compute():
    source = {name:json.loads((DATA / name).read_text())
              for name in ("receipt.json", "row-audit.json", "pool-audit.json")}
    rows, pool, orbit = source["row-audit.json"], source["pool-audit.json"], source["receipt.json"]
    primes = rows["primes"]
    height = isqrt(pool["N"] // 3)
    certificate = next(c for c in orbit["certificates"] if c["N"] == pool["N"])
    base = lzma.decompress((ROOT / certificate["flags_file"]).read_bytes())
    results = []
    for case in rows["results"]:
        p, K, h = case["prime"], case["K"], case["period"]
        i = primes.index(p)
        jump = pool["derived_generator_indices"][str(p)]
        flags = bytes(base[(k * jump) % h] for k in range(h))
        assert K == least_power_above(p, height)
        # These prescribed prefixes are below sqrt(|Delta|)/2, where
        # the primitive prime-power ideal is already reduced.
        assert 4 * p ** (2 * (K - 1)) < pool["N"]
        assert all(flags[k] for k in range(K))
        words = {k:[k if j == i else 0 for j in range(len(primes))] for k in range(K)}
        for row in case["rows"]:
            k = row["k"]
            words[k] = [(k if j == i else 0) - a for j, a in enumerate(row["row"])]
            assert prod(q ** abs(e) for q, e in zip(primes, words[k])) == row["leading"] <= height
        assert set(words) == {k for k in range(h) if flags[k]}
        events = event_flags(flags, words, i)
        pure_entries = bytes(flags[k] * (1 - flags[k - 1]) for k in range(h))
        lengths = segments(flags, events)
        bound = 2 * K - 1
        assert max(lengths) <= bound and sum(flags) <= bound * sum(events)
        assert events == pure_entries
        results.append(dict(prime=p, period=h, K=K, multiplicity_bound=bound,
                            hits=sum(flags), events=sum(events),
                            maximum_identical_row_segment=max(lengths),
                            event_flags_sha256=sha256(events).hexdigest(),
                            equals_previous_pure_entry_process=True))
    # A sharp all-good periodic-word model. The group is Z/3, its
    # canonical words are 0,1,-1, and every norm is at most H=2.
    toy_flags = bytes([1, 1, 1])
    toy_words = {0:[0], 1:[1], 2:[-1]}
    toy_events = event_flags(toy_flags, toy_words, 0)
    assert toy_events == bytes([0, 0, 1])
    assert segments(toy_flags, toy_events) == [3]
    assert all(2 ** abs(toy_words[k][0]) <= 2 and toy_words[k][0] % 3 == k for k in range(3))
    # Factor-producing phases count as events without needing a word.
    factor_events = event_flags(bytes([1, 0, 0]), {0:[0]}, 0, bytes([0, 1, 0]))
    assert factor_events == bytes([1, 1, 0])
    return dict(scope="Periodic relation-change density bridge; integer postprocessing only, no new arithmetic cycle.",
                input_semantic_sha256={name:sha256(json.dumps({k:v for k,v in data.items() if k != "seconds"},
                                                             sort_keys=True, separators=(",", ":")).encode()).hexdigest()
                                       for name, data in source.items()},
                results=results,
                all_good_control=dict(group_order=3, prime=2, height=2, K=2,
                                      hits=3, events=1, pure_entries=0, bound=3,
                                      density_bound_attained=True),
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
