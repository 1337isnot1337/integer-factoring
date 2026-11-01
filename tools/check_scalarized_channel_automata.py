#!/usr/bin/env python3
"""Exact product-space automaton check over Mat_2(F_2)."""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import time


def matrix_product(a, b):
    answer = 0
    for i, j in product(range(2), repeat=2):
        value = sum(((a >> (2*i+k)) & 1)*((b >> (2*k+j)) & 1) for k in range(2)) & 1
        answer |= value << (2*i+j)
    return answer


def basis(vectors):
    rows = [value for value in vectors if value]
    rank = 0
    for bit in range(3, -1, -1):
        pivot = next((i for i in range(rank, len(rows)) if rows[i] >> bit & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> bit & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return tuple(rows[:rank])


def space_product(a, b):
    return basis(matrix_product(x, y) for x, y in product(a, b))


def value(functional, matrix):
    return (functional & matrix).bit_count() & 1


def bilinear_rank(functional, left, right):
    rows = [sum(value(functional, matrix_product(a, b)) << j for j, b in enumerate(right))
            for a in left]
    return len(basis(rows))


def automaton(channel_spaces):
    states = [((9,), (9,))]
    indices = {states[0]: 0}
    transitions = []
    for left, right in states:
        targets = []
        for u, v in channel_spaces:
            target = space_product(left, u), space_product(right, v)
            if target not in indices:
                indices[target] = len(states)
                states.append(target)
            targets.append(indices[target])
        transitions.append(targets)
    return states, transitions


def coreachable(transitions, terminal):
    reverse = [set() for _ in transitions]
    for q, row in enumerate(transitions):
        for target in row:
            reverse[target].add(q)
    found = {q for q, weight in enumerate(terminal) if weight}
    queue = deque(found)
    while queue:
        q = queue.popleft()
        for previous in reverse[q]:
            if previous not in found:
                found.add(previous)
                queue.append(previous)
    return found


def type_counts(transitions, terminal, labels, k):
    dynamic = {(0, (0, 0, 0)): 1}
    for _ in range(k):
        following = defaultdict(int)
        for (state, kind), count in dynamic.items():
            for channel, target in enumerate(transitions[state]):
                next_kind = list(kind)
                next_kind[labels[channel]] += 1
                following[target, tuple(next_kind)] += count
        dynamic = following
    counts = defaultdict(int)
    for (state, kind), count in dynamic.items():
        counts[kind] += count*terminal[state]
    return {kind: count for kind, count in counts.items() if count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    nilpotents = (2, 4, 15)
    labels = [i for i in range(3) for _ in range(4)]
    channel_spaces = [((nilpotents[i],), (9,)) for i in labels]
    states, transitions = automaton(channel_spaces)
    assert len(states) == 11
    transition_checks = 0
    for q, (left, right) in enumerate(states):
        for channel, (u, v) in enumerate(channel_spaces):
            target = states[transitions[q][channel]]
            assert target == (space_product(left, u), space_product(right, v))
            transition_checks += 1
    functionals = [ell for ell in range(16) if value(ell, 9) == 1]
    assert len(functionals) == 8
    # A zero terminal at the current prefix need not remain zero after
    # another channel.  E01 is invisible to ell11, but E01 E10 = E00.
    ell11_terminal = [bilinear_rank(1, left, right) for left, right in states]
    prefix = transitions[0][0]
    continuation = transitions[prefix][4]
    assert states[prefix] == ((2,), (9,))
    assert ell11_terminal[prefix] == 0
    assert ell11_terminal[continuation] == 1
    assert prefix in coreachable(transitions, ell11_terminal)
    direct_word_checks = rank_checks = supervector_checks = 0
    summaries = []
    for ell in functionals:
        terminal = [bilinear_rank(ell, left, right) for left, right in states]
        assert all(rank in (0, 1) for rank in terminal)
        rank_checks += len(terminal)
        alive = coreachable(transitions, terminal)
        assert len(alive) in (7, 10) and 0 in alive
        h = [Fraction(3, 2) if q == 0 else Fraction(1) for q in range(len(states))]
        for q in alive:
            actual = sum(h[target] for target in transitions[q] if target in alive)
            assert actual <= 8*h[q]
            supervector_checks += 1
        totals = []
        for k in range(1, 7):
            counts = type_counts(transitions, terminal, labels, k)
            total = sum(counts.values())
            assert total <= h[0]*8**k
            if ell == 1:
                assert total == 4**k*(2**(k+1)+(-1)**k)//3
            if k <= 3:
                direct = defaultdict(int)
                for word in product(range(3), repeat=k):
                    coefficient = 9
                    for letter in word:
                        coefficient = matrix_product(coefficient, nilpotents[letter])
                    if value(ell, coefficient):
                        kind = tuple(word.count(i) for i in range(3))
                        direct[kind] += 4**k
                    direct_word_checks += 1
                assert counts == dict(direct)
            totals.append(total)
        summaries.append({"functional_mask": ell,
                          "live_states": len(alive),
                          "nonzero_terminal_states": sum(bool(t) for t in terminal),
                          "counts_k1_through_k6": totals})
    output = {"status": "passed", "reachable_product_space_states": len(states),
              "normalized_functionals_checked": len(functionals),
              "transition_space_identities": transition_checks,
              "terminal_bilinear_rank_checks": rank_checks,
              "direct_word_checks": direct_word_checks,
              "rational_supervector_checks": supervector_checks,
              "zero_terminal_prefix_requires_retention": True,
              "functional_summaries": summaries,
              "scope": "Exact scalarized automaton and terminal-rank checks; this toy still has base8 and no score improvement.",
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "seconds": time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
