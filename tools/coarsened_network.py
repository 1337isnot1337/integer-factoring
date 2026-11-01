"""Tiny exact reference for network coarsening and padded clique tables.

This enumerates states and labels explicitly. It is a finite integration
checker, not the uniform factoring evaluator or a runtime benchmark.
"""

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, product
from math import prod
from typing import Callable


@dataclass(frozen=True)
class Edge:
    first: int
    second: int
    bits: int
    left: Callable[[object, int], int]
    right: Callable[[object, int], int]


@dataclass(frozen=True)
class Network:
    states: tuple[tuple, ...]
    weights: tuple[Callable[[object], int], ...]
    edges: tuple[Edge, ...]

    def __post_init__(self):
        if len(self.states) != len(self.weights) or any(not s for s in self.states):
            raise ValueError("each node must have weights and a nonempty state set")
        pairs = [(e.first, e.second) for e in self.edges]
        if len(set(pairs)) != len(pairs) or any(
            not 0 <= e.first < e.second < len(self.states) or e.bits < 0
            for e in self.edges
        ):
            raise ValueError("edges must have distinct ordered endpoints and nonnegative widths")


def network_value(network: Network, modulus: int | None = None) -> int:
    """Contract the original network by enumerating settings and channels."""
    if modulus is not None and modulus < 2:
        raise ValueError("modulus must be at least 2")
    total = 0
    for state in product(*network.states):
        value = prod(weight(z) for weight, z in zip(network.weights, state))
        for e in network.edges:
            value *= sum(e.left(state[e.first], x) * e.right(state[e.second], x)
                         for x in range(1 << e.bits))
        total += value
    return total if modulus is None else total % modulus


def clique_value(network: Network, grouping: tuple[int, ...], groups: int,
                 width: int, modulus: int | None = None) -> int:
    """Build group tables and contract a common-width padded clique.

    Each table independently sums the node states and the channels of its
    internal edges. Cross-group channels are packed in original edge order;
    higher padding bits are forced to zero. Empty groups are permitted.
    """
    if (groups < 1 or width < 0 or len(grouping) != len(network.states)
            or any(not 0 <= group < groups for group in grouping)):
        raise ValueError("invalid grouping or clique width")
    if modulus is not None and modulus < 2:
        raise ValueError("modulus must be at least 2")

    pairs = tuple(combinations(range(groups), 2))
    nodes = tuple(tuple(i for i, g in enumerate(grouping) if g == group)
                  for group in range(groups))
    cross = {pair: tuple(e for e in network.edges
                         if {grouping[e.first], grouping[e.second]} == set(pair))
             for pair in pairs}
    cross_bits = {pair: sum(e.bits for e in cross[pair]) for pair in pairs}
    if any(bits > width for bits in cross_bits.values()):
        raise ValueError("common width cannot hold all crossing channels")
    internal = tuple(tuple(e for e in network.edges
                           if grouping[e.first] == group == grouping[e.second])
                     for group in range(groups))

    @lru_cache(maxsize=None)
    def table(group: int, labels: tuple[int, ...]) -> int:
        incident = tuple(pair for pair in pairs if group in pair)
        label_by_pair = dict(zip(incident, labels))
        if any(label_by_pair[pair] >= 1 << cross_bits[pair] for pair in incident):
            return 0
        total = 0
        for settings in product(*(network.states[n] for n in nodes[group])):
            z = dict(zip(nodes[group], settings))
            value = prod(network.weights[n](z[n]) for n in nodes[group])
            for edge in internal[group]:
                value *= sum(edge.left(z[edge.first], x) * edge.right(z[edge.second], x)
                             for x in range(1 << edge.bits))
            for pair in incident:
                packed = label_by_pair[pair]
                for edge in cross[pair]:
                    channel = packed & ((1 << edge.bits) - 1)
                    packed >>= edge.bits
                    if grouping[edge.first] == group:
                        value *= edge.left(z[edge.first], channel)
                    else:
                        value *= edge.right(z[edge.second], channel)
            total += value
        return total if modulus is None else total % modulus

    total = 0
    for labels in product(range(1 << width), repeat=len(pairs)):
        by_pair = dict(zip(pairs, labels))
        total += prod(table(group, tuple(by_pair[pair] for pair in pairs if group in pair))
                      for group in range(groups))
        if modulus is not None:
            total %= modulus
    return total
