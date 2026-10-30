"""Finite checks for the written-table adapter used by the half-rate proof.

This is a small exhaustive reference, not the large factoring algorithm or
the Appendix B evaluator. Bits and off-diagonal blocks are lowest first,
as in both pinned source projects. Integer table entries are kept exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from math import prod


def to_word(length: int, value: int) -> tuple[int, ...]:
    if length < 0 or not 0 <= value < 1 << length:
        raise ValueError("value does not fit the word")
    return tuple((value >> i) & 1 for i in range(length))


def of_word(word: tuple[int, ...]) -> int:
    if any(bit not in (0, 1) for bit in word):
        raise ValueError("word contains a non-bit")
    return sum(bit << i for i, bit in enumerate(word))


def unflat(k: int, width: int, vertex: int, index: int) -> tuple[tuple[int, ...], ...]:
    """Appendix B's unflat: the self block is empty, not width zero bits."""
    if k < 1 or width < 1 or not 0 <= vertex < k:
        raise ValueError("invalid clique shape")
    word = to_word((k - 1) * width, index)
    return tuple(() if other == vertex else
                 word[(other if other < vertex else other - 1) * width:
                      (other if other < vertex else other - 1) * width + width]
                 for other in range(k))


def source_tables_to_pinst(rows: tuple[tuple[int, ...], ...], width: int
                           ) -> tuple[tuple[int, ...], ...]:
    """Remove a source table's self block, evaluating it at all-zero bits."""
    k = len(rows)
    if k < 1 or width < 1 or any(len(row) != 1 << (k * width) for row in rows):
        raise ValueError("invalid source clique tables")
    result = []
    for vertex, row in enumerate(rows):
        cut = vertex * width
        mask = (1 << cut) - 1
        result.append(tuple(row[(index & mask) | ((index >> cut) << (cut + width))]
                            for index in range(1 << ((k - 1) * width))))
    return tuple(result)


def source_clique_value(rows: tuple[tuple[int, ...], ...], width: int) -> int:
    """Literal zero-diagonal/symmetric-array definition, for tiny inputs."""
    k = len(rows)
    if k < 1 or width < 1 or any(len(row) != 1 << (k * width) for row in rows):
        raise ValueError("invalid source clique tables")
    total = 0
    for array in product(range(1 << width), repeat=k * k):
        if any(array[v * k + v] for v in range(k)):
            continue
        if any(array[v * k + u] != array[u * k + v]
               for v, u in combinations(range(k), 2)):
            continue
        total += prod(rows[v][sum(array[v * k + u] << (u * width)
                                  for u in range(k))] for v in range(k))
    return total


def pinst_value(rows: tuple[tuple[int, ...], ...], width: int,
                modulus: int | None = None) -> int:
    """Sum over one word per unordered edge, the target's zval definition."""
    k = len(rows)
    if k < 1 or width < 1 or any(len(row) != 1 << ((k - 1) * width) for row in rows):
        raise ValueError("invalid written instance")
    if modulus is not None and modulus < 2:
        raise ValueError("modulus must be at least two")
    pairs = tuple(combinations(range(k), 2))
    total = 0
    for values in product(range(1 << width), repeat=len(pairs)):
        edges = dict(zip(pairs, values))
        factor = 1
        for vertex in range(k):
            blocks = [to_word(width, edges[min(vertex, other), max(vertex, other)])
                      for other in range(k) if other != vertex]
            index = of_word(tuple(bit for block in blocks for bit in block))
            factor *= rows[vertex][index]
        total += factor
        if modulus is not None:
            total %= modulus
    return total


@dataclass(frozen=True)
class MessageEdge:
    first: int
    second: int
    width: int
    first_messages: tuple[int, ...]
    second_messages: tuple[int, ...]


@dataclass(frozen=True)
class EqualityNetwork:
    """A small 0/1 gate network with deterministic outgoing messages."""

    weights: tuple[tuple[int, ...], ...]
    edges: tuple[MessageEdge, ...]

    def __post_init__(self) -> None:
        if not self.weights or any(not row or any(w not in (0, 1) for w in row)
                                   for row in self.weights):
            raise ValueError("nonempty 0/1 gate weights required")
        seen = set()
        for edge in self.edges:
            if not 0 <= edge.first < edge.second < len(self.weights) or edge.width < 0:
                raise ValueError("invalid edge")
            pair = edge.first, edge.second
            if pair in seen:
                raise ValueError("duplicate edge")
            seen.add(pair)
            for vertex, messages in ((edge.first, edge.first_messages),
                                     (edge.second, edge.second_messages)):
                if len(messages) != len(self.weights[vertex]):
                    raise ValueError("one message per state required")
                if any(not 0 <= value < 1 << edge.width for value in messages):
                    raise ValueError("message does not fit its width")

    def value(self) -> int:
        total = 0
        for states in product(*(range(len(row)) for row in self.weights)):
            if all(self.weights[v][states[v]] for v in range(len(states))) and all(
                    edge.first_messages[states[edge.first]] ==
                    edge.second_messages[states[edge.second]] for edge in self.edges):
                total += 1
        return total

    def group_tables(self, groups: tuple[int, ...], k: int, width: int
                     ) -> tuple[tuple[int, ...], ...]:
        """Scatter each internally valid local setting into its unique row entry.

        Each group pair uses the same source-edge order at both ends. Unused
        high channel bits stay zero. No internal channel assignments are read.
        """
        if len(groups) != len(self.weights) or k < 1 or width < 1:
            raise ValueError("invalid grouping")
        if any(not 0 <= group < k for group in groups):
            raise ValueError("group outside range")
        edge_order = sorted(self.edges, key=lambda edge: (edge.first, edge.second))
        cross = {pair: tuple(edge for edge in edge_order
                             if tuple(sorted((groups[edge.first], groups[edge.second]))) == pair)
                 for pair in combinations(range(k), 2)}
        if any(sum(edge.width for edge in edges) > width for edges in cross.values()):
            raise ValueError("grouped edge does not fit the common width")
        result = []
        for vertex in range(k):
            members = [v for v, group in enumerate(groups) if group == vertex]
            internal = [edge for edge in edge_order
                        if groups[edge.first] == groups[edge.second] == vertex]
            row = [0] * (1 << ((k - 1) * width))
            for choices in product(*(range(len(self.weights[v])) for v in members)):
                state = dict(zip(members, choices))
                if not all(self.weights[v][state[v]] for v in members):
                    continue
                if any(edge.first_messages[state[edge.first]] !=
                       edge.second_messages[state[edge.second]] for edge in internal):
                    continue
                index = 0
                block_index = 0
                for other in range(k):
                    if other == vertex:
                        continue
                    label = offset = 0
                    for edge in cross[min(vertex, other), max(vertex, other)]:
                        if groups[edge.first] == vertex:
                            message = edge.first_messages[state[edge.first]]
                        else:
                            message = edge.second_messages[state[edge.second]]
                        label |= message << offset
                        offset += edge.width
                    index |= label << (block_index * width)
                    block_index += 1
                row[index] += 1
            result.append(tuple(row))
        return tuple(result)
