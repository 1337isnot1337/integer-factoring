"""Deterministic bounded-width bidirectional word search in a form class group.

Pruning makes failure inconclusive.  Every positive answer is checked as
an exact bounded word; a negative answer is exact only when neither side
ever pruned a class.  No general success guarantee at small width is used.
Sorted arrays give deterministic lookup and avoid a hash-table assumption.
"""

from bisect import bisect_left
from heapq import merge

from quadratic_form_arithmetic import (
    compose_primitive, form_discriminant, inverse_primitive,
    power_primitive, reduce_primitive,
)


def choose_frontier(entries, width, selection):
    """Choose at most width entries; entries are sorted by canonical form."""
    if selection == "height":
        return entries[:width]
    if selection != "stratified":
        raise ValueError("selection must be height or stratified")
    buckets = []
    for entry in entries:
        height_bits = entry[0][0].bit_length()
        if not buckets or buckets[-1][0] != height_bits:
            buckets.append((height_bits, []))
        buckets[-1][1].append(entry)
    chosen, level = [], 0
    while len(chosen) < width:
        row = [entries[level] for _, entries in buckets if len(entries) > level]
        if not row:
            break
        chosen.extend(row[:width-len(chosen)])
        level += 1
    return chosen


class BoundedClassWordSearch:
    """Reuse an identity-side frontier for several target classes."""

    def __init__(self, discriminant, generators, width, left_depth, right_depth,
                 selection="stratified"):
        if width < 1 or min(left_depth, right_depth) < 0:
            raise ValueError("positive width and nonnegative depths are required")
        if discriminant >= 0 or discriminant % 4 not in (0, 1):
            raise ValueError("a negative order discriminant is required")
        if selection not in ("height", "stratified"):
            raise ValueError("unknown frontier selection rule")
        self.discriminant = discriminant
        self.generators = [reduce_primitive(generator) for generator in generators]
        if any(form_discriminant(generator) != discriminant for generator in self.generators):
            raise ValueError("every generator must have the requested discriminant")
        parity = discriminant % 2
        self.identity = (1, parity, (parity-discriminant)//4)
        self.width, self.left_depth, self.right_depth = width, left_depth, right_depth
        self.selection = selection
        self.moves = [(index, sign,
                       generator if sign == 1 else inverse_primitive(generator))
                      for index, generator in enumerate(self.generators)
                      for sign in (-1, 1)]
        self.left_entries, self.left_statistics, _ = self._walk(self.identity, left_depth)
        self.left_keys = [entry[0] for entry in self.left_entries]

    def _walk(self, start, budget, opposite=None, opposite_keys=None):
        zero = (0,)*len(self.generators)
        seen, seen_keys = [(start, zero, 0)], [start]
        frontier = seen.copy()
        stats = {"group_compositions": 0, "pruned_class_occurrences": 0,
                 "maximum_frontier": 1, "layers": []}

        def meeting(form, word, depth):
            if opposite is None:
                return None
            position = bisect_left(opposite_keys, form)
            if position == len(opposite_keys) or opposite_keys[position] != form:
                return None
            left_form, left_word, left_depth = opposite[position]
            return {"meeting_class": left_form, "left_word": left_word,
                    "left_path_depth": left_depth, "right_word": word,
                    "right_path_depth": depth,
                    "word": tuple(a-b for a, b in zip(left_word, word))}

        answer = meeting(start, zero, 0)
        if answer is not None:
            stats["retained_classes"] = 1
            return seen, stats, answer
        for depth in range(1, budget+1):
            candidates = []
            for form, word, _ in frontier:
                for index, sign, generator in self.moves:
                    neighbor = compose_primitive(form, generator)
                    stats["group_compositions"] += 1
                    new_word = list(word)
                    new_word[index] += sign
                    new_word = tuple(new_word)
                    answer = meeting(neighbor, new_word, depth)
                    if answer is not None:
                        stats["retained_classes"] = len(seen)
                        stats["meeting_in_partial_layer"] = depth
                        return seen, stats, answer
                    position = bisect_left(seen_keys, neighbor)
                    if position < len(seen_keys) and seen_keys[position] == neighbor:
                        continue
                    candidates.append((neighbor, new_word, depth))
            candidates.sort()
            distinct = []
            for entry in candidates:
                if not distinct or entry[0] != distinct[-1][0]:
                    distinct.append(entry)
            chosen = choose_frontier(distinct, self.width, self.selection)
            pruned = len(distinct)-len(chosen)
            stats["pruned_class_occurrences"] += pruned
            stats["layers"].append({"depth": depth, "input_frontier": len(frontier),
                                    "distinct_new_classes": len(distinct),
                                    "retained_frontier": len(chosen),
                                    "pruned_new_classes": pruned})
            stats["maximum_frontier"] = max(stats["maximum_frontier"], len(chosen))
            if not chosen:
                break
            seen = list(merge(seen, sorted(chosen)))
            seen_keys = [entry[0] for entry in seen]
            frontier = chosen
        stats["retained_classes"] = len(seen)
        return seen, stats, None

    def statistics(self):
        return {"width": self.width, "left_depth": self.left_depth,
                "right_depth": self.right_depth, "selection": self.selection,
                "word_length_bound": self.left_depth+self.right_depth,
                "left": self.left_statistics}

    def solve(self, form):
        if form_discriminant(form) != self.discriminant:
            raise ValueError("the target has the wrong discriminant")
        target = reduce_primitive(form)
        _, stats, meeting = self._walk(target, self.right_depth,
                                      self.left_entries, self.left_keys)
        if meeting is None:
            complete = (self.left_statistics["pruned_class_occurrences"] == 0
                        and stats["pruned_class_occurrences"] == 0)
            return {"found": False,
                    "status": "no_bounded_word" if complete else "incomplete_search",
                    "negative_answer_is_exact": complete, "right_statistics": stats}
        word = meeting["word"]
        assert sum(abs(value) for value in word) <= self.left_depth+self.right_depth
        value = self.identity
        for generator, exponent in zip(self.generators, word):
            value = compose_primitive(value, power_primitive(generator, exponent))
        assert value == target
        return {"found": True, "status": "found", **meeting,
                "word_l1": sum(abs(value) for value in word),
                "right_statistics": stats}
