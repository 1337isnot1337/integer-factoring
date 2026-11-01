"""Select independent original integer rows using primitive temporary residuals."""

from functools import reduce
from math import gcd


class PrimitiveIntegerRowBasis:
    """The same original-row subsequence as rational Gaussian elimination.

Content division applies only to an internal rank-test residual. The
returned rows are untouched input relations; no half-relation is created.
    """

    def __init__(self, dimension):
        if not isinstance(dimension, int) or dimension < 0:
            raise ValueError('Require a nonnegative integer row dimension.')
        self.dimension = dimension
        self.pivots = {}
        self.rows = []
        self.statistics = {'row_attempts': 0, 'elimination_steps': 0,
                           'content_normalizations': 0, 'maximum_temporary_entry_bits': 0}

    def _primitive(self, row):
        if not row:
            return row
        self.statistics['maximum_temporary_entry_bits'] = max(
            self.statistics['maximum_temporary_entry_bits'],
            max(abs(value).bit_length() for value in row.values()))
        content = reduce(gcd, row.values(), 0)
        if content > 1:
            row = {j: value // content for j, value in row.items()}
            self.statistics['content_normalizations'] += 1
        if row[min(row)] < 0:
            row = {j: -value for j, value in row.items()}
        return row

    def add(self, row):
        if len(row) != self.dimension or any(not isinstance(value, int) for value in row):
            raise ValueError('One integer entry per generator is required.')
        self.statistics['row_attempts'] += 1
        residual = self._primitive({j: value for j, value in enumerate(row) if value})
        for pivot in sorted(self.pivots):
            coefficient = residual.get(pivot, 0)
            if not coefficient:
                continue
            previous = self.pivots[pivot]
            multiplier = previous[pivot]
            updated = {j: multiplier * value for j, value in residual.items()}
            for j, value in previous.items():
                value = updated.get(j, 0) - coefficient * value
                if value:
                    updated[j] = value
                else:
                    updated.pop(j, None)
            residual = self._primitive(updated)
            self.statistics['elimination_steps'] += 1
        if not residual:
            return False
        self.pivots[min(residual)] = residual
        self.rows.append(list(row))
        return True
