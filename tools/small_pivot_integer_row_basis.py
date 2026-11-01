"""Select original relations with exact primitive residuals and small pivots."""

from primitive_integer_row_basis import PrimitiveIntegerRowBasis


class SmallPivotIntegerRowBasis(PrimitiveIntegerRowBasis):
    """Keep the original-row subsequence while changing only internal pivot order.

    Stored residuals vanish at every earlier pivot. Dictionary insertion
    order therefore gives an exact forward elimination even when pivot
    columns are not numerically ordered. Content division never changes
    an emitted relation.
    """

    def add(self, row):
        if len(row) != self.dimension or any(not isinstance(value, int) for value in row):
            raise ValueError('One integer entry per generator is required.')
        self.statistics['row_attempts'] += 1
        residual = self._primitive({j: value for j, value in enumerate(row) if value})
        for pivot, previous in self.pivots.items():
            if not residual:
                return False
            coefficient = residual.get(pivot, 0)
            if not coefficient:
                continue
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
        pivot = min(residual, key=lambda j: (abs(residual[j]), j))
        if residual[pivot] < 0:
            residual = {j: -value for j, value in residual.items()}
        assert pivot not in self.pivots
        self.pivots[pivot] = residual
        self.rows.append(list(row))
        return True
