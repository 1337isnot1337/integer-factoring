"""Exact small-leaf quotient updates for the oriented atom product tree."""

from incremental_oriented_coprime_basis import _ProductTree, IncrementalOrientedBasis


class _QuotientUpdateTree(_ProductTree):
    """Replace each ancestor P by (P / old_leaf) * new_leaf exactly."""

    def __init__(self):
        super().__init__()
        self.exact_divisions = 0
        self.division_operand_bits = 0
        self.maximum_divisor_bits = 0
        self.maximum_update_multiplier_bits = 0

    def _set(self, slot, norm):
        previous_leaf = self.leaves[slot]
        assert norm >= 1 and previous_leaf >= 1
        self.leaves[slot] = norm
        node = self.capacity + slot
        self.products[node] = norm
        self.leaf_updates += 1
        while node > 1:
            node //= 2
            value = self.products[node]
            if previous_leaf != 1:
                self.exact_divisions += 1
                self.division_operand_bits += value.bit_length() + previous_leaf.bit_length()
                self.maximum_divisor_bits = max(self.maximum_divisor_bits,
                                                previous_leaf.bit_length())
                value, remainder = divmod(value, previous_leaf)
                assert remainder == 0
            if norm != 1:
                self.maximum_update_multiplier_bits = max(self.maximum_update_multiplier_bits,
                                                          norm.bit_length())
                value = self._multiply(value, norm)
            self.products[node] = value


class QuotientUpdateOrientedBasis(IncrementalOrientedBasis):
    """Identical refinement and overlap queries with cheaper tree updates."""

    def __init__(self, discriminant):
        super().__init__(discriminant)
        self.tree = _QuotientUpdateTree()

    def result(self, *, dense=False, verify=True):
        record = super().result(dense=dense, verify=verify)
        record['tree_update_method'] = 'exact leaf quotient'
        record['statistics'].update(
            tree_exact_divisions=self.tree.exact_divisions,
            tree_division_operand_bits=self.tree.division_operand_bits,
            tree_maximum_divisor_bits=self.tree.maximum_divisor_bits,
            tree_maximum_update_multiplier_bits=self.tree.maximum_update_multiplier_bits)
        return record


def quotient_update_oriented_coprime_basis(discriminant, keys, *, dense=False, verify=True):
    state = QuotientUpdateOrientedBasis(discriminant)
    for norm, root in keys:
        state.add(norm, root)
    return state.result(dense=dense, verify=verify)
