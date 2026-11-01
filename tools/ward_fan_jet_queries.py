#!/usr/bin/env python3
"""Exact truncated-eta queries for the fixed 84-label duplicated fan map."""

from collections import Counter

from check_ward_fan_leading_block import FORMS, local_query_table
from complement_shell_tensor import multiply
from ward_fan_entry_queries import fan_entry_query


STRIDE = 16


class FanJetOracle:
    def __init__(self, code, order=1, endpoint=False):
        if not 0 <= order <= 7:
            raise ValueError('The requested quotient jet order must lie between zero and seven.')
        self.code = code
        self.order = order
        self.endpoint = endpoint
        self.mask = sum(((1 << (8 - b)) - 1) << (STRIDE * b) for b in range(order + 1))
        y = 1 << STRIDE
        local = (local_query_table(0, y, 1, 1 ^ y) if endpoint
                 else local_query_table(3, 1 ^ y, 1, y))
        self.local = [(a & self.mask, b & self.mask) for a, b in local]
        self.pair_cache = {}
        self.physical_cache = {}
        self.logical_cache = {}
        self.counts = Counter()

    def mul(self, first, second):
        # Every surviving x exponent is <=7, so no x carry reaches stride16.
        return multiply(first, second) & self.mask

    def valuation(self, value):
        if not value:
            return self.order + 1
        return ((value & -value).bit_length() - 1) // STRIDE

    def physical(self, x, y, u, v):
        self.counts['physical_query_requests'] += 1
        keys = []
        lower_order = 0
        for edge in range(2 * FORMS - 1):
            entries = []
            for copy in range(2):
                site = 2 * edge + copy
                index = (((x >> (2 * site) & 3) * 4 + (y >> (2 * site) & 3)) * 4
                         + (u >> (2 * site) & 3)) * 4 + (v >> (2 * site) & 3)
                entries.append(self.local[index])
            key = tuple(sorted(entries))
            if key not in self.pair_cache:
                (m1, b1), (m2, b2) = key
                self.pair_cache[key] = (self.mul(m1, m2),
                                        self.mul(b1, m2) ^ self.mul(m1, b2))
            value = self.pair_cache[key]
            lower_order += min(self.valuation(value[0]), self.valuation(value[1]))
            if lower_order > self.order:
                self.counts['local_valuation_zero_queries'] += 1
                return 0
            keys.append(value)
        key = tuple(keys)
        if key not in self.physical_cache:
            sites = [value for pair in key for value in (pair, (1, 0))]
            self.physical_cache[key] = fan_entry_query(FORMS, iter(sites), self.mul)
            self.counts['distinct_physical_jets_computed'] += 1
        return self.physical_cache[key]

    def polynomial(self, x, y, u, v):
        key = x, y, u, v
        if key not in self.logical_cache:
            result = 0
            for yy in self.code[y]['preimage']:
                for vv in self.code[v]['input']:
                    result ^= self.physical(self.code[x]['output'], yy,
                                            self.code[u]['output'], vv)
            if not self.endpoint:
                # Independently certified complete lowest boundary.
                assert result & 255 == 0
                for b in range(self.order + 1):
                    forbidden = (1 << max(0, 3 - b)) - 1
                    assert not (result >> (STRIDE * b)) & forbidden
            self.logical_cache[key] = result
        return self.logical_cache[key]


def original_jet(code, labels, local, y_order):
    """Literal two-input contraction in F2[x,y]/(x^8,y^(y_order+1))."""
    x, y, u, v = labels
    mask = sum(255 << (STRIDE * b) for b in range(y_order + 1))
    xx, uu = code[x]['output'], code[u]['output']
    result = 0
    for yy in code[y]['input']:
        for vv in code[v]['input']:
            value = 1
            for site in range(10):
                index = (((xx >> (2 * site) & 3) * 4 + (yy >> (2 * site) & 3)) * 4
                         + (uu >> (2 * site) & 3)) * 4 + (vv >> (2 * site) & 3)
                value = multiply(value, local[index]) & mask
                if not value:
                    break
            result ^= value
    return result
