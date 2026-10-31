"""Exact variable-cap smooth-square search with fully charged box cost.

The search enumerates two balanced rectangles in the exponent box.  It is
deterministic and costs O(sqrt(product(2*t_i+1))) times polynomial factors;
that product can be exponential in the arithmetic input parameters.  It
does not implement a polynomial-time oracle for the universal-square gate.
"""

from bisect import bisect_left
from math import gcd, isqrt, prod

from quadratic_form_arithmetic import (
    compose_primitive, form_discriminant, inverse_primitive,
    power_primitive, reduce_primitive,
)
from universal_smooth_square import representation_from_square_word


def balanced_rectangles(radices):
    """Partition a mixed-radix box into at most two balanced rectangles.

    A variable (index, step, count) contributes step*v to its original
    coordinate, for 0<=v<count.  The two halves may share the split
    coordinate, but their choices are independent inside each rectangle.
    """
    radices = tuple(radices)
    if any(radix < 1 for radix in radices):
        raise ValueError("positive radices are required")
    total = prod(radices)
    root = isqrt(total)
    prefix, pivot = 1, None
    for index, radix in enumerate(radices):
        if prefix*radix > root:
            pivot = index
            break
        prefix *= radix
    if pivot is None:
        assert total == 1
        return [{"left": [], "right": [], "offset": [0]*len(radices),
                 "left_domain": 1, "right_domain": 1}]
    low_count = max(1, root//prefix)
    full_count, remainder = divmod(radices[pivot], low_count)
    left_prefix = [(index, 1, radices[index]) for index in range(pivot) if radices[index] > 1]
    right_suffix = [(index, 1, radices[index]) for index in range(pivot+1, len(radices))
                    if radices[index] > 1]
    rectangles = []
    for high_count, low, offset_value in ((full_count, low_count, 0),
                                           (1, remainder, full_count*low_count)):
        if not low:
            continue
        left = left_prefix+([(pivot, 1, low)] if low > 1 else [])
        right = ([(pivot, low_count, high_count)] if high_count > 1 else [])+right_suffix
        offset = [0]*len(radices)
        offset[pivot] = offset_value
        left_domain = prod(variable[2] for variable in left)
        right_domain = prod(variable[2] for variable in right)
        assert left_domain <= root and right_domain <= 2*(root+1)
        rectangles.append({"left": left, "right": right, "offset": offset,
                           "left_domain": left_domain, "right_domain": right_domain})
    assert sum(rectangle["left_domain"]*rectangle["right_domain"]
               for rectangle in rectangles) == total
    return rectangles


def _deduplicate(entries):
    entries.sort()
    result = []
    for entry in entries:
        if not result or entry[0] != result[-1][0]:
            result.append(entry)
    return result


class SquareBoxSearch:
    """Reusable exact search index for one discriminant and cap vector."""

    def __init__(self, discriminant, primes, roots, caps):
        self.primes, self.caps = tuple(primes), tuple(caps)
        self.roots = dict(roots)
        self.discriminant = discriminant
        if (discriminant >= 0 or discriminant % 4 not in (0, 1)
                or len(self.primes) != len(self.caps)
                or len(set(self.primes)) != len(self.primes)
                or any(cap < 0 for cap in self.caps)
                or any(prime < 2 or gcd(prime, discriminant) != 1 for prime in self.primes)):
            raise ValueError("a negative order discriminant and distinct unramified primes with nonnegative caps are required")
        parity = discriminant % 2
        self.identity = (1, parity, (parity-discriminant)//4)
        self.generators = []
        for prime in self.primes:
            if prime not in self.roots or (self.roots[prime]**2-discriminant) % (4*prime):
                raise ValueError("a valid splitting root is required for each prime")
            middle = self.roots[prime]
            self.generators.append(reduce_primitive((prime, middle,
                                                     (middle*middle-discriminant)//(4*prime))))
        self.squares = [power_primitive(generator, 2) for generator in self.generators]
        self.box_size = prod(2*cap+1 for cap in self.caps)
        self.index_group_multiplications = 0
        self.shift = self.identity
        for square, cap in zip(self.squares, self.caps):
            self.shift = self._index_multiply(self.shift, power_primitive(square, cap))
        self.rectangles = []
        for specification in balanced_rectangles([2*cap+1 for cap in self.caps]):
            left = self._enumerate_half(specification["left"])
            right = self._enumerate_half(specification["right"])
            offset_class = self.identity
            for square, exponent in zip(self.squares, specification["offset"]):
                offset_class = self._index_multiply(offset_class, power_primitive(square, exponent))
            self.rectangles.append({**specification, "left_entries": left,
                                    "left_keys": [entry[0] for entry in left],
                                    "right_entries": [(inverse_primitive(value), word)
                                                      for value, word in right],
                                    "inverse_offset_class": inverse_primitive(offset_class)})

    def _index_multiply(self, left, right):
        if left == self.identity:
            return right
        if right == self.identity:
            return left
        self.index_group_multiplications += 1
        return compose_primitive(left, right)

    def _enumerate_half(self, variables):
        states = [(self.identity, (0,)*len(self.primes))]
        for index, step, count in variables:
            base = power_primitive(self.squares[index], step)
            powers, value = [], self.identity
            for exponent in range(count):
                powers.append((value, step*exponent))
                if exponent+1 < count:
                    value = self._index_multiply(value, base)
            enlarged = []
            for state, word in states:
                for value, coefficient in powers:
                    new_word = list(word)
                    new_word[index] += coefficient
                    enlarged.append((self._index_multiply(state, value), tuple(new_word)))
            states = _deduplicate(enlarged)
        return states

    def statistics(self):
        return {"box_size": self.box_size,
                "multiplication_count_scope": "Index/lookup loops; polynomial binary-power setup and coordinate construction are excluded",
                "index_group_multiplications": self.index_group_multiplications,
                "rectangles": [{"left_domain": item["left_domain"],
                                "right_domain": item["right_domain"],
                                "left_classes": len(item["left_entries"]),
                                "right_classes": len(item["right_entries"]),
                                "offset": item["offset"]}
                               for item in self.rectangles]}

    def solve(self, form, with_representation=True):
        """Return an exact witness or NO; all rectangle searches are complete."""
        if form_discriminant(form) != self.discriminant:
            raise ValueError("the target has the wrong discriminant")
        target = reduce_primitive(form)
        multiplications = 0

        def multiply(left, right):
            nonlocal multiplications
            if left == self.identity:
                return right
            if right == self.identity:
                return left
            multiplications += 1
            return compose_primitive(left, right)

        shifted = multiply(target, self.shift)
        for rectangle_index, rectangle in enumerate(self.rectangles):
            required = multiply(shifted, rectangle["inverse_offset_class"])
            for right_inverse, right_word in rectangle["right_entries"]:
                needed = multiply(required, right_inverse)
                position = bisect_left(rectangle["left_keys"], needed)
                if position == len(rectangle["left_keys"]) or rectangle["left_keys"][position] != needed:
                    continue
                left_word = rectangle["left_entries"][position][1]
                word = [left+right+offset-cap for left, right, offset, cap
                        in zip(left_word, right_word, rectangle["offset"], self.caps)]
                assert all(abs(exponent) <= cap for exponent, cap in zip(word, self.caps))
                value = self.identity
                for square, exponent in zip(self.squares, word):
                    value = multiply(value, power_primitive(square, exponent))
                assert value == target
                answer = {"represented": True, "half_word": word,
                          "rectangle_index": rectangle_index,
                          "query_group_multiplications": multiplications}
                if with_representation:
                    base = representation_from_square_word(form, self.primes, self.roots,
                                                           word, max((1,)+self.caps))
                    content = prod(prime**(cap-abs(exponent))
                                   for prime, cap, exponent in zip(self.primes, self.caps, word))
                    x, y = (content*coordinate for coordinate in base["primitive_vector"])
                    root_norm = prod(prime**cap for prime, cap in zip(self.primes, self.caps))
                    a, b, c = form
                    assert a*x*x+b*x*y+c*y*y == root_norm*root_norm
                    assert gcd(x, y) == content
                    answer["representation"] = {"vector": (x, y), "content": content,
                                                "primitive_vector": base["primitive_vector"],
                                                "primitive_norm": base["primitive_norm"],
                                                "root_norm": root_norm,
                                                "fixed_norm": root_norm*root_norm}
                return answer
        return {"represented": False, "half_word": None,
                "query_group_multiplications": multiplications}
