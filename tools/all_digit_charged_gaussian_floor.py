#!/usr/bin/env python3
"""Exact interval bounds for the derivative of Gaussian pressure in digit count.

The kernel is h(x)=x/(exp(x)-1), extended continuously by h(0)=1.
The underlying exponential bounds are the repository's integer-outward
Taylor implementation. No floating-point operation proves a sign.
"""

from functools import lru_cache

from charged_gaussian_floor import Dyadic


class DigitDerivativeBounds(Dyadic):
    """Use h in the existing finite and periodic energy expressions."""

    @lru_cache(maxsize=262144)
    def softplus_point(self, numerator):
        """Bound h at one exactly represented dyadic argument.

        The inherited ``softplus`` method uses monotonicity to evaluate
        interval endpoints. This override deliberately changes its kernel.
        """
        if numerator == 0:
            return self.one
        absolute = (abs(numerator), abs(numerator))
        exponential = self.exp_negative(absolute)
        denominator = self.subtract(self.one, exponential)
        assert denominator[0] > 0, "Increase precision for this near-zero point."
        thermal = self.divide(self.multiply(absolute, exponential), denominator)
        # h(-x)=h(x)+x, including the complete linear part.
        return (thermal[0] - min(numerator, 0),
                thermal[1] - min(numerator, 0))

    def continuation(self, scale, charge, terms=12):
        """A rigorous lower bound for A(J,t)+A(0,t)-2sqrt(2)A(J/2,sqrt(2)t)."""
        assert 0 <= charge[0] <= charge[1] <= 2 * self.scale
        positive = self.add(self.paired_sum(scale, charge, terms),
                            self.paired_sum(scale, self.zero, terms))
        physical_scale = self.multiply(self.root_two, scale)
        half_charge = self.divide_integer(charge, 2)
        negative = self.multiply(
            self.multiply(self.rational(2), self.root_two),
            self.paired_sum(physical_scale, half_charge, terms))

        # Both missing physical tails start at least at ``first`` and
        # have spacing ``step``. For x>=0, h(x)<=(x+1)exp(-x).
        first = self.multiply(physical_scale,
                              self.subtract(self.rational(2 * terms + 1), half_charge))
        step = self.multiply(self.rational(2), physical_scale)
        ratio = self.exp_negative(step)
        denominator = self.subtract(self.one, ratio)
        geometric = self.add(
            self.divide(self.add(first, self.one), denominator),
            self.divide(self.multiply(step, ratio),
                        self.multiply(denominator, denominator)))
        tail = self.multiply(
            self.multiply(self.rational(4), self.root_two),
            self.multiply(self.exp_negative(first), geometric))
        # Missing positive tails are discarded only in this lower bound.
        return self.subtract(positive, negative)[0] - tail[1]
