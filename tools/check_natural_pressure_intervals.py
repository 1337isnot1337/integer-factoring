#!/usr/bin/env python3
"""Independent rational-power and alternating-series pressure controls."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from certify_natural_pressure_interval import certify, load_profile, log_interval, power_interval


def alternating_log(value, terms=300):
    """Independent log(1+u) expansion for 0 <= u <= 1/2."""
    u = value - 1
    assert 0 <= u <= Fraction(1, 2) and terms % 2 == 0
    total, power = Fraction(0), Fraction(1)
    for j in range(1, terms + 1):
        power *= u
        total += power / j if j % 2 else -power / j
    return total, total + power * u / (terms + 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    counts = {'exact_power_enclosures': 0, 'independent_log_enclosures': 0,
              'exact_pressure_signs': 0, 'published_integer_bound_controls': 0,
              'rejection_controls': 0}
    for denominator in range(2, 12):
        for numerator in range(1, denominator + 1):
            base = Fraction(numerator, denominator)
            for exponent in [0, 1, 2, 3, 7, 16, 31]:
                exact = base ** exponent
                for bits in [8, 20]:
                    lo, hi = power_interval(base, exponent, bits)
                    assert Fraction(lo, 1 << bits) <= exact <= Fraction(hi, 1 << bits)
                    counts['exact_power_enclosures'] += 1
    assert power_interval(Fraction(1, 2), 1 << 80, 64) == (0, 1)
    counts['exact_power_enclosures'] += 1
    factors = [Fraction(3, 2), Fraction(4, 3), Fraction(5, 4)]
    reference = [alternating_log(x) for x in factors]
    for powers in product(range(-2, 3), repeat=3):
        value = Fraction(1)
        reference_lo = reference_hi = Fraction(0)
        for factor, (lo, hi), exponent in zip(factors, reference, powers):
            value *= factor ** exponent
            reference_lo += exponent * (lo if exponent >= 0 else hi)
            reference_hi += exponent * (hi if exponent >= 0 else lo)
        lo, hi = log_interval(value, 96, 64)
        if value == 1:
            assert lo == hi == 0 and reference_lo <= 0 <= reference_hi
        else:
            assert Fraction(lo, 1 << 96) <= reference_lo <= reference_hi <= Fraction(hi, 1 << 96)
        counts['independent_log_enclosures'] += 1
    examples = [([0, 1, 2, 4], [0, 1]),
                ([0, 2, 3, 4, 7, 9, 10, 11, 13], [0, 2, 4])]
    for (orders, energies), base, rho in product(
            examples, [Fraction(1, 2), Fraction(2, 3), Fraction(4, 5)],
            [Fraction(1, 3), Fraction(2, 3), Fraction(3, 4), Fraction(1)]):
        support = 2 * rho
        a, b = support.numerator, support.denominator
        p = sum((base ** (a * x) for x in orders), Fraction(0))
        q = sum((base ** (b * x) for x in energies), Fraction(0))
        exact = q ** a - p ** b
        result = certify(orders, energies, rho, base, bits=128, terms=72)
        lower, upper = result['a_log_q_minus_b_log_p_dyadic_interval']
        assert result['passed'] == (exact > 0)
        assert (upper < 0 if exact < 0 else (lower <= 0 <= upper if exact == 0 else lower > 0))
        counts['exact_pressure_signs'] += 1
    source_inputs = {}
    for path, rho, base in [
        ('evidence/hot64-charge-rational-pressure.json', Fraction(32, 45), Fraction(459, 460)),
        ('evidence/hot96-charge-profile.json', Fraction(49, 69), Fraction(768, 769)),
    ]:
        data = load_profile(Path(path))
        assert certify(data['orders_w'], data['energies_x'], rho, base)['passed']
        source_inputs[path] = sha256(Path(path).read_bytes()).hexdigest()
        counts['published_integer_bound_controls'] += 1
    bad_calls = [lambda: power_interval(Fraction(0), 2, 32),
                 lambda: power_interval(Fraction(2), 2, 32),
                 lambda: power_interval(Fraction(1, 2), -1, 32),
                 lambda: power_interval(Fraction(1, 2), True, 32),
                 lambda: log_interval(Fraction(0), 32, 20),
                 lambda: log_interval(Fraction(2), 32, 0),
                 lambda: certify([0, -1], [0, 1], Fraction(1), Fraction(1, 2))]
    for call in bad_calls:
        try:
            call()
        except ValueError:
            counts['rejection_controls'] += 1
        else:
            raise AssertionError('Invalid input was accepted.')
    paths = [Path(__file__), Path(__file__).with_name('certify_natural_pressure_interval.py')]
    record = {'passed': True, 'counts': counts, 'input_sha256': source_inputs,
              'source_sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
              'scope': 'Exact powers and independent alternating-series controls; '
                       'large polynomial profiles are loaded, not reconstructed.'}
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record))


if __name__ == '__main__':
    main()
