#!/usr/bin/env python3
"""Exact rational pressure certificate for the q115 height-gauge upper profile."""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from certify_coherent_hot_f4_hyperoval import certify, certify_small_type

Q115 = Path('evidence/f8-q115-height-gauge-crossing.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
RHO = Fraction(706733, 10**6)
BASE_RATIONAL = Fraction(2651072, 2651073)
TAIL_ORDER = 48


def main() -> None:
    q = json.loads(Q115.read_text())
    assert q['cutoff'] == 47 and q['tails'] == [12290, 12290]
    assert q['upper_below_critical'] is True
    zero = q['zero_orders_t'] + [TAIL_ORDER] * q['tails'][0]
    infinity = q['infinity_orders_t'] + [TAIL_ORDER] * q['tails'][1]
    profile = {
        'zero_orders_t': zero,
        'infinity_orders_t': infinity,
        'heights': q['custom_heights'],
    }
    assert len(zero) == len(infinity) == len(q['custom_heights'])**2 == 115**2
    assert 0 in zero and 0 in infinity and 0 in q['custom_heights']

    certificate = certify(profile, RHO, BASE_RATIONAL, bits=192, terms=96)
    higher = certify(profile, RHO, BASE_RATIONAL, bits=224, terms=112)
    negative = certify(profile, Fraction(706732, 10**6), BASE_RATIONAL, bits=192, terms=96)
    assert certificate['passed'] and higher['passed']
    assert not negative['passed']
    rational_type = certify_small_type(profile, certificate)
    assert rational_type['passed']

    sources = [
        'certify_f8_q115_height_gauge.py',
        'certify_coherent_hot_f4_hyperoval.py',
        'certify_natural_pressure_interval.py',
    ]
    out = {
        'completed': True,
        'q': 115,
        'cutoff': 47,
        'tail_orders_assigned': TAIL_ORDER,
        'unresolved_tails_per_cut': q['tails'],
        'strict_rho_bound': [RHO.numerator, RHO.denominator],
        'strict_rho_bound_decimal': float(RHO),
        'critical_rho': 1 / (2 ** 0.5),
        'margin': 1 / (2 ** 0.5) - float(RHO),
        'base_sha256': sha256(BASE.read_bytes()).hexdigest(),
        'q115_profile_sha256': sha256(Q115.read_bytes()).hexdigest(),
        'certificate': certificate,
        'higher_precision_certificate': higher,
        'negative_nearby_bound_control': negative,
        'small_rational_type': rational_type,
        'tail_augmented_profile_is_upper_bound_for_complete_orders': True,
        'monotonicity_used': 'Unresolved tail orders are at least 48. Replacing them by 48 makes every channel sum Z no smaller, so a positive exact gap for this profile certifies the complete q115 profile is below 706733/1000000.',
        'duplicate_height_caveat': 'The crossing uses w11=11, duplicating mode 5. The next odd distinct height 13 fails. This certificate does not establish that the weighted physical gauge or the uniform transfer is admissible.',
        'full_q115_profile_recomputed_this_run': False,
        'uniform_transfer_proved': False,
        'factoring_exponent_improved': False,
        'source_sha256': {name: sha256(Path('tools', name).read_bytes()).hexdigest() for name in sources},
        'scope': 'Exact arithmetic for the q115 height-gauge cutoff-47 tail-augmented two-cut pressure. The weighted map/receipt is verified by check_f8_q115_height_gauge_crossing.py. No complete profile, transfer, score, or factoring exponent is claimed.',
    }
    Path('evidence/f8-q115-height-gauge-pressure.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({
        'status': 'PASS',
        'strict_rho_bound': [RHO.numerator, RHO.denominator],
        'decimal_bound': float(RHO),
        'margin': 1 / (2 ** 0.5) - float(RHO),
        'gap_lower_integer': certificate['gap_4a_logQ_minus_b_logZ0Zinf'][0],
        'negative_control_failed': True,
        'type_denominator': rational_type['denominator'],
        'minimum_type_count': rational_type['minimum_count'],
        'duplicate_height_caveat': True,
        'factoring_exponent_improved': False,
    }))


if __name__ == '__main__':
    main()
