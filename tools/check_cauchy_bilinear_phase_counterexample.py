#!/usr/bin/env python3
"""Exact 4-by-5 counterexample to omitting the Cauchy degree condition."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_frobenius_rational_trace import field_of_degree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    field = field_of_degree(4)
    omega = next(x for x in range(2, field.size) if field.mul(x,x)^x^1 == 0)
    omega2 = field.mul(omega, omega)
    subfield = [0, 1, omega, omega2]
    torus = [v for v in range(1, field.size) if field.power(v, 5) == 1]
    assert len(torus) == 5
    iota = lambda x: field.inverse(x) if x else 0

    def value(coefficients, u):
        total, power = 0, 1
        for coefficient in coefficients:
            total ^= field.mul(coefficient, power)
            power = field.mul(power, u)
        return total

    actual_total = naive_total = defect_total = 0
    rows, pole_points, phase_controls = [], 0, 0
    for v in torus:
        first = [1, field.mul(1^v, omega2), field.mul(1^v, omega)]
        second = [1, 1^field.mul(v,v), 1^field.mul(v,v)]
        coefficients = [0]*5
        for i, a in enumerate(first):
            for j, b in enumerate(second):
                coefficients[i+j] ^= field.mul(a,b)
        leading = field.mul(omega, field.mul(1^v, 1^field.mul(v,v)))
        assert coefficients[4] == leading
        reduced = coefficients[:4]
        reduced[1] ^= leading
        actual = 0
        for encoded_integer, u in enumerate(subfield):
            phase = field.power(v, encoded_integer)
            assert value(coefficients, u) == phase == value(reduced, u)
            denominator = 1 ^ field.mul(u,v)
            pole_points += int(denominator == 0)
            actual ^= field.mul(phase, iota(denominator))
            phase_controls += 1
        inverse_v = field.inverse(v)
        inverse_denominator = iota(1^field.mul(v,v))
        naive = field.mul(value(coefficients, inverse_v), inverse_denominator)
        corrected = field.mul(value(reduced, inverse_v), inverse_denominator)
        if v == 1:
            derivative = reduced[1] ^ reduced[3]
            corrected ^= derivative
            assert coefficients == [1,0,0,0,0] and derivative == 0
        defect = field.mul(leading, inverse_v)
        assert corrected == actual and actual^naive == defect
        assert bool(leading) == (v != 1)
        rows.append({'v':v,'phase_coefficients_in_U':coefficients,
            'reduced_coefficients_in_U':reduced,'actual_inner_sum':actual,
            'naive_inner_sum':naive,'corrected_inner_sum':corrected,'missing_term':defect})
        actual_total ^= actual
        naive_total ^= naive
        defect_total ^= defect
    assert phase_controls == 20 and pole_points == 1
    assert actual_total == 1 and naive_total == omega2 and defect_total == omega
    assert actual_total != naive_total
    report = {'schema_version':1,
        'scope':'Exact failure of the bounded-degree interpolation formula when an encoded bilinear Fourier phase violates its formal degree hypothesis. The proved bounded-degree algorithm is unaffected.',
        'field_polynomial':field.polynomial,'omega':omega,'omega_squared':omega2,
        'subfield_by_encoded_integer':subfield,'torus':torus,'rows':rows,
        'phase_encoding_controls':phase_controls,'zero_denominator_points':pole_points,
        'actual_joint_sum':actual_total,'naive_joint_sum':naive_total,'missing_joint_sum':defect_total,
        'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                         for name in ['check_cauchy_bilinear_phase_counterexample.py',
                                      'check_frobenius_rational_trace.py','check_carry_fourier_elimination.py']}}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ['phase_encoding_controls','zero_denominator_points',
        'actual_joint_sum','naive_joint_sum','missing_joint_sum']}))


if __name__ == '__main__':
    main()
