#!/usr/bin/env python3
"""Global modular product circuits and unit traces for bilinear interval images.

The represented polynomial is prod_(u in I,v in J) (Z-q(u,v)), where
q(u,v)=A*u*v+B*u+C*v+D. Generation does not enumerate either interval.
The unit-trace routine requires its denominators Z and Z-1 to be units
for every nonempty parity class. It does not evaluate a normalized pole
at an integer center.
"""

from __future__ import annotations

from math import comb


def prefix_power_sums(stop: int, degree: int) -> tuple[int, ...]:
    """Exact sums of t^j, 0<=t<stop, for 0<=j<=degree."""
    if stop < 0 or degree < 0:
        raise ValueError("Require nonnegative endpoint and degree")
    result = [stop]
    for j in range(1, degree + 1):
        numerator = stop ** (j + 1) - sum(
            comb(j + 1, i) * result[i] for i in range(j)
        )
        assert numerator % (j + 1) == 0
        result.append(numerator // (j + 1))
    return tuple(result)


def interval_power_sums(interval: tuple[int, int], degree: int) -> tuple[int, ...]:
    lower, upper = interval
    if lower < 0 or upper < lower:
        raise ValueError("Require a nonnegative half-open interval")
    low = prefix_power_sums(lower, degree)
    high = prefix_power_sums(upper, degree)
    return tuple(b - a for a, b in zip(low, high))


def _parity_section(interval: tuple[int, int], parity: int) -> tuple[int, int]:
    """u=parity+2*x maps the given interval to this ordinary x interval."""
    lower, upper = interval
    return ((lower - parity + 1) // 2, (upper - parity + 1) // 2)


def _box_moments(coefficients: tuple[int, int, int, int],
                 interval_u: tuple[int, int], interval_v: tuple[int, int],
                 degree: int) -> tuple[int, ...]:
    """Exact power sums of a bilinear polynomial on one box."""
    a, b, c, d = coefficients
    sums_u = interval_power_sums(interval_u, degree)
    sums_v = interval_power_sums(interval_v, degree)
    polynomial = {(0, 0): 1}
    monomials = ((1, 1, a), (1, 0, b), (0, 1, c), (0, 0, d))
    result = []
    for j in range(degree + 1):
        result.append(sum(value * sums_u[i] * sums_v[k]
                          for (i, k), value in polynomial.items()))
        if j == degree:
            break
        following: dict[tuple[int, int], int] = {}
        for (i, k), value in polynomial.items():
            for du, dv, coefficient in monomials:
                if coefficient:
                    key = (i + du, k + dv)
                    following[key] = following.get(key, 0) + value * coefficient
        polynomial = {key: value for key, value in following.items() if value}
    return tuple(result)


def bilinear_parity_moments(coefficients: tuple[int, int, int, int],
                            interval_u: tuple[int, int],
                            interval_v: tuple[int, int],
                            degree: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """For c=0,1 return sums of ((q(u,v)-c)/2)^j over q(u,v)=c mod 2.

    The zeroth moment is the cardinality of that parity class. The four
    input parity boxes are summed by polynomial power sums, never scanned.
    Signed polynomial coefficients and zero/empty intervals are allowed.
    """
    if len(coefficients) != 4 or degree < 0:
        raise ValueError("Require four coefficients and nonnegative degree")
    for lower, upper in (interval_u, interval_v):
        if lower < 0 or upper < lower:
            raise ValueError("Require nonnegative half-open intervals")
    a, b, c, d = coefficients
    moments = [[0] * (degree + 1), [0] * (degree + 1)]
    for e in (0, 1):
        section_u = _parity_section(interval_u, e)
        if section_u[0] == section_u[1]:
            continue
        for f in (0, 1):
            section_v = _parity_section(interval_v, f)
            if section_v[0] == section_v[1]:
                continue
            corner = a * e * f + b * e + c * f + d
            parity = corner % 2
            reduced = (2 * a, b + a * f, c + a * e, (corner - parity) // 2)
            box = _box_moments(reduced, section_u, section_v, degree)
            moments[parity] = [x + y for x, y in zip(moments[parity], box)]
    expected = (interval_u[1] - interval_u[0]) * (interval_v[1] - interval_v[0])
    assert moments[0][0] + moments[1][0] == expected
    return tuple(tuple(row) for row in moments)


def _elementaries(moments: tuple[int, ...], degree: int) -> tuple[int, ...]:
    result = [1]
    for k in range(1, degree + 1):
        numerator = sum((-1) ** (j - 1) * result[k - j] * moments[j]
                        for j in range(1, k + 1))
        assert numerator % k == 0
        result.append(numerator // k)
    return tuple(result)


def bilinear_product_circuit(coefficients: tuple[int, int, int, int],
                             interval_u: tuple[int, int],
                             interval_v: tuple[int, int],
                             precision: int) -> dict:
    """Write a +,* circuit for the whole product modulo 2^precision.

    Every constant in the node list is reduced modulo 2^precision.
    Exponents are expanded into binary-powering nodes. The node format
    matches interval_packet_circuits.evaluate_circuit.
    """
    if precision < 1:
        raise ValueError("Require positive precision")
    moments = bilinear_parity_moments(coefficients, interval_u, interval_v,
                                     precision - 1)
    modulus = 1 << precision
    nodes: list[list] = [["x"]]

    def constant(value):
        nodes.append(["const", value % modulus])
        return len(nodes) - 1

    def binary(operation, left, right):
        nodes.append([operation, left, right])
        return len(nodes) - 1

    def power(base, exponent):
        result = constant(1)
        while exponent:
            if exponent & 1:
                result = binary("mul", result, base)
            exponent >>= 1
            if exponent:
                base = binary("mul", base, base)
        return result

    result = constant(1)
    packets = []
    for parity in (0, 1):
        count = moments[parity][0]
        degree = min(count, precision - 1)
        elementary = _elementaries(moments[parity], degree)
        packet = [((-2) ** j * value) % modulus
                  for j, value in enumerate(elementary)]
        packets.append({"parity": parity, "count": count,
                        "retained_degree": degree, "coefficients": packet})
        if not count:
            continue
        base = 0 if parity == 0 else binary("add", 0, constant(-1))
        series = constant(1)
        for coefficient in packet[1:]:
            series = binary("add", binary("mul", series, base), constant(coefficient))
        factor = binary("mul", power(base, count - degree), series)
        result = binary("mul", result, factor)
    return {"coefficients": list(coefficients), "interval_u": list(interval_u),
            "interval_v": list(interval_v), "precision": precision,
            "input_degree": moments[0][0] + moments[1][0],
            "packets": packets, "nodes": nodes, "output": result}


def unit_pole_trace(coefficients, interval_u, interval_v, ring, argument):
    """Sum 1/(argument-q(u,v)) using only p parity moments.

    ``ring`` supplies precision, constant, add, scale, multiply and inverse.
    Every nonempty parity class requires argument-parity to be a unit.
    In particular argument=w+T is valid in quadratic 2-adic truncated jets.
    """
    moments = bilinear_parity_moments(coefficients, interval_u, interval_v,
                                     ring.precision - 1)
    result = ring.constant(0)
    for parity in (0, 1):
        if not moments[parity][0]:
            continue
        inverse = ring.inverse(ring.add(argument, ring.constant(-parity)))
        inverse_power = inverse
        for j, moment in enumerate(moments[parity]):
            result = ring.add(result, ring.scale(inverse_power, (1 << j) * moment))
            if j + 1 < ring.precision:
                inverse_power = ring.multiply(inverse_power, inverse)
    return result


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficients", nargs=4, type=int, default=[1, 0, 0, 0])
    parser.add_argument("--u", nargs=2, type=int, required=True)
    parser.add_argument("--v", nargs=2, type=int, required=True)
    parser.add_argument("--precision", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    circuit = bilinear_product_circuit(tuple(args.coefficients), tuple(args.u),
                                       tuple(args.v), args.precision)
    Path(args.output).write_text(json.dumps(circuit, indent=2) + "\n")
    print(f"Wrote {len(circuit['nodes'])} nodes for degree {circuit['input_degree']}.")
