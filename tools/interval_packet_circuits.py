#!/usr/bin/env python3
"""Small global circuits for consecutive products modulo a power of two."""

from interval_unit_jets import prefix_power_sums


def paired_elementaries(lower, pairs, degree):
    """Elementary functions of (lower+t)*(lower+2*pairs-1-t)/2."""
    sums = prefix_power_sums(pairs, 2 * degree)
    quadratic = (lower * (lower + 2 * pairs - 1), 2 * pairs - 1, -1)
    power = (1,)
    moments = [pairs]
    for j in range(1, degree + 1):
        result = [0] * (len(power) + 2)
        for i, a in enumerate(power):
            for k, b in enumerate(quadratic):
                result[i + k] += a * b
        power = tuple(result)
        numerator = sum(a * sums[i] for i, a in enumerate(power))
        assert numerator % (1 << j) == 0
        moments.append(numerator >> j)
    elementary = [1]
    for k in range(1, degree + 1):
        numerator = sum((-1) ** (j - 1) * elementary[k - j] * moments[j]
                        for j in range(1, k + 1))
        assert numerator % k == 0
        elementary.append(numerator // k)
    return tuple(elementary)


def interval_product_circuit(interval, precision, number=None):
    """Emit an ordinary +,* circuit for P_I or Q_(number,I) mod 2^precision."""
    lower, upper = interval
    assert 1 <= lower <= upper and precision >= 1
    modulus = 1 << precision
    nodes = [['x']]

    def constant(value):
        nodes.append(['const', value % modulus])
        return len(nodes) - 1

    def binary(operation, left, right):
        nodes.append([operation, left, right])
        return len(nodes) - 1

    def scale(node, value):
        return binary('mul', node, constant(value))

    def power(node, exponent):
        result = constant(1)
        while exponent:
            if exponent % 2:
                result = binary('mul', result, node)
            exponent //= 2
            if exponent:
                node = binary('mul', node, node)
        return result

    length = upper - lower
    pairs = length // 2
    result = constant(1)
    if pairs:
        center = 2 * lower + 2 * pairs - 1
        degree = min(pairs, precision - 1)
        elementary = paired_elementaries(lower, pairs, degree)
        x_squared = binary('mul', 0, 0)
        if number is None:
            base = binary('add', x_squared, scale(0, -center))
        else:
            base = binary('add', constant(number * number), scale(0, -center * number))
        series, x_power = constant(1), constant(1)
        for j in range(1, degree + 1):
            coefficient = (1 << j) * elementary[j]
            if number is None:
                term = constant(coefficient)
            else:
                x_power = binary('mul', x_power, x_squared)
                term = scale(x_power, coefficient)
            series = binary('add', binary('mul', series, base), term)
        result = binary('mul', power(base, pairs - degree), series)
    if length % 2:
        last = upper - 1
        factor = (binary('add', 0, constant(-last)) if number is None
                  else binary('add', constant(number), scale(0, -last)))
        result = binary('mul', result, factor)
    return {'interval': list(interval), 'number': number, 'precision': precision,
            'input_degree': length, 'nodes': nodes, 'output': result}


def evaluate_circuit(circuit, ring, argument):
    values = []
    for instruction in circuit['nodes']:
        operation = instruction[0]
        if operation == 'x':
            value = argument
        elif operation == 'const':
            value = ring.constant(instruction[1])
        elif operation == 'add':
            value = ring.add(values[instruction[1]], values[instruction[2]])
        else:
            assert operation == 'mul'
            value = ring.multiply(values[instruction[1]], values[instruction[2]])
        values.append(value)
    return values[circuit['output']]


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lower', type=int, required=True)
    parser.add_argument('--upper', type=int, required=True)
    parser.add_argument('--precision', type=int, required=True)
    parser.add_argument('--number', type=int)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    circuit = interval_product_circuit((args.lower, args.upper), args.precision, args.number)
    Path(args.output).write_text(json.dumps(circuit, indent=2) + '\n')
    print(f"Wrote {len(circuit['nodes'])} circuit nodes for degree {circuit['input_degree']}.")
