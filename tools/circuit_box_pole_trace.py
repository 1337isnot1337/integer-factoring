#!/usr/bin/env python3
"""Complete normalized pole traces for one-variable regular circuits.

Compute 2^K sum 1/(N+2^K*w-q(i)), where w^2+w+1=0. A finite-precision
function quotient removes the input's formal degree; a degree-budgeted
dyadic tree removes the singularity. The input interval is never listed.
This theorem is univariate and does not supply the bivariate factor-box
oracle for q(a,b)=ab.
"""

from __future__ import annotations

from math import comb

from circuit_box_averages import normalized_box_average, validate_regular_circuit
from circuit_box_bezout import Builder, valuation_two
from circuit_box_products import IntegerFunctionQuotient


def _trim(coefficients):
    result = list(coefficients)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def pointwise_circuit_polynomial(circuit, target, precision):
    """A degree below 2p representative of q(X)-target at every integer."""
    numerator_degree, denominator_degree = validate_regular_circuit(circuit, 1)
    ring = IntegerFunctionQuotient(1, precision)
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = ring.variable(0)
        elif operation == "const":
            value = ring.constant(instruction[1])
        elif operation == "add":
            value = ring.add(values[instruction[1]], values[instruction[2]])
        elif operation == "mul":
            value = ring.multiply(values[instruction[1]], values[instruction[2]])
        else:
            assert operation == "inv"
            # Euler powering is used only for pointwise residues here.
            # Every inverse input is odd at both integer parities.
            base = values[instruction[1]]
            exponent = (1 << (precision - 1)) - 1
            value = ring.constant(1)
            while exponent:
                if exponent & 1:
                    value = ring.multiply(value, base)
                exponent >>= 1
                if exponent:
                    base = ring.multiply(base, base)
        values.append(value)
    reduced = ring.add(values[circuit["output"]], ring.constant(-target))
    coefficients = [reduced.get((degree,), 0) for degree in range(ring.degree)]
    return {"coefficients": _trim(coefficients), "quotient_degree": ring.degree,
            "formal_numerator_degree_bound": numerator_degree,
            "formal_denominator_degree_bound": denominator_degree}


def _parity_substitute(coefficients, parity, modulus):
    """Coefficients of f(parity+2X), reduced before coefficient growth."""
    if parity == 0:
        return _trim([(value << degree) % modulus
                      for degree, value in enumerate(coefficients)])
    return _trim([(sum(coefficients[index] * comb(index, degree)
                       for index in range(degree, len(coefficients))) << degree) % modulus
                  for degree in range(len(coefficients))])


def _residue_interval(interval, residue, exponent):
    step = 1 << exponent
    lower, upper = interval
    return (lower - residue + step - 1) // step, (upper - residue + step - 1) // step


def _root_multiplicity_mod_two(coefficients, parity):
    polynomial = sum((value & 1) << degree for degree, value in enumerate(coefficients))
    assert polynomial
    if parity == 0:
        return valuation_two(polynomial)
    multiplicity = 0
    while polynomial.bit_count() % 2 == 0:
        degree = polynomial.bit_length() - 1
        quotient = 0
        remainder = polynomial
        for index in range(degree - 1, -1, -1):
            if remainder & (1 << (index + 1)):
                quotient ^= 1 << index
                remainder ^= 3 << index
        assert remainder == 0 and quotient
        polynomial = quotient
        multiplicity += 1
    return multiplicity


def dyadic_valuation_partition(coefficients, cutoff, precision, interval):
    """Partition an interval into constant valuations below K and a >=K part.

    The represented polynomial is known modulo 2^precision, with
    precision>cutoff. A leaf has i=a+2^t*y and q(i)-N=2^v*u(y), at the
    available precision. High leaves use v=cutoff, with no unit promise.
    """
    if cutoff < 0 or precision <= cutoff:
        raise ValueError("Require 0 <= cutoff < coefficient precision")
    if len(interval) != 2 or interval[0] > interval[1]:
        raise ValueError("Require an ordered half-open interval")
    coefficients = _trim([value % (1 << precision) for value in coefficients])
    initial_budget = max(1, len(coefficients) - 1)
    pending = [(0, 0, 0, coefficients, initial_budget)]
    active_budget = initial_budget
    maximum_budget = initial_budget
    processed = 0
    leaves = []

    def append_leaf(residue, exponent, valuation, polynomial, high):
        local_interval = _residue_interval(interval, residue, exponent)
        count = local_interval[1] - local_interval[0]
        if count:
            leaves.append({"residue": residue, "residue_exponent": exponent,
                           "valuation": valuation, "high_valuation": high,
                           "coefficients": _trim(polynomial),
                           "local_interval": local_interval, "count": count})

    while pending:
        residue, exponent, offset, polynomial, budget = pending.pop()
        active_budget -= budget
        local_interval = _residue_interval(interval, residue, exponent)
        if local_interval[0] == local_interval[1]:
            continue
        processed += 1
        remaining = precision - offset
        content = min((valuation_two(value) for value in polynomial if value), default=remaining)
        if offset + content >= cutoff:
            division = cutoff - offset
            assert division >= 0
            append_leaf(residue, exponent, cutoff,
                        [value >> division for value in polynomial], True)
            continue
        polynomial = [value >> content for value in polynomial]
        offset += content
        remaining -= content
        degree_mod_two = max(index for index, value in enumerate(polynomial) if value & 1)
        assert degree_mod_two <= budget
        roots = [parity for parity in (0, 1)
                 if (polynomial[0] if parity == 0 else sum(polynomial)) % 2 == 0]
        if not roots:
            append_leaf(residue, exponent, offset, polynomial, False)
            continue
        children_budget = 0
        for parity in (0, 1):
            child = _parity_substitute(polynomial, parity, 1 << remaining)
            child_residue = residue + (parity << exponent)
            if parity in roots:
                child_budget = _root_multiplicity_mod_two(polynomial, parity)
                assert child_budget >= 1
                children_budget += child_budget
                pending.append((child_residue, exponent + 1, offset, child, child_budget))
                active_budget += child_budget
            else:
                append_leaf(child_residue, exponent + 1, offset, child, False)
        assert children_budget <= degree_mod_two
        maximum_budget = max(maximum_budget, active_budget)
        assert maximum_budget <= initial_budget
    assert processed <= initial_budget * (cutoff + 1)
    assert sum(leaf["count"] for leaf in leaves) == interval[1] - interval[0]
    assert sum(leaf["high_valuation"] for leaf in leaves) <= initial_budget
    return {"leaves": leaves, "processed_states": processed,
            "maximum_active_degree_budget": maximum_budget,
            "initial_degree_budget": initial_budget,
            "processed_state_bound": initial_budget * (cutoff + 1)}


def _leaf_resolvent_circuits(coefficients, scale):
    builder = Builder()
    polynomial = builder.constant(coefficients[-1])
    for coefficient in reversed(coefficients[:-1]):
        polynomial = builder.operation("add", builder.operation("mul", polynomial, 0),
                                       builder.constant(coefficient))
    shifted = builder.operation("add", polynomial, builder.constant(scale))
    norm = builder.operation("add", builder.operation("mul", polynomial, shifted),
                             builder.constant(scale * scale))
    builder.nodes.append(["inv", norm])
    inverse = len(builder.nodes) - 1
    first = builder.scale(builder.operation("mul", shifted, inverse), -scale)
    second = builder.scale(inverse, -scale * scale)
    return ({"nodes": builder.nodes, "output": first},
            {"nodes": builder.nodes, "output": second})


def _sum_regular(circuit, interval, precision):
    count = interval[1] - interval[0]
    if not count:
        return 0
    average = normalized_box_average(circuit, (interval,), precision + 1)
    doubled = count * average["scaled_average_residue"] % (1 << (precision + 1))
    assert doubled % 2 == 0
    return doubled >> 1


def normalized_circuit_pole_trace(circuit, interval, target, scale_exponent, precision):
    """The complete 2^K sum 1/(N+2^K*w-q(i)) modulo 2^p, in basis 1,w."""
    if scale_exponent < 0 or precision < 1:
        raise ValueError("Require nonnegative pole scale and positive precision")
    interval = tuple(interval)
    if len(interval) != 2 or any(not isinstance(value, int) for value in interval) or interval[0] > interval[1]:
        raise ValueError("Require an integer half-open interval")
    construction_precision = scale_exponent + precision
    compressed = pointwise_circuit_polynomial(circuit, target, construction_precision)
    partition = dyadic_valuation_partition(compressed["coefficients"], scale_exponent,
                                            construction_precision, interval)
    modulus = 1 << precision
    answer = [0, 0]
    valuation_counts = {}
    roots = []
    evaluated = 0
    for leaf in partition["leaves"]:
        valuation = leaf["valuation"]
        valuation_counts[str(valuation)] = valuation_counts.get(str(valuation), 0) + leaf["count"]
        if leaf["high_valuation"]:
            first = leaf["residue"] + (leaf["local_interval"][0] << leaf["residue_exponent"])
            roots.append({"residue": leaf["residue"], "residue_exponent": leaf["residue_exponent"],
                          "count": leaf["count"], "first_index": first})
        shift = scale_exponent - valuation
        if shift >= precision:
            continue
        scale = 1 << shift
        circuits = _leaf_resolvent_circuits(leaf["coefficients"], scale)
        for component in (0, 1):
            answer[component] = (answer[component]
                                 + _sum_regular(circuits[component], leaf["local_interval"], precision)) % modulus
        evaluated += 1
    return {"interval": list(interval), "target": target,
            "scale_exponent": scale_exponent, "precision": precision,
            "normalized_trace": answer,
            "congruence_root_count": sum(root["count"] for root in roots),
            "root_balls": roots, "valuation_counts_capped_at_scale": valuation_counts,
            "input_circuit_nodes": len(circuit["nodes"]),
            "formal_numerator_degree_bound": compressed["formal_numerator_degree_bound"],
            "formal_denominator_degree_bound": compressed["formal_denominator_degree_bound"],
            "quotient_degree": compressed["quotient_degree"],
            "compressed_polynomial_degree": len(compressed["coefficients"]) - 1,
            "construction_precision": construction_precision,
            "maximum_average_precision": precision + 1 if evaluated else 0,
            "processed_states": partition["processed_states"],
            "processed_state_bound": partition["processed_state_bound"],
            "maximum_active_degree_budget": partition["maximum_active_degree_budget"],
            "valuation_leaves": len(partition["leaves"]),
            "evaluated_unit_norm_leaves": evaluated,
            "scope": "Complete singular normalized pole trace for one-variable regular circuits; no bivariate factor-box oracle."}
