"""Deterministic basis extraction in a black-box finite abelian 2-group.

Elements must have canonical hashable identifiers. Only the supplied group
operations are used; no coordinates, discrete-log oracle, or random oracle
are assumed. See research/arithmetic-torsion-extraction.md for proof/cost.
"""


def power(value, exponent, identity, multiply):
    result = identity
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


def order_log_two(value, identity, multiply, bound):
    current = value
    for exponent in range(bound+1):
        if current == identity:
            return exponent
        current = multiply(current, current)
    raise ValueError("element order exceeds the supplied 2-power bound")


def subset_products(values, identity, multiply):
    result = [identity]
    for value in values:
        result += [multiply(item, value) for item in result]
    return result


def order_two_coordinates(target, basis, identity, multiply):
    """Meet in the middle; the basis consists of independent order-two elements."""
    split = len(basis)//2
    left = {value: mask for mask, value in enumerate(subset_products(basis[:split], identity, multiply))}
    for mask, value in enumerate(subset_products(basis[split:], identity, multiply)):
        needed = multiply(target, value)  # order-two elements are self-inverse
        if needed in left:
            return left[needed] | mask << split
    return None


def two_group_basis(generators, identity, multiply, inverse, order_bit_bound):
    """Return (cyclic generator, log2(order)) pairs for the generated 2-group."""
    basis = []
    for value in generators:
        candidate = value
        candidate_order = order_log_two(candidate, identity, multiply, order_bit_bound)
        while candidate_order:
            tops = [power(element, 1 << (order-1), identity, multiply) for element, order in basis]
            top = power(candidate, 1 << (candidate_order-1), identity, multiply)
            coordinates = order_two_coordinates(top, tops, identity, multiply)
            if coordinates is None:
                basis.append((candidate, candidate_order))
                break
            involved = [j for j in range(len(basis)) if coordinates >> j & 1]
            assert involved  # a nonidentity order-two element cannot have empty coordinates
            pivot = min(involved, key=lambda j: basis[j][1])
            if basis[pivot][1] < candidate_order:
                # Keep the larger-order generator; the replaced one becomes the
                # candidate. The same binary dependency holds after this swap.
                old = basis[pivot]
                basis[pivot] = (candidate, candidate_order)
                candidate, candidate_order = old
            old_order = candidate_order
            for j in involved:
                element, order = basis[j]
                assert order >= candidate_order
                correction = power(element, 1 << (order-candidate_order), identity, multiply)
                candidate = multiply(candidate, inverse(correction))
            candidate_order = order_log_two(candidate, identity, multiply, order_bit_bound)
            assert candidate_order < old_order
    return basis


def project_to_two_part(generators, annihilator, identity, multiply):
    """An annihilator is enough; its factorization and exact group order are unused."""
    annihilator = abs(annihilator)
    if not annihilator:
        raise ValueError("annihilator must be nonzero")
    odd = annihilator
    while odd % 2 == 0:
        odd //= 2
    return [power(value, odd, identity, multiply) for value in generators]
