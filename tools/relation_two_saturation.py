"""Find a useful involution directly from independent integer relations.

The algorithm halves the relation-lattice index only when the proposed
half-relation still evaluates to the identity. It never computes a
determinant, a group order, or an annihilator. The generic contract
requires linearly independent genuine relations in an abelian group.
The quadratic-form wrapper checks strict diagonal dominance to establish
nonsingularity without an integer determinant computation.
"""

from math import gcd

from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from quadratic_ideal_hnf import compose_forms_hnf


def first_binary_dependencies(rows, limit=2):
    """Return the first independent dependencies among rows modulo two.

    Fewer than ``limit`` results means the entire kernel was exhausted.
    Each returned bitmask records the selected original row indices.
    """
    pivots, answer = {}, []
    for index, row in enumerate(rows):
        value = sum((entry & 1) << column for column, entry in enumerate(row))
        coefficients = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value, coefficients
                break
            old_value, old_coefficients = pivots[pivot]
            value ^= old_value
            coefficients ^= old_coefficients
        if not value:
            answer.append(coefficients)
            if len(answer) == limit:
                break
    return tuple(answer)


def half_relation(rows, dependency):
    """The integral vector one half of the selected relation-row sum."""
    if dependency <= 0 or dependency >= 1 << len(rows):
        raise ValueError("a nonzero dependency on the supplied rows is required")
    totals = [0] * len(rows[0])
    for index, row in enumerate(rows):
        if dependency >> index & 1:
            totals = [left + right for left, right in zip(totals, row)]
    if any(value & 1 for value in totals):
        raise ValueError("the row selection is not a binary dependency")
    return tuple(value // 2 for value in totals)


def first_useful_relation_top(generators, relation_matrix, identity,
                              multiply, inverse, inspect_top, *,
                              verify_relations=True, keep_updates=False):
    """Extract a useful order-two element from independent genuine relations.

    ``relation_matrix`` has at most as many rows as generators, and its
    rows are independent over the rationals. Every row evaluates to the
    identity on the supplied generators.
    The group is abelian. Among its nonidentity involutions,
    ``inspect_top`` may reject at most one by returning ``None``.

    The return value is ``(top, certificate, stats)``. A missing top
    certifies that the image of the presented quotient's torsion subgroup
    has no useful involution. For a square matrix this is the whole
    generated group. Generation of a larger group is not tested.
    Verification of the initial relation identities can be omitted only
    when those identities have already been established by the caller.
    """
    generators = tuple(generators)
    size = len(generators)
    rows = [list(row) for row in relation_matrix]
    if len(rows) > size or any(len(row) != size for row in rows):
        raise ValueError("independent relation rows cannot outnumber generators")
    if any(not isinstance(value, int) for row in rows for value in row):
        raise ValueError("relation entries must be integers")
    row_norms = [sum(abs(value) for value in row) for row in rows]
    if any(value == 0 for value in row_norms):
        raise ValueError("independent relation rows cannot include a zero row")
    # The positive gcd of maximal minors is at most their Hadamard bound,
    # the product of row l1 norms. Its binary valuation has this cap.
    saturation_bound = sum((value - 1).bit_length() for value in row_norms)
    stats = {
        "relation_dimension": size,
        "relation_rows": len(rows),
        "relation_deficiency": size - len(rows),
        "initial_relations_verified": 0,
        "determinant_computed": False,
        "annihilator_supplied": False,
        "saturation_bound": saturation_bound,
        "saturations": 0,
        "kernel_passes": 0,
        "tops_evaluated": 0,
        "identity_tops": 0,
        "rejected_nonidentity_tops": 0,
        "paired_bad_top_updates": 0,
        "group_multiplications": 0,
        "group_inversions": 0,
        "maximum_entry_bits": max((abs(value).bit_length() for row in rows
                                    for value in row), default=0),
    }
    if keep_updates:
        stats["updates"] = []

    def mul(left, right):
        stats["group_multiplications"] += 1
        return multiply(left, right)

    positive_powers = [[value] for value in generators]
    negative_powers = [None] * size

    def word_value(exponents):
        result = identity
        for index, exponent in enumerate(exponents):
            if not exponent:
                continue
            if exponent < 0:
                if negative_powers[index] is None:
                    stats["group_inversions"] += 1
                    negative_powers[index] = [inverse(generators[index])]
                powers = negative_powers[index]
                exponent = -exponent
            else:
                powers = positive_powers[index]
            while len(powers) < exponent.bit_length():
                powers.append(mul(powers[-1], powers[-1]))
            bit = 0
            while exponent:
                if exponent & 1:
                    result = mul(result, powers[bit])
                bit += 1
                exponent >>= 1
        return result

    if verify_relations:
        for row in rows:
            if word_value(row) != identity:
                raise ValueError("a supplied row is not a group relation")
            stats["initial_relations_verified"] += 1

    rejected_top = None
    while True:
        dependencies = first_binary_dependencies(rows)
        stats["kernel_passes"] += 1
        if not dependencies:
            stats["termination"] = "odd presented torsion"
            return None, None, stats
        chosen, new_row = None, None
        rejected_dependencies = []
        for dependency in dependencies:
            vector = half_relation(rows, dependency)
            top = word_value(vector)
            stats["tops_evaluated"] += 1
            if mul(top, top) != identity:
                raise ValueError("a half-relation did not give an involution")
            if top == identity:
                stats["identity_tops"] += 1
                chosen, new_row = dependency, vector
                break
            certificate = inspect_top(top)
            if certificate is not None:
                stats["termination"] = "useful involution"
                return top, certificate, stats
            stats["rejected_nonidentity_tops"] += 1
            if rejected_top is None:
                rejected_top = top
            elif top != rejected_top:
                raise ValueError("the predicate rejected two distinct nonidentity involutions")
            rejected_dependencies.append(dependency)
        if chosen is None:
            if len(dependencies) == 1:
                stats["termination"] = "cyclic two-part with its sole top rejected"
                return None, None, stats
            assert len(rejected_dependencies) == 2
            chosen = rejected_dependencies[0] ^ rejected_dependencies[1]
            new_row = half_relation(rows, chosen)
            # XOR differs from the sum by twice a genuine relation.
            # Its half therefore has image rejected_top squared = 1.
            stats["paired_bad_top_updates"] += 1
        if stats["saturations"] >= saturation_bound:
            raise ValueError("the independent-relation index bound was exceeded")
        pivot = (chosen & -chosen).bit_length() - 1
        rows[pivot] = list(new_row)
        stats["saturations"] += 1
        stats["maximum_entry_bits"] = max(
            stats["maximum_entry_bits"],
            max((abs(value).bit_length() for value in new_row), default=0),
        )
        if keep_updates:
            stats["updates"].append({
                "dependency": chosen, "replaced_row": pivot,
                "new_relation": list(new_row),
            })


def factor_from_dominant_relations(integer, generators, relation_matrix, *,
                                   verify_relations=True, keep_updates=False):
    """Split odd N from a strictly dominant matrix of class-group relations.

    A proper factor is guaranteed if N has at least two distinct prime
    factors and the generators span C_{-4N}. The supplied relation rows
    must be genuine; their identities are verified unless the caller
    explicitly records that it has already established them.
    """
    size = len(generators)
    if len(relation_matrix) != size or any(len(row) != size for row in relation_matrix):
        raise ValueError("one square relation matrix is required")
    return factor_from_partial_relations(
        integer, generators, relation_matrix, tuple(range(size)),
        verify_relations=verify_relations, keep_updates=keep_updates,
    )


def factor_from_partial_relations(integer, generators, relation_matrix,
                                  pivot_columns=None, *,
                                  verify_relations=True, keep_updates=False):
    """Try independent partial relations, checking a strictly dominant minor.

    ``pivot_columns`` identifies one distinct generator column per row;
    by default these are the first columns. The corresponding square
    minor must be strictly diagonally dominant. A returned factor is exact.
    Failure concerns only the torsion image of this partial presentation.
    With full class-group generation, failure implies N has at most
    ``len(generators)-len(relation_matrix)+1`` distinct prime factors.
    """
    if integer <= 1 or integer % 2 == 0:
        raise ValueError("the target must be odd and greater than one")
    size, rank = len(generators), len(relation_matrix)
    if rank > size or any(len(row) != size for row in relation_matrix):
        raise ValueError("one generator column per supplied element is required")
    if pivot_columns is None:
        pivot_columns = tuple(range(rank))
    if (len(pivot_columns) != rank or len(set(pivot_columns)) != rank
            or any(not isinstance(column, int) or not 0 <= column < size
                   for column in pivot_columns)):
        raise ValueError("one distinct valid pivot column per row is required")
    if any(abs(row[pivot_columns[index]]) <=
           sum(abs(row[column]) for other, column in enumerate(pivot_columns) if other != index)
           for index, row in enumerate(relation_matrix)):
        raise ValueError("strict diagonal dominance of the designated minor is required")
    discriminant = -4 * integer
    canonical = []
    for form in generators:
        a, middle, c = form
        if (a <= 0 or middle * middle - 4 * a * c != discriminant
                or gcd(gcd(a, middle), c) != 1):
            raise ValueError("a generator is not a primitive positive form of the required discriminant")
        canonical.append(reduce_form(form))
    top, factor, stats = first_useful_relation_top(
        canonical, relation_matrix, (1, 0, integer),
        compose_forms_hnf, inverse_form,
        lambda form: ambiguous_factor(form, integer),
        verify_relations=verify_relations, keep_updates=keep_updates,
    )
    return {"factor": factor, "ambiguous_form": top, "stats": stats}
