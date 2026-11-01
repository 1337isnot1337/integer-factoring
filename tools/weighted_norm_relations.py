"""Exact logarithm-weighted independence checks for bounded norm relations.

The column weights are log(weight_bases[j]); no real logarithm is evaluated.
Before comparing integer products, the wrapper bounds their written sizes
by the supplied norm height. The relation-halving engine is unchanged.
"""

from math import gcd

from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top


def least_power_above(base, height_bound):
    """Least positive k with base**k > height_bound, using exact integers."""
    if not isinstance(base, int) or base < 2:
        raise ValueError("The base must be an integer at least two.")
    if not isinstance(height_bound, int) or height_bound < 1:
        raise ValueError("The height bound must be a positive integer.")
    exponent, value = 0, 1
    while value <= height_bound:
        exponent += 1
        value *= base
    return exponent


def check_weighted_dominance(relation_matrix, weight_bases, pivot_columns,
                             height_bound):
    """Certify a dominant designated minor by bounded integer products.

    For rows K_i e_i-v_i constructed from n_i<=H and the least K_i
    with p_i**K_i>H, the checked written-product bound always holds.
    Arbitrarily large binary exponents are rejected before exponentiation.
    """
    bases = tuple(weight_bases)
    size, rank = len(bases), len(relation_matrix)
    if not isinstance(height_bound, int) or height_bound < 1:
        raise ValueError("The height bound must be a positive integer.")
    if any(not isinstance(base, int) or base < 2 for base in bases):
        raise ValueError("Every weight base must be an integer at least two.")
    if rank > size or any(len(row) != size for row in relation_matrix):
        raise ValueError("Independent relation rows cannot outnumber columns.")
    if any(not isinstance(entry, int) for row in relation_matrix for entry in row):
        raise ValueError("Relation entries must be integers.")
    if (len(pivot_columns) != rank or len(set(pivot_columns)) != rank
            or any(not isinstance(column, int) or not 0 <= column < size
                   for column in pivot_columns)):
        raise ValueError("One distinct valid pivot column per row is required.")
    bit_weights = [base.bit_length() for base in bases]
    product_bit_budget = 4 * height_bound.bit_length() + 2 * max(bit_weights, default=0)
    row_budgets = [sum(abs(entry) * bits for entry, bits in zip(row, bit_weights))
                   for row in relation_matrix]
    if any(value > product_bit_budget for value in row_budgets):
        raise ValueError("A row exceeds the bounded norm-product size budget.")
    comparisons = []
    for index, row in enumerate(relation_matrix):
        pivot = pivot_columns[index]
        left = bases[pivot] ** abs(row[pivot])
        right = 1
        for other, column in enumerate(pivot_columns):
            if other != index:
                right *= bases[column] ** abs(row[column])
        if left <= right:
            raise ValueError("The designated minor is not strictly logarithm-weighted dominant.")
        comparisons.append((left.bit_length(), right.bit_length()))
    return {
        "weight_rule": "logarithm of the supplied positive integer bases",
        "floating_logarithms_used": False,
        "norm_height_bits": height_bound.bit_length(),
        "product_bit_budget": product_bit_budget,
        "maximum_row_product_bit_bound": max(row_budgets, default=0),
        "dominant_minor_rows": rank,
        "comparison_bit_lengths": comparisons,
    }


def factor_from_weighted_relations(integer, generators, relation_matrix,
                                    weight_bases, pivot_columns=None, *,
                                    height_bound, verify_relations=True,
                                    keep_updates=False):
    """Try bounded rows whose designated minor is logarithm-weighted dominant.

    The weights establish rational independence only. Every row must still
    be a genuine relation in the primitive form class group of -4N.
    Verification can be skipped only after the caller has established it.
    Failure has the same partial torsion-image meaning as the unweighted
    relation-halving wrapper. No full-generation assertion is made here.
    """
    if not isinstance(integer, int) or integer <= 1 or integer % 2 == 0:
        raise ValueError("The target must be an odd integer greater than one.")
    if len(weight_bases) != len(generators):
        raise ValueError("One positive weight base per generator is required.")
    if pivot_columns is None:
        pivot_columns = tuple(range(len(relation_matrix)))
    dominance = check_weighted_dominance(
        relation_matrix, weight_bases, pivot_columns, height_bound)
    discriminant = -4 * integer
    canonical = []
    for form in generators:
        a, middle, c = form
        if (a <= 0 or middle * middle - 4 * a * c != discriminant
                or gcd(gcd(a, middle), c) != 1):
            raise ValueError("A generator is not a primitive positive form of the required discriminant.")
        canonical.append(reduce_form(form))
    top, factor, stats = first_useful_relation_top(
        canonical, relation_matrix, (1, 0, integer),
        compose_forms_hnf, inverse_form,
        lambda form: ambiguous_factor(form, integer),
        verify_relations=verify_relations, keep_updates=keep_updates,
    )
    stats["weighted_dominance"] = dominance
    return {"factor": factor, "ambiguous_form": top, "stats": stats}
