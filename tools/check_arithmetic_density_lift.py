#!/usr/bin/env python3
"""Dense smooth subset -> global class-number multiple: finite exact controls."""
from functools import lru_cache
from math import gcd, lcm, prod
import json

from arithmetic_two_group import power
from binary_quadratic_forms import (
    ambiguous_factor, compose_forms, reduced_forms, reduce_form, smooth_prime_word,
)


def check(discriminant):
    forms = reduced_forms(discriminant)
    identity = next(form for form in forms if form[0] == 1)
    multiply = lru_cache(None)(compose_forms)
    generators = []
    for prime in (2, 3, 5, 7):
        if gcd(prime, discriminant) != 1:
            continue
        roots = [middle for middle in range(prime+1)
                 if (middle*middle-discriminant) % (4*prime) == 0]
        if roots:
            middle = roots[0]
            generators.append(reduce_form((prime, middle, (middle*middle-discriminant)//(4*prime))))
    subgroup, queue = {identity}, [identity]
    while queue:
        value = queue.pop()
        for generator in generators:
            new = multiply(value, generator)
            if new not in subgroup:
                subgroup.add(new)
                queue.append(new)
    smooth = {form for form in forms if smooth_prime_word(form, 7) is not None}
    assert smooth <= subgroup
    index_bound = (len(forms)+len(smooth)-1)//len(smooth)
    assert len(forms) % len(subgroup) == 0
    index = len(forms)//len(subgroup)
    assert index <= index_bound
    orders = []
    for generator in generators:
        value, order = generator, 1
        while value != identity:
            value, order = multiply(value, generator), order+1
        orders.append(order)
    determinant = prod(orders)
    assert determinant % len(subgroup) == 0
    correction = lcm(*range(1, index_bound+1))
    global_multiple = determinant*correction
    assert global_multiple % len(forms) == 0
    assert all(power(form, global_multiple, identity, multiply) == identity for form in forms)
    global_ambiguous = {form for form in forms if multiply(form, form) == identity}
    local_ambiguous = global_ambiguous & subgroup
    return {
        'discriminant': discriminant, 'class_number': len(forms),
        'factor_base_subgroup_order': len(subgroup), 'smooth_form_count': len(smooth),
        'actual_index': index, 'density_index_bound': index_bound,
        'subgroup_relation_determinant': determinant, 'index_correction': correction,
        'global_class_number_multiple': global_multiple,
        'subgroup_determinant_already_annihilates_global_group': all(
            power(form, determinant, identity, multiply) == identity for form in forms),
        'global_two_torsion_size': len(global_ambiguous),
        'subgroup_two_torsion_size': len(local_ambiguous),
        'local_ambiguous_forms': sorted(local_ambiguous),
        'global_ambiguous_forms': sorted(global_ambiguous),
    }


def main():
    records = [check(-value) for value in range(3, 513) if (-value) % 4 in (0, 1)]
    control = next(row for row in records if row['discriminant'] == -348)
    assert (control['class_number'], control['factor_base_subgroup_order'],
            control['smooth_form_count'], control['density_index_bound']) == (6, 3, 3, 2)
    assert control['subgroup_relation_determinant'] == 3
    assert control['global_class_number_multiple'] == 6
    assert not control['subgroup_determinant_already_annihilates_global_group']
    assert control['subgroup_two_torsion_size'] == 1 < control['global_two_torsion_size']
    assert all(ambiguous_factor(form, 87) is None for form in control['local_ambiguous_forms'])
    assert any(ambiguous_factor(form, 87) == 3 for form in control['global_ambiguous_forms'])
    print(json.dumps({
        'scope': 'Exact density/index transfer; classes and orders enumerated only for finite verification; no fast arithmetic witness extractor',
        'discriminants_checked': len(records),
        'proper_factor_base_subgroups': sum(row['actual_index'] > 1 for row in records),
        'index_correction_needed_for_annihilation': sum(not row['subgroup_determinant_already_annihilates_global_group'] for row in records),
        'subgroups_missing_global_two_torsion': sum(row['subgroup_two_torsion_size'] < row['global_two_torsion_size'] for row in records),
        'maximum_density_index_bound': max(row['density_index_bound'] for row in records),
        'missing_factoring_direction_control': control,
    }, indent=2))


if __name__ == '__main__':
    main()
