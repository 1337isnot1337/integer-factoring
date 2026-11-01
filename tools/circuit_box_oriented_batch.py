"""Two coprime-basis passes for supplied unramified quadratic ideals."""

from math import gcd

from binary_quadratic_forms import reduce_form
from classical_coprime_basis import factor_into_coprimes
from oriented_coprime_basis import _form, _key
from oriented_norm_basis_lift import verify_lifted_oriented_basis


def _product_tree(values, statistics):
    levels = [list(values)]
    while len(levels[-1]) > 1:
        previous, following = levels[-1], []
        for index in range(0, len(previous), 2):
            if index + 1 == len(previous):
                following.append(previous[index])
            else:
                first, second = previous[index], previous[index + 1]
                statistics["normalization_multiplications"] += 1
                statistics["normalization_multiplication_operand_bits"] += first.bit_length() + second.bit_length()
                following.append(first * second)
        levels.append(following)
    return levels


def _remainder_tree(value, moduli, statistics):
    if not moduli:
        return []
    levels = _product_tree(moduli, statistics)

    def remainder(first, second):
        statistics["normalization_remainders"] += 1
        statistics["normalization_remainder_operand_bits"] += first.bit_length() + second.bit_length()
        return first % second

    residues = [remainder(value, levels[-1][0])]
    for layer in reversed(levels[:-1]):
        residues = [remainder(residues[index // 2], modulus) for index, modulus in enumerate(layer)]
    return residues


def two_pass_oriented_basis(discriminant, keys, *, unitary_preconditioning=True, verify=True):
    """Construct norms, batch every orientation split, and emit exact words."""
    if not isinstance(discriminant, int) or discriminant >= 0 or discriminant % 4 not in (0, 1):
        raise ValueError("a negative quadratic-order discriminant is required")
    keys = [_key(discriminant, *key) for key in keys]
    statistics = {
        "normalization_multiplications": 0,
        "normalization_multiplication_operand_bits": 0,
        "normalization_remainders": 0,
        "normalization_remainder_operand_bits": 0,
        "normalized_root_incidences": 0,
        "unitary_preconditioning_gcds": 0,
        "unitary_preconditioning_gcd_operand_bits": 0,
        "orientation_gcd_calls_after_basis_passes": 0,
    }
    first_pass = factor_into_coprimes([norm for norm, _ in keys], verify=verify)
    parents = first_pass["basis"]
    parent_words = first_pass["words"]
    incidence = [[] for _ in parents]
    normalized = {}
    for index, ((_, root), word) in enumerate(zip(keys, parent_words)):
        residues = _remainder_tree(root, [2 * parents[column] for column, _ in word], statistics)
        for (column, exponent), residue in zip(word, residues):
            incidence[column].append((index, exponent))
            normalized[index, column] = residue
            statistics["normalized_root_incidences"] += 1
    references = [normalized[rows[0][0], column] for column, rows in enumerate(incidence)]
    second_inputs = list(parents)
    differences = {}
    difference_rows = {}
    for column, rows in enumerate(incidence):
        for index, _ in rows:
            difference = ((normalized[index, column] - references[column]) // 2) % parents[column]
            if unitary_preconditioning and difference:
                statistics["unitary_preconditioning_gcds"] += 1
                statistics["unitary_preconditioning_gcd_operand_bits"] += parents[column].bit_length() + difference.bit_length()
                difference = gcd(parents[column], difference)
            differences[index, column] = difference
            if difference:
                difference_rows[index, column] = len(second_inputs)
                second_inputs.append(difference)
    second_pass = factor_into_coprimes(second_inputs, verify=verify)
    second_basis, second_words = second_pass["basis"], second_pass["words"]
    children = second_words[:len(parents)]
    parent_of = {}
    for parent, word in enumerate(children):
        for column, _ in word:
            if column in parent_of:
                raise AssertionError("a retained component belongs to two coprime parents")
            parent_of[column] = parent
    retained = sorted(parent_of)
    if unitary_preconditioning:
        if len(retained) != len(second_basis):
            raise AssertionError("a unitary input produced an auxiliary component")
        if any(exponent != 1 for word in children for _, exponent in word):
            raise AssertionError("a natural unitary component is not a whole parent divisor")
    output_column = {column: index for index, column in enumerate(retained)}
    reference_roots = {}
    for parent, word in enumerate(children):
        residues = _remainder_tree(references[parent], [2 * second_basis[column] for column, _ in word], statistics)
        for (column, _), residue in zip(word, residues):
            reference_roots[column] = residue
    atoms, column_signs = [], {}
    for column in retained:
        norm, reference = second_basis[column], reference_roots[column]
        root = min(reference, (-reference) % (2 * norm))
        column_signs[column] = 1 if root == reference else -1
        atoms.append({"norm": norm, "root": root, "form": reduce_form(_form(discriminant, norm, root))})
    words = [[] for _ in keys]
    for index, parent_word in enumerate(parent_words):
        for parent, exponent in parent_word:
            difference = differences[index, parent]
            powers = dict(second_words[difference_rows[index, parent]]) if difference else {}
            for column, parent_power in children[parent]:
                difference_power = powers.get(column, 0)
                if difference and 0 < difference_power < parent_power:
                    raise AssertionError("an unramified orientation split cut a prime-power component")
                sign = 1 if not difference or difference_power else -1
                coefficient = column_signs[column] * sign * exponent * parent_power
                words[index].append((output_column[column], coefficient))
        words[index].sort()
    input_bits = sum(norm.bit_length() for norm, _ in keys)
    second_bits = sum(value.bit_length() for value in second_inputs)
    floor_log_sum = sum(norm.bit_length() - 1 for norm, _ in keys)
    if second_bits > 4 * floor_log_sum:
        raise AssertionError("second-pass normalized data exceeded its input-size bound")
    if statistics["unitary_preconditioning_gcd_operand_bits"] > 4 * floor_log_sum:
        raise AssertionError("unitary gcd operands exceeded the input-size bound")
    usage = sum(map(len, words))
    if usage > floor_log_sum:
        raise AssertionError("sparse output exceeded the norm-product bound")
    statistics.update({
        "first_pass_basis_size": len(parents),
        "second_pass_input_count": len(second_inputs),
        "second_pass_input_bits": second_bits,
        "second_pass_basis_size": len(second_basis),
        "retained_atoms": len(atoms),
        "discarded_auxiliary_atoms": len(second_basis) - len(atoms),
        "zero_difference_incidences": sum(value == 0 for value in differences.values()),
        "sparse_usage_entries": usage,
        "first_coprime_pass": first_pass["statistics"],
        "second_coprime_pass": second_pass["statistics"],
    })
    result = {
        "discriminant": discriminant, "keys": keys, "atoms": atoms,
        "sparse_words": words, "input_bit_bound": input_bits,
        "second_pass_data_bit_bound": 4 * floor_log_sum,
        "statistics": statistics, "factoring_oracle_calls": 0,
        "residual_primality_tests": 0,
        "unitary_preconditioning": bool(unitary_preconditioning),
        "scope": "Complete two-pass supplied-key construction; no relation-supply bound.",
    }
    if verify:
        result["verification"] = verify_lifted_oriented_basis(result)
    return result
