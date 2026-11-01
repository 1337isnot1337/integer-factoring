"""Compact exact oriented ideals with batched discriminant validation.

Atoms are norm/root pairs. Reduced forms and expanded HNF coordinates are
not constructed. Fast-arithmetic complexity and Python timings are separate.
"""
from math import gcd
from classical_coprime_basis import factor_into_coprimes
from circuit_box_oriented_batch import _remainder_tree
from oriented_norm_basis_lift import verify_coprime_norms


def _statistics():
    return {'normalization_multiplications':0,'normalization_multiplication_operand_bits':0,
            'normalization_remainders':0,'normalization_remainder_operand_bits':0,
            'discriminant_batches':0,'discriminant_operand_bits':0,
            'local_congruence_checks':0,'unramified_gcds':0,'unramified_gcd_operand_bits':0}


def validate_compact_keys(discriminant, keys, *, statistics=None):
    """Normalize and validate all keys with one discriminant remainder tree."""
    if not isinstance(discriminant,int) or discriminant>=0 or discriminant%4 not in (0,1):
        raise ValueError('Require a negative quadratic-order discriminant.')
    stats = _statistics() if statistics is None else statistics
    normalized = []
    input_bits = discriminant.bit_length()
    for key in keys:
        if len(key)!=2: raise ValueError('Each key consists of its norm and root.')
        norm,root = key
        if not isinstance(norm,int) or not isinstance(root,int) or norm<1 or norm%2!=1:
            raise ValueError('Require a positive odd norm and an integer root.')
        input_bits += norm.bit_length()+max(1,root.bit_length())
        normalized.append((norm,root%(2*norm)))
    if normalized:
        stats['discriminant_batches'] += 1
        stats['discriminant_operand_bits'] += discriminant.bit_length()
        residues = _remainder_tree(discriminant,[4*norm for norm,_ in normalized],stats)
        for (norm,root),residue in zip(normalized,residues):
            stats['local_congruence_checks'] += 1
            if (root*root-residue)%(4*norm): raise ValueError('A root does not define the stated ideal.')
            stats['unramified_gcds'] += 1
            stats['unramified_gcd_operand_bits'] += norm.bit_length()+residue.bit_length()
            if gcd(norm,residue)!=1: raise ValueError('An ideal norm is ramified.')
    return normalized,{'input_bits':input_bits,'norm_bits':sum(norm.bit_length() for norm,_ in normalized),
                       'statistics':stats}


def _balanced_product(values):
    layer = list(values)
    if not layer: return 1
    while len(layer)>1:
        layer = [layer[i]*layer[i+1] if i+1<len(layer) else layer[i] for i in range(0,len(layer),2)]
    return layer[0]


def verify_compact_oriented_basis(record, *, expected_discriminant=None, expected_keys=None):
    """Verify literal ideal words by norms and unique unramified root lifting."""
    discriminant = record['discriminant']
    keys = record['keys']; atoms = record['atoms']; words = record['sparse_words']
    if expected_discriminant is not None and discriminant!=expected_discriminant:
        raise ValueError('The certificate discriminant differs from its input.')
    if expected_keys is not None:
        normalized_input,_ = validate_compact_keys(discriminant,expected_keys)
        if normalized_input!=[tuple(k) for k in keys]:
            raise ValueError('The certificate keys differ from their input.')
    if len(words)!=len(keys): raise ValueError('One sparse word per supplied key is required.')
    if any(set(atom)!={'norm','root'} for atom in atoms): raise ValueError('Atoms must be compact norm/root pairs.')
    norm_bits = sum(norm.bit_length() for norm,_ in keys)
    if any(not isinstance(a['norm'],int) or a['norm']<=1 for a in atoms):
        raise ValueError('Every atom norm must exceed one.')
    if sum(a['norm'].bit_length() for a in atoms)>2*norm_bits:
        raise ValueError('The atom data exceed their norm-product input bound.')
    validated,validation = validate_compact_keys(discriminant,list(keys)+[(a['norm'],a['root']) for a in atoms])
    normalized_keys = validated[:len(keys)]; normalized_atoms = validated[len(keys):]
    if normalized_keys!=[tuple(k) for k in keys]: raise ValueError('Input keys are not normalized.')
    if normalized_atoms!=[(a['norm'],a['root']) for a in atoms]: raise ValueError('Atom roots are not normalized.')
    for norm,root in normalized_atoms:
        if root!=min(root,(-root)%(2*norm)): raise ValueError('Atom roots are not canonical.')
    coprime = verify_coprime_norms([a['norm'] for a in atoms])
    stats = _statistics(); used = set(); incidences = 0
    for (norm,root),word in zip(normalized_keys,words):
        seen = set(); weighted = 0
        for column,exponent in word:
            if (not isinstance(column,int) or not 0<=column<len(atoms) or column in seen
                    or not isinstance(exponent,int) or exponent==0):
                raise ValueError('Invalid sparse exponent word.')
            seen.add(column); used.add(column)
            weighted += (atoms[column]['norm'].bit_length()-1)*abs(exponent)
            if weighted>norm.bit_length()-1: raise ValueError('An exponent word exceeds its claimed norm.')
        terms = [atoms[column]['norm']**abs(exponent) for column,exponent in word]
        if _balanced_product(terms)!=norm: raise ValueError('An exponent word has the wrong norm.')
        residues = _remainder_tree(root,[2*atoms[column]['norm'] for column,_ in word],stats)
        for (column,exponent),residue in zip(word,residues):
            atom = atoms[column]
            expected = (atom['root'] if exponent>0 else -atom['root'])%(2*atom['norm'])
            if residue!=expected: raise ValueError('A root orientation is inconsistent.')
            incidences += 1
    if used!=set(range(len(atoms))): raise ValueError('An unused atom was supplied.')
    return {'arithmetic_literal_ideal_identities':len(keys),'sparse_orientation_incidences':incidences,
            'compact_atoms':len(atoms),'coprimality':coprime,'key_validation':validation,
            'root_remainder_statistics':stats,'expanded_hnf_products':0,'reduced_forms_constructed':0}


# The two-pass core is adapted without reduced-form construction from
# circuit_box_oriented_batch.py SHA256 9dfee586c29373d8769992a17f84e6f8d87c4c4c713baa3d6f683cbe0c06cba4.
def compact_oriented_coprime_basis(discriminant, keys, *, unitary_preconditioning=True, verify=True):
    """Construct norms, batch every orientation split, and emit exact words."""
    if not isinstance(discriminant, int) or discriminant >= 0 or discriminant % 4 not in (0, 1):
        raise ValueError("a negative quadratic-order discriminant is required")
    keys, input_validation = validate_compact_keys(discriminant, keys)
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
    first_pass = factor_into_coprimes([norm for norm, _ in keys], verify=False)
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
    second_pass = factor_into_coprimes(second_inputs, verify=False)
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
        atoms.append({"norm": norm, "root": root})
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
        "scope": "Compact supplied-key construction with batched discriminant checks; no relation-supply bound.",
        "input_validation": input_validation, "total_input_bits": input_validation["input_bits"],
        "reduced_forms_constructed": 0,
    }
    if verify:
        result["verification"] = verify_compact_oriented_basis(result)
    return result
