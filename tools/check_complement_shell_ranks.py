#!/usr/bin/env python3
"""Complete ranks of the exact complementary-shell crossing tensor.

The parity of single occupation at each complementary mode pair is conserved.
Both tensor flattenings therefore split into exact XOR-syndrome blocks.  The
checker retains all blocks and all polynomial coefficients; it does not use
sampled minors or omit intermediate high modes.
"""

from collections import Counter, defaultdict
from itertools import product
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from check_complement_shell_crossings import shell
from check_modular_crossing import column_bits, rank_two
from complement_shell_tensor import local_tensor, multiply, pair_occupations


PAIR_TABLE = tuple(tuple(local_tensor(x, y, u, v)
                        for v, u in product(range(4), repeat=2))
                   for x, y in product(range(4), repeat=2))


def syndrome(occupations):
    return sum((part.bit_count() & 1) << i for i, part in enumerate(occupations))


def pair_word(first, second):
    return sum((4 * a + b) << (4 * i) for i, (a, b) in enumerate(zip(first, second)))


def entry(first, second):
    """first encodes (x,y), second encodes (v,u), in four bits per site."""
    result = 1
    occupied = first | second
    while occupied:
        shift = ((occupied & -occupied).bit_length() - 1) // 4 * 4
        value = PAIR_TABLE[(first >> shift) & 15][(second >> shift) & 15]
        if not value:
            return 0
        result = multiply(result, value)
        occupied &= ~(15 << shift)
    return result


def verify_local_symmetry():
    checked = nonzero = 0
    for x, y, u, v in product(range(4), repeat=4):
        value = local_tensor(x, y, u, v)
        checked += 1
        if value:
            assert ((x ^ y ^ u ^ v).bit_count() & 1) == 0
            nonzero += 1
    singleton_pairs = ((1, 0), (2, 0), (0, 1), (0, 2))
    minor_checks = 0
    for (x, y), (v, u) in product(singleton_pairs, repeat=2):
        mode_equal = (x | y) == (v | u)
        expected = (2 if v else 6) if x else 1
        assert local_tensor(x, y, u, v) == (expected if mode_equal else 0)
        minor_checks += 1
    for (x, u), (v, y) in product(singleton_pairs, repeat=2):
        # The column invariant complements the state mode on both sides.
        first_mode = x if x else u ^ 3
        second_mode = v if v else y ^ 3
        expected = (2 if v else 3) if x else 1
        assert local_tensor(x, y, u, v) == (expected if first_mode == second_mode else 0)
        minor_checks += 1
    return {"local_entries_checked": checked, "local_nonzero_entries": nonzero,
            "protected_minor_local_identities_checked": minor_checks}


def protected_central_minors(modulus, parameters, states=None):
    """Exact identity-minor sizes from complement-pair-disjoint labels.

    Row-cut keys are a union b; column-cut keys are a union (M-b).
    A label containing both modes of any complementary pair is excluded.
    Every distinct union gives a different diagonal pivot at degree L.
    """
    states = parameters if states is None else states

    def simple(parts):
        values = pair_occupations(parts, modulus)
        if 3 in values:
            return None
        support = sum(bool(value) << i for i, value in enumerate(values))
        mask = sum(1 << ((part - 1) // 2) for part in parts)
        opposite = sum(1 << ((modulus - part - 1) // 2) for part in parts)
        return support, mask, opposite

    first = [value for parts in parameters if (value := simple(parts)) is not None]
    second = [value for parts in states if (value := simple(parts)) is not None]
    row_unions, column_unions = set(), set()
    disjoint_pairs = 0
    for first_support, first_mask, _ in first:
        for second_support, second_mask, second_complement in second:
            if not first_support & second_support:
                disjoint_pairs += 1
                row_unions.add(first_mask | second_mask)
                column_unions.add(first_mask | second_complement)
    r, c = len(row_unions), len(column_unions)
    return {
        "simple_parameter_labels": len(first), "simple_state_labels": len(second),
        "complement_pair_disjoint_ordered_pairs": disjoint_pairs,
        "row_identity_minor": r, "column_identity_minor": c,
        "rank_sum_lower_bound": math.sqrt(r * c),
    }


def rank_labels(modulus, parameters, states=None, compare_dense=False):
    """Complete ranks on equally sized parameter and state shell subsets.

    Parameters and states can differ, but every label must have the same
    length and energy.  The partner labels are mode complements of the
    parameters.  Both cuts are (x,y)|(v,u) and (x,u)|(v,y).
    """
    started = time.monotonic()
    parameters = tuple(tuple(parts) for parts in parameters)
    states = parameters if states is None else tuple(tuple(parts) for parts in states)
    if not parameters or len(parameters) != len(states):
        raise ValueError("parameter and state domains must have equal positive size")
    length, weight = len(parameters[0]), sum(parameters[0])
    if weight > length * modulus // 2:
        raise ValueError("a nonempty shell of weight at most L*M/2 is required")
    if (len(set(parameters)) != len(parameters) or len(set(states)) != len(states)
            or any(len(parts) != length or sum(parts) != weight
                   for parts in parameters + states)):
        raise ValueError("domains must be distinct labels from one length/energy shell")
    d = len(parameters)
    parameter_occupations = [pair_occupations(parts, modulus) for parts in parameters]
    state_occupations = [pair_occupations(parts, modulus) for parts in states]
    parameter_sigmas = [syndrome(parts) for parts in parameter_occupations]
    state_sigmas = [syndrome(parts) for parts in state_occupations]
    blocks = defaultdict(list)
    for x, y in product(range(d), repeat=2):
        blocks[parameter_sigmas[x] ^ state_sigmas[y]].append((x, y))
    indices = {pair: (sigma, i) for sigma, pairs in blocks.items()
               for i, pair in enumerate(pairs)}
    row = {sigma: [[0] * len(pairs) for _ in range(2 * length + 1)]
           for sigma, pairs in blocks.items()}
    col = {sigma: [[0] * len(pairs) for _ in range(2 * length + 1)]
           for sigma, pairs in blocks.items()}
    digest = hashlib.sha256()
    nonzero_entries = 0
    for sigma in sorted(blocks):
        pairs = blocks[sigma]
        words = [pair_word(parameter_occupations[x], state_occupations[y]) for x, y in pairs]
        for i, (x, y) in enumerate(pairs):
            for k, (v, u) in enumerate(pairs):
                value = entry(words[i], words[k])
                assert value.bit_length() <= 2 * length + 1
                digest.update(value.to_bytes((2 * length + 8) // 8, 'little'))
                assert value.bit_count() % 2 == int(x == v and y == u)
                if not value:
                    continue
                nonzero_entries += 1
                c_sigma, c_i = indices[x, u]
                c_sigma_check, c_k = indices[v, y]
                assert c_sigma == c_sigma_check
                while value:
                    bit = value & -value
                    degree = bit.bit_length() - 1
                    row[sigma][degree][i] ^= 1 << k
                    col[c_sigma][degree][c_i] ^= 1 << c_k
                    value ^= bit
    profile = []
    rank_distributions = []
    for degree in range(2 * length + 1):
        row_ranks = [rank_two(row[sigma][degree]) for sigma in sorted(blocks)]
        col_ranks = [rank_two(col[sigma][degree]) for sigma in sorted(blocks)]
        r, c = sum(row_ranks), sum(col_ranks)
        if not r:
            assert not c
            continue
        profile.append([(degree - length) * modulus, r, c])
        rank_distributions.append({
            "degree": degree,
            "row_block_rank_counts": sorted(Counter(row_ranks).items()),
            "column_block_rank_counts": sorted(Counter(col_ranks).items()),
        })
    if compare_dense:
        for degree in range(2 * length + 1):
            dense = [0] * (d * d)
            for x, y, u, v in product(range(d), repeat=4):
                value = entry(pair_word(parameter_occupations[x], state_occupations[y]),
                              pair_word(parameter_occupations[v], state_occupations[u]))
                if value >> degree & 1:
                    dense[x * d + y] ^= 1 << (v * d + u)
            r, c = rank_two(dense), rank_two(column_bits(dense, d))
            assert r == sum(rank_two(row[sigma][degree]) for sigma in blocks)
            assert c == sum(rank_two(col[sigma][degree]) for sigma in blocks)
    threshold = d ** math.sqrt(2)
    rank_sum = sum(math.sqrt(r * c) for _, r, c in profile)
    minors = protected_central_minors(modulus, parameters, states)
    central = next(((r, c) for shift, r, c in profile if shift == 0), (0, 0))
    assert central[0] >= minors['row_identity_minor']
    assert central[1] >= minors['column_identity_minor']
    return {
        "mode_cutoff": modulus, "length": length, "operator_weight": weight,
        "partner_weight": length * modulus - weight, "labels": d,
        "nonzero_profile": profile,
        "xor_syndrome_blocks": len(blocks),
        "largest_block": max(map(len, blocks.values())),
        "evaluated_entries": sum(len(pairs) ** 2 for pairs in blocks.values()),
        "dense_entry_count": d ** 4, "nonzero_entries": nonzero_entries,
        "block_order_profile_sha256": digest.hexdigest(),
        "rank_sum": rank_sum, "finite_seed_threshold": threshold,
        "seed_margin": threshold / rank_sum,
        "maximum_rank_product": max(r * c for _, r, c in profile),
        "protected_central_minors": minors,
        "rank_distributions": rank_distributions,
        "independent_dense_rank_comparison": compare_dense,
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def run_case(modulus, length, weight, compare_dense=False):
    return rank_labels(modulus, shell(modulus, length, weight),
                       compare_dense=compare_dense)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="16:4:32,32:3:47,32:4:64,64:3:95")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dense-comparison-limit", type=int, default=8)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "exact complete-profile ranks using conserved pair-occupation parity; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "local_tensor_sha256": hashlib.sha256(Path(__file__).with_name(
            'complement_shell_tensor.py').read_bytes()).hexdigest(),
        "local_symmetry": verify_local_symmetry(), "cases": [],
    }
    for description in filter(None, args.cases.split(',')):
        modulus, length, weight = map(int, description.split(':'))
        dense = len(shell(modulus, length, weight)) <= args.dense_comparison_limit
        result = run_case(modulus, length, weight, dense)
        report["cases"].append(result)
        print(json.dumps({key: result[key] for key in (
            'mode_cutoff', 'length', 'labels', 'nonzero_profile', 'rank_sum',
            'finite_seed_threshold', 'seed_margin', 'elapsed_seconds')}), flush=True)
        if args.output:
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
