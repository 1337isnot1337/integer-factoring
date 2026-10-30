#!/usr/bin/env python3
"""Bounded Frobenius-module certificate and honest P10 shortfall.

This audits an improved common-functional compression without constructing a
full ten-site matrix. It also verifies explicit constant unipotent tensor
splittings, rather than treating character identities as module isomorphisms.
"""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import hashlib
import json
from math import comb
from pathlib import Path

from assemble_two_even_site_cores import monomials
from check_common_character_compressions import low_basis, matvec
from check_folded_gram_channels import verify_elementary_split
from check_modular_crossing import rank_two
from check_two_even_site_P8_compressions import audit_witness
from check_two_even_site_jordan_compressions import build_cores, embedding_controls, tensor
from check_two_even_site_row_corners import physical_certificate, P7_CERTIFICATE, P8_CERTIFICATE
from find_tagged_common_character_compressions import search_block
from folded_gram_ranks import frobenius_weights, tensor_power_blocks


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "evidence/two-even-site-P10-frobenius-certificate.json"
CUBE_COLUMNS = [1, 22, 104, 128, 6, 96, 20, 40]
PAIRS, LABELS = 10, 11520
CHARACTER = [1] + [0] * 7 + [1] + [0] * 9


def weighted_free(weights):
    result = {0: [1]}
    for weight in weights:
        old = result
        size = 2 * len(old[0])
        result = {degree: [0] * size for a in old for degree in (a, a + weight)}
        for degree, matrix in old.items():
            for shift, local in ((0, [1, 2]), (weight, [2, 0])):
                result[degree + shift] = [a ^ b for a, b in zip(result[degree + shift], tensor(matrix, local))]
    return result


def permute_columns(weights, columns, order):
    count = len(weights)
    reordered = []
    for index in range(1 << count):
        old = sum(((index >> (count - 1 - j)) & 1) << (count - 1 - position)
                  for j, position in enumerate(order))
        reordered.append(columns[old])
    return [weights[i] for i in order], reordered


def explicit_free_split(power):
    pending = [([1] * power, [1 << i for i in range(1 << power)])]
    result = []
    while pending:
        weights, columns = pending.pop()
        counts = Counter(weights)
        repeated = next((weight for weight in sorted(counts) if counts[weight] >= 3), None)
        if repeated is None:
            order = sorted(range(len(weights)), key=lambda i: weights[i])
            weights, columns = permute_columns(weights, columns, order)
            assert tuple(weights) == frobenius_weights(sum(weights))
            result.append((weights, columns))
            continue
        selected = [i for i, weight in enumerate(weights) if weight == repeated][:3]
        order = selected + [i for i in range(len(weights)) if i not in selected]
        weights, columns = permute_columns(weights, columns, order)
        tail = weights[3:]
        tail_size = 1 << len(tail)
        for new_weights, local_columns in (([2 * repeated, repeated], CUBE_COLUMNS[:4]),
                                            ([repeated], CUBE_COLUMNS[4:6]),
                                            ([repeated], CUBE_COLUMNS[6:8])):
            new_columns = []
            for column in local_columns:
                for index in range(tail_size):
                    vector = 0
                    for row in range(8):
                        if column >> row & 1:
                            vector ^= columns[row * tail_size + index]
                    new_columns.append(vector)
            pending.append((new_weights + tail, new_columns))
    return result


def verify_free_splits():
    reports = []
    for power in range(7):
        source = weighted_free([1] * power)
        blocks = explicit_free_split(power)
        columns = [vector for _, vectors in blocks for vector in vectors]
        assert len(columns) == rank_two(columns) == 1 << power
        multiplicities = Counter((sum(weights), (power - sum(weights)) // 2) for weights, _ in blocks)
        assert multiplicities == Counter({(k, d): count for k, d, count in tensor_power_blocks(power)})
        checks = 0
        for weights, vectors in blocks:
            target = weighted_free(weights)
            for degree in range(power + 1):
                first = source.get(degree, [0] * (1 << power))
                second = target.get(degree, [0] * len(vectors))
                for j, vector in enumerate(vectors):
                    expected = 0
                    for i, image in enumerate(vectors):
                        if second[i] >> j & 1:
                            expected ^= image
                    assert matvec(first, vector) == expected
                    checks += 1 << power
        reports.append({"natural_factors": power, "ambient_dimension": 1 << power,
                        "complete_polynomial_intertwining_entries_checked": checks,
                        "blocks": [{"weights": weights, "basis_columns": [hex(vector) for vector in vectors]}
                                   for weights, vectors in blocks]})
    return reports


def column_module(cores, h, syndrome, weights):
    data = cores[h]
    block = data["syndromes"].index(syndrome)
    free = weighted_free(weights)
    size = len(data["cuts"][1][0][block]) * (1 << len(weights))
    positions = {term: i for i, term in enumerate(monomials(PAIRS))}
    result = [[0] * size for _ in positions]
    for (degree, charge), matrix in zip(data["monomials"], data["cuts"][1]):
        for extra, other in free.items():
            index = positions[degree + extra, charge]
            result[index] = [a ^ b for a, b in zip(result[index], tensor(matrix[block], other))]
    return result


def free_multiplicities(copies):
    result = defaultdict(int)
    for natural in range(copies + 1):
        for highest, twist, multiplicity in tensor_power_blocks(natural):
            result[highest] += comb(copies, natural) * 2 ** (copies - natural) * multiplicity
    assert sum(count * (1 << len(frobenius_weights(highest))) for highest, count in result.items()) == 4 ** copies
    return result


def physical_total(values):
    total = comb(PAIRS, 4) * values[4, 15] + values[5, 0]
    restrictions = []
    for syndrome in combinations(range(PAIRS), 2):
        choices = []
        for anchors in combinations([site for site in range(PAIRS) if site not in syndrome], 3):
            sites = sorted(syndrome + anchors)
            relative = sum(1 << sites.index(site) for site in syndrome)
            choices.append((relative, anchors))
        relative, anchors = max(choices, key=lambda item: values[5, item[0]])
        amount = values[5, relative]
        total += amount
        restrictions.append({"syndrome_sites": list(syndrome), "anchors": list(anchors),
                             "core_syndrome": relative, "compression": amount})
    return total, restrictions


def audit(cores, certificate, discover=False):
    assert certificate["pairs"] == PAIRS and certificate["labels"] == LABELS
    assert certificate["character"] == CHARACTER
    cube = verify_elementary_split()
    splits = verify_free_splits()
    expected_types = {(4, 15)} | {(5, s) for s in cores[5]["syndromes"] if s.bit_count() in (0, 2)}
    assert {(item["core_pairs"], item["syndrome"]) for item in certificate["types"]} == expected_types
    values, leading, uniform, reports = {}, {}, {}, []
    uniform_mask = (1 << 8) | (1 << 16)  # A8 + B6; chi is one.
    for supplied in certificate["types"]:
        h, syndrome = supplied["core_pairs"], supplied["syndrome"]
        multiplicities = free_multiplicities(PAIRS - h)
        minimum_highest = 8 - (2 if h == 4 else 4)
        # Both A8 and B6 vanish on every discarded summand. This makes
        # the same coefficient-rank ceiling valid before those summands
        # are discarded, including compressions that mix the summands.
        core_block = cores[h]["syndromes"].index(syndrome)
        maximum_degrees = (2, 0) if h == 4 else (4, 2)
        for (degree, charge), generator in zip(cores[h]["monomials"], cores[h]["cuts"][1]):
            if degree > maximum_degrees[charge]:
                assert not any(generator[core_block])
        assert 8 - maximum_degrees[0] == 6 - maximum_degrees[1] == minimum_highest
        expected_highest = {k for k in multiplicities if k >= minimum_highest}
        assert {item["highest"] for item in supplied["modules"]} == expected_highest
        total = lead = ceiling = 0
        modules = []
        for item in supplied["modules"]:
            highest = item["highest"]
            assert tuple(item["weights"]) == frobenius_weights(highest)
            assert item["multiplicity"] == multiplicities[highest]
            matrices = column_module(cores, h, syndrome, item["weights"])
            if discover:
                assert search_block(matrices, CHARACTER, seed=0) == item["witness"]
            report = audit_witness(matrices, CHARACTER, item["witness"])
            rank = rank_two(matrices[8])
            assert rank == item["leading_rank"]
            upper = rank_two([a ^ b for a, b in zip(matrices[8], matrices[16])])
            total += item["multiplicity"] * item["witness"]["compression"]
            lead += item["multiplicity"] * rank
            ceiling += item["multiplicity"] * upper
            modules.append(dict(report, highest=highest, weights=item["weights"], multiplicity=item["multiplicity"],
                                leading_coefficient_rank=rank, uniform_coefficient_rank=upper))
        assert total == supplied["weighted_compression"] and lead == supplied["leading_rank_ceiling"]
        values[h, syndrome], leading[h, syndrome], uniform[h, syndrome] = total, lead, ceiling
        reports.append({"core_pairs": h, "core_syndrome": syndrome, "weighted_compression": total,
                        "leading_coefficient_rank_ceiling": lead, "uniform_coefficient_rank_ceiling": ceiling,
                        "modules": modules})
    total, restrictions = physical_total(values)
    leading_total, _ = physical_total(leading)
    uniform_total, _ = physical_total(uniform)
    assert total == certificate["column_total"] == 5738
    assert restrictions == certificate["restrictions"]
    assert leading_total == certificate["leading_rank_ceiling"] == 10638
    assert uniform_total == 9810
    row = (2048 * comb(PAIRS, 4) + 6144 * comb(PAIRS, 2) - 9077) * 4 ** (PAIRS - 7)
    assert row == certificate["row_bound"] == 44638912
    product = row * total
    assert product == certificate["product"] == 256138077056
    assert product ** 16 < LABELS ** 45 and 45 ** 2 < 2 * 32 ** 2
    prior7, prior8 = [json.loads(path.read_text()) for path in (P7_CERTIFICATE, P8_CERTIFICATE)]
    corner = json.loads((ROOT / "evidence/two-even-site-row-corner-certificate.json").read_text())
    corner_values = {item["core_syndrome"]: item["rank_dimension"] for item in corner["corners"]}
    physical = physical_certificate(PAIRS, corner_values, prior7, prior8)
    for record, column in zip(physical["three_anchor_restrictions"], restrictions):
        assert record["syndrome_sites"] == column["syndrome_sites"]
        record["column"] = {key: column[key] for key in ("anchors", "core_syndrome", "compression")}
    controls = embedding_controls(cores, physical)
    return {"status": "pass", "pairs": PAIRS, "labels": LABELS, "character": CHARACTER,
            "character_scalar_polynomial": "(q+r)^2 (p+r)^8",
            "generic_natural_cube_intertwiner": cube, "explicit_unipotent_direct_sums": splits,
            "modules": reports, "rank_certificate_dimensions": [row, total], "rank_product": product,
            "scalar_compression_entries_checked": sum(module["all_generator_compression_entries_checked"]
                                                       for record in reports for module in record["modules"]),
            "literal_ten_site_embedding_controls": controls,
            "largest_matrix_dimension": max(item["witness"]["size"] for record in certificate["types"] for item in record["modules"]),
            "leading_coefficient_column_ceiling": leading_total,
            "retained_module_uniform_column_ceiling": uniform_total,
            "uniform_ceiling_coefficient_mask": uniform_mask,
            "uniform_ceiling_scope": "The full selected core principal restrictions, including all free-factor summands, allowing the best physical anchor choice for each syndrome. This does not cover enlarging the selected cores or changing the functional.",
            "exact_shortfall": {"rank_product_power": 16, "label_power": 45,
                                "difference": LABELS ** 45 - product ** 16},
            "strict_seed_excluded": False,
            "scope": "An all-field, all-power rank-cost lower certificate for this seed; its value does not reach the target exponent."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--basis5", type=Path)
    parser.add_argument("--discover", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    cores = build_cores(args.basis5)
    result = audit(cores, certificate, args.discover)
    result["certificate_sha256"] = hashlib.sha256(args.certificate.read_bytes()).hexdigest()
    names = (Path(__file__).name, "check_folded_gram_channels.py", "folded_gram_ranks.py",
             "check_two_even_site_jordan_compressions.py", "check_two_even_site_P8_compressions.py",
             "check_two_even_site_row_corners.py", "assemble_two_even_site_cores.py",
             "check_common_character_compressions.py", "find_tagged_common_character_compressions.py")
    result["source_sha256"] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}
    result["source_certificate_sha256"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (P7_CERTIFICATE, P8_CERTIFICATE, ROOT / "evidence/two-even-site-row-corner-certificate.json")}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
