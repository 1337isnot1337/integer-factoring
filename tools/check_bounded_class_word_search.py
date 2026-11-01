#!/usr/bin/env python3
"""Replay bounded-width paths, exact unpruned controls, and factoring rows."""

from pathlib import Path
import json

from binary_quadratic_forms import reduced_forms
from bounded_class_word_search import BoundedClassWordSearch
from check_square_relation_factoring import rational_determinant, trial_prime
from early_stop_class_factoring import factor_from_annihilated_generators
from quadratic_form_arithmetic import (
    compose_primitive, inverse_primitive, power_primitive, reduce_primitive,
)
from run_square_relation_factoring import first_split_primes
from run_bounded_class_factoring import run


def evaluate(generators, word, identity):
    value = identity
    for generator, exponent in zip(generators, word):
        value = compose_primitive(value, power_primitive(generator, exponent))
    return value


def check_meeting(answer, target, generators, identity, left_budget, right_budget):
    assert answer["found"] and answer["status"] == "found"
    left, right, word = answer["left_word"], answer["right_word"], answer["word"]
    assert len(left) == len(right) == len(word) == len(generators)
    assert sum(abs(value) for value in left) <= answer["left_path_depth"] <= left_budget
    assert sum(abs(value) for value in right) <= answer["right_path_depth"] <= right_budget
    assert list(word) == [a-b for a, b in zip(left, right)]
    assert sum(abs(value) for value in word) == answer["word_l1"] <= left_budget+right_budget
    meeting = tuple(answer["meeting_class"])
    assert evaluate(generators, left, identity) == meeting
    assert compose_primitive(target, evaluate(generators, right, identity)) == meeting
    assert evaluate(generators, word, identity) == target


def check_certificate(path):
    record = json.loads(path.read_text())
    odd_n = record["odd_n"]
    identity = (1, 0, odd_n)
    primes, roots = first_split_primes(-4*odd_n, len(record["factor_base_primes"]))
    assert primes == record["factor_base_primes"]
    assert {str(prime): root for prime, root in roots.items()} == record["roots"]
    generators = [reduce_primitive((prime, roots[prime],
                                    (roots[prime]**2+4*odd_n)//(4*prime)))
                  for prime in primes]
    left, right = record["search"]["left_depth"], record["search"]["right_depth"]
    diagonal = left+right+1
    assert record["diagonal_fixed_before_targets"] == diagonal
    assert record["status"] == "factored" and not record["missing_relation_indices"]
    assert not any(record[name] for name in ("class_number_supplied", "annihilator_supplied",
                                            "input_factors_supplied", "complete_generation_claimed",
                                            "uniform_success_guarantee_claimed"))
    assert len(record["targets"]) == len(generators)
    rows = []
    for index, entry in enumerate(record["targets"]):
        target = power_primitive(generators[index], diagonal)
        assert entry["index"] == index and tuple(entry["target"]) == target
        check_meeting(entry, target, generators, identity, left, right)
        row = [-value for value in entry["word"]]
        row[index] += diagonal
        assert row == entry["relation_row"]
        assert row[index] > sum(abs(value) for j, value in enumerate(row) if j != index)
        assert evaluate(generators, row, identity) == identity
        rows.append(row)
    assert rows == record["relation_matrix"]
    annihilator = abs(rational_determinant(rows))
    assert annihilator == record["annihilator"] > 0
    assert all(power_primitive(generator, annihilator) == identity for generator in generators)
    extraction = factor_from_annihilated_generators(odd_n, generators, annihilator, with_trace=True)
    assert json.loads(json.dumps(extraction)) == record["extraction"]
    factor = extraction["factor"]
    assert factor and 1 < factor < odd_n and odd_n % factor == 0
    factors = sorted((factor, odd_n//factor))
    assert record["factorization"] == factors
    assert max(factors) < 1 << 31 and all(trial_prime(value) for value in factors)
    return {"certificate": path.name, "input_bits": odd_n.bit_length(),
            "rows_checked": len(rows), "maximum_returned_word_l1":
                max(entry["word_l1"] for entry in record["targets"]),
            "rational_elimination_determinant": annihilator,
            "factorization": factors, "returned_factors_checked_prime": True,
            "search_rerun": False}


def unpruned_controls():
    tested = positive = negative = 0
    for absolute in range(3, 129):
        discriminant = -absolute
        if discriminant % 4 not in (0, 1):
            continue
        forms = reduced_forms(discriminant)
        parity = discriminant % 2
        identity = (1, parity, (parity-discriminant)//4)
        generators = [form for form in forms if form[0] in (2, 3, 5, 7)]
        moves = generators+[inverse_primitive(form) for form in generators]
        reachable = {identity}
        for depth in range(5):
            if depth:
                reachable |= {compose_primitive(form, move)
                              for form in reachable for move in moves}
            search = BoundedClassWordSearch(discriminant, generators, len(forms),
                                            depth//2, depth-depth//2)
            assert search.left_statistics["pruned_class_occurrences"] == 0
            for form in forms:
                answer = search.solve(form)
                assert answer["found"] == (form in reachable)
                if answer["found"]:
                    check_meeting(answer, form, generators, identity,
                                  depth//2, depth-depth//2)
                    positive += 1
                else:
                    assert answer["status"] == "no_bounded_word"
                    assert answer["negative_answer_is_exact"]
                    negative += 1
                tested += 1
    return {"decisions_checked": tested, "positive_decisions": positive,
            "exact_negative_decisions": negative}


def trap_navigation_control():
    odd_n = 1022117
    primes, roots = first_split_primes(-4*odd_n, 3)
    generators = [reduce_primitive((prime, roots[prime],
                                    (roots[prime]**2+4*odd_n)//(4*prime)))
                  for prime in primes]
    target, identity = (17, 10, 60126), (1, 0, odd_n)
    results = []
    for width in (1, 16):
        search = BoundedClassWordSearch(-4*odd_n, generators, width, 6, 6)
        answer = search.solve(target)
        if width == 1:
            assert not answer["found"] and answer["status"] == "incomplete_search"
            assert not answer["negative_answer_is_exact"]
        else:
            check_meeting(answer, target, generators, identity, 6, 6)
            path, form = [target], target
            for generator, exponent in zip(generators, answer["word"]):
                move = inverse_primitive(generator) if exponent > 0 else generator
                for _ in range(abs(exponent)):
                    form = compose_primitive(form, move)
                    path.append(form)
            assert form == identity and path[1][0] > target[0]
            answer["verified_path_to_identity"] = path
        results.append({"width": width, "left_statistics": search.left_statistics,
                        "answer": answer})
    return {"odd_n": odd_n, "target": target, "factor_base_primes": primes,
            "word_length_bound": 12, "runs": results}


def pruning_calibration():
    records = []
    for selection in ("height", "stratified"):
        result = run(1621982935400141, 9, 256, 64, selection)
        missing = result["missing_relation_indices"]
        assert missing
        assert all(not result["targets"][index]["negative_answer_is_exact"] for index in missing)
        records.append({"selection": selection, "width": 256, "word_length_bound": 64,
                        "found_targets": [entry["index"] for entry in result["targets"] if entry["found"]],
                        "missing_relation_indices": missing,
                        "left_group_compositions": result["search"]["left"]["group_compositions"],
                        "right_group_compositions": sum(entry["right_statistics"]["group_compositions"]
                                                        for entry in result["targets"]),
                        "all_misses_reported_inconclusive": True})
    assert len(records[1]["found_targets"]) > len(records[0]["found_targets"])
    return {"odd_n": 1621982935400141, "runs": records}


def main():
    root = Path(__file__).resolve().parent.parent
    small = unpruned_controls()
    trap = trap_navigation_control()
    calibration = pruning_calibration()
    certificates = [check_certificate(root/f"evidence/bounded-class-factoring-{bits}bit.json")
                    for bits in (51, 60)]
    print(json.dumps({"scope": "Exact positive paths and unpruned decisions; bounded-width misses remain inconclusive",
                      "small_unpruned_controls": small, "height_trap_navigation": trap,
                      "pruning_calibration": calibration, "factoring_certificates": certificates,
                      "uniform_success_guarantee_claimed": False}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
