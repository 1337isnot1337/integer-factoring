#!/usr/bin/env python3
"""Bounded non-diagonal F8 probes of the five-direction cubic family."""

import argparse
import hashlib
import heapq
from itertools import product
import json
from math import isqrt
from pathlib import Path
import random
import subprocess

from explore_cubic_complement_outputs import ComplementOutputs


def multiply(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 11
    return value


TABLE = tuple(tuple(multiply(a, b) for b in range(8)) for a in range(8))
INVERSE = (0,) + tuple(next(b for b in range(1, 8) if TABLE[a][b] == 1) for a in range(1, 8))


def scale(value, scalar, width):
    mask = ((1 << (3 * width)) - 1) // 7
    result = 0
    for i in range(3):
        plane = (value >> i) & mask
        coefficient = TABLE[scalar][1 << i]
        for j in range(3):
            if coefficient >> j & 1:
                result ^= plane << j
    return result


def rank8(rows, width):
    pivots = {}
    for value in rows:
        while value:
            pivot = ((value & -value).bit_length() - 1) // 3
            scalar = (value >> (3 * pivot)) & 7
            if pivot not in pivots:
                pivots[pivot] = scale(value, INVERSE[scalar], width)
                break
            value ^= scale(pivots[pivot], scalar, width)
    return len(pivots)


def matrices(family, complement):
    d = family.pairs
    output = [[(i, 1)] + [(d + j, value) for j, value in enumerate(complement[i]) if value]
              for i in range(d)]
    planes = [[sum(((value >> k) & 1) << j for j, value in row) for k in range(3)]
              for row in output]
    rows = [[0] * (d * d) for _ in range(5)]
    columns = [[0] * (d * d) for _ in range(5)]
    for k, y, v in product(range(5), range(d), range(d)):
        parent = family.parent[k][y][v]
        for x in range(d):
            acted = [0, 0, 0]
            for index, scalar in output[x]:
                for b in range(3):
                    if scalar >> b & 1:
                        acted[b] ^= parent[index]
            for u in range(d):
                value = 0
                for a, b in product(range(3), repeat=2):
                    if (acted[a] & planes[u][b]).bit_count() & 1:
                        value ^= TABLE[1 << a][1 << b]
                rows[k][d * x + y] ^= value << (3 * (d * v + u))
                columns[k][d * x + u] ^= value << (3 * (d * v + y))
    assert rows[0] == columns[0] == [1 << (3 * i) for i in range(d * d)]
    gram = [[int(x == u) for u in range(d)] for x in range(d)]
    for x, u, j in product(range(d), repeat=3):
        gram[x][u] ^= TABLE[complement[x][j]][complement[u][j]]
    for x, y, v, u in product(range(d), repeat=4):
        high = TABLE[complement[x][y]][complement[u][v]]
        assert ((rows[4][d * x + y] >> (3 * (d * v + u))) & 7) == high
        assert (((rows[1][d * x + y] ^ rows[3][d * x + y]) >> (3 * (d * v + u))) & 7) == (high ^ (gram[x][u] if y == v else 0))
    return rows, columns


def native_input(rows, columns):
    lines = ["7 5 16"]
    for cut in (rows, columns):
        for matrix in cut:
            for row in matrix:
                planes = [sum(((row >> (3 * j + b)) & 1) << j for j in range(49)) for b in range(3)]
                lines.append(" ".join(format(plane, "x") for plane in planes))
    return "\n".join(lines) + "\n"


def affordable_decision(records):
    pivots = {}
    for record in records:
        value = record["code"]
        while value:
            pivot = (value.bit_length() - 1) // 3
            scalar = (value >> (3 * pivot)) & 7
            if pivot not in pivots:
                pivots[pivot] = scale(value, INVERSE[scalar], 5)
                break
            value ^= scale(pivots[pivot], scalar, 5)
    residual = 1
    while residual:
        pivot = (residual.bit_length() - 1) // 3
        if pivot not in pivots:
            break
        residual ^= scale(pivots[pivot], (residual >> (3 * pivot)) & 7, 5)
    result = {"affordable_span_dimension": len(pivots), "endpoint_remainder": residual}
    if residual:
        return dict(result, status="excluded: affordable span misses the crossing")
    denominator = 10 ** 12
    weights, products = {}, {}
    for record in records:
        for scalar in range(1, 8):
            value = scale(record["code"], scalar, 5)
            weights[value] = isqrt(record["product"] * denominator * denominator)
            products[value] = record["product"]
    distances, previous = {0: 0}, {}
    queue = [(0, 0)]
    cap = 16 * denominator
    while queue:
        cost, state = heapq.heappop(queue)
        if cost != distances[state]:
            continue
        if cost >= cap:
            return dict(result, status="excluded: integer cost gate", cost_lower_numerator=cap,
                        cost_denominator=denominator)
        if state == 1:
            witness = []
            while state:
                old, value = previous[state]
                witness.append(value)
                state = old
            upper = sum(weights[value] + int(weights[value] ** 2 != products[value] * denominator ** 2)
                        for value in witness)
            if upper ** 985 < 7 ** 1393 * denominator ** 985:
                status = "strict seed candidate: independent full tensor audit required"
            elif cost ** 408 > 7 ** 577 * denominator ** 408:
                status = "excluded: shortest cost exceeds the threshold"
            else:
                status = "unresolved exact threshold interval"
            return dict(result, status=status, witness=witness, cost_lower_numerator=cost,
                        witness_upper_numerator=upper, cost_denominator=denominator)
        for value, weight in weights.items():
            successor, candidate = state ^ value, cost + weight
            if candidate < distances.get(successor, candidate + 1):
                distances[successor] = candidate
                previous[successor] = state, value
                heapq.heappush(queue, (candidate, successor))
    raise AssertionError("Span membership and shortest-path reachability disagree.")


def verify_screen(rows, columns, text, seed):
    events = [json.loads(line) for line in text.splitlines()]
    start, complete = events[0], events[-1]
    assert start == {"event": "start", "labels": 7, "dimension": 5, "field_size": 8,
                     "field_polynomial": 11, "integer_cost_gate": 16, "product_cap": 255}
    assert complete["event"] == "complete" and complete["projective_tested"] == 4681
    records = [event for event in events if event["event"] == "eligible"]
    by_code = {record["code"]: record for record in records}
    assert len(records) == len(by_code) == complete["eligible"]
    for record in records:
        code = record["code"]
        assert code >> (3 * ((code.bit_length() - 1) // 3)) == 1
        assert record["product"] == record["row_rank"] * record["column_rank"] < 256
    selected = set(by_code)
    rng = random.Random(seed)
    for _ in range(3):
        code = rng.randrange(1, 8 ** 5)
        pivot = (code.bit_length() - 1) // 3
        selected.add(scale(code, INVERSE[(code >> (3 * pivot)) & 7], 5))
    for code in selected:
        ranks = []
        for cut in (rows, columns):
            matrix = [0] * 49
            for i, source in enumerate(cut):
                scalar = (code >> (3 * i)) & 7
                if scalar:
                    matrix = [a ^ scale(b, scalar, 49) for a, b in zip(matrix, source)]
            ranks.append(rank8(matrix, 49))
        assert (code in by_code) == (ranks[0] * ranks[1] < 256)
        if code in by_code:
            assert ranks == [by_code[code]["row_rank"], by_code[code]["column_rank"]]
    return dict(affordable_decision(records), projective_vectors_tested=4681,
                affordable_count=len(records), independently_checked_vectors=len(selected),
                native_screen_sha256=hashlib.sha256(text.encode()).hexdigest())


def candidates(seed, scalars):
    d = 7
    powers = [1]
    for _ in range(6):
        powers.append(TABLE[powers[-1]][2])
    assert len(set(powers)) == 7 and TABLE[powers[-1]][2] == 1
    yield {"family": "zero_control"}, [[0] * d for _ in range(d)]
    yield {"family": "identity_control"}, [[int(i == j) for j in range(d)] for i in range(d)]
    for scalar in scalars:
        for subset in range(1, 127):
            coefficients = []
            for j in range(7):
                value = 0
                for k in range(7):
                    if subset >> k & 1:
                        value ^= powers[(-k * j) % 7]
                coefficients.append(TABLE[scalar][value])
            matrix = [[coefficients[(j - i) % 7] for j in range(7)] for i in range(7)]
            yield {"family": "scaled_order7_spectral_projector", "eigenvalue_subset": subset,
                   "scalar": scalar}, matrix
    for scalar in range(8):
        yield {"family": "order7_shift_plus_scalar", "scalar": scalar}, [
            [int(j == (i + 1) % 7) ^ (scalar if i == j else 0) for j in range(7)] for i in range(7)]
        yield {"family": "jordan_plus_scalar", "scalar": scalar}, [
            [int(j == i + 1) ^ (scalar if i == j else 0) for j in range(7)] for i in range(7)]
    rng = random.Random(seed)
    for rank in (1, 2):
        for case in range(16):
            matrix = [[0] * d for _ in range(d)]
            for _ in range(rank):
                first, second = [[rng.randrange(8) for _ in range(d)] for _ in range(2)]
                for i, j in product(range(d), repeat=2):
                    matrix[i][j] ^= TABLE[first[i]][second[j]]
            yield {"family": "random_rank_at_most_" + str(rank), "case": case}, matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--spectral-scalars", default="1,2,4")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.replay and not args.native:
        parser.error("--native is required for a new enumeration")
    family = ComplementOutputs(7)
    source = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in (Path(__file__).name, "screen_cubic_gf8.cpp", "explore_cubic_complement_outputs.py",
                           "explore_cubic_multirow_outputs.py")}
    result = {"field": "F2[z]/(z^3+z+1)", "labels": 7, "formal_coefficient_coordinates": 5,
              "seed": args.seed, "source_sha256": source, "cases": [], "complete": False}
    if args.replay:
        saved = json.loads(args.replay.read_text())
        assert saved["complete"] and saved["source_sha256"] == source
        jobs = [(dict(record), record["complement_matrix"]) for record in saved["cases"]]
    else:
        jobs = candidates(args.seed, list(map(int, args.spectral_scalars.split(","))))
    seen = set()
    for descriptor, complement in jobs:
        key = tuple(map(tuple, complement))
        if key in seen:
            continue
        seen.add(key)
        rows, columns = matrices(family, complement)
        payload = native_input(rows, columns)
        if args.replay:
            screen = descriptor["native_screen"]
        else:
            screen = subprocess.run([str(args.native)], input=payload, text=True,
                                    stdout=subprocess.PIPE, check=True).stdout
        decision = verify_screen(rows, columns, screen, args.seed + len(result["cases"]))
        if args.replay:
            assert all(descriptor[k] == value for k, value in decision.items())
            assert descriptor["native_input_sha256"] == hashlib.sha256(payload.encode()).hexdigest()
        record = dict(descriptor, complement_matrix=complement, native_screen=screen,
                      native_input_sha256=hashlib.sha256(payload.encode()).hexdigest(), **decision)
        result["cases"].append(record)
        if len(result["cases"]) % 20 == 0:
            args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
            print(json.dumps({"cases": len(result["cases"]), "latest_status": decision["status"]}), flush=True)
    result["complete"] = True
    result["all_recorded_spaces_excluded_over_F8"] = all(r["status"].startswith("excluded:") for r in result["cases"])
    result["projective_vectors_tested"] = sum(r["projective_vectors_tested"] for r in result["cases"])
    result["independently_checked_vectors"] = sum(r["independently_checked_vectors"] for r in result["cases"])
    result["scope"] = "The recorded complementary matrices over F8; all formal five-coordinate coefficient vectors are enumerated for each matrix."
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("cases", "source_sha256")}), flush=True)


if __name__ == "__main__":
    main()
