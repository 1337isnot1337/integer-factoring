#!/usr/bin/env python3
"""Exact F4 boundary retractions and a checkpointed q32 orbit catalogue.

Complete physical Smith profiles use exact F4 polynomial arithmetic.
Floating scalar minimization is exploratory. All-scale comparisons use
integer Bernstein certificates. No stronger factoring score is assumed.
"""

from collections import Counter
import ctypes
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import comb
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile
import threading
import time

from coherent_hot_compounds import height


BOUNDARY = (24, 36, 66, 129)
ROOT = Path(__file__).resolve().parents[1]
LOCK = threading.RLock()
SOURCES = ("coherent_f4_boundary.py", "coherent_f4_boundary.cpp",
           "coherent_hot_compounds_fast.cpp", "selected_gaussian_coefficients.cpp",
           "native_rectangular_smith_f4.cpp")


def multiply(a, b):
    a0, a1, b0, b1 = a & 1, a >> 1, b & 1, b >> 1
    return (a0 * b0 ^ a1 * b1) | ((a0 * b1 ^ a1 * b0 ^ a1 * b1) << 1)


def inverse(a):
    if not a:
        raise ZeroDivisionError("zero field element")
    return multiply(a, a)


def dot(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= multiply(x, y)
    return result


def boundary_maps(kernel, covector):
    kernel, covector = tuple(kernel), tuple(covector)
    if len(kernel) != 4 or len(covector) != 4 or any(x not in range(4) for x in kernel + covector):
        raise ValueError("use two vectors in F4^4")
    if dot(kernel, covector) != 1:
        raise ValueError("the boundary kernel and covector must pair to one")
    lower = sorted((w for w in range(256) if height(w) < 16), key=lambda w: (height(w), w))
    assert len(lower) == 28
    embedding = [{w: 1} for w in lower + [15]]
    retraction = [{w: 1} for w in lower + [15]]
    pivot = next(i for i, value in enumerate(covector) if value)
    for i in range(4):
        if i == pivot:
            continue
        column = {BOUNDARY[i]: 1}
        coefficient = multiply(covector[i], inverse(covector[pivot]))
        if coefficient:
            column[BOUNDARY[pivot]] = coefficient
        row = {BOUNDARY[i]: 1}
        for j in range(4):
            value = row.get(BOUNDARY[j], 0) ^ multiply(kernel[i], covector[j])
            if value:
                row[BOUNDARY[j]] = value
            else:
                row.pop(BOUNDARY[j], None)
        embedding.append(column)
        retraction.append(row)
    return embedding, retraction


def _subtract_multiple(value, basis, coefficient):
    for index, entry in basis.items():
        updated = value.get(index, 0) ^ multiply(coefficient, entry)
        if updated:
            value[index] = updated
        else:
            value.pop(index, None)


def _basis(vectors):
    pivots = {}
    for vector in vectors:
        value = dict(vector)
        while value:
            pivot = max(value)
            if pivot in pivots:
                _subtract_multiple(value, pivots[pivot], value[pivot])
            else:
                scale = inverse(value[pivot])
                pivots[pivot] = {index: multiply(scale, entry) for index, entry in value.items()}
                break
    return pivots


def validate_maps(embedding, retraction, modes):
    q = len(embedding)
    if not 1 <= modes <= 16 or q < 1 or len(retraction) != q:
        raise ValueError("incompatible physical dimensions")
    for vector in list(embedding) + list(retraction):
        if not vector or any(not 0 <= word < 1 << modes or coefficient not in (1, 2, 3)
                             for word, coefficient in vector.items()):
            raise ValueError("physical vectors must have valid nonzero coefficients")
    heights, charges = [], []
    for left, right in zip(retraction, embedding):
        words = set(left) | set(right)
        hs, ns = {height(word) for word in words}, {word.bit_count() for word in words}
        if len(hs) != 1 or len(ns) != 1:
            raise ValueError("each paired label must have one height and occupation")
        heights.append(hs.pop())
        charges.append(ns.pop())
    for i, left in enumerate(retraction):
        for j, right in enumerate(embedding):
            words = tuple(sorted(set(left) | set(right)))
            if dot([left.get(w, 0) for w in words], [right.get(w, 0) for w in words]) != int(i == j):
                raise ValueError("physical LP=I failed")
    checks = 0
    for vectors in (embedding, retraction):
        pivots = _basis(vectors)
        if len(pivots) != q:
            raise ValueError("dependent physical vectors")
        for vector in vectors:
            for bit in range(modes):
                value = {word ^ (1 << bit): coefficient for word, coefficient in vector.items()
                         if word >> bit & 1}
                while value:
                    pivot = max(value)
                    if pivot not in pivots:
                        raise ValueError("physical map space is not deletion closed")
                    _subtract_multiple(value, pivots[pivot], value[pivot])
                checks += 1
    return heights, charges, checks


@lru_cache(None)
def _library(words):
    source = Path(__file__).with_suffix(".cpp")
    dependencies = [source.with_name(name) for name in SOURCES[1:]]
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError("the exact F4 boundary engine requires C++17")
    flags = ["-O3", "-std=c++17", "-shared", "-fPIC", f"-DF4_SERIES_WORDS={words}"]
    try:
        cpu = Path("/proc/cpuinfo").read_text()
    except OSError:
        cpu = ""
    if platform.machine().lower() in ("x86_64", "amd64", "i386", "i686") and "pclmulqdq" in cpu:
        flags += ["-mpclmul", "-msse2"]
    digest = sha256(b"".join(path.read_bytes() for path in dependencies)
                    + repr((compiler, flags, platform.machine())).encode()).hexdigest()[:24]
    cache = Path(tempfile.gettempdir()) / ("integer-factoring-coherent-boundary-f4-" + str(os.getuid()))
    cache.mkdir(mode=0o700, parents=True, exist_ok=True)
    target = cache / (digest + ".so")
    if not target.exists():
        descriptor, temporary = tempfile.mkstemp(prefix=digest + "-", suffix=".so", dir=cache)
        os.close(descriptor)
        try:
            result = subprocess.run(compiler + flags + [str(source), "-o", temporary],
                                    capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr[-12000:])
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    library = ctypes.CDLL(str(target))
    u32 = ctypes.POINTER(ctypes.c_uint32)
    library.coherent_boundary_f4_block.argtypes = [ctypes.c_int, ctypes.c_int,
        u32, u32, u32, u32, u32, u32, ctypes.POINTER(ctypes.c_int64)]
    library.coherent_boundary_f4_block.restype = ctypes.c_int
    return library


def _pack(vectors):
    offsets, masks, coefficients = [0], [], []
    for vector in vectors:
        for mask, coefficient in sorted(vector.items()):
            masks.append(mask)
            coefficients.append(coefficient)
        offsets.append(len(masks))
    return ((ctypes.c_uint32 * len(offsets))(*offsets),
            (ctypes.c_uint32 * len(masks))(*masks),
            (ctypes.c_uint32 * len(coefficients))(*coefficients))


def _tensor(first, second, modes):
    return {a | (b << modes): multiply(ca, cb)
            for a, ca in first.items() for b, cb in second.items()}


def profiles(embedding, retraction, modes=8):
    heights, charges, checks = validate_maps(embedding, retraction, modes)
    q, output, metadata = len(embedding), [], []
    for infinity in (False, True):
        orders, blocks = [], []
        for charge in sorted({a + b for a in charges for b in charges}):
            pairs = [(x, y) for x in range(q) for y in range(q) if charges[x] + charges[y] == charge]
            rows = [_tensor(retraction[x], embedding[y] if infinity else retraction[y], modes)
                    for x, y in pairs]
            columns = [_tensor(retraction[y] if infinity else embedding[y], embedding[x], modes)
                       for x, y in pairs]
            count = len(pairs)
            coefficient_bound = max(heights[x] + heights[y] for x, y in pairs)
            words = 2
            while 64 * words <= coefficient_bound:
                words *= 2
            attempts = []
            while words <= 128:
                result = (ctypes.c_int64 * count)()
                with LOCK:
                    status = _library(words).coherent_boundary_f4_block(
                        modes, count, *_pack(rows), *_pack(columns), result)
                attempts.append([64 * words, status])
                if status == 0:
                    break
                if status != 6:
                    raise ArithmeticError(f"F4 block failed with status {status}")
                words *= 2
            else:
                raise ArithmeticError("complete F4 profile exceeded the series budget")
            orders.extend(map(int, result))
            blocks.append({"charge": charge, "size": count, "precision_attempts": attempts,
                           "endpoint_entries": count * count})
        assert len(orders) == q * q
        output.append(sorted(orders))
        metadata.append(blocks)
    return {"heights": heights, "charges": charges, "zero_orders_t": output[0],
            "infinity_orders_t": output[1], "deletion_checks": checks,
            "retraction_entries": q * q, "blocks": metadata}


def frobenius(key):
    return tuple(multiply(x, x) for x in key)


def transpose(key):
    kernel, covector = key[:4], key[4:]
    first = next(value for value in covector if value)
    return (tuple(multiply(value, inverse(first)) for value in covector)
            + tuple(multiply(value, first) for value in kernel))


def orbit(key):
    return {tuple(key), frobenius(key), transpose(key), frobenius(transpose(key))}


def catalogue():
    keys = []
    for pivot in range(4):
        for tail in product(range(4), repeat=3 - pivot):
            kernel = (0,) * pivot + (1,) + tail
            other = [i for i in range(4) if i != pivot]
            for values in product(range(4), repeat=3):
                covector = [0] * 4
                for i, value in zip(other, values):
                    covector[i] = value
                covector[pivot] = 1 ^ dot(kernel, covector)
                keys.append(kernel + tuple(covector))
    assert len(set(keys)) == 5440
    representatives = {min(orbit(key)) for key in keys}
    assert len(representatives) == 1416
    best = min(orbit((1, 1, 1, 0, 1, 1, 1, 1)))
    return sorted(representatives, key=lambda key: (
        min(sum(a != b for a, b in zip(value, best)) for value in orbit(key)), key))


def counting_record(record):
    return {"heights": record["heights"],
            "zero_counts": sorted(Counter(record["zero_orders_t"]).items()),
            "infinity_counts": sorted(Counter(record["infinity_orders_t"]).items())}


def expand_counts(counts):
    return [order for order, count in counts for _ in range(count)]


def pressure(record):
    import numpy as np
    from scipy.optimize import brentq, minimize_scalar
    from scipy.special import logsumexp
    zero, infinity, heights = (np.asarray(record[name]) for name in
                               ("zero_orders_t", "infinity_orders_t", "heights"))
    def at_scale(scale):
        log_pressure = float(logsumexp(-scale * zero) + logsumexp(-scale * infinity))
        return brentq(lambda rho: log_pressure - 4 * rho * float(logsumexp(-scale * heights / rho)),
                      1e-7, 2, xtol=1e-13)
    result = minimize_scalar(at_scale, bounds=(.05, .8), method="bounded",
                             options={"xatol": 1e-12})
    return {"rho": float(result.fun), "log_s": float(result.x),
            "scope": "exploratory minimization on 0.05 <= log(s) <= 0.8"}


def _positive_bernstein(difference):
    nonzero = [index for index, value in enumerate(difference) if value]
    if not nonzero:
        return {"comparison": "equal"}
    shift = nonzero[0]
    polynomial = list(difference[shift:nonzero[-1] + 1])
    multiplicity = 0
    while sum(polynomial) == 0:
        quotient, running = [], 0
        for value in polynomial[:-1]:
            running += value
            quotient.append(running)
        assert running + polynomial[-1] == 0
        polynomial = quotient
        multiplicity += 1
    original_degree = len(polynomial) - 1
    for degree in sorted({original_degree, original_degree + 16, 2 * original_degree}):
        bernstein = [sum(polynomial[i] * comb(degree - i, k - i)
                         for i in range(min(k, original_degree) + 1))
                     for k in range(degree + 1)]
        if not all(value >= 0 for value in bernstein):
            continue
        assert any(bernstein)
        rebuilt = [sum(bernstein[k] * comb(degree - k, j - k) * (-1) ** (j - k)
                       for k in range(j + 1)) for j in range(degree + 1)]
        assert rebuilt == polynomial + [0] * (degree - original_degree)
        return {"comparison": "strictly positive on 0 < r < 1", "initial_power": shift,
                "one_minus_r_multiplicity": multiplicity, "quotient_coefficients": polynomial,
                "bernstein_degree": degree}
    return {"comparison": "not certified by the retained Bernstein degrees"}


def compare(reference, candidate):
    def counting_product(record):
        answer = Counter()
        for a, ca in Counter(record["zero_orders_t"]).items():
            for b, cb in Counter(record["infinity_orders_t"]).items():
                answer[a + b] += ca * cb
        return answer
    a, b = counting_product(reference), counting_product(candidate)
    assert min(a) >= 0 and min(b) >= 0
    degree = max(max(a), max(b))
    return _positive_bernstein([b[i] - a[i] for i in range(degree + 1)])


def run(output, maximum=1416, seconds=7200, resume=False):
    output = Path(output)
    representatives = catalogue()
    catalogue_hash = sha256(json.dumps(representatives, separators=(",", ":")).encode()).hexdigest()
    source_hashes = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest() for name in SOURCES}
    reference_key = min(orbit((1, 1, 1, 0, 1, 1, 1, 1)))
    reference = profiles(*boundary_maps(reference_key[:4], reference_key[4:]))
    reference_pressure = pressure(reference)
    assert abs(reference_pressure["rho"] - .712115012500715) < 1e-12
    if resume:
        report = json.loads(output.read_text())
        assert report["catalogue_sha256"] == catalogue_hash and report["source_sha256"] == source_hashes
        assert [tuple(record["key"]) for record in report["records"]] == representatives[:len(report["records"])]
    else:
        report = {"schema": 1, "complete": False, "projectors": 5440, "orbits": 1416,
                  "catalogue_sha256": catalogue_hash, "source_sha256": source_hashes,
                  "reference_key": reference_key, "reference": counting_record(reference),
                  "reference_pressure": reference_pressure, "profiles": {}, "records": [],
                  "in_progress": None, "factoring_exponent_improved": False}
    started = time.monotonic()
    previous_seconds = report.get("elapsed_seconds", 0)
    best = min([reference_pressure["rho"]] + [p["pressure"]["rho"] for p in report["profiles"].values()])
    def save():
        report["elapsed_seconds"] = previous_seconds + time.monotonic() - started
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(output.suffix + ".tmp")
        temporary.write_text(json.dumps(report, separators=(",", ":")) + "\n")
        os.replace(temporary, output)
    save()
    for index in range(len(report["records"]), min(maximum, len(representatives))):
        if time.monotonic() - started >= seconds:
            break
        key = representatives[index]
        report["in_progress"] = key
        save()
        result = reference if key == reference_key else profiles(*boundary_maps(key[:4], key[4:]))
        encoded = counting_record(result)
        digest = sha256(json.dumps(encoded, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        if digest not in report["profiles"]:
            encoded["pressure"] = pressure(result)
            encoded["comparison_to_reference"] = compare(reference, result)
            report["profiles"][digest] = encoded
        value = report["profiles"][digest]["pressure"]["rho"]
        report["records"].append({"key": key, "orbit_size": len(orbit(key)), "profile": digest,
                                   "deletion_checks": result["deletion_checks"], "blocks": result["blocks"]})
        report["in_progress"] = None
        if value < best - 1e-13:
            best = value
            report["best_key"] = key
            print(json.dumps({"new_best": best, "key": key, "completed": index + 1}), flush=True)
        if (index + 1) % 20 == 0:
            print(json.dumps({"completed": index + 1, "orbits": 1416, "best": best,
                              "distinct_profiles": len(report["profiles"]),
                              "seconds": time.monotonic() - started}), flush=True)
        save()
    report["complete"] = len(report["records"]) == len(representatives)
    report["completed_orbits"] = len(report["records"])
    report["covered_projectors"] = sum(record["orbit_size"] for record in report["records"])
    save()
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-orbits", type=int, default=1416)
    parser.add_argument("--seconds", type=float, default=7200)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.max_orbits < 1 or args.seconds <= 0:
        raise ValueError("use positive catalogue and time budgets")
    report = run(args.output, args.max_orbits, args.seconds, args.resume)
    print(json.dumps({key: report[key] for key in ("complete", "completed_orbits", "covered_projectors")}, indent=2))


if __name__ == "__main__":
    main()
