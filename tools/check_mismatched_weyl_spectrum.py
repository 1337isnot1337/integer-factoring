#!/usr/bin/env python3
"""Replay the all-polynomial obstruction for the F2, radix-three h=63 curve.

Only Python's standard library and a C++17 compiler are needed.  Native block
dumps are rebuilt by default; --prefix accepts existing dumps for a fast replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile
from time import monotonic


def rank(rows):
    pivots = {}
    for value in rows:
        while value:
            p = value.bit_length()-1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return len(pivots)


def product(A, B):
    answer = []
    for row in A:
        value = 0
        while row:
            bit = row & -row
            value ^= B[bit.bit_length()-1]
            row ^= bit
        answer.append(value)
    return answer


def read_dump(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"WBIN001\0"
        h,d,nb,pw,nw,cut = struct.unpack("<6I", f.read(24))
        assert pw == 1 and cut == 1
        weights = struct.unpack("<%dI" % nw, f.read(4*nw))
        labels = [struct.unpack("<%dI" % nw, f.read(4*nw)) for _ in range(d)]
        blocks = []
        for _ in range(nb):
            n = struct.unpack("<I", f.read(4))[0]
            pairs = struct.unpack("<%dI" % n, f.read(4*n))
            polys = struct.unpack("<%dQ" % (n*n), f.read(8*n*n))
            blocks.append((pairs, polys))
        assert not f.read()
    return h,d,weights,labels,blocks


def certify(prefix):
    h,d,weights,labels,blocks = read_dump(str(prefix)+".row.bin")
    assert (h,d,weights) == (63,132,(1,3,9,27))
    assert all(sum(w*a for w,a in zip(weights,label)) == h for label in labels)
    assert len(set(labels)) == d
    pair_codes = [code for pairs,_ in blocks for code in pairs]
    assert sorted(pair_codes) == list(range(d*d))
    for pairs,_ in blocks:
        differences = {tuple(a-b for a,b in zip(labels[code//d],labels[code%d]))
                       for code in pairs}
        assert len(differences) == 1
    ranks = {1:0, 32:0, 64:0}
    zero_rows = 0
    witnesses = {}
    endpoint_entries = 0
    for pairs,polys in blocks:
        n = len(pairs)
        A = []
        for i,code in enumerate(pairs):
            row = 0
            for j in range(n):
                polynomial = polys[i*n+j]
                assert polynomial.bit_count() % 2 == (i == j)
                row |= (polynomial & 1) << j
                endpoint_entries += 1
            A.append(row)
            if not row:
                zero_rows += 1
                alpha,y = divmod(code,d)
                sigma = tuple(a+b for a,b in zip(labels[alpha],labels[y]))
                witnesses.setdefault(sigma,(alpha,y))
        ranks[1] += rank(A)
        power = A
        for _ in range(5):
            power = product(power,power)
        ranks[32] += rank(power)
        ranks[64] += rank(product(power,power))
    assert len(blocks) == 1341 and endpoint_entries == 632732
    assert ranks == {1:11924,32:7280,64:7280}
    assert zero_rows == 2531 and len(witnesses) == 496
    row_floor = d*d-ranks[32]
    col_floor = len(witnesses)
    assert row_floor*col_floor > d**3
    payload = json.dumps(sorted((sigma,pair) for sigma,pair in witnesses.items()),separators=(",",":"))
    return {
        "energy":h,"labels":d,"weights":weights,
        "row_blocks":len(blocks),"endpoint_entries_checked":endpoint_entries,
        "power_ranks":ranks,"zero_rows":zero_rows,
        "protected_column_sum_blocks":col_floor,
        "nilpotent_kernel_dimension":row_floor,
        "same_piece_rank_product_lower":row_floor*col_floor,
        "dimension_cubed":d**3,
        "witness_sha256":hashlib.sha256(payload.encode()).hexdigest(),
        "dump_sha256":hashlib.sha256(Path(str(prefix)+".row.bin").read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix")
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    here = Path(__file__).resolve().parent
    source = here/"explore_mismatched_weyl.cpp"
    with tempfile.TemporaryDirectory(prefix="mismatch-weyl-") as temporary:
        temporary = Path(temporary)
        if args.prefix:
            prefix = Path(args.prefix)
        else:
            executable = temporary/"source"
            subprocess.run(["g++","-std=c++17","-O2",str(source),"-o",str(executable)],check=True)
            prefix = temporary/"h63"
            native = subprocess.run([str(executable),"63","--radix","3","--dump-prefix",str(prefix)],check=True,capture_output=True,text=True)
            raw = json.loads(native.stdout)
            assert raw["labels"] == 132
        result = certify(prefix)
    result["source_sha256"] = {
        path.name:hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__),source)
    }
    result["seconds"] = round(monotonic()-started,3)
    result["scope"] = "All decompositions of I inside K[A0], over every characteristic-two extension K, at the one fixed shell h=63; not the full coefficient or word algebra."
    text = json.dumps(result,indent=2)+"\n"
    if args.output:
        Path(args.output).write_text(text)
    print(text,end="")


if __name__ == "__main__":
    main()
