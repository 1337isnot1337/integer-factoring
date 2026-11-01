#!/usr/bin/env python3
"""Independent four-leg checks for syndrome ranks, with an optional finite screen."""
import argparse
import hashlib
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from check_trace_bilinear_kernel import Field,rank


def source_path():
    return Path(__file__).resolve().with_name("trace_weight_restriction.cpp")


def direct_restricted(t,m,k,polynomial):
    field=Field(polynomial)
    labels=[x for x in range(1<<t) if x.bit_count()==m]
    index={x:i for i,x in enumerate(labels)}
    n=len(labels)
    values=[]
    alpha=[field.power(field.alpha,j) for j in range(k)]
    for z in range(field.order):
        value=1
        for a in alpha:
            value&=field.trace(field.mul(a,z))
        values.append(value)
    rows=[0]*(n*n)
    columns=[0]*(n*n)
    for i,x in enumerate(labels):
        for j,y in enumerate(labels):
            s=x^y
            for l,u in enumerate(labels):
                v=x^y^u
                if v not in index or not values[field.mul(s,x^u)]:
                    continue
                q=index[v]
                rows[i*n+j]|=1<<(l*n+q)
                columns[i*n+l]|=1<<(j*n+q)
    return rank(rows),rank(columns)


def cap(t,m,k):
    total=0
    for s in range(t//2+1):
        if 0<=m-s<=t-2*s:
            block=math.comb(2*s,s)*math.comb(t-2*s,m-s)
            total+=math.comb(t,2*s)*min(block,1<<k)
    return total


def binary_product(a,b):
    out=[]
    for row in a:
        value=0
        while row:
            bit=row&-row
            row^=bit
            value^=b[bit.bit_length()-1]
        out.append(value)
    return out


def kernel(field,k):
    alpha=[field.power(field.alpha,j) for j in range(k)]
    values=[]
    for z in range(field.order):
        v=1
        for a in alpha:
            v&=field.trace(field.mul(a,z))
        values.append(v)
    return [sum(values[field.mul(s,w)]<<w for w in range(field.order))
            for s in range(field.order)]


def composition_checks():
    checked=0
    for t,p in [(3,11),(4,19),(5,37),(6,67)]:
        field=Field(p)
        matrices={k:kernel(field,k) for k in range(1,t)}
        for k in matrices:
            for l in matrices:
                if k+l>=t:
                    continue
                assert not any(binary_product(matrices[k],matrices[l]))
                checked+=1
                for a in range(t-k-l):
                    mask=sum(1<<w for w in range(field.order)
                             if w&((1<<a)-1)==((1<<a)-1))
                    controlled=[row&mask for row in matrices[k]]
                    assert not any(binary_product(controlled,matrices[l]))
                    checked+=1
    return checked


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full",action="store_true")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    cpp=source_path()
    small=[]
    full=[]
    with tempfile.TemporaryDirectory(prefix="trace-weight-check-") as temp:
        folder=Path(temp)
        executable=folder/"trace_weight_restriction"
        subprocess.run(["g++","-O3","-std=c++17",str(cpp),"-o",str(executable)],check=True)
        for t,m,kmax in [(4,2,3),(7,2,6)]:
            path=folder/f"small-{t}-{m}.json"
            subprocess.run([str(executable),str(t),str(m),str(kmax),str(path)],check=True,capture_output=True,text=True)
            result=json.loads(path.read_text())
            assert result["label_count"]==math.comb(t,m)
            for profile in result["profiles"]:
                actual=direct_restricted(t,m,profile["k"],result["field_polynomial"])
                assert actual==(profile["row_rank"],profile["column_rank"])
                assert profile["block_size_cap"]==cap(t,m,profile["k"])
            small.append(result)
        if args.full:
            for t,m in [(12,3),(12,4),(12,5),(14,4),(14,5),(14,6)]:
                path=folder/f"full-{t}-{m}.json"
                subprocess.run([str(executable),str(t),str(m),"6",str(path)],check=True,capture_output=True,text=True)
                result=json.loads(path.read_text())
                assert result["label_count"]==math.comb(t,m)
                assert len(result["profiles"])==6
                for profile in result["profiles"]:
                    assert profile["block_size_cap"]==cap(t,m,profile["k"])
                    assert profile["row_rank"]==profile["column_rank"]<=profile["block_size_cap"]
                full.append(result)
    report={"status":"passed","native_source_sha256":hashlib.sha256(cpp.read_bytes()).hexdigest(),
            "independent_restricted_profiles":sum(len(r["profiles"]) for r in small),
            "small_results":small,"ordinary_composition_checks":composition_checks(),
            "full_screen":args.full,"screen_profiles":sum(len(r["profiles"]) for r in full),
            "screen_results":full,
            "scope":"Exact ranks of individual restricted pieces; no crossing synthesis."}
    temporary=args.output.with_suffix(args.output.suffix+".tmp")
    with temporary.open("w") as out:
        json.dump(report,out,indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    temporary.replace(args.output)
    print(json.dumps({k:v for k,v in report.items() if not k.endswith("results")},indent=2))


if __name__=="__main__":
    main()
