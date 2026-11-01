#!/usr/bin/env python3
"""Independent dense F2 elimination of two five-label twist controls."""
import argparse
import hashlib
import json
from pathlib import Path


def dense_rank(matrix):
    a=[row.copy() for row in matrix]
    r=0
    for col in range(len(a[0]) if a else 0):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:
            continue
        a[r],a[pivot]=a[pivot],a[r]
        for i in range(r+1,len(a)):
            if a[i][col]:
                for j in range(col,len(a[i])):
                    a[i][j]=(a[i][j]+a[r][j])%2
        r+=1
    return r


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    d=5
    ident=[[int(i==j) for j in range(d)] for i in range(d)]
    j01=[[int((i,j) in ((0,1),(1,0))) for j in range(d)] for i in range(d)]
    j23=[[int((i,j) in ((2,3),(3,2))) for j in range(d)] for i in range(d)]
    cases=[("catalogue_attaining_E_plus_R_plus_C",[(ident,ident)],[11,11]),
           ("outside_catalogue_disjoint_pair_row_updates",[(j01,j01),(j23,j23)],[8,18])]
    results=[]
    for name,row_terms,expected in cases:
        def tensor(x,y,u,v):
            value=ident[x][v]*ident[y][u]+ident[x][u]*ident[y][v]
            value+=sum(a[x][y]*b[u][v] for a,b in row_terms)
            return value%2
        rows=[[tensor(x,y,u,v) for u in range(d) for v in range(d)]
              for x in range(d) for y in range(d)]
        cols=[[tensor(x,y,u,v) for y in range(d) for v in range(d)]
              for x in range(d) for u in range(d)]
        ranks=[dense_rank(rows),dense_rank(cols)]
        assert ranks==expected
        results.append({"case":name,"cut_ranks":ranks,"rank_product":ranks[0]*ranks[1],
                        "row_matrices":[{"A":a,"B":b} for a,b in row_terms],
                        "column_matrices":{"C":ident,"D":ident}})
    evidence={"field":2,"external_labels":d,"method":"independent dense Gaussian elimination",
              "results":results,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(evidence,indent=2)+"\n")
    print(json.dumps({"cases":[{k:v for k,v in item.items() if k in ("case","cut_ranks","rank_product")}
                              for item in results],"source_sha256":evidence["source_sha256"]},indent=2))


if __name__=="__main__":
    main()
