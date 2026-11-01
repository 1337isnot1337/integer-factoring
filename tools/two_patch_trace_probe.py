#!/usr/bin/env python3
"""Exact polynomial crossings from two bilinear trace patches."""
import argparse,json,math,os,time
from pathlib import Path
from check_trace_bilinear_kernel import Field,rank,irreducibility_by_division


class CirculantRank:
    def __init__(self,t):
        self.t=t
        self.d=1<<t
        self.lower=[sum(1<<x for x in range(self.d) if not(x>>i&1)) for i in range(t)]
        self.cache={0:0}
    def translate(self,word,x):
        for i in range(self.t):
            if x>>i&1:
                delta=1<<i
                mask=self.lower[i]
                word=((word&mask)<<delta)|((word&~mask)>>delta)
        return word
    def __call__(self,kernel):
        if kernel in self.cache:
            return self.cache[kernel]
        if kernel.bit_count()%2:
            answer=self.d
        elif any((kernel&~m).bit_count()%2 for m in self.lower):
            answer=self.d//2
        else:
            answer=rank([self.translate(kernel,x) for x in range(self.d)])
        self.cache[kernel]=answer
        return answer


def profile(t,polynomial):
    field=Field(polynomial)
    d=1<<t
    assert field.degree==t
    powers=[field.power(field.alpha,j) for j in range(t)]
    syndromes=[sum(field.trace(field.mul(a,z))<<j for j,a in enumerate(powers)) for z in range(d)]
    assert len(set(syndromes))==d
    squares=[field.mul(s,s) for s in range(d)]
    rows=[[0]*d for _ in range(2*t+1)]
    columns=[[0]*d for _ in range(2*t+1)]
    for s in range(d):
        for w in range(d):
            b0=field.mul(s,w)^squares[s]
            b1=b0^s^w
            weight=syndromes[b0].bit_count()+syndromes[b1].bit_count()
            assert (weight==0)==(s==w)
            k=weight
            while True:
                rows[k][s]|=1<<w
                columns[k][w]|=1<<s
                if k==0:
                    break
                k=(k-1)&weight
    ranker=CirculantRank(t)
    r=[sum(ranker(v) for v in block) for block in rows]
    c=[sum(ranker(v) for v in block) for block in columns]
    products=[a*b for a,b in zip(r,c)]
    assert all(bool(a)==bool(b) for a,b in zip(r,c))
    rank_sum=sum(math.sqrt(p) for p in products)
    return {"t":t,"field_polynomial":polynomial,"row_ranks":r,"column_ranks":c,
            "rank_sum_floor":sum(math.isqrt(p) for p in products),
            "rank_sum_estimate":rank_sum,"maximum_rank_product":max(products),
            "rank_sum_rho":math.log2(rank_sum)/(2*t),
            "maximum_product_rho":math.log2(max(products))/(4*t),
            "cached_circulant_kernels":len(ranker.cache)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths",nargs="+",type=int,default=[2,3,4,5,6])
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    results=[]
    for t in args.lengths:
        polynomial=(1<<t)|1
        while True:
            try:
                irreducibility_by_division(polynomial)
                break
            except AssertionError:
                polynomial+=2
        start=time.monotonic()
        result=profile(t,polynomial)
        result["elapsed_seconds"]=time.monotonic()-start
        results.append(result)
        print(json.dumps(result),flush=True)
        temporary=args.output.with_suffix(args.output.suffix+".tmp")
        with temporary.open("w") as out:
            json.dump(results,out,indent=2);out.write("\n");out.flush();os.fsync(out.fileno())
        temporary.replace(args.output)
    assert len(json.loads(args.output.read_text()))==len(args.lengths)


if __name__=="__main__":
    main()
