#!/usr/bin/env python3
"""Exact integer interpolation and literal quadratic Galois-ring controls."""
import argparse
from hashlib import sha256
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/precision-residue-indicator-controls.json'


def add(p,q):
    out = [0]*max(len(p),len(q))
    for i,c in enumerate(p): out[i] += c
    for i,c in enumerate(q): out[i] += c
    return out


def mul(p,q):
    out = [0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): out[i+j] += a*b
    return out


def polynomial_power(p,n):
    out = [1]
    for _ in range(n): out = mul(out,p)
    return out


def smoothstep(L):
    return [0]*L + [(-1)**(k-L)*comb(2*L-1,k)*comb(k-1,L-1)
                     for k in range(L,2*L)]


def evaluate(p,x,modulus):
    out = 0
    for c in reversed(p): out = (out*x+c)%modulus
    return out


def integer_controls():
    records = []
    scalar_count = 0
    for L in range(1,25):
        p = smoothstep(L)
        expanded = [0]*(2*L)
        for j in range(L,2*L):
            term = [0]*j + [comb(2*L-1,j)*c
                            for c in polynomial_power([1,-1],2*L-1-j)]
            expanded = add(expanded,term)
        assert expanded == p
        complemented = [0]*(2*L)
        for j,c in enumerate(p):
            complemented = add(complemented,[c*v for v in polynomial_power([1,-1],j)])
        expected = [1]+[-c for c in p[1:]]
        assert complemented == expected
        assert all(c==0 for c in p[:L])
        assert max(abs(c).bit_length() for c in p) <= 4*L
        inv = smoothstep(L+1)[1:]
        assert len(inv)-1 == 2*L and not any(inv[:L])
        modulus = 1 << min(L,10)
        for x in range(modulus):
            assert evaluate(p,x,modulus) == x%2
            value = evaluate(inv,x,modulus)
            assert value == (pow(x,-1,modulus) if x%2 else 0)
            scalar_count += 1
        records.append(dict(precision=L,indicator_coefficients=p,
                            inverse_coefficients=inv,checked_modulus=modulus))
    return dict(exact_polynomial_identities=len(records),scalar_residues=scalar_count,
                coefficients=records)


class QuadraticRing:
    """Pairs a+b*omega modulo 2^L, omega^2=-1-omega."""
    def __init__(self,L): self.L,self.modulus = L,1 << L
    def add(self,x,y): return tuple((a+b)%self.modulus for a,b in zip(x,y))
    def mul(self,x,y):
        a,b = x; c,d = y; m = self.modulus
        return ((a*c-b*d)%m,(a*d+b*c-b*d)%m)
    def power(self,x,e):
        value = (1,0)
        for _ in range(e): value = self.mul(value,x)
        return value
    def sigma(self,x):
        a,b=x
        return ((a-b)%self.modulus,(-b)%self.modulus)
    def norm(self,x):
        value=self.mul(x,self.sigma(x))
        assert value[1]==0
        return value[0]
    def scale(self,x,c): return tuple((v*c)%self.modulus for v in x)
    def indicator(self,x): return (1-evaluate(smoothstep(self.L),self.norm(x),self.modulus))%self.modulus
    def inverse(self,x):
        return self.scale(self.sigma(x),evaluate(smoothstep(self.L+1)[1:],self.norm(x),self.modulus))
    def teichmuller(self):
        return [(0,0),(1,0),(0,1),(self.modulus-1,self.modulus-1)]


def polynomial_value(ring,terms,values,conjugate=False):
    total=(0,0)
    for coefficient,exponents in terms:
        value=ring.sigma(coefficient) if conjugate else coefficient
        for x,e in zip(values,exponents): value=ring.mul(value,ring.power(x,e))
        total=ring.add(total,value)
    return total


def ring_controls():
    cases=[]; elements=0; point_values=0; inverse_pairs=0
    for L in range(1,7):
        ring=QuadraticRing(L);m=ring.modulus
        norm_histogram={};unit_count=zero_residue_count=0;inverses={}
        for a in range(m):
            for b in range(m):
                x=(a,b);norm=ring.norm(x)
                assert norm == (a*a-a*b+b*b)%m
                assert ring.sigma(ring.sigma(x))==x
                unit = bool(a%2 or b%2)
                assert bool(norm%2)==unit
                assert ring.indicator(x)==int(not unit)
                inverse=ring.inverse(x)
                inverses[x]=inverse
                if unit:
                    assert ring.mul(x,inverse)==(1,0)
                    unit_count+=1
                else:
                    assert inverse==(0,0)
                    zero_residue_count+=1
                norm_histogram[str(norm)]=norm_histogram.get(str(norm),0)+1
                elements+=1
        if L<=3:
            for x,ix in inverses.items():
                for y,iy in inverses.items():
                    assert inverses[ring.mul(x,y)]==ring.mul(ix,iy)
                    assert ring.sigma(ring.mul(x,y))==ring.mul(ring.sigma(x),ring.sigma(y))
                    inverse_pairs+=1
        teich=ring.teichmuller()
        assert len(set(teich))==4 and all(ring.power(x,4)==x for x in teich)
        assert all(ring.sigma(x)==ring.power(x,2) for x in teich)
        terms=[((1,0),(1,0)),((1,0),(0,1)),((0,1),(1,1)),
               ((2%m,0),(2,0)),((0,2%m),(0,1))]
        count=0;lifted_count=0;values=[];rational_values=[];rational_sum=(0,0)
        for x in teich:
            for y in teich:
                value=polynomial_value(ring,terms,[x,y])
                conjugate=polynomial_value(ring,terms,[ring.power(x,2),ring.power(y,2)],True)
                assert ring.sigma(value)==conjugate
                is_root=int(not(value[0]%2 or value[1]%2))
                indicator=ring.indicator(value)
                assert indicator==is_root
                count+=is_root;lifted_count=(lifted_count+indicator)%m
                point_values+=1
                values.append(dict(inputs=[x,y],value=value,norm=ring.norm(value),indicator=indicator))
                first=ring.add(ring.add(x,y),(1,0))
                second=ring.add(x,ring.mul((0,1),y))
                numerator=ring.add(x,ring.scale(y,2))
                direct=ring.mul(numerator,ring.mul(ring.inverse(first),ring.inverse(second)))
                collected=ring.mul(numerator,ring.inverse(ring.mul(first,second)))
                assert direct==collected
                rational_sum=ring.add(rational_sum,collected)
                rational_values.append(dict(inputs=[x,y],numerator=numerator,
                                            denominators=[first,second],value=collected))
        assert lifted_count==count%m
        if L>4: assert lifted_count==count
        even_constant=(2%m,0)
        assert ring.indicator(even_constant)==1
        square_counterexample=ring.power((1,1),2)!=ring.sigma((1,1))
        assert square_counterexample == (L>1)
        packed_count=0;packed_values=[]
        for x in (0,1):
            for y in (0,1):
                f=x+y+2*x*y;g=x*y+2*x
                packed=(f%m,g%m)
                expected=int(not(f%2 or g%2))
                assert ring.indicator(packed)==expected
                packed_count+=expected
                packed_values.append(dict(inputs=[x,y],polynomials=[f,g],packed=packed,indicator=expected))
        assert packed_count==1
        cases.append(dict(precision=L,ring_elements=m*m,units=unit_count,
                          zero_residue_elements=zero_residue_count,norm_histogram=norm_histogram,
                          direct_residue_root_count=count,lifted_count=lifted_count,
                          exact_count_recovered=L>4,polynomial_terms=terms,teichmuller_values=values,
                          even_constant=list(even_constant),even_constant_nonzero=L>1,
                          even_constant_norm=ring.norm(even_constant),
                          power_is_not_frobenius=square_counterexample,
                          collected_inverse_values=rational_values,collected_inverse_sum=rational_sum,
                          packed_binary_equation_values=packed_values,packed_binary_root_count=packed_count))
    return dict(exhaustive_ring_elements=elements,teichmuller_polynomial_values=point_values,
                inverse_multiplicativity_pairs=inverse_pairs,cases=cases)


def report():
    invoices=[]
    for w,v,delta in ((1,2,2),(2,2,2),(8,4,3),(16,8,2)):
        L=v*w+1
        invoices.append(dict(w=w,v=v,delta=delta,exact_count_precision=L,
                             indicator_queries=L+1,unit_inverse_queries=L+1,
                             indicator_S=(2*L-1)*w*delta,
                             unit_inverse_S=((2*L+1)*w-1)*delta,
                             collected_denominator_S_example=5+((2*L+1)*w-1)*7,
                             packed_three_equation_S=(2*L-1)*3*w*delta))
    result=dict(scope='Linear-precision residue-root indicator and total unit inverse; no precision-free contraction or stronger factoring bound.',
                integer=integer_controls(),quadratic_ring=ring_controls(),invoices=invoices,
                source_sha256={'tools/check_precision_residue_indicators.py':sha256(Path(__file__).read_bytes()).hexdigest()})
    return json.loads(json.dumps(result))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write: RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else: assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],integer_identities=result['integer']['exact_polynomial_identities'],
                          scalar_residues=result['integer']['scalar_residues'],
                          exhaustive_ring_elements=result['quadratic_ring']['exhaustive_ring_elements'],
                          inverse_multiplicativity_pairs=result['quadratic_ring']['inverse_multiplicativity_pairs'],
                          teichmuller_values=result['quadratic_ring']['teichmuller_polynomial_values'])))


if __name__=='__main__': main()
