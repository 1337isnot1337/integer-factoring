"""Exact sparse F2[a,b,z] arithmetic for determinantal certificates."""
from dataclasses import dataclass
MASK=(1<<20)-1
@dataclass(frozen=True)
class Poly:
    terms:frozenset=frozenset()
    def __bool__(self):return bool(self.terms)
    def __len__(self):return len(self.terms)
    def __add__(self,b):return Poly(self.terms^b.terms)
    def __mul__(self,b):
        # Guard packed exponent addition: no carry between variables.
        if self and b:
            for shift in (0,20,40):
                if max((m>>shift)&MASK for m in self.terms)+max((m>>shift)&MASK for m in b.terms)>MASK:
                    raise OverflowError('Packed polynomial exponent limit exceeded')
        if self==b:return Poly(frozenset(2*x for x in self.terms))
        out=set()
        for x in self.terms:
            for y in b.terms:
                z=x+y
                if z in out:out.remove(z)
                else:out.add(z)
        return Poly(frozenset(out))
    def __pow__(self,n):
        r=ONE;a=self
        while n:
            if n&1:r=r*a
            n//=2
            if n:a=a*a
        return r
    def exquo(self,b):
        q,r=self.div(b)
        if r:raise ArithmeticError('Nonexact division')
        return q
    def div(self,b):
        if not b:raise ZeroDivisionError()
        lead=max(b.terms);a=set(self.terms);q=set();rem=set()
        while a:
            m=max(a)
            if all(((m>>s)&MASK)>=((lead>>s)&MASK) for s in (0,20,40)):
                v=m-lead;q.add(v);a.symmetric_difference_update(v+t for t in b.terms)
            else:rem.add(m);a.remove(m)
        return Poly(frozenset(q)),Poly(frozenset(rem))
    def compose(self,x,y):
        pos={A:40,B:20,Z:0}[x];out=ZERO
        for m in self.terms:
            exp=(m>>pos)&MASK;out=out+Poly(frozenset({m-(exp<<pos)}))*y**exp
        return out
    def __str__(self):return str([[(m>>s)&MASK for s in (40,20,0)] for m in sorted(self.terms)])
ZERO=Poly();ONE=Poly(frozenset({0}));A=Poly(frozenset({1<<40}));B=Poly(frozenset({1<<20}));Z=Poly(frozenset({1}))
