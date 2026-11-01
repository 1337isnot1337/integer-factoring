"""Exact factorization in F8[t], using x^3+x+1 for the coefficient field."""
import random
from three_row_hecke_paths import M,F

def trim(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return a

def monic(a):
    a=trim(a)
    return [M[F.power(a[-1],6)][v] for v in a] if a else []

def add(a,b):
    out=list(a)+[0]*max(0,len(b)-len(a))
    for i,v in enumerate(b):out[i]^=v
    return trim(out)

def multiply(a,b):
    if not a or not b:return []
    out=[0]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        if c:
            row=M[c]
            for j,d in enumerate(b):
                if d:out[i+j]^=row[d]
    return trim(out)

def divrem(a,b):
    a=trim(a);b=trim(b)
    if not b:raise ZeroDivisionError()
    out=[0]*max(0,len(a)-len(b)+1);inv=F.power(b[-1],6)
    while len(a)>=len(b):
        k=len(a)-len(b);c=M[a[-1]][inv];out[k]=c;row=M[c]
        for i,v in enumerate(b):
            if v:a[k+i]^=row[v]
        while a and not a[-1]:a.pop()
    return trim(out),a

def quotient(a,b):
    q,r=divrem(a,b)
    if r:raise ArithmeticError('Nonexact polynomial division')
    return q

def gcd(a,b):
    while b:a,b=b,divrem(a,b)[1]
    return monic(a)

def derivative(a):return trim([a[i+1] if (i+1)&1 else 0 for i in range(len(a)-1)])

def square_mod(a,f):
    out=[0]*max(0,2*len(a)-1)
    for i,c in enumerate(a):out[2*i]=M[c][c]
    return divrem(out,f)[1]

def qpower_mod(a,f):
    for _ in range(3):a=square_mod(a,f)
    return a

def sqrt_poly(a):
    if any(a[1::2]):raise ValueError('Not a polynomial square')
    return trim([F.power(c,4) for c in a[::2]])

def distinct_degree(f):
    f=monic(f);x=[0,1];h=divrem(x,f)[1];d=1;out=[]
    while 2*d<=len(f)-1:
        h=qpower_mod(h,f);g=gcd(add(h,x),f)
        if len(g)>1:
            out.append((g,d));f=quotient(f,g)
            if f==[1]:break
            h=divrem(h,f)[1]
        d+=1
    if len(f)>1:out.append((f,len(f)-1))
    return out

def equal_degree(f,d,rng):
    if len(f)-1==d:return [f]
    assert (len(f)-1)%d==0
    for _ in range(256):
        a=[rng.randrange(8) for _ in range(len(f)-1)];trace=[];h=a
        for _ in range(3*d):trace=add(trace,h);h=square_mod(h,f)
        g=gcd(trace,f)
        if 1<len(g)<len(f):return equal_degree(g,d,rng)+equal_degree(quotient(f,g),d,rng)
    raise ArithmeticError('Equal-degree split attempts exhausted')

def factor(f,seed=803771):
    rng=random.Random(seed);out=[]
    def visit(f,mult):
        f=monic(f)
        if len(f)<=1:return
        der=derivative(f)
        if not der:return visit(sqrt_poly(f),2*mult)
        c=gcd(f,der);w=quotient(f,c);i=1
        while w!=[1]:
            y=gcd(w,c);z=quotient(w,y)
            if z!=[1]:
                for block,d in distinct_degree(z):
                    for g in equal_degree(block,d,rng):out.append((g,mult*i))
            w=y;c=quotient(c,y);i+=1
        if c!=[1]:visit(sqrt_poly(c),2*mult)
    visit(f,1);return sorted(out,key=lambda r:(len(r[0]),r[0],r[1]))

def irreducible(f):
    f=monic(f);d=len(f)-1
    if d<1:return False
    primes=[];v=d;p=2
    while p*p<=v:
        if v%p==0:
            primes.append(p)
            while v%p==0:v//=p
        p+=1
    if v>1:primes.append(v)
    tests={d//p for p in primes};x=divrem([0,1],f)[1];h=x
    for i in range(1,d+1):
        h=qpower_mod(h,f)
        if i in tests and gcd(add(h,x),f)!=[1]:return False
    return h==x
