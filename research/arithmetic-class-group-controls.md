# Actual quadratic-form controls for the arithmetic framework

The abstract torsion routine is now exercised on 40 quadratic-form class
groups, including nonfundamental discriminants. The exact checker validates
28,421 pair products and 2,960,531 associative triples. Diagonal relations
on supplied class generators produce an unsaturated relation lattice; its
determinant projects the generators to their two-parts. The resulting
two-torsion basis yields proper factors for all 18 composite controls,
including 10403=101*103, 10807=101*107 and 16637=127*131.

These are finite arithmetic controls. The code enumerates all reduced forms
and finds each input generator's order by cycling through its powers. These
searches have exponential input-bit cost. They are not implementations of
the fast generator/relation primitives required for a stronger factoring
bound. The strongest complete score remains 1/2.

## Exact composition used by the checker

The implementation is [binary_quadratic_forms.py](../tools/binary_quadratic_forms.py).
It uses primitive positive definite forms of discriminant Delta<0 with the
canonical convention -a<b<=a<=c, and b>=0 if a=c. This agrees with
Lenstra--Pomerance Section 2. Integer shears and the proper swap
(a,b,c)->(c,-b,a) give a direct reduction routine.

For composition, first replace (a,b,c) by an equivalent form with leading
coefficient coprime to the second form's coefficient A. The search uses
Q(1,k)=a+bk+ck^2 for 0<=k<A. For odd p dividing A, this polynomial cannot
vanish at every residue, since it has degree at most two and primitiveness
excludes all coefficients being divisible by p. The only obstruction at
p=2 is a even, b and c odd; the proper swap fixes it. CRT therefore
guarantees that the explicit finite search terminates.

The shear changes the middle coefficient to b+2ck. Once gcd(a,A)=1,
choose t satisfying

```math
t\equiv ((B-b)/2)a^{-1}\pmod A,
```

where B is the second form's middle coefficient. The composed form before
reduction has coefficients

```math
(aA,\ b+2at,\ ((b+2at)^2-\Delta)/(4aA)).
```

This is coprime-leading Gauss composition. The congruences make the last
coefficient integral, and reduction gives its canonical class identifier.
The coprime-leading search is deliberately simple; it is not the efficient
composition algorithm quoted in the asymptotic theorem.

## What the checker independently compares

For every listed class group it checks closure, identity/inverse behavior,
commutativity, all associative triples, and canonical invariance under
several integer shears and proper swaps. It also checks that forms whose
first coefficient is 7-smooth and coprime to Delta equal the products of
their signed prime-form factors from Lenstra--Pomerance Section 2.8. There
are 110 such smooth-form controls across the listed instances.

Every nonidentity form is supplied as a generator, so generation is explicit
and exhaustive. The diagonal matrix of their individual orders is a
full-rank relation matrix; its determinant generally exceeds the class
number. The generic algorithm uses that determinant without factoring it.
The independently enumerated set of self-inverse classes equals the span
of the extracted basis tops. This verifies complete two-torsion recovery,
rather than merely finding one successful factor.

For Delta=-4N, an ambiguous form has b=0, b=a or a=c. Respectively its
discriminant factorization uses (a,c), (a,4c-a), or (2a-b,2a+b), after
discarding powers of two. Gcds with N recover a proper factor from at least
one basis top for each tested N. The target list includes both squarefree
integers and integers with repeated prime factors.

Run `python3 tools/check_arithmetic_class_groups.py` or
`make verify-arithmetic-framework`. The
[evidence record](../evidence/arithmetic-class-groups.json) lists each
discriminant, class number, two-rank, relation-determinant bit length and
factor. It deliberately reports the exhaustive search scope. This work
uses substantial OpenAI Codex assistance and adds no Lean formalization.
