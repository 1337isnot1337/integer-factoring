# Polynomial factor extraction from class-group generators and an annihilator

Let N>1 be odd. Suppose a list of m forms generates the full proper class
group of discriminant -4N, and a supplied nonzero integer E annihilates
those generators. If N has at least two distinct prime factors, a
deterministic algorithm finds a proper divisor of N using

```math
O\!\left(m\,\operatorname{bits}(E)+m b^2\right)
```

class-group operations, where b=bits(4N). The forms have polynomial bit
cost per group operation. The remaining integer arithmetic is polynomial
in the input bit lengths. No factorization of E, exact class number,
squarefreeness promise, or random class is needed.

This improves the postprocessing bound in
[two-part extraction](arithmetic-torsion-extraction.md) when the goal is
to find a factor. The earlier routine computes a complete Sylow-two basis
and can pay an exponential factor in its rank. Here the algorithm stops
at the first useful ambiguous form. Before that event, it retains at most
one cyclic basis element. Construction of the complete generators and the
annihilator remains a separate obligation; this note does not supply
either input or improve the unconditional factoring score.

## 1. Exactly one possible unhelpful nonidentity class

Use the canonical reduced representative (A,B,C), with

```math
-A<B\le A\le C,
\qquad B\ge0\text{ if }A=C.
```

The inverse class is represented by (A,-B,C). Uniqueness of the reduced
representative shows that the class is its own inverse precisely at one
of the three boundaries

```math
B=0,\qquad B=A,\qquad A=C.                               \tag{1}
```

These are the usual ambiguous-form cases. The classical association with
coprime factorizations is discussed by Lenstra and Pomerance,
[*A rigorous time bound for factoring integers*, Section 2, especially
Theorem 2.5, pages 487--488](https://community.ams.org/journals/jams/1992-05-03/S0894-0347-1992-1137100-0/S0894-0347-1992-1137100-0.pdf).
For the particular discriminant -4N, all formulas and exceptional cases
needed here can be checked directly.

If B=0, the discriminant identity gives N=AC. Primitivity gives
gcd(A,C)=1. This factorization is trivial only when A=1, because A<=C.
The corresponding class is the identity

```math
I=(1,0,N).
```

If B=A, then A(4C-A)=4N. Since N is odd, A is even; writing A=2a gives

```math
N=a(2C-a).                                                \tag{2}
```

Both factors are odd. Primitivity gives gcd(2a,C)=1, which implies that
the two factors in (2) are coprime. Reduction gives C>=2a, so 2C-a>=3a.
The factorization is trivial exactly when a=1. Its reduced form is

```math
U=(2,2,(N+1)/2),                                         \tag{3}
```

which is primitive exactly when N=1 mod 4. For N>1 in that residue class,
it is reduced and differs from I.

Finally, if A=C, then B is even and

```math
N=(A-B/2)(A+B/2).                                        \tag{4}
```

Primitivity gives gcd(A,B)=1. It forces A odd, and because N is odd the
two factors in (4) are odd. Their gcd divides gcd(2A,B), hence is one.
Both factors are positive. If the smaller factor were one, put
v=B/2>=0. The equations A-v=1 and 2v=B<=A imply v<=1. Primitivity and
oddness force v even, so v=0 and N=1, contrary to the hypothesis. Thus
this case always gives a proper factor for N>1.

Consequently the only ambiguous classes whose displayed factorization
is trivial are I and, when N=1 mod 4, U. There is at most **one
nonidentity unhelpful class**. This conclusion holds for nonfundamental
discriminants and nonsquarefree N. No prime factor of N was used to test
any of the formulas.

## 2. Projecting the supplied generators to the two-part

Write E=2^s u with u odd. Let g_1,...,g_m be the supplied generators and
put h_i=g_i^u. Since E annihilates the whole group, u kills its odd-order
part and acts invertibly on its Sylow-two subgroup P. Therefore the h_i
generate P.

Computing u only strips powers of two. Raising the generators to u takes
O(m bits(E)) group operations by binary exponentiation. The implementation
also checks g_i^E=I for each generator it processes, within the same
order of cost. Completeness of the generator list remains a promise.

The actual order of an element of P is at most the class number. An
elementary bound sufficient here is |G|<=|Delta|=4N: for a reduced form,
A<=sqrt(|Delta|/3), and summing the at most 2A+1 possible middle
coefficients over A gives at most |Delta| possibilities. Thus at most
b=bits(4N) squarings find the order of any element known to lie in P.
The bound uses the actual class size, not the possibly much larger
two-adic valuation s of an oversaturated annihilator.

## 3. Early-stop insertion with one cyclic generator

Process the h_i in order. Maintain either an empty basis or a single
element z of order 2^ell. Its order-two top, z^(2^(ell-1)), is U. The
processed inputs generate the same subgroup as z; this is the invariant
after each input finishes. A candidate h of order 2^k is handled as
follows.

Compute its nonidentity order-two top

```math
t=h^{2^{k-1}}.
```

Apply the formulas of Section 1 to t. If they give a proper factor,
return that factor immediately. Otherwise t=U.

If the basis is empty, retain h as z and finish this input. If the basis
contains z and k<=ell, replace the candidate by

```math
h\,z^{-2^{\ell-k}}.                                     \tag{5}
```

The two factors have the same top U, so the new candidate has order at
most 2^(k-1). The subgroup generated by the basis and candidate is
unchanged. If k>ell, first swap the basis and the candidate, together
with their order exponents, and then perform (5). This also preserves
the generated subgroup and strictly lowers the order of the displaced
candidate. Repeat until the candidate is the identity or a factor is
found.

This is the ordinary prime-power basis insertion used in the existing
[generic two-part routine](arithmetic-torsion-extraction.md). For the
present goal, a second independent top would already factor N.
Every nonidentity top seen before success is the single explicitly known
class U, so there is no subset-product or discrete-log search among a
growing list of tops.

The candidate order exponent decreases at each cancellation. Each input
therefore uses at most b cancellations, plus at most one initial basis
insertion. Computing orders, tops, and correction powers takes O(b)
group operations per cancellation. The total after projection is
O(m b^2). It retains O(1) group elements, apart from the input and optional
diagnostic trace.

## 4. Why complete generators force a factor

Suppose the routine finishes without finding a factor. Its subgroup
invariant says that P is either trivial or cyclic, and that its only
possible nonidentity element of order two is U.

If N has at least two distinct prime factors, there exists a coprime
factorization N=ac with 1<a<=c. The form

```math
(a,0,c)
```

is primitive and already reduced, with discriminant -4N. It is
ambiguous, and its displayed factorization is proper. Hence it differs
from both I and U. It lies in P, contradicting the preceding conclusion.

This argument uses the existence of a nontrivial coprime factorization
only to prove termination with a factor. The algorithm is not supplied
that factorization. It also does not require the ambiguous-form
factorization map to be a homomorphism or require a genus-theory rank
formula.

If N is a prime power, a nontrivial coprime factorization need not exist,
and the routine can finish without a factor. Standard exact perfect-power
handling is separate. With incomplete supplied generators it can also
finish without a factor for a composite N. Its negative return value is
therefore labeled `no_factor_in_supplied_subgroup`, rather than a
primality or prime-power certificate.

## 5. Interaction with the arithmetic program

A full-rank relation matrix on complete generators gives an annihilator
by taking the absolute determinant, as proved in
[two-part extraction](arithmetic-torsion-extraction.md). This routine
then gives polynomial postprocessing in the encoded relation data and b.
The earlier roughness trial division and generic rank-exponential tail
are unnecessary for this specific discriminant and early-stop goal.

There is no change to the unresolved arithmetic search requirement.
The [universal smooth-square gate](universal-smooth-square.md) constructs
the relation matrix only if its bounded representations can be found.
The [decision reduction](smooth-square-decision-reduction.md) replaces
that search requirement with an exact variable-cap decision oracle.
The [small-generator theorem](zero-free-class-generators.md) is conditional
on its explicit analytic inputs. None of these inputs is created by
the present postprocessing result.

## 6. Implementation and checked evidence

The implementation is
[early_stop_class_factoring.py](../tools/early_stop_class_factoring.py).
The checker is
[check_early_stop_class_factoring.py](../tools/check_early_stop_class_factoring.py).
Run

```console
python3 tools/check_early_stop_class_factoring.py > evidence/early-stop-class-factoring.json
```

The exact [evidence](../evidence/early-stop-class-factoring.json) covers
all 511 odd integers from 3 through 1023. It enumerates each complete
small class group solely to supply independently known generators and
an annihilator, then runs two orderings of the complete generator list.
Every annihilator is intentionally oversaturated. The controls include
97 nonsquarefree inputs and 323 inputs with at least two distinct prime
factors. All 1,022 complete-generator runs return a proper factor exactly
when such a coprime factorization exists.

The checker verifies 1,423 ambiguous forms. Across those runs it exercises
2,896 dependent-top cancellations and 64 order swaps. The maximum
retained basis rank is one. A proper-subgroup control at N=21 returns
no factor, and an invalid annihilator is rejected.

The additional mixed-top control has N=161, discriminant -644, and
generators

```math
(2,2,81),\quad(3,2,54),\quad(5,4,33).
```

They generate all 16 classes. With E=240 their initial projected tops
are all the unhelpful class U=(2,2,81). The algorithm performs one order
swap and three dependent-top cancellations, then exposes (7,0,23) and
returns 7. The evidence records the full trace. Simply taking the final
order-two power of each original generator would miss this factor.

The finite computations check the implementation, the exceptional-form
classification, and the cancellation cases. The proof above supplies the
uniform statement. Neither part constructs the missing generators or
annihilator for a large factoring input.
