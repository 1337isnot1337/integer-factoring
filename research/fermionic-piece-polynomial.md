# Reflected monomial duals and a four-piece crossing certificate

## Result and scope

There is an explicit characteristic-two dual family consisting of one
monomial per label. It uses bit-complemented mode indices. The dual
monomial may have greater weight than the label; the source locality
argument extends to that case with a correspondingly larger finite
piece range.

On the three labels

```math
(1,11),\quad(3,9),\quad(5,7),
```

the reflected duals are `(5,15),(7,13),(9,11)`. They give exactly four
nonzero crossing pieces over `F_2`:

| Shift | Row rank | Column rank | Rank product |
| ---: | ---: | ---: | ---: |
| -32 | 1 | 9 | 9 |
| -16 | 2 | 4 | 8 |
| 0 | 4 | 5 | 20 |
| 16 | 6 | 2 | 12 |

The pieces sum to the crossing with multiplier one. The maximum rank
product is **20**, improving the finite value 27 obtained on the same
three-label domain by the role-permutation certificate in
[PR22's modular construction](https://github.com/1337isnot1337/integer-factoring/blob/c12da897fcc6081685e8c663217ce9083db87b4a/research/past-half-pieces.md).

This is a finite tensor improvement. Three labels do not encode an
arbitrary number of bits, and no extrapolation to a uniform coefficient
below `1/sqrt(2)` is made. Tensoring the four terms individually would
create exponentially many pieces. Their aggregated coefficient tensors
require their own rank analysis.

The underlying operators and crossing argument are Anthropic's Appendix B,
pinned to
[`16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The concurrent
[PR22 Pfaffian note](https://github.com/1337isnot1337/integer-factoring/blob/c12da897fcc6081685e8c663217ce9083db87b4a/research/modular-pfaffian.md)
already supplies a folded, polynomial-time Pfaffian evaluator for the
original inhomogeneous duals. The additions here are the arbitrary-weight
monomial formulation, reflected dual family, and four-piece certificate.
The accompanying matching implementation is an independent evaluator
for these generalized monomial inputs. No new Lean formalization is claimed.

## 1. Generalize the locality parameter

Let `a` be homogeneous of weight `A`, let `b` have weight `h`, and let
the input and output states both have weight `h`. There is no need to
assume `A<=h`. Define

```math
Q_s(a,b)=\sum_{i=0}^{A+h}(-1)^i\binom{A+h}{i}
                 a[s-h+i]b[h-s-i].
\tag{1}
```

Then the source's locality proof gives

```math
Q_s(a,b)=\sum_{i=0}^{A+h}(-1)^i\binom{A+h}{i}
                 b[h-s-i]a[s-h+i],
\tag{2}
```

and, on the weight-`h` state space,

```math
Q_s(a,b)=0\quad\text{outside }-h-A\le s\le2h,
\qquad
\sum_sQ_s(a,b)=\bigl(a[-h]b\bigr)\,\mathrm{Id}.
\tag{3}
```

Here is the change to the source proof. Its common Laurent polynomial is
obtained by multiplying the two formal operator orders by
`(omega_1-omega_2)^(A+h)`. The matching expansion has poles of order at
most `A+h` for arbitrary nonnegative `A,h`. Thus the two expansions still
give the same Laurent polynomial. Evaluating that polynomial at `(1,1)`
leaves the complete contraction, namely `a[-h]b`.

For the upper support bound, use the original order: its inner operator
would send a weight-`h` input to weight `2h-s-i`, which is negative when
`s>2h`. For the lower support bound, use the exchanged order: its inner
operator sends that input to weight `s+i`, which is negative for every
`i<=A+h` when `s<-A-h`. The assumption `A<=h` in the original statement
was used to put this asymmetric interval inside `[-2h,2h]`.

The proof is an identity of integral matrix coefficients before reducing
modulo two. The stable state quotient in
[`modular-crossing.md`](modular-crossing.md) therefore preserves it.
For a polynomial `a` of several weights, apply (1)--(3) to its homogeneous
parts and add. The largest occurring weight gives a common finite range.

## 2. Bit-complemented monomials are exact duals

Fix a power of two `M`. Choose any collection of distinct-odd partitions
`mu` of one weight `h`, with every part less than `2M`. Write their parts
as `2a+1`, where `0<=a<M`. Define the reflected partition

```math
\mu^*=\{2M-m:m\in\mu\},\qquad a^\mu=t^{\mu^*},
\qquad b^\mu=t^\mu.
\tag{4}
```

The particle counts of the labels may differ. If a label has `ell` parts,
its reflected monomial has weight `A=2M ell-h`; its locality parameter
is `A+h=2M ell`.

**Reflected duality.** For all chosen labels `mu,nu`,

```math
[1](t^{\mu^*})[-h]t^\nu=[\mu=\nu]\quad\text{in }\mathbf F_2.
\tag{5}
```

Only annihilations contribute. Different particle counts give zero.
For equal counts, pair an operator part `2(M-1-a)+1` with an input part
`2b+1`. Its coefficient is one precisely when

```math
(M-1-a)\mathbin{\mathrm{AND}}b=0,
\quad\text{equivalently }b\subseteq a\text{ bitwise}.
```

In particular `b<=a`. The input and original labels have equal total
weight and equal particle count, so their sums of indices are equal.
A matching with all `b<=a` must therefore have equality in every pair.
It exists exactly when `nu=mu`, and is then unique because all parts
are distinct. This proves (5).

Combining (3) and (5) gives a scale-one crossing family on every such
label set. The two intermediate-state factorizations remain exact.
They use weights `s+i` and `2h-s-i`; the allowed range of `i` may now
depend on the first label, and can be padded to the common maximum.

For a common particle count `ell`, all dual weights coincide. If `ell`
is also a power of two, `A+h=2M ell` is a power of two. In characteristic
two only the endpoints `i=0,A+h` in (1) survive. This explains the small
number of terms in the certificate; it does not itself bound matrix ranks.

## 3. An unfurled polynomial matching formula

This section gives an exact evaluator for an arbitrary regular monomial
`a`, rather than folding the original `aFam` expansion. Let the parts of
`a,b,u,y` be distinct and odd, with `|b|=|u|=|y|=h` and `|a|=A`.
Put `t=1+z` for the formulas below, and define

```math
C(m,k)=\binom{k-1}{m-1}\bmod2,
\qquad D(m,k)=\binom{k+m-1}{m-1}\bmod2.
```

For odd `m,k`, these are the binary tests

```math
C(m,k)=[(m-1)\subseteq(k-1)],\qquad
D(m,k)=[(m-1)\mathbin{\mathrm{AND}}k=0].
```

Form an alternating polynomial matrix `N(z)` whose four vertex groups
are the output parts `U`, the operator parts `B`, the operator parts `A`,
and the input parts `Y`. All within-group entries are zero. Its six
cross-group entries are:

| Vertices | Entry |
| --- | --- |
| `U_u,B_b` | `C(b,u) t^b` |
| `U_u,A_a` | `C(a,u) z^u t^a` |
| `U_u,Y_y` | `[u=y] z^y` |
| `B_b,A_a` | `D(b,a) z^a` |
| `B_b,Y_y` | `D(b,y) z^y t^b` |
| `A_a,Y_y` | `D(a,y) t^a` |

Then

```math
\boxed{\operatorname{Pf}N(z)
       =z^{A+h}\sum_s[t^u]Q_s(t^a,t^b)t^y\,z^s.}
\tag{6}
```

Every displayed exponent is nonnegative. Every matching has degree at
most `A+3h`, so the formula is a finite polynomial identity.

### Derivation

On the characteristic-two state quotient, creation and annihilation
satisfy the square-zero Clifford relations. Commuting derivatives to the
right gives a matching expansion. Contractions within either already
normally ordered operator block are forbidden. In characteristic two,
the usual matching signs all equal one.

The only infinite contraction is from a part `b` of the left operator
to a part `a` of the right operator. Its formal generating function is

```math
\sum_{k\ge1,\ k\text{ odd}}D(b,k)C(a,k)z^k
=D(b,a)\frac{z^a}{(1+z)^{a+b}}.
\tag{7}
```

To verify this directly, write `a=2alpha+1`, `b=2beta+1`,
`k=2kappa+1`. The summand is nonzero exactly when
`alpha subseteq kappa` and `beta AND kappa=0`. If `alpha AND beta`
is nonzero there are no terms. Otherwise summing the free binary digits
of `kappa` gives the rational expression in (7). It is interpreted as
a formal power series at zero.

Write the field of `a` with a variable `z` and retain the opposite shift
of the other operator by equality of input and output weights. The
definition (1) multiplies that series by
`z^(-A)(1+z)^(A+h)`.
Distribute the factor `(1+z)^a` to each `A` vertex, `(1+z)^b` to each
`B` vertex, and `z^y` to each `Y` vertex. Their total is
`z^h(1+z)^(A+h)`.
Equation (7) then becomes the polynomial edge `D(b,a)z^a`.
The other five contractions give the remaining table entries directly.
Dividing the complete matching polynomial by `z^(A+h)` proves (6).

The matching recursion in the implementation is independent of the
source's normal-order recursion and its intermediate-partition sums.
It also differs from PR22's folded original-dual matrix. Both formulations
can use a fraction-free polynomial Pfaffian algorithm; the implementation
here retains subset matching as a small independent verification path.

## 4. Effective search at large mode indices

Suppose every monomial has at most `L` particles and every part is at
most `B`. The matrix in (6) has at most `4L` vertices and degree
`O(LB)`. The subset matching algorithm computes its Pfaffian in
`2^(O(L)) poly(LB)` bit operations. Its coefficients are bits; it does
not hide rational coefficient heights or an ambient partition enumeration.

In particular, `L=O(lambda)` and `B<=2^(O(lambda))` suffice for a
`2^(O(lambda))` entry algorithm. A polynomial containing `2^(O(lambda))`
explicit monomials can be evaluated term by term within the same form
of bound. The reflected family needs only one monomial per label.

For `2^lambda` labels, all final piece matrices and exact Gaussian rank
factorizations can also be computed in `2^(O(lambda))` time under these
conditions. This includes computing the width metadata; it uses the
broader interface in
[`subexponential-piece-transfer.md`](subexponential-piece-transfer.md).
For example, if `B=2^(O(lambda^beta))` with fixed `beta<1`, the complete
shift interval has logarithmic cardinality `O(lambda^beta+log lambda)`.
That piece entropy is admissible in the broader transfer.

These cost bounds enlarge the space of explicit families one can search.
They do not give a uniform rank upper bound. A family at large energy
must still prove its actual two-flattening width inequality.

## 5. The exact three-label certificate

Use `M=8`, `h=12`, and the three labels in the opening table.
They have two particles each, so their reflected weight is `A=20`
and `A+h=32`. Thus only `i=0,32` occur in (1).

Rows are ordered `(x,y)` and columns `(v,u)`, with each label index in
`0,1,2`. Encode a row as an integer whose bit `3v+u` is that entry.
The four matrices are exactly:

```text
s=-32: [273, 0, 0, 0, 273, 0, 0, 0, 273]
s=-16: [256, 0, 0, 0, 256, 0, 0, 0, 17]
s=  0: [272, 2, 0, 8, 257, 0, 0, 0, 17]
s= 16: [256, 0, 4, 0, 256, 32, 64, 128, 17]
```

All other shifts in the support interval `[-32,24]` are zero.
XORing the four matrices gives `I_9`, which is the desired crossing in
the row flattening. Exact elimination gives the ranks stated at the start.

The checker independently evaluates both operator orders using the
literal source normal-order recursion, retaining the generalized angle
parameter 32. It checks every shift in `[-32,24]`, including all zeros.
Thus the finite certificate is not validated only by the new Pfaffian
formula.

Run:

```sh
python3 tools/check_fermionic_piece_polynomial.py \
  --output evidence/fermionic-piece-polynomial.json
```

The report records 550 original-family source comparisons, nine reflected
duality entries, 4,617 generalized entries checked in both operator orders,
and all 81 entries of the summed crossing. Dropping a specified term is
a negative control. The literal matrices and rank profile are included
in the evidence file.

The next target is a growing reflected or related family with a proved
small width sum. The finite four-piece improvement is a useful primitive,
but neither its small rank product nor the new evaluator completes that
asymptotic step.

The starting operator construction is prior work by Anthropic; the folded
original-family Pfaffian is concurrent PR22 work cited above. The
generalized reflected-dual analysis and certificate used substantial
OpenAI assistant support.
