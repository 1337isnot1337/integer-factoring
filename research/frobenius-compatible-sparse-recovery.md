# Sparse recovery compatible with weighted Frobenius digit boxes

## Result and purpose

The weighted factoring evaluator need not work over every commutative
characteristic-two algebra. A smaller, explicitly constructible family
of tensor-product algebras is sufficient for deterministic recovery.
In that family, both operand digit boxes have **invertible** Fourier
weights, and full-width boxes are exact Frobenius coboundaries.

This removes a concrete obstacle to using the arithmetic identities in
[carry-fourier-elimination.md](carry-fourier-elimination.md): unweighted
boxes vanish at zero frequency, while unrestricted polynomial weights
are not fixed by the digit-width Frobenius. The construction below
avoids both problems without changing the factor search being performed.

For any fixed `0<alpha<beta<1`, and a top group count growing at least
as a fixed positive power of `n` under `P^40<=n`, a uniform evaluator
of the resulting weighted queries with cost

```math
\exp\!\left(O\!\left(n/(\log n)^\alpha\right)\right)
```

would give complete deterministic factoring with the same bound.
This is a transfer theorem and an admissible-algebra construction; it
does not supply the improved evaluator. The established score remains
`1/2`.

The sparse recovery itself uses odd-prime extension degrees at most
`n^(1/2+o(1))` after an explicitly charged roughness pretest. Rounding
the adaptive network's digit width to those degrees adds only
`n^(19/20+o(1))` capacity bits under `P^40<=n`. The construction does
not build a large compositum field or assume a free primitive root.

## 1. A punctured-field hitting theorem for sparse polynomials

Let `A` be any commutative unital `F_2`-algebra. It may have zero
divisors. Let

```math
F(X,Y)=\sum_{i=1}^t a_iX^{u_i}Y^{v_i},
\qquad 0<t\le K,\qquad 0\le u_i,v_i<2^n,
\tag{1}
```

where the exponent pairs are distinct and every coefficient `a_i` is
nonzero. The same argument works with any fixed number of variables.

For an odd prime degree `ell`, write

```math
E_\ell=\mathbb F_{2^\ell},\qquad r_\ell=2^\ell-1,
\qquad E_\ell^\circ=E_\ell\setminus\{0,1\}.
\tag{2}
```

Evaluation takes place in `A tensor_(F_2) E_ell`. The natural map
`A -> A tensor E_ell` is injective: extend `1` to an `F_2`-basis of
`E_ell` and view the tensor product as a direct sum of copies of `A`.

### Punctured interpolation

Suppose `r_ell>K` and (1) vanishes at every pair in
`(E_ell^circ)^2`. Reduce its exponents modulo `r_ell` in each variable,
obtaining `G` of degree below `r_ell` in each variable and with at most
`K` nonzero monomials.

Fix a permitted value of `Y`. A polynomial of degree below `r_ell`
that vanishes at the `r_ell-1` values in `E_ell^circ` must be a
multiple of

```math
S_{r_\ell}(X)=1+X+\cdots+X^{r_\ell-1}.             \tag{3}
```

This remains true over the tensor-product algebra. Differences between
distinct elements of `E_ell` are units after embedding, so ordinary
Vandermonde interpolation is valid. Equivalently, divide by the monic
polynomial (3); its remainder of degree below `r_ell-1` vanishes at
that many distinct points and must be zero.

The polynomial `G(X,Y)` has at most `K<r_ell` distinct `X` exponents.
At least one coefficient is identically zero, which forces the scalar
multiple in (3) to be zero. Every coefficient as a polynomial in `Y`
therefore vanishes on `E_ell^circ`. Apply the same argument to each of
those coefficients, using its at most `K` terms. The result is

```math
F=0\pmod{X^{r_\ell}-1,\ Y^{r_\ell}-1}.          \tag{4}
```

The strict condition `K<r_ell` is essential. Polynomial (3) itself
vanishes on the punctured field but is nonzero at one.

### The anchored difference product

Choose one support term in (1), with exponent pair `(u_0,v_0)`.
For each other term, choose one coordinate where its exponent differs
from that anchor and let `Delta_i` be the positive absolute difference
in that coordinate. Put

```math
\Pi=\prod_{i\ne0}\Delta_i,
\qquad1\le\Pi\le2^{n(t-1)}\le2^{n(K-1)}.
\tag{5}
```

The empty product is one. If (4) holds, another support term must have
the same two residues as the anchor. Otherwise its residue coefficient
would remain the nonzero coefficient `a_0`. For that colliding term,
`r_ell` divides every coordinate difference, and hence divides the
selected `Delta_i`. Thus every degree whose entire punctured grid is
bad satisfies

```math
r_\ell\mid\Pi.                                  \tag{6}
```

This argument does not assume that coefficients are binary, that
cancellation is absent, or that `A` is a field. It uses only the
nonzero anchor and the injective scalar extension.

For distinct prime degrees `ell,ell'`,

```math
\gcd(2^\ell-1,2^{\ell'}-1)
=2^{\gcd(\ell,\ell')}-1=1.                       \tag{7}
```

Consequently, any collection of allowed odd prime degrees with

```math
\prod_\ell(2^\ell-1)>2^{n(K-1)}                  \tag{8}
```

is a nonzeroness family: every nonzero `F` in (1) is nonzero at some
pair in one of its punctured field grids.

### The coefficients may vary with the degree

There is a useful stronger form. Fix the distinct exponent pairs in
(1), but for each selected degree allow a different commutative
`E_ell`-algebra `B_ell` and different nonzero coefficients
`a_(i,ell) in B_ell`. Keep those coefficients fixed while evaluating
that degree's whole punctured grid. Then at least one of these
degree-dependent polynomials is nonzero at a permitted grid point.

The proof is unchanged. Distinct field elements have unit differences
in every `B_ell`, so punctured interpolation still gives (4). A nonzero
anchor coefficient still requires an exponent collision. The common
integer `Pi` in (5) depends only on the exponent support, not on any
coefficient or coefficient algebra. Thus an entirely bad family would
make the product in (8) divide the same smaller positive integer `Pi`.

For example, a nonzero weight attached to each genuine divisor may
depend on the field degree and the radix used by that query. Such
weights preserve recovery even when the corresponding weighted
polynomials are not scalar extensions of one fixed polynomial.

## 2. An explicit degree bound and deterministic implementation

Set

```math
H=n(K-1)+1,\quad a=\lceil\log_2(H+1)\rceil,
\quad d_0=\lceil\log_2(K+1)\rceil,
\quad L=64(a+1)2^{\lceil a/2\rceil}.              \tag{9}
```

Use all odd prime degrees `d_0<ell<=L`, or stop as soon as the
integer product in (8) exceeds its target. Every chosen degree has
`r_ell>K`.

Here is an elementary proof that the stopping point occurs before `L`.
The central-binomial estimate used in
[sparse-divisor-nonzeroness.md](sparse-divisor-nonzeroness.md), Section 3,
gives, for even `v`,

```math
\pi(v)\ge\frac{v-\log_2(v+1)}{\log_2v}.
```

For `v>=32`, this implies at least `v/(2log2 v)` odd primes at most
`v`. One can check the numerical implication without prime asymptotics:
if `2^j<=v<2^(j+1)` and `j>=5`, then
`log2 v+log2(v+1)<=2(j+1)<=v/2`.

Write `R=2^ceil(a/2)`. We have `H>=K`, `d_0<=a<=R`, and
`log2 L<=8a`. There are therefore at least `4R` odd primes at most
`L`; discarding those at most `d_0` leaves at least `3R` primes.
If they are ordered increasingly as `ell_i`, then `ell_i-1>=2i`.
Hence

```math
\sum_i(\ell_i-1)\ge(3R)^2>H,
\qquad
\prod_i(2^{\ell_i}-1)>2^H>2^{n(K-1)}.           \tag{10}
```

The bound is deliberately conservative. It gives

```math
L=O\!\left(\sqrt{nK}\,\log(nK+2)
ight).         \tag{11}
```

Enumerate prime degrees by trial division. For each degree, construct
an irreducible binary polynomial by exhaustive search with the exact
irreducibility test described in
[carry-fourier-elimination.md](carry-fourier-elimination.md), Section 6.
This costs `2^O(ell)`. Enumerate all nonzero field elements except one,
and all ordered pairs of them. The entire family, including field
construction and the grid sizes, costs `2^O(L)` times a polynomial
factor. No unproved supply of special primes or irreducible polynomials
is assumed.

## 3. Roughness makes the degree bound sublinear

Fix a rational `0<beta<1`. Before using the weighted evaluator, test
divisors in increasing order up to `2^u`, stopping if trial division
already determines the smallest prime factor. Choose

```math
u\asymp n/(\log n)^\beta.                       \tag{12}
```

All values can be computed with integer comparisons. For example, if
`beta=a/b`, the least integer `u` with
`u^b floor(log2 n)^a>=n^b` is an exact choice. The optional larger
choice

```math
u=\max\left(
 \left\lceil n/\lfloor\log_2n\rfloor^\beta\right\rceil,
 \left\lceil100n/\lfloor\log_2n\rfloor\right\rceil
 \right)                                         \tag{13}
```

has the same asymptotic size and provides the convenient uniform bound
`K<=n^(1/100)` below. Finitely many small lengths may use trial division
through the square root.

If no divisor up to `2^u` exists, every prime factor exceeds `2^u`.
Counting multiplicities gives `Omega(N)u<n`, so

```math
\tau(N)\le2^{\Omega(N)}
\le K:=2^{\lfloor(n-1)/u\rfloor}
=n^{o(1)}.                                      \tag{14}
```

The pretest costs `2^O(u)` bit operations. Equations (11)--(14) give

```math
L=n^{1/2+o(1)}.                                 \tag{15}
```

With choice (13), `n/u<=floor(log2 n)/100`, which proves the stated
uniform `K` bound. The additional pretest term remains smaller than
`n/(log n)^alpha` for every fixed `alpha<beta<1`.

## 4. Round the digit width and keep the arithmetic budget

Use the adaptive network dimensions

```math
h=\text{least power of two at least }P^9,
\qquad D=2h^2,
\qquad P^{40}\le n.
```

For a chosen degree `ell`, replace its preliminary digit width `w_0`
by

```math
w_\ell=\ell\lceil w_0/\ell\rceil,
\qquad
m_\ell=D\left\lceil
       \frac{2w_\ell+\log_2D}{D}\right\rceil.
\tag{16}
```

Here `w_0=ceil(n/D)` for the first adaptive constructor. The
operand-length constructor may instead use its preliminary width from
`T-A-1`; rounding upward only decreases its required support lengths.
All exact convolution, carry, and factor-assignment proofs remain valid
with the new `w_ell,m_ell`.

The extra capacity is less than `D ell` per `D`-digit block and is
`O(D ell)` across the two operands, including the operand-length
constructor's at most `T=2D` positions. Since
`D<8P^18=O(n^(9/20))`,

```math
D\ell\le DL=n^{19/20+o(1)}=o(n).                \tag{17}
```

Thus all leading wire and setting-volume bounds of the two adaptive
constructors remain unchanged. Their fixed explicit numerical constants
are not asserted unchanged after this rounding. Their `O(n)` volume,
`O(n/P)` local table exponent, and `2^O(n/P)` written construction cost
continue to hold. The algorithm can check the integer condition
`DL<=n/16` before using this branch; it eventually holds uniformly under
the cutoff, and failure can use the finite trial-division fallback.

Moreover `m_ell=O(n/P^18+L)=O(n/P^18)` eventually. The field construction
below is still absorbed by the existing written-table budget.

## 5. Tensor-product algebras and invertible box weights

For each degree, construct separately

```math
F_\ell=\mathbb F_{2^{2m_\ell}},\qquad
E_\ell=\mathbb F_{2^\ell},\qquad
\mathcal A_\ell=F_\ell\otimes_{\mathbb F_2}E_\ell.
\tag{18}
```

Construct a primitive `M_ell=2^(m_ell)+1`-st root `zeta` in `F_ell`
by the deterministic `2^O(m_ell)` procedure already proved in the carry
Fourier note. The two fields cost `2^O(m_ell+ell)` to construct.
The tensor algebra has `F_2`-dimension `2m_ell ell`; using the two
defining polynomials gives polynomial bit cost for its arithmetic.
There is no search for a primitive element of a compositum of degree
`lcm(2m_ell,ell)`. The tensor algebra may have zero divisors, which the
following argument permits.

For odd `ell`,

```math
\gcd(2^{m_\ell}+1,2^\ell-1)=1.                  \tag{19}
```

Indeed, a common odd prime divisor would give an odd multiplicative
order for two, dividing `ell`, while `2^(m_ell)=-1` would give an element
of order two in that odd-order group.

Take any `alpha in E_ell^circ` and any frequency `t`. Write
`Y=alpha zeta^t` in the tensor algebra. Its two factors are units, and
`alpha^(M_ell)!=1` by (19). Therefore

```math
(1+Y)\,S_{M_\ell}(Y)=1+Y^{M_\ell}
 =1+\alpha^{M_\ell},
\qquad
(1+Y)^{-1}=
 \frac{S_{M_\ell}(Y)}{1+\alpha^{M_\ell}}.        \tag{20}
```

The denominator belongs to `E_ell` and is nonzero. Equation (20) is an
explicit inverse in the possibly non-field tensor algebra. It can be
computed with `O(m_ell)` algebra operations using division-free
geometric sums and an ordinary inverse in `E_ell`.

Let `b=2^(w_ell)`. Because `ell|w_ell`, the Frobenius
`sigma:z->z^b` fixes every element of `E_ell`. It acts only on the
`F_ell` coefficient factor in (18). The full digit-box weight is

```math
S_b(Y)=\frac{1+Y^b}{1+Y}
      =\frac{\sigma(1+Y)}{1+Y}.                 \tag{21}
```

Both numerator and denominator are units by (20). Hence this weight
is a unit at **every** frequency, including zero.

The same is true for every dyadic digit interval, not just a full box.
For an interval `[a,a+2^r-1]`, its weight is

```math
Y^a S_{2^r}(Y)=Y^a\frac{1+Y^{2^r}}{1+Y}.        \tag{22}
```

Frobenius sends nontrivial `alpha` to another nontrivial field element,
so the numerator is also a unit. This covers partial top digits,
leading-bit constraints, fixed digits, and the partial digit in a prefix
query. Formula (21), with a common Frobenius fixing the recovery field,
applies to the full-width blocks; (22) handles the remaining blocks
exactly.

Finally, `alpha^(b^i)=alpha` for every digit position `i`. The
multiplicative weight `alpha^p` of an entire first operand is therefore
the product of the same digit weight `alpha^(p_i)` at every position.
Choose a second nontrivial `beta in E_ell^circ` for the other operand.
Its box weights satisfy all the same unit statements. Neither operand
has the zero-frequency exception of an unweighted characteristic-two
box.

These identities make additional algebraic transformations legal. They
do not by themselves reduce a tensor rank: multiplying by units or
renaming field elements preserves the relevant linear information.
The evaluator still needs a useful joint contraction identity.

## 6. Complete deterministic recovery and the permitted evaluator

For a prefix query, use the bivariate polynomial

```math
\mathcal D_{N,k,x}(X,Y)
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}
      X^dY^{N/d}.                              \tag{23}
```

Its coefficients are one, its monomials are distinct, both exponents
are below `2^n`, and its support has at most `K` terms after the
roughness pretest. Add monomial weights to **both** input families of
the exact multiplication network. The proved unique-assignment
bijection gives (23) in characteristic two as well as over the integers.

Assume a uniform procedure computes the value of (23) in the algebras
(18) for every selected degree and every pair
`(alpha,beta) in (E_ell^circ)^2`, with the rounded network parameters.
Its correctness is required for these algebras and these inputs.
It is **not** required to compute a formal polynomial circuit valid
over every commutative `F_2`-algebra. This allows a prospective evaluator
to use the particular Frobenius, unit, and root identities proved above.

Section 1 makes these evaluations a deterministic nonzeroness test.
For smallest-prime-factor recovery, remove the known contribution
`X Y^N` when the queried prefix contains divisor one. Find the smallest
remaining divisor by successive prefix queries. There is always such
a divisor for `N>=2`, including `N` itself, and the smallest one is
prime. All subtractions, monomial evaluations, prefix comparisons,
integer divisions, and output costs are polynomial in the relevant
bit sizes. Repeating the procedure on the quotient yields complete
prime factorization in at most `n` factor-recovery calls.

If the specialized evaluation procedure has cost at most `A(n)` per
query, excluding the reusable field setup and written input construction,
use a nondecreasing upper
envelope for `A` over smaller input lengths. The complete bit cost is

```math
\exp\!\left(O\!\left(
 \log(A(n)+2)+n/P+u+L+\max_\ell m_\ell+\log n
 \right)\right).                                \tag{24}
```

This charges the `2^O(n/P)` written weighted tables, all field points,
both coordinates, every chosen degree, all field and primitive-root
construction, and every prefix query.
When the top group count grows at least as a fixed positive power of
`n`, as in the maximal-level schedule of the completed composition,
`max m_ell` is `n^(1-delta)` for some fixed `delta>0`, up to the smaller
rounding term. Equations (12), (15), and (24) therefore preserve every
target exponent `alpha<beta<1`.

The gain in this result is the admissible coefficient family and its
everywhere invertible digit weights. It does not assert that the
remaining coupled interval/Hankel sums can yet be evaluated at an
exponent better than one half.

## 7. Reproducible exact checks

Run:

```sh
python3 tools/check_frobenius_sparse_recovery.py \
  --output evidence/frobenius-sparse-recovery.json
```

The checker constructs the small fields, verifies sparse recovery with
nonbinary and degree-dependent coefficients and deliberate cancellation examples, checks the
strict puncture-size condition, and compares tensor-algebra box sums
with their explicit inverses and Frobenius formulas. It includes an
actual tensor algebra with zero divisors and reconstructs small complete
factorizations through the bivariate prefix interface. Exact integer parameter
checks cover the prime-degree product target, width rounding, and
capacity costs; no floating-point threshold is used.
