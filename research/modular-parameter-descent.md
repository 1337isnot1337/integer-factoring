# The modular quotient also descends in the operator parameter

## Result and scope

The quotient in [modular-crossing.md](modular-crossing.md) can be applied to
both arguments of the source operator.  Fix a prime `p` and write

\[
 V_p=\mathbf F_p[t_1,t_2,\ldots],\qquad
 I_p=(t_{pk}:k\geq1)+(t_k^p:k\geq1),\qquad W_p=V_p/I_p.
\]

**Parameter descent.** For every `a` in `I_p`, every mode `r`, and every
state `v` in `W_p`, the source operator satisfies `a[r]v=0`.  Consequently
`a[r]v` depends only on the classes of both `a` and `v` in `W_p`.

This supplies a missing justification for deleting nonregular monomials
from an operator parameter.  The original note correctly observes that
preservation of the state ideal alone would not justify that deletion;
the separate argument below does.  It permits a smaller representation
of the inhomogeneous dual polynomials without changing any crossing
piece, rank, channel width, or factoring exponent.

A second statement identifies the limit of further quotienting:
`W_p` has no proper nonzero subspace invariant under every mode of `t_1`.
Thus an additional quotient that preserves all those modes cannot shrink
the intermediate space.  This is a representation statement, not a
lower bound on factoring or on arbitrary crossing decompositions.

The source operators are those of Claude's Appendix B in
[`anthropics/formal-math` at `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5),
as implemented by the existing
[source recursion](../tools/check_modular_crossing.py).
The quotient and normalized crossing are the previous repository
construction in [modular-crossing.md](modular-crossing.md).
The new claims here are ordinary proofs with independent finite checks;
there is no new Lean formalization.

## 1. A normal-order symbol calculation

Use commuting formal symbols `X_j` and `D_j` for multiplication by `t_j`
and differentiation by `t_j`.  A monomial in these symbols is interpreted
in normal order: all differentiations act on the original input before
all multiplications.  Multiplication of symbols is commutative, although
composition of the corresponding operators need not be commutative.

For a positive integer `m`, the symbol of the source field is

\[
 L_m(z)=
 \sum_{j\geq1}\binom{j-1}{m-1}X_jz^{j-m}
 +(-1)^{m-1}\sum_{j\geq1}
 j\binom{j+m-1}{m-1}D_jz^{-j-m}.
 \tag{1}
\]

The binomial coefficient in the first sum is zero when `j<m`.
For a monomial `a=t_{m_1}\cdots t_{m_k}`, the symbol is the product
`L_{m_1}(z)\cdots L_{m_k}(z)`.  Its coefficient of
`z^(r-wt(a))`, interpreted in normal order, is `a[r]`.
This is the defining normal-order expansion; it does not replace
normal ordering by ordinary operator composition.

The two sums can be manipulated coefficientwise on polynomial states.
For a fixed state and output monomial, only finitely many differentiated
and created variables can contribute.  Frobenius below is therefore an
identity of finite coefficients and does not require convergence of a
two-sided Laurent series.

In `W_p`, every normal-order symbol divisible by one of

\[
 X_{pj},\quad X_j^p,\quad D_j^p
 \tag{2}
\]

acts as zero.  The first two claims follow from the state ideal.
For the last, every surviving exponent is less than `p`, so a `p`th
derivative is zero.  Normal multiplication by additional symbols
preserves this zero action: a created forbidden factor remains
forbidden, and an additional derivative does not remove a factor that
will only be created afterward.

## 2. Multiples of `t_(pk)` act as zero

Set `m=pk`.  The units digit of `m-1` in base `p` is `p-1`.
Lucas's identity shows that

\[
 \binom{j-1}{pk-1}\not\equiv0\pmod p
 \quad\Longrightarrow\quad p\mid j.
\]

Every surviving creation term in (1) is therefore divisible by
`X_(pj)`.  For the annihilation coefficient, if `p` does not divide `j`,
the units digit of `j+pk-1` is smaller than `p-1`, so its binomial
coefficient vanishes.  If `p` divides `j`, the explicit factor `j`
vanishes.  Thus the entire annihilation sum is zero.

Every normal monomial in `L_(pk)` acts as zero for the reason in (2).
Multiplying its symbol by the symbol of any other parameter retains
that property.  Hence every parameter divisible by `t_(pk)` acts as
zero on `W_p`.

## 3. Multiples of `t_k^p` act as zero

In the commutative normal-order symbol algebra, Frobenius gives

\[
 L_k(z)^p=
 \sum_{j\geq1}\binom{j-1}{k-1}^{p}X_j^p z^{p(j-k)}
 +(-1)^{p(k-1)}\sum_{j\geq1}
 \left(j\binom{j+k-1}{k-1}\right)^pD_j^p z^{-p(j+k)}.
 \tag{3}
\]

Every term acts as zero by (2), including after normal multiplication
by the symbol of an arbitrary other parameter.  Therefore every
parameter divisible by `t_k^p` acts as zero.  The generators in
Sections 2 and 3 generate `I_p`; linearity proves parameter descent.

The argument works after scalar extension to every commutative
`F_p`-algebra.  It uses no inverse or nonzero test on that coefficient
ring.

## 4. Reduce the dual polynomials before evaluating them

For a regular partition `mu`, the source dual parameter is

\[
 a^\mu=\prod_{m\in\mu}
 \left((-1)^{m-1}\sum_{j=1}^{m}\binom mjt_j\right).
 \tag{4}
\]

Let `red_p` delete monomials containing a part divisible by `p` or a
multiplicity at least `p`, and reduce the remaining coefficients
modulo `p`.  The theorem gives

\[
 a^\mu[r]v=\mathrm{red}_p(a^\mu)[r]v
 \qquad(v\in W_p).
 \tag{5}
\]

The source angle operation multiplies the weight-`w` component by
`binom(h+w,i)`.  The ideal is homogeneous, so this operation preserves
`I_p`; (5) also holds after every angle operation.  Each normalized
crossing piece can therefore be computed with the reduced dual.
The resulting tensor is identical, entry by entry, to the existing
full-parameter implementation.  An operator parameter must still
retain all its surviving weights; quotienting does not homogenize it.

At `p=2`, regular monomials are products of distinct odd-index variables.
Writing `m=2a+1` and `j=2b+1`, Lucas's identity gives the especially simple
one-part formula

\[
 \mathrm{red}_2(a^{(2a+1)})
   =\sum_{b=0}^{a}\binom ab\,t_{2b+1}\pmod2.
 \tag{6}
\]

Products in (6) are formed in the square-zero algebra.  A repeated odd
variable contributes zero.  This is a direct exact compiler for the
reduced dual, rather than a claim that all source field identities can
be inferred from ordinary first derivatives in characteristic two.

## 5. No smaller quotient preserves every `t_1` mode

The modes of `t_1` include multiplication by each `t_j` and the
derivation `j partial/partial t_j`.  For the surviving variables
`p` does not divide `j`, so the latter differs from the ordinary
partial derivative by a unit.

Let `U` be a nonzero subspace of `W_p` invariant under all those modes,
and choose a nonzero polynomial `f` in `U`.  Among its monomials choose
one with exponent vector `alpha` maximal in the componentwise partial
order.  Apply `partial^alpha` to `f`.  Every other monomial vanishes:
if it survived, its exponent vector would be at least `alpha`, contrary
to maximality.  The chosen monomial becomes its coefficient times
`prod alpha_j!`, a nonzero scalar since each `alpha_j<p`.

Thus `1` lies in `U`.  Closure under multiplication produces every
regular monomial, so `U=W_p`.  This proves the statement.

An infinite quotient is not being treated as a finite-dimensional
matrix here.  Each chosen polynomial has finite support, and the
derivative argument uses finitely many modes.  The conclusion does not
cover restrictions preserving only selected label operators, joint
refactorizations of complete pieces, or a different operator algebra.

## Reproduction and open question

Run

```sh
python3 tools/check_modular_parameter_descent.py
```

The checker compares the unreduced source recursion with zero on
nonregular operator parameters, compares full and directly reduced
dual angles, checks the complete normalized crossing using the new
compiler, and includes an incorrect characteristic-zero deletion as a
negative control.  All comparisons are exact finite-field arithmetic.
The emitted record gives the tested parameter and state ranges and
the number of entries actually checked.

The next question is whether selecting only a proper family of label
operators permits a smaller invariant channel space while preserving
nondegenerate duality and both crossing factorizations.  Section 5
excludes the unrestricted `t_1`-mode quotient, and the existing
declared channel entropy still yields the same half-rate exponent.

This note and checker used substantial OpenAI assistant support.
