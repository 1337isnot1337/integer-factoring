# Exact separation for continuously gauged Smith profiles

## Result and scope

A weighted Smith profile can be computed without listing every minor.
Minimum-valuation pivoting, implemented with fraction-free polynomial
elimination, produces a nested sequence of minors attaining all
determinantal minima. This gives an exact separation oracle for the
linear constraints used in continuous Gaussian gauge searches.

The polynomial arithmetic and rational-weight comparisons are exact.
An optimization driver may still use floating-point heights; such
search output needs an exact rational replay before it can certify a
rank-width bound. This note establishes the oracle and the associated
convex formulation. It does not establish a new factoring exponent.

The fraction-free recurrence is classical: see Erwin H. Bareiss,
[“Sylvester's Identity and Multistep Integer-Preserving Gaussian
Elimination,” Mathematics of Computation 22 (1968), 565–578](https://www.ams.org/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf).
The weighted pivot rule and its use in this repository's pressure
searches are explained below. The implementation and checks were
developed with substantial OpenAI assistant support.

## 1. Weighted determinantal minima

Let \(A(z)\) be an invertible \(n\)-by-\(n\) matrix over
\(\mathbb F_2(z)\), represented initially by polynomial entries.
Use either the zero valuation
\(\nu(f)=\operatorname{ord}_0f\), or the infinity valuation
\(\nu(f)=-\deg f\) for nonzero polynomials.

For row weights \(r_i\), column weights \(c_j\), and a scale
\(\ell\ge0\), define

~~~math
\Delta_k(r,c,\ell)
 =\min_{\substack{|I|=|J|=k\\\det A_{I,J}\ne0}}
 \left(\ell\,\nu(\det A_{I,J})
       +\sum_{i\in I}r_i+\sum_{j\in J}c_j\right).        \tag{1}
~~~

Set \(\Delta_0=0\). The sorted weighted Smith orders are

~~~math
\kappa_k=\Delta_k-\Delta_{k-1}.                          \tag{2}
~~~

For rational weights and scale, a common ramification clears their
denominators. Real weights can be interpreted in an ordered value
group; the optimization formulation below also follows by continuity
from rational weights.

## 2. A minimum pivot finds every prefix

Suppose row and column sets \(I,J\) have already been selected,
with nonsingular pivot minor \(d=\det A_{I,J}\). The remaining
Schur-complement entry at \(i,j\) is

~~~math
\frac{\det A_{I\cup\{i\},\,J\cup\{j\}}}{d}.             \tag{3}
~~~

Signs disappear in characteristic two. After applying the row and
column weights, its valuation is

~~~math
\ell\left[
\nu(\det A_{I\cup\{i\},\,J\cup\{j\}})-\nu(d)
\right]+r_i+c_j.                                       \tag{4}
~~~

Choose a nonzero entry minimizing (4). If its valuation is \(m\),
every entry of the remaining matrix has valuation at least \(m\).
Dividing any entry in its pivot row or column by the pivot therefore
gives a regular element of the valuation ring. The corresponding row
and column eliminations are unimodular.

The next Schur complement has entries of valuation at least \(m\):
both the old entry and the eliminated product have that property.
Repeating this argument produces a diagonal Smith form with
nondecreasing pivot valuations.

At step \(k\), the determinant of the selected \(k\)-minor is the
product of the first \(k\) pivots, including their row and column
weights. Consequently the selected sets attain (1) simultaneously
for every \(k\). The nested prefix minors are a certificate of
which linear constraints are tight.

## 3. Polynomial arithmetic throughout

The implementation stores the bordered minors themselves. If
\(B_{a,b}\) denotes the current bordered-minor matrix, \(p=B_{i,j}\)
is the selected pivot, and \(d\) is the preceding selected minor,
the next bordered minors are

~~~math
B'_{a,b}
  =\frac{pB_{a,b}+B_{a,j}B_{i,b}}{d}.                    \tag{5}
~~~

Sylvester's determinant identity makes every division exact. Initially
\(B=A\) and \(d=1\). The arithmetic therefore remains in
\(\mathbb F_2[z]\), even though (3) is a rational Schur complement.

One query uses a polynomial number of polynomial multiplications and
exact divisions. If the entries have degree at most \(D\), every
stored minor has degree at most \(nD\). Cached states depend on the
selected row and column sets, so repeated gauge queries can reuse
them. The cache can be cleared when memory use matters.

## 4. Convex pressure optimization with a separation oracle

For fixed row and column heights, let
\(\boldsymbol\kappa\) be the true sorted Smith vector. Introduce
a sorted candidate vector \(\boldsymbol x\) with

~~~math
\sum_{j=1}^{k}x_j
 \le \ell\,\nu(\det A_{I,J})
       +\sum_{i\in I}r_i+\sum_{j\in J}c_j
\quad
\text{for every nonzero \(k\)-minor},                    \tag{6}
~~~

and impose equality at \(k=n\). These are linear constraints in
the heights, scale, and candidate orders.

The ascending partial sums in (6) say that
\(\boldsymbol x\) is at least as spread out as
\(\boldsymbol\kappa\), with the same total sum. Convexity gives

~~~math
\sum_j\log(1+e^{-x_j})
 \ge\sum_j\log(1+e^{-\kappa_j}).                         \tag{7}
~~~

The true vector is feasible. Thus minimizing the left side over all
candidate vectors in (6) recovers the true Smith pressure. For fixed
physical type, the gauge mean shift adds only a linear term.
Regularity constraints for the internal Gaussian factors are linear
as well. This is a convex optimization problem.

To separate a proposed point, run the oracle at its heights and
scale. At each \(k\), compare its proposed partial sum with the
attaining prefix minor from Section 2. A violated inequality supplies
one explicit row of (6); if all prefixes pass, every minor inequality
passes. It is therefore unnecessary to materialize the
\(\binom{2n}{n}-1\) possible nonempty minors of a generic square
matrix.

Optimizing the physical type as well is a further step. Alternating a
fixed-type convex solve with a Gibbs type update is a search method;
this note makes no global optimality claim for that alternating
procedure.

## 5. Exact API and reproduction

The reusable implementation is
[weighted_smith_oracle.py](../tools/weighted_smith_oracle.py).
Its matrix entries are binary polynomial bit masks. Integer and
fractional inputs give exact output:

~~~python
from fractions import Fraction
from weighted_smith_oracle import WeightedSmithOracle

oracle = WeightedSmithOracle([[1, 2], [2, 1]], infinity=False)
orders, prefixes = oracle.profile(
    [Fraction(0), Fraction(1, 2)],
    [Fraction(-1, 3), Fraction(2, 3)],
    scale=Fraction(1, 2),
)
~~~

Each prefix records its row mask, column mask, polynomial determinant,
and unweighted valuation. Passing floating-point inputs is useful for
search, but floating comparisons are not exact certificates.

Run

~~~bash
python tools/check_weighted_smith_oracle.py \
  --output evidence/weighted-smith-oracle.json
~~~

The independent checker enumerates minors by determinant expansion.
It compares every prefix minimum, including negative rational weights,
zero and positive scales, both places, random polynomial matrices, and
the nonscalar four-mode Schur core. It also checks the selected minor
polynomials themselves and rejects a singular negative control.
