# Complete binary quadratic coefficient catalogue

## Result and scope

For every scalar tensor A with sixteen arbitrary entries over F2, form

```math
K_A(z)=I+zA+z^2(E+I+A),\qquad
I(x,y,u,v)=[x=u][y=v],\quad E(x,y,u,v)=[x=v][y=u].
```

The four labels are bits. Then K_A(0)=I and K_A(1)=E. This includes
non-Gaussian tensors and makes no prescribed choice of regular Smith
frames. All 65,536 tensors produce 162 distinct paired profiles. For the
row cut xy|uv at infinity and the column cut xu|yv at zero, their bare
Smith pressure satisfies, at every real s>=1,

```math
B_A(s)\ge\frac{256}{27}.
```

This is an exact finite computational theorem about the specified
coefficient space and certificate. It neither excludes higher degree,
different constant terms, auxiliary states, physical gauges, larger fields,
or rank improvements from combining coefficients, nor establishes a
factoring lower bound. It supplies no factoring exponent above 1/2.

The polynomial/Smith method and the degree-one baseline are developed in
[PR #30](https://github.com/1337isnot1337/integer-factoring/pull/30) and
[PR #39](https://github.com/1337isnot1337/integer-factoring/pull/39), respectively.
The coefficient catalogue here is separate from PR #39's prescribed
affine-target and polynomial-frame searches. Base main: de7bedb.
Substantial OpenAI Codex assistance was used. No Lean files change.

## Profile convention and coefficient cost

Let mu be the four Smith orders of t^2 K_A(t^-1) in the row cut,
and nu those of K_A(t) in the column cut. Both matrices have full generic
rank: at t=1 their flattenings are permutation matrices. Set x=1/s. Then

```math
B_A(s)=s^2\left(\sum_i s^{-\mu_i}\right)
                 \left(\sum_j s^{-\nu_j}\right)
       =x^{-2}P_{\mu,\nu}(x).
```

The degree regularization remains two even when the actual leading
coefficient vanishes. The infinity orders then increase accordingly, so
the expression is unchanged by this bookkeeping choice.

For clarity, the usual coefficient bound follows directly from a local
Smith factorization over F2[[t]]. At zero the coefficient-k matrix of a
tensor power of length m has rank at most (k+1)s^k(sum s^-nu)^m.
At infinity, extract its degree-(2m-k) coefficient from the reversed
polynomial: its rank is at most (2m-k+1)s^(2m-k)(sum s^-mu)^m.
Indeed each tensor-product Smith channel of order a contributes at most
(k-a+1) to the coefficient rank at zero, and only channels with a<=k
contribute. Bound their number by s^k(sum s^-nu)^m; the infinity
argument is identical. Multiplying eliminates k and gives the displayed
pressure, up to factors polynomial in m. This is a sufficient rank
certificate, so a lower bound on it is not a lower bound on actual ranks.

## Exhaustive calculation and independent replay

The [native enumerator](../tools/enumerate_quadratic_crossings.cpp) visits
every 16-bit coefficient mask in increasing order. Bit 8x+4y+2u+v is
A(x,y,u,v). It computes the valuation of every nonzero minor of sizes
one through four in both polynomial matrices. Successive differences
of the minimum valuations give the Smith orders. The determinant degree
is at most eight. No numerical rank threshold is used.

The [independent checker](../tools/check_quadratic_crossings.py) instead
uses ordinary F2 ranks of block Toeplitz matrices modulo t^k. If R_k is
this rank, R_k-R_(k-1) counts orders below k. Precision nine recovers all
orders. A full replay checks every tensor, profile multiplicity, and first
representative against the [complete catalogue](../evidence/quadratic-crossing-catalogue.json).
The [full replay record](../evidence/quadratic-crossing-full-replay.json)
confirms all 65,536 masks; it is separate from the default sampled run.

Eighteen profiles are maximal under coordinatewise order dominance.
Every other profile has pressure at least one of these eighteen at every
s>=1. For each maximal profile the checker proves

```math
27P_{\mu,\nu}(x)-256x^2\ge0\qquad(0\le x\le1).
```

The proof converts the polynomial to its rational Bernstein coefficients
and subdivides intervals at exact midpoints until all coefficients are
nonnegative. Bernstein basis functions are nonnegative and sum to one,
so each accepted interval certifies the inequality throughout that interval.
The saved [certificate summary](../evidence/quadratic-crossing-certificate.json)
records all eighteen profiles and subdivision counts; the checker
reconstructs every coefficient rather than trusting those counts.

For the two equality profiles it first divides out an exact square:
(1-3x)^2 for ((1,2,2,3),(0,1,1,2)), and (1-3x^2)^2 for
((0,2,2,4),(0,2,2,4)). It then certifies the quotient's nonnegativity
in the same way. The bound is sharp: mask 1632 is the elementary
pencil (1+z)I+zE, with B=(1+x)^4/x and equality at x=1/3;
mask zero replaces z by z^2, with equality at x=1/sqrt(3).

## Reproduction and next step

```sh
c++ -O3 -std=c++17 -Wall -Wextra -Werror tools/enumerate_quadratic_crossings.cpp -o /tmp/enumerate-quadratic
/tmp/enumerate-quadratic > /tmp/quadratic-catalogue.json
cmp /tmp/quadratic-catalogue.json evidence/quadratic-crossing-catalogue.json
python3 tools/check_quadratic_crossings.py
python3 tools/check_quadratic_crossings.py --full
make verify
```

The default checker verifies every saved representative and a deterministic
sample of other tensors. The full replay is separate and covers all masks.
The next constructive question is whether cubic or higher interpolations
with these endpoints, or a changed constant tensor, beat the elementary
pressure. The quadratic result does not answer that question.
