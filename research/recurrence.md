# Quantitative recurrence for the new parameters

This file proves the arithmetic recurrence estimates. Applying them to a
factoring machine also requires the modified crossing factors, generalized
grouping theorem, and implementation/cost interfaces stated in the integration
note. None of the tests below substitutes for those interfaces.

## An integer schedule

For `d>=0`, use

\[
 K_0=48,\quad D_d=(d+48)^2,\quad
 B_d=\lfloor\log_2K_d\rfloor+1,
 \quad G_d=12D_d^2K_d^2B_d,\quad K_{d+1}=G_d^2.
\]

These are exact integers. Since `B_d>=ln K_d`, the grouping theorem's
condition `G_d>=12D_d^2K_d^2 ln K_d` holds without evaluating a real logarithm.

The elementary bounds

\[
 K_d\ge48\cdot2^d\ge d+48,\qquad
 D_d\le K_d^2,\qquad B_d\le K_d,\qquad
 K_d^2\le G_d\le K_d^8
\]

hold by induction. The last upper bound uses `12<=K_d`. Consequently

\[
 K_d^4\le K_{d+1}\le K_d^{16}.
\]

The upper exponent 16 in this coarse inequality is useful only for bounding
overhead. It is not the sharp growth rate of the schedule.

## The actual growth rate is four

Put `a_d=ln K_d`. The exact recurrence is

\[
 a_{d+1}=4a_d+2\ln12+8\ln(d+48)+2\ln B_d.
\]

The coarse upper bound gives `a_d<=16^d ln48` and therefore
`B_d<=7*16^d`. Also `ln(d+48)<=ln48+d/48`. Thus

\[
 4a_d\le a_{d+1}\le4a_d+C_0+C_1d,
\]

where the fixed constants
`C_0=2ln12+8ln48+2ln7` and `C_1=1/6+2ln16` suffice.
Dividing by `4^(d+1)` and summing proves

\[
 (\ln48)4^d\le\ln K_d\le
 \left(\ln48+\frac{C_0}{3}+\frac{C_1}{9}\right)4^d.
\]

In particular `ln K_d=Theta(4^d)`, with fixed constants.

For sufficiently large input bit length `n`, choose the largest `ell` with
`K_ell^40<=n`, and write `P=K_ell`. Then

\[
 \ell=\log_4\ln n+O(1),\qquad P\le n^{1/40}.
\]

The finitely many smaller input lengths can use ordinary trial division; this
does not change the asymptotic bound. An implementation can find `ell` by
successive integer powers and comparisons.

## Summable contraction errors

Let

\[
 \rho=\frac{8+\log_2 7}{12},\qquad
 q_d=\rho+\frac1{G_d}+\frac1{D_d}.
\]

The rational envelope `rho=281/312` may be used throughout instead. Both
choices are greater than `9/10` and less than one.

For every interval of levels,

\[
 \prod_{j=d}^{s-1}q_j\le\frac{11}{10}\rho^{s-d}.
\]

To prove this, use `1+x<=exp(x)` and the convergent sums

\[
 \sum_{j\ge0}\frac1{D_j}\le\frac1{47},\qquad
 \sum_{j\ge0}\frac1{G_j}\le
 \sum_{j\ge0}\frac1{48^2 4^j}=\frac1{1728}.
\]

Their sum divided by `rho` is smaller than `1/11`. For `0<=x<1`,
`exp(x)<=1/(1-x)` follows by comparing power series. The claimed factor
`11/10` follows. This also explains why a summable `1/D_d` was selected:
a fixed positive slack would worsen the leading contraction constant.

## Every level's cost

Suppose the clique at height `d` has volume `N_d`, and the new crossing and
grouping construction satisfies

\[
 N_d\le q_dN_{d+1}+3G_d^4.
\]

Let the top height be `ell` and `N_ell=N_top`. Iterating gives, with the
fixed constant `C=11/10`,

\[
 N_d\le C\rho^{\ell-d}N_{\rm top}
       +3C\sum_{j=d}^{\ell-1}\rho^{j-d}G_j^4
 \le C\rho^{\ell-d}N_{\rm top}+3C\ell P^2.
\]

The last inequality uses `G_j^4=K_(j+1)^2<=P^2`. In particular the brute
force at the bottom needs exponent
`O(rho^ell N_top+ell P^2)`.

The table construction at height `d` has its principal exponent proportional
to `N_(d+1)/K_d`. This is controlled at every height, not only at the bottom:

\[
 \frac{\rho^{\ell-d-1}}{K_d}
 \le\frac{\rho^\ell}{48\rho}(2\rho)^{-d}
 \le\frac{\rho^\ell}{48\rho}.
\]

Hence all principal table exponents are
`O(rho^ell N_top+ell P^2)` as well. This bound prevents an intermediate
level from silently dominating the bottom-level estimate.

## Overhead requirements and their scale

The following are the required implementation interfaces. They are stated
explicitly so that a future formalization can verify each one.

* Keep the original top multiplication-network grouping at
  `g_top=3P^8`, separately from the new recursive grid sizes `G_d`. Its
  bound is `N_top=O(n+P^8 n^(3/4))=O(n)` under `P^40<=n`.
* Exhaustive search over groupings costs `K_d^G_d` candidates. Its binary
  logarithm is `G_d log2 K_d`. Even using the top parameter, this is
  `O(P^8 log P)`, and it fits below the final main exponent.
* Choosing one coefficient per rectangle creates at most
  `(O(n)+1)^(O(ell P^2))` paths for the actual recursive levels. The
  coarser bound `exp(O(ell G_ell^4 log n))` also suffices.
* Computing one compressed crossing channel, including its binomial
  coefficient and the new integer scale, takes polynomial time in the
  input and label lengths. The exact channel functions are given in the
  crossing-compression note.
* Table access, arithmetic, and compilation to a deterministic finite
  machine must introduce a fixed polynomial overhead, with degree
  independent of `ell`. A separate machine for each fixed depth would
  not establish a uniform bound.

For coarse budgets comparable to those used upstream, `g_top^4=81P^32`,
`ell<=P`, and `log n=O(log P)` under the maximality of `ell`. Thus even
an overhead exponent `O(P^38)` is `O(n^(19/20))`. The top boundary
term `g_top n^(3/4)` is at most `3n^(19/20)`. Both are
`o(n/(log n)^alpha)` for the constants below. A fixed polynomial
simulation changes only the constant multiplying the final exponent.

With these interfaces established, the recurrence yields

\[
 T(n)=\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right),
 \qquad
 \alpha=\frac{\log(1/\rho)}{\log4}.
\]

The limiting value is `alpha=0.07551046810477115...`; the rational envelope
gives `alpha=log(312/281)/log4=0.07548794932066111...`. Here `n` is the
bit length of the integer to factor, not the integer itself.
