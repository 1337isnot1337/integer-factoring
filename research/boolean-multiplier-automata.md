# An exact finite description of arbitrary one-support multipliers

## Status

This note reduces an unbounded polynomial-filter family to a finite set of
explicit cost functions. The reduction is constructive. It includes arbitrary
polynomial degree, arbitrary characteristic-two coefficient extensions, and
arbitrary allowed ramification. It does not prove a factoring improvement or a
uniform pressure floor for the resulting family.

The family has the form

\[
 f_k(x)=1+R(x)g(k),                                      \tag{1}
\]

where \(g\) is a fixed nonconstant Boolean polynomial and every endpoint
\(f_k(1)\) is nonzero. Thus \(R(1)\ne1\). The multiplier \(R\) is an ordinary
polynomial, not a Laurent polynomial. General coefficient vectors that do not
lie in the two-dimensional span of \(1\) and \(g\) are outside this result.
The choice \(R=0\) is the always-available baseline: every filter coefficient
is 1, its infinity degree is 0, and its relative cost is 1. No ordinary
polynomial degree is assigned to the zero multiplier itself.

Write \(m=nP\), \(X=x^m\), \(\tau=nt\), and decompose the multilinear Boolean
polynomial into homogeneous parts \(g=g_0+\cdots+g_d\). Here and below, sums of
Boolean values are in \(\mathbf F_2\).

## 1. Only the homogeneous-layer pattern matters

For a physical occupation state \(k\), put

\[
 H_k(Y)=\sum_{j=0}^{d}g_j(k)Y^j.
\]

The exact finite-difference transform of (1) is

\[
 \beta_k=X^{|k|}+R(x)T_k(X),\qquad
 T_k(X)=\sum_jg_j(k)X^{|k|-j}.                          \tag{2}
\]

If \(H_k=0\), then \(\beta_k=X^{|k|}\). Otherwise let
\(a=\deg H_k\), and write

\[
 T_k(X)=X^{|k|-a}A(X),\qquad
 A(X)=X^aH_k(X^{-1}),\quad A(0)=1,\quad\deg A\le a.
\]

Define the binary rational series

\[
 S_H(X)=\frac{X^a}{A(X)}.                              \tag{3}
\]

If \(h\) is the first exponent at which \(R\) differs from \(S_H\), measured
in units of \(m\), then

\[
 \frac{\operatorname{ord}_0\beta_k}{m}-|k|=h-a.         \tag{4}
\]

Thus the zero-place relative exponent depends only on the \((d+1)\)-bit
pattern \((g_0(k),\ldots,g_d(k))\), not on \(|k|\). There are at most
\(2^{d+1}\) patterns, even when the filter uses many physical modes.

For a binary integral multiplier of degree \(D\), meaning \(R\in\mathbf F_2[X]\),
endpoint admissibility is exactly \(R(1)=0\). Its degree at infinity is
\(Dg(k)\) in units of \(m\). If \(q_H\) is the probability of layer pattern
\(H\) in the charged physical measure and \(p=\Pr(g=1)\) in the uncharged
measure, the full relative cost is

\[
 \left(\sum_Hq_H e^{-2\tau\Delta_H}\right)
       \left((1-p)+p e^{2\tau D}\right),               \tag{5}
\]

where \(\Delta_H\) is the excess in (4), and \(\Delta_0=0\). No equality
between the two physical measures is assumed.

## 2. The first exceptional coefficient gives an exact continuous family

Call a term exceptional if its exponent in units of \(m\) is not an integer,
or if its coefficient is outside \(\{0,1\}\). Suppose \(\ell\) is the first
exceptional exponent of a multiplier. All earlier terms form a binary
polynomial prefix \(c(X)\).

For each target (3), either that prefix already differs from the target at an
integer exponent below \(\ell\), or the first difference occurs at \(\ell\).
In the latter case, no later term can cancel the first exceptional term.
Consequently every zero order is unchanged if the original multiplier is
replaced by

\[
 R_\omega(x)=c(X)+\omega x^{m\ell},\qquad
       \omega\in\mathbf F_4\setminus\mathbf F_2.        \tag{6}
\]

Terms of \(c\) at the same exponent as \(\ell\), when \(\ell\) is an integer,
may be included or omitted: adding \(\omega\) still gives a nonbinary
coefficient. The degree of (6) is \(\ell\), which is no larger than the
original degree. Its endpoints are nonzero because
\(1+c(1)+\omega\in\mathbf F_4^*\). Therefore (6) weakly improves the full cost.

This argument also covers a multiplier with binary coefficients whose first
exception is only a nonintegral exponent. It is not necessary to repair an
endpoint with a second nearby exponent. The nonbinary coefficient in (6)
already supplies a valid endpoint color.

Let a binary prefix agree with a set of active targets. The inactive targets
have fixed first differences. On an interval of possible cutoffs, their total
zero-place contribution is a constant \(A\), while the active contribution is
\(B e^{-2\tau\ell}\). Thus (5) becomes

\[
 (A+B e^{-2\tau\ell})\bigl((1-p)+p e^{2\tau\ell}\bigr). \tag{7}
\]

For \(A,B>0\) and \(0<p<1\), its optimum is obtained by clipping

\[
 \ell_*=
 \frac{\log B+\log(1-p)-\log A-\log p}{4\tau}           \tag{8}
\]

to the permitted interval. Zero-probability terms have the corresponding
one-sided limits. Every rational cutoff allowed by \(m=nP\) is realized by an
actual polynomial. Such cutoffs become dense when odd \(P\) is allowed to grow,
so the real optimization gives the exact infimum over ramifications. A strict
gain at a real cutoff can be preserved by a sufficiently close rational cutoff.

## 3. A finite binary prefix tree

For a nonzero layer pattern of degree at most \(d\), the target (3) is already
in lowest terms: the denominator has constant coefficient 1. Two distinct
targets \(X^a/A\) and \(X^b/B\) differ in an exponent at most \(2d\), since
the nonzero numerator \(X^aB+X^bA\) has degree at most \(2d\).

Build the binary prefix tree of all these rational series. At a node of depth
\(s\), a fixed prefix has coefficients through exponent \(s-1\). Cutoffs in
\([s-1,s]\) give a continuous family (7). If a node has just one active
rational target, every later cutoff has the same inactive contribution and
the same active coefficient. That node supplies the entire ray
\([s-1,\infty)\), so its descendants need not be stored.

If a binary branch leaves every target, all beta orders are fixed. Among
endpoint-trivial binary completions, choose one of least degree: keep the
prefix if it has even parity, or add a single coefficient at the next exponent
if it has odd parity. Higher coefficients cannot improve any zero order and
can only increase the infinity degree. This handles every such branch by one
finite profile.

It remains to describe a binary branch that follows one target for an
arbitrarily long time. That is the only source of unbounded discrete choices.

## 4. Binary tails are finite arithmetic progressions

Fix a singleton target \(S=X^a/A\), and suppose its first difference from a
binary multiplier is at exponent \(h\). The least-degree endpoint-trivial
multiplier with that first difference is obtained as follows:

1. Copy the coefficients of \(S\) below \(h\).
2. Flip the coefficient at \(h\).
3. If the resulting polynomial has odd parity, add \(X^{h+1}\).

Let its degree be \(D=h+\delta\). Every other multiplier with that same first
difference has the same beta orders and degree at least \(D\).

The coefficient sequence of \(S\), its prefix parity, and the gap since its
last nonzero coefficient form a finite eventually periodic state. An exact
way to generate it starts with remainder \(u_0=X^a\). At each position, set

\[
 s_h=u_h(0),\qquad u_{h+1}=(u_h+s_hA)/X.                \tag{9}
\]

After the finite initial segment, \(\deg u_h<\deg A\). If \(A\ne1\), the
remainder is nonzero and the recurrence is invertible on this finite set.
A run of \(\deg A\) zero coefficients would force the remainder to be zero,
so the gap is bounded. If \(A=1\), then \(S=X^a\), and the tail is constant
after its single nonzero coefficient. The implementation records the exact
cycle rather than relying on an estimated period bound.

It follows that every tail is a finite startup followed by finitely many
arithmetic progressions

\[
 h=h_0+Lz,\qquad D=h+\delta,\qquad z=0,1,2,\ldots,      \tag{10}
\]

with constant \(\delta\) on each progression. If \(\delta\ge0\), the continuous
cutoff (6) at \(\ell=h\) has the same zero orders and no larger infinity
degree, so that progression is dominated. Only negative degree gaps need to
remain in the exact lower envelope.

For a retained progression the cost is

\[
 (A+B e^{-2\tau h})\bigl((1-p)+p e^{2\tau(h+\delta)}\bigr). \tag{11}
\]

Its unconstrained optimum is

\[
 h_*=
 \frac{\log B+\log(1-p)-\log A-\log p-2\tau\delta}{4\tau}. \tag{12}
\]

Convexity as a function of \(h\) shows that only the two nearest allowed
positions in (10), clipped at \(h_0\), need to be evaluated. Hence the
unbounded polynomial degree is eliminated exactly.

For universal layer patterns of degrees 1 through 5, the implementation
produces respectively 5, 13, 29, 61, and 125 continuous pieces, together with
1, 4, 13, 42, and 164 negative-gap periodic rows. No multiplier degree cutoff
appears in those counts.

## 5. Code-coset filters on arbitrarily many physical modes

Let \(g(k)\) indicate the binary syndrome equation
\(\sum_i h_i k_i=b\), where \(h_i\in(\mathbf Z/2)^c\). Write \(e_h\) for a
group-algebra basis element. Its homogeneous-layer generating polynomial is

\[
 H_k(Y)=[e_b]\prod_{i:k_i=1}\bigl(1+Y(1+e_{h_i})\bigr). \tag{13}
\]

Indeed, multiplying by \(X^{|k|}\) and substituting \(Y=X^{-1}\) turns the
product into \(\prod_i(X+1+e_{h_i})\), whose coefficient of \(e_b\) is the
original transform \(T_k(X)\). The augmentation ideal of this group algebra
has nilpotence degree \(c+1\), so (13) has degree at most \(c\).

Moreover, \((1+e_h)^2=0\), whence

\[
 (1+Y(1+e_h))^2=1.                                    \tag{14}
\]

Thus only the parity of occupations in each group of equal nonzero columns
matters. Distinct column groups use disjoint physical modes. Under independent
occupation probabilities \(p_i\), their odd-parity probabilities are exactly

\[
 \Pr(z_h=1)=\frac{1-\prod_{i:h_i=h}(1-2p_i)}2.           \tag{15}
\]

For the physical Gibbs measures, \(1-2p_i=\tanh(t\lambda_i/2)\), or its charged
version with \(\lambda_i+J\). Formula (15) retains every actual mode energy.
Zero columns are spectators. A fixed set of column types therefore gives an
exact computation for arbitrarily many physical modes, using only one parity
bit per distinct nonzero type and at most \(2^{c+1}\) layer patterns.

## Implementation

[`boolean_multiplier_profiles.py`](../tools/boolean_multiplier_profiles.py)
contains exact binary arithmetic, finite-degree prefix covers, the complete
unbounded automaton, and the group-algebra signatures.
[`boolean_multiplier_pressure.py`](../tools/boolean_multiplier_pressure.py)
evaluates its explicit continuous and arithmetic-progression minima.

Numerical optimization of scale, charge, or the code-column assignment remains
an exploratory search. The finite automaton removes the multiplier degree and
coefficient-field cutoffs; it does not remove those outer optimization tasks.

The independent replay in
[`check_boolean_multiplier_automata.py`](../tools/check_boolean_multiplier_automata.py)
passed 36 complete finite-degree cover comparisons, 20,475 exact integral
multiplier routes through the unbounded automaton, and 9,915 explicit periodic
tail controls. It also checked 411 literal F4 cutoff witnesses on 3,948 physical
states, 70,656 independent code homogeneous-layer values, and 172 exact
heterogeneous pattern distributions. Forty floating cost comparisons serve as
additional implementation checks; they are not used to prove the reductions.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_boolean_multiplier_automata.py
```

The record is
[`boolean-multiplier-automata.json`](../evidence/boolean-multiplier-automata.json).
