# A certified limitation of the stabilized centered-bound computation

## Scope

This note concerns the **computed upper bound** from
[the stable Pfaffian module](stable-pfaffian-modules.md).
It is not a lower bound on actual coefficient ranks, on the best
decomposition of the crossing tensor, or on a restricted code.

For a full grade \(L\) on \(N=2P\) physical modes, let
\(d_0=\binom NL\). The centered construction uses
\(p=1+x,q=1+y,r=1\), and takes the minimum of its coefficient-rank upper
bounds over six constant ancillary entries.

For every fixed ratio \(0<\alpha\leq1/2\), with \(L/N\to\alpha\), this
particular computed cost satisfies

\[
 \boxed{
 \liminf_{N\to\infty}
 \frac{\log S_{\rm computed}}{\log d_0}
 \geq\frac{29}{20}>\sqrt2.
 }
 \tag{1}
\]

Thus the computed \(\rho\) has lower limit at least \(29/40=0.725\).
The sufficient finite-seed target is \(1/\sqrt2\).

The proof is elementary. A specific block coefficient is always an
all-ones matrix. A subset of paths in the already-proved two-dimensional
fusion recurrence gives enough multiplicity to establish (1).
No asymptotic theorem about tilting modules is assumed.

This explains why changing grade ratios or binary digit patterns inside
the existing full-grade bound is unlikely to be productive. A new code
can still remove the corresponding incoming directions. Such a code
needs its own rank analysis; substituting its smaller dimension into the
unchanged upper bound cannot remove the block identified here.

## 1. A block that survives every ancillary choice

At jet \(x^a y^b\), put

\[
 c=L-a-b,\qquad
 m=2c+1,\qquad
 n=N+1-2L+2b,\qquad d=a.
 \tag{2}
\]

Then \(m+n+2d=N+2\), as required by the stable decomposition.
For \(a,b\geq0\), \(a+b\leq L\leq N/2\), all relevant weight slices exist.

Every one of the six allowed ancillary rows has grade sum one.
Consequently the same \(m,n,d\) works for all six choices.
In the notation of the stable module, the \(h\) row weight is \(v_r=b\).
After removing \(x^d\), its coefficient is

\[
 [x^0y^b]
 (1+x)^{u_c-w}(1+y)^{b-z}(1+x+y)^z
 =[y^b]\,(1+y)^b=1.
 \tag{3}
\]

Thus the entire rectangular coefficient block is all ones and has rank
one. This holds for every weighted intersection, independently of binary
digit carries.

Let \(M_P(m,n,d)\) be the exact stable multiplicity. Each ancillary
upper bound contains this contribution, so their minimum does too:

\[
 U_{a,b}\geq M_P(m,n,d).
 \tag{4}
\]

For a symmetric jet \(a=b\), both cut bounds contain the same term.
After the usual dimension cap, that jet alone contributes at least

\[
 \min\{d_0^2,M_P(m,n,d)\}
 \tag{5}
\]

to the computed rank sum.

## 2. A direct lower bound on fusion multiplicities

Write \(\mu_s(m)\) for the multiplicity of \(B_m\) in the decomposition of
\(V^{\otimes s}\), suppressing its forced determinant twist.

The elementary fusion recurrence always has an upward edge
\(m\mapsto m+1\). From a positive even level \(2j\), it also has a
downward edge \(2j\mapsto2j-1\) with multiplicity two. The following odd
step returns upward.

Start at level zero and consider pairs of steps. Retain only:

- an advance from \(2j\) to \(2j+2\);
- a down-and-up loop at a positive even level, with weight two.

To reach \(2k\) after \(2s\) steps, the first advance is forced.
Choose the other \(k-1\) advances among the remaining \(s-1\) pairs;
all other pairs are loops. Appending one upward step gives

\[
 \boxed{
 \mu_{2s+1}(2k+1)
 \geq 2^{s-k}\binom{s-1}{k-1},
 \qquad 1\leq k\leq s.
 }
 \tag{6}
\]

This counts actual direct summands produced by the constructive
recurrence. All omitted paths can only increase the multiplicity.

In the stable sum, choose an even integer \(j\) for the number of \(g\)
factors before the ancillary factor. Equations (2), (6) give

\[
 M_P(m,n,d)\geq
 \binom Nj\,
 2^d
 \binom{j/2-1}{(m-3)/2}
 \binom{(N-j)/2-1}{(n-3)/2},
 \tag{7}
\]

whenever the binomial arguments are valid.

## 3. Entropy form of the bound

Let \(H\) denote binary entropy, with base-two logarithms. Set

\[
 \frac LN\to\alpha,\qquad
 \frac aN=\frac bN\to x,\qquad
 \frac mN\to u=2\alpha-4x,\qquad
 \frac nN\to v=1-2\alpha+2x,\qquad
 \frac jN\to t.
\]

For \(u<t<1-v\), the elementary binomial entropy estimates applied to
(7) give the exponential lower rate

\[
 \mathcal E(\alpha,x,t)=
 H(t)+\frac12\left[
 1-u-v+tH(u/t)+(1-t)H(v/(1-t))
 \right].
 \tag{8}
\]

The label entropy is \(H(\alpha)\). Integer rounding and parity choices
change these expressions only by \(o(N)\) in their logarithms.

We now choose \(x,t\) explicitly; no numerical optimization is needed for
the constant in (1).

## 4. Ratios \(0<\alpha\leq2/5\)

Choose

\[
 x=\alpha/4,\qquad t=7\alpha/5.
\]

Then \(u=\alpha\), \(v=1-3\alpha/2\), and both binomial intervals are
strictly valid.
Let \(H_3\) be the entropy of a three-part probability vector. The chain
rule rewrites (8) as

\[
 \mathcal E=
 \frac12H(7\alpha/5)+\frac12H(3\alpha/2)
 +\left[\frac34H_3(2/3,4/15,1/15)+\frac14\right]\alpha.
 \tag{9}
\]

Define \(D(\alpha)=\mathcal E-(29/20)H(\alpha)\).
It has \(D(0)=0\) by continuity and is concave on \([0,2/5]\).
Indeed, its second derivative is negative because

\[
 \frac{7/10}{1-7\alpha/5}
 +\frac{3/4}{1-3\alpha/2}
 >
 \frac{29/20}{1-\alpha}
 \qquad(\alpha>0).
\]

At the other endpoint,

\[
 100D(2/5)=
 \log_2
 \frac{3^{87}5^{15}}{2^{16}7^{28}11^{22}}>0.
 \tag{10}
\]

The integer numerator exceeds the integer denominator. Concavity and the
two endpoint values give \(D(\alpha)\geq0\) throughout this interval.

## 5. Ratios \(2/5\leq\alpha\leq1/2\)

Choose

\[
 x=\alpha/3,\qquad t=1/2.
\]

Then \(u=2\alpha/3,\ v=1-4\alpha/3\), and

\[
 \mathcal E=
 1+\frac{\alpha}{3}
 +\frac14\left[H(4\alpha/3)+H(2-8\alpha/3)\right].
 \tag{11}
\]

On the stated interval,

\[
 \mathcal E\geq
 1+\frac{2}{15}
 +\frac14\left[H(2/3)+H(14/15)\right]
 >\frac{29}{20}.
 \tag{12}
\]

The last strict inequality is equivalent to the exact integer comparison

\[
 45^{15}>2^{29}14^{14}.
 \tag{13}
\]

Since \(H(\alpha)\leq1\), (12) also gives
\(\mathcal E/H(\alpha)>29/20\).
Together with Section 4, this proves (1).

The dimension cap in (5) does not affect the conclusion: its alternative
exponent is \(2\), which is larger than \(29/20\).

## 6. Finite exact diagnostics

The checker verifies the fusion-path inequality against the exact
recurrence, verifies the rank-one coefficient for all six ancillary
choices, and checks both large-integer inequalities.

It also computes exact stable multiplicities at the numerically
interesting boundary near \(\alpha=0.2255,\ x=0.05657\).
These are lower bounds only on the stabilized upper-bound computation:

| \(P\) | \(L\) | Symmetric jet \(a=b\) | Lower \(\rho\) from this one block |
|---:|---:|---:|---:|
| 16 | 7 | 2 | 0.7439715726 |
| 32 | 14 | 4 | 0.7482075212 |
| 64 | 29 | 7 | 0.7548905194 |
| 128 | 58 | 14 | 0.7568398933 |
| 256 | 115 | 29 | 0.7583194416 |
| 512 | 231 | 58 | 0.7591326331 |

These values are consistent with the larger entropy envelope seen in
finite screens, but no limiting equality is asserted.

## Reproduce

~~~bash
python3 tools/check_stable_boundary_envelope.py --maximum-pairs 512 \
  --output evidence/stable-boundary-envelope.json
~~~

The relevant next target is a restricted code whose incoming coordinates
annihilate the forms (3). The full-grade stable decomposition alone cannot
decide whether such a code retains enough label dimension.
