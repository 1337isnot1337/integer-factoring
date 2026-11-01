# An asymptotic floor for the computed Ward rectangular bound

## Scope and conclusion

The [dual-incidence theorem](noncommutative-crossing-ward-charge-filtration.md)
removes the column-channel multiplier from the Ward descent. Its use
with unrestricted rectangular Pfaffian bounds still has a structural
limitation.

Let \(P\to\infty\), and suppose binary forms and a homogeneous input
code have parameters
\[
 k/P\to\kappa,\qquad L/P\to\lambda,
 \qquad 0\le2\kappa<\lambda\le1.
\]
Let \(2\le d\le D_L\), where \(D_L\) is the full homogeneous
product-image dimension. Let \(S_{\rm computed}\) be the bound
implemented by
[ward_channel_free_bounds.py](../tools/ward_channel_free_bounds.py),
after restricting its column descent to the exact shift \(k-j\).
This includes its minima with the original full-grade bound, the
earlier Ward bounds, physical rectangle dimensions, and \(d^2\).
It also permits the binary joint-order restriction.

Then
\[
 \boxed{\qquad
 \liminf_{P\to\infty}
 \frac{\log_2 S_{\rm computed}}{2\log_2d}
 \ge\frac{29}{40}=0.725>\frac1{\sqrt2}.
 \qquad}
 \tag{1}
\]
This is a lower bound on a particular **computed upper bound**. It
is not a lower bound on the actual crossing ranks or on factoring.
It does not cover \(L/P>1\), a vanishing excess ratio
\((L-2k)/P\to0\), or an additional minimum with a different rank
factorization. In particular, the global-form Euler parity bound is
a separate computation. Finite seeds with better actual ranks are
not ruled out.

## 1. Four protected blocks cover all implemented minima

Put \(\Delta=L-2k\). Choose an integer \(a\) with
\(0\le2a\le\Delta\), and set
\[
 b=a+2k,\qquad
 m=2\Delta+1-4a,\qquad
 n=2P+1-2\Delta+2a.
 \tag{2}
\]
Let \(M_P(m,n,t)\) denote the stable multiplicity of the block with
two highest weights \(m,n\) and determinant twist \(t\). The
following four blocks occur in the bounds for the coefficient
\(G[a,b]\). The residual jet is measured after removing the twist.

| Bound component | Block \((m',n',t)\) | Rectangle grades \((R_A,R_Y;C_V,C_U)\) | Residual jet |
|---|---|---|---|
| Ward row, term \(j=k\) | \((m,n,a)\) | \((\Delta,\Delta;L,L)\) | \((0,a)\) |
| Original full row | \((m,n,a)\) | \((L,L;L,L)\) | \((0,b)\) |
| Original full column | \((m,n-4k,b)\) | \((L,L;L,L)\) | \((0,a)\) |
| New and old Ward column, term \(j=k\) | \((m+2k,n-2k,a)\) | \((\Delta,L;L,\Delta)\) | \((0,a)\) |

Every one of the six ancillary choices in
[the rectangular engine](../tools/ward_rectangular_bounds.py)
retains each listed block. All its row and column weight slices are
nonempty. This follows directly from \(0\le2a\le\Delta\) and
\(L\le P\); in particular, \(m\ge1\) and \(n-4k\ge1\).

The residual coefficient matrix is all ones. To see this, its row
weight in the second factor is precisely the displayed residual
\(y\)-degree \(v_r\). At residual \(x\)-degree zero, the local
intersection polynomial becomes
\[
 (1+y)^{v_r-z}(1+y)^z=(1+y)^{v_r}.
\]
Its top coefficient is one, independently of the intersection
parameters. Thus each nonempty block has rank one, for every
ancillary choice. Taking the minimum of the six ancillary bounds
does not remove its multiplicity.

The old column channel count is one in the term \(j=k\), so its
minimum also retains the last block. The source coefficient indices
are nonnegative. Moreover \(a+b=2k+2a\ge2k\), so the joint-order
restriction does not remove the chosen coefficient.

### Physical dimension caps

The rectangular engine also caps its rank bound by the physical
rectangle dimension. This changes the protected multiplicity by at
most a polynomial factor. Each augmented weight space for the
four-dimensional ancillary module is a sum of four physical weight
spaces whose two grades differ from the target grades by at most
two. For valid binomial coefficients,
\[
 \binom N{r+h}\le(N+1)^2\binom Nr,
 \qquad |h|\le2.
\]
Consequently the augmented row and column dimensions are at most
\(4(2P+1)^4\) times their physical counterparts. A protected block
multiplicity cannot exceed either augmented dimension. Therefore
the physical cap preserves its exponential rate. The checker also
tests the deliberately looser factor \(16(2P+1)^4\).

The separate logical cap \(d^2\) is retained explicitly in Section 4.

## 2. A fusion entropy comparison

The constructive two-dimensional fusion recurrence in
[the stable-boundary envelope](stable-boundary-envelope.md)
gives a multiplicity lower bound with the following rate. For
\(u=m/(2P)\), \(v=n/(2P)\), and \(u+v\le1\), define
\[
 \mathcal E(u,v)=\max_{u\le t\le1-v}
 \left\{
 H(t)+\frac12\left[
 1-u-v+tH(u/t)+(1-t)H(v/(1-t))
 \right]\right\},
 \tag{3}
\]
where \(H\) is binary entropy. Integer rounding changes logarithms
by \(o(P)\). The multiplicity is at least
\(2^{2P\mathcal E(u,v)-o(P)}\).

For completeness, the underlying finite lower bound selects an even
number \(j\) of first-factor sites and uses the actual fusion paths
\[
 M_P(m,n,t_0)\ge
 \binom{2P}{j}\,2^{t_0}
 \binom{j/2-1}{(m-3)/2}
 \binom{(2P-j)/2-1}{(n-3)/2}.
\]
It counts a subset of direct summands. Maximizing its entropy over
\(j/(2P)\) gives (3).

The maximum in (3) occurs at
\[
 t=\frac{1-v}{2-u-v}.
\]
The function \(\mathcal E\) is symmetric and jointly concave: the
expression in braces is jointly concave in \(u,v,t\), and partial
maximization preserves concavity on this convex domain.

There is also a useful explicit form. Put \(d_*=1-u-v\). Then
\[
 2\mathcal E(u,v)
 =H_3(u,v,d_*)
 +(1+d_*)H\!\left(\frac{u+d_*}{1+d_*}\right)+d_*.
 \tag{4}
\]
For fixed \(u\), the last two terms form an increasing function of
\(d_*\), since its derivative is
\[
 1+\log_2\frac{1+d_*}{u+d_*}>0.
 \tag{5}
\]

Let \(\epsilon=\lambda-2\kappa>0\), and choose
\[
 x=\begin{cases}
 \epsilon/4,&\epsilon\le2/3,\\
 \epsilon/3,&\epsilon\ge2/3.
 \end{cases}
 \qquad
 a/P\to x.
 \tag{6}
\]
The row block in Section 1 has normalized weights
\[
 u=\epsilon-2x,\qquad v=1-\epsilon+x,
 \qquad d_*=x.
\]
The new Ward column has weights \((u+\kappa,v-\kappa)\); the
original full column has weights \((u,v-2\kappa)\).

Since \(\lambda\le1\), one has
\(\kappa\le(1-\epsilon)/2\). The choice (6) ensures
\(\kappa\le v-u\). Symmetry and concavity along a line of constant
\(u+v\) therefore give
\[
 \mathcal E(u+\kappa,v-\kappa)\ge\mathcal E(u,v).
 \tag{7}
\]
For the other column, \(v-2\kappa\in[x,v]\). Concavity reduces its
minimum to the two endpoints. At \(v\) the value is unchanged; at
\(x\), the variables \(v\) and \(d_*=x\) in (4) are exchanged.
The entropy \(H_3\) is unchanged, and (5), together with \(v\ge x\),
shows that the value can only increase. Hence
\[
 \mathcal E(u,v-2\kappa)\ge\mathcal E(u,v).
 \tag{8}
\]
All four protected components thus have the same lower rate
\(\mathcal E_0:=\mathcal E(u,v)\).

## 3. An explicit uniform constant

Write \(\alpha=\epsilon/2\). When \(\epsilon\le2/3\), the first
entropy lemma in the stable-boundary envelope applies with
\(x=\epsilon/4\) and \(t=7\alpha/5\), and gives
\[
 \mathcal E_0\ge\frac{29}{20}H(\alpha).
 \tag{9}
\]

For \(\epsilon\ge2/3\), choose \(t=\alpha\) in (3). The entropy
chain rule simplifies the result to
\[
 \mathcal E_0\ge
 \frac12H(\alpha)+\frac12H(4\alpha/3)+\frac43\alpha.
 \tag{10}
\]
This is at least \((29/20)H(\alpha)\) throughout
\(1/3\le\alpha\le1/2\). Here is an elementary exact verification.
Let
\[
 D(\alpha)=\frac12H(4\alpha/3)-\frac{19}{20}H(\alpha)
             +\frac43\alpha.
\]
At the first endpoint,
\[
 180D(1/3)=\log_2\frac{2^{114}3^9}{5^{50}}>0.
\]
Its derivative is
\[
 D'(\alpha)=
 \frac23\log_2\frac{3-4\alpha}{4\alpha}
 -\frac{19}{20}\log_2\frac{1-\alpha}{\alpha}
 +\frac43.
\]
For \(1/3\le\alpha\le2/5\), the first logarithm is greater than
\(-1/4\) and the second is at most one, giving
\(D'>13/60\). For \(2/5\le\alpha\le1/2\), they are bounded below
by \(-1\) and above by \(2/3\), respectively, giving
\(D'>1/30\). The logarithmic estimates follow from
\((7/8)^4>1/2\) and \((3/2)^3<4\).
This proves (9) in both ranges.

## 4. Label loss and the logical cap

The [binary charge entropy theorem](ward-code-charge-entropy.md)
gives, for every positive real \(z\),
\[
 d\le D_L\le
 z^{-(L-2k)}\frac{(1+z)^{2(P-k)}}{(1+z^2)^k}
 \le z^{-\Delta}(1+z)^{2P}.
\]
Optimizing the last expression yields
\[
 \limsup\frac{\log_2d}{P}\le2H(\alpha).
 \tag{11}
\]

Every implemented row and column minimum retains a protected rate
of at least \(\mathcal E_0\), up to a polynomial factor, until the
logical cap is applied. The selected coefficient therefore contributes
at least
\[
 \min\{d^2,\ 2^{2P\mathcal E_0-o(P)}\}
\]
to the computed cost. If \(\log d/P\) has a positive limit
\(\delta\), this gives
\[
 \liminf\rho_{\rm computed}
 \ge\min\{1,\mathcal E_0/\delta\}
 \ge\frac{29}{40},
\]
by (9) and (11). If \(\log d/P\to0\), the cap gives a lower
limit of one instead. The same reasoning on convergent subsequences
proves (1) without assuming that \(\log d/P\) converges.

## 5. Verification and next structural requirement

The exact checker
[check_ward_rectangular_envelope.py](../tools/check_ward_rectangular_envelope.py)
checks 3,000 rank-one blocks across all ancillary choices in 125
finite parameter cases, plus 500 physical-cap comparisons. It also
checks the exact integer inequalities used above. Its 8,200 floating
entropy controls illustrate the analytic argument; they are not
substitutes for that proof. The recorded output is
[ward-rectangular-envelope.json](../evidence/ward-rectangular-envelope.json).

The obstruction is now specific: the complete rectangular factors
retain large incoming stable modules even after the Ward operator
would annihilate or identify some of them. Improving this method
requires a rank calculation that keeps those constant maps or the
actual code correlations. Increasing the number of binary forms
inside the present unrestricted-rectangle computation cannot remove
the protected blocks in the fixed-ratio regime stated above.

