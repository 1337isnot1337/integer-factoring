# A minimum-order tradeoff for the one-hole operator block

## Result and scope

This is an obstruction for a specific route to improving the
[occupation-transfer kernels](coupled-identity-operator-filters.md).
It applies over every fixed coefficient extension of characteristic
\(p\), for every prime-power interval length \(n=p^k\), and for
arbitrary distinct selected identity indices. It permits all
off-diagonal entries of the operator kernel.

Assume the coefficient on the fully occupied selected state is a
nonzero monomial. If the infinity minimum in the next occupation
block improves by \(2g\), with **\(g>0\)**, then the zero minimum
in that block loses at least \(nP+g\) powers of \(q\), or twice
that number in the zero uniformizer \(z\).

The conclusion is independent of constant row flags. It does not
exclude improvements in later Smith orders while both minima stay
fixed, changes to other occupation blocks, a nonmonomial top
coefficient, or additional physical state gauges.

## 1. The highest two Weyl blocks

Let the selected indices be
\(j_1<\cdots<j_r\), put

\[
 t=(p-1,\ldots,p-1),\qquad T=r(p-1),\qquad
 m=nP,\qquad c=1-q^m.
 \tag{1}
\]

The [prime-power reflected change](prime-power-joint-identity-filters.md)
gives the scalar contraction matrix \(C=cI\). Write
\(d=F_{tt}(q)\), and let \(B(q)\) be the submatrix of \(F\)
whose row and column labels are the one-hole states

\[
 b_i=t-e_i,\qquad 1\le i\le r.
 \tag{2}
\]

In the Weyl representation, the occupation-\(T\) block is \(d\).
The occupation-\(T-1\) block is exactly

\[
 dI-c^{-1}B.
 \tag{3}
\]

Here the row basis is the monomials \(Y^{b_i}\). Multiplying (3)
by the zero unit \(c\) gives \(cdI-B\).

To check the normalization, the universal Weyl formula sends a
matrix unit \(E_{ab}\) to

\[
 (-1)^r a!\,\partial^{t-a}(-c^{-1}Y)^{t-b}.
 \tag{4}
\]

On a one-hole input monomial, multiplication in (4) vanishes when
\(|t-b|\ge2\), since that input has only one available occupation.
The top matrix unit acts as the identity by
\((p-1)!=-1\) in the prime field. For
\(a=t-e_i,b=t-e_j\), one has
\(a!=(-1)^{r-1}\), using \((p-2)!=1\). The remaining derivative
has coefficient \(p-1=-1\), so the resulting matrix unit is
\(-c^{-1}E_{ij}\). This proves (3), including odd characteristic.

At zero, the weighted block is therefore regularly equivalent to

\[
 \operatorname{diag}(q^{w_i})(cdI-B),
 \qquad w_i=P\sum_s j_s(b_i)_s.
 \tag{5}
\]

At infinity it is, up to constant factorial units,

\[
 B(z^{-2})\operatorname{diag}(z^{\kappa_i}),
 \qquad
 \kappa_i=P\sum_s(2j_s+1-2n)(b_i)_s.
 \tag{6}
\]

Since every state in (2) has occupation \(T-1\), the orderings of
\(w_i\) and \(\kappa_i\) agree. Their common smallest label is
\(b_*=t-e_r\), which removes one occupation from the largest
selected index. Write \(w_*=w_r\) and \(\kappa_*=\kappa_r\).

For the Gaussian kernel \(F=I\), equation (5) has diagonal
\((c-1)q^{w_i}=-q^{m+w_i}\). Its zero minimum is therefore
\(w_*+m\), measured in powers of \(q\). The infinity minimum
is \(\kappa_*\), measured in powers of \(z\).

## 2. The tradeoff

Suppose \(d=a q^h\), with \(a\ne0\). Multiplying the whole
operator kernel by \(a^{-1}q^{-h}\) makes \(d=1\). This common
scalar shifts the zero and infinity profiles oppositely, so it
does not change their product cost. The filtered endpoint remains
invertible and has the same regular physical repair.

Let \(\mu_0\) denote the smallest \(q\)-order of (5), and let
\(\mu_\infty\) be the smallest \(z\)-order of (6), after this
normalization. For every positive integer \(g\),

\[
 \boxed{\quad
 \mu_\infty\ge\kappa_*+2g
 \quad\Longrightarrow\quad
 \mu_0\le w_*-g.
 \quad}
 \tag{7}
\]

Indeed, the minimum Smith order of a matrix is the smallest order
of any entry. Its \((* ,*)\) entry in (6) consequently gives

\[
 \deg_q B_{**}\le-g.
 \tag{8}
\]

This diagonal polynomial is nonzero. The selected indices are
distinct, so the one-hole endpoint energies are distinct. A
filtration-preserving invertible endpoint is triangular in those
energies and has every diagonal entry nonzero. In particular,
\(B_{**}(1)\ne0\).

Since \(g>0\), equations (8) and
\(\operatorname{val}_q c=0\) imply

\[
 \operatorname{val}_q(c-B_{**})
 =\operatorname{val}_q B_{**}\le-g.
\]

The corresponding diagonal entry of (5) has order at most
\(w_*-g\), proving (7). Its loss from the Gaussian minimum is
at least

\[
 (w_*+m)-(w_*-g)=m+g.
 \tag{9}
\]

The strict inequality \(g>0\) is necessary. At \(g=0\), the
Gaussian diagonal \(B_{**}=1\) cancels the constant term of
\(c\), and its zero order is \(m\), not zero.

For a complete active zero profile at aggregate charge \(J=1\),
use the fully occupied auxiliary state in the Weyl reconstruction.
Equation (7) then supplies an active zero order at most

\[
 -2Pj_r-2g,
 \tag{10}
\]

whereas the corresponding Gaussian order is
\(2P(n-j_r)\). An aggregate charge change gives the ordinary
shift for charge \(-1\); its physical type cost still has to be
included.

## 3. Consequences for search design

The theorem rules out any proposed row-lattice target that improves
this infinity minimum while retaining the zero block within
\(m+g-1\) powers of its Gaussian minimum. Enumerating additional
constant row flags or enlarging the fixed coefficient field cannot
rescue such a target.

This explained a failed three-direction target for
\(n=32\), \(P=1\), and selected indices \((0,30,31)\).
That target tried to improve the one-hole infinity minimum while
preserving its zero minimum. A larger constant-flag screen was
stopped once (7) supplied the uniform obstruction.

A different target is to change a deeper occupation block while
keeping the top and one-hole blocks Gaussian. For example, four
binary modes with an additive relation between their selected
indices have equal-energy two-hole states. Their endpoint admits
an invertible block within that energy level, and (7) places no
restriction on a perturbation confined to that deeper block.

## 4. Independent exact controls

The [checker](../tools/check_identity_operator_minimum_tradeoff.py)
applies the universal differential operator directly to the top
and one-hole monomials in characteristics two, three, five, seven,
and eleven. It checks 550 matrix-unit actions and 2,281 vanishing
actions from lower occupation sectors. Twelve additional complete
polynomial Weyl calculations verify (7) for positive gains and
several selected geometries and ramifications.

The [recorded evidence](../evidence/identity-operator-minimum-tradeoff.json)
is reproduced by:

```sh
python3 tools/check_identity_operator_minimum_tradeoff.py \
  --output /tmp/operator-minimum-tradeoff-replay.json
```
