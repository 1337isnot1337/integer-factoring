# Spherical rank formula for the full Catalan crossing

## Result and limits

For the full spin-zero path family of length `2r`, the two cut ranks of
every six-vertex coefficient tensor, and of every linear recombination of
those coefficients, have an exact characteristic-zero formula in terms of
Jack spherical polynomials. This replaces dense rank computation by a
polynomial-block support calculation.

The formula is a structural tool, not a post-half bound. For the natural
`(r+1)(r+2)/2` coefficient pieces, exact rational calculations give total
rank cost `1608.377200168389` already at `r=4` (14 full-path labels), and
`1034173106197.841373251718` at `r=12` (208,012 labels). These pieces do
not meet the target. The formula does not settle arbitrary decompositions
of the coefficient span, power-of-two restricted label sets, or a fixed
characteristic-five transfer for unbounded heights. The score remains
`1/2`.

## 1. Multiplicity-free block decomposition

Let `V=C^2` and `W=V tensor V`, with the two binary registers as the two
factors of `W`. The full Dyck-path encoder space of length `2r` is the
Specht module `[r,r]` for simultaneous position permutations. Therefore
the path-pair space is `[r,r] tensor [r,r]`.

Schur-Weyl duality on `W^(tensor 2r)` writes this pair space as a sum of
Specht blocks `[mu]`. Manivel's Theorem 1 gives the multiplicity-free
decomposition

```math
[r,r]\otimes[r,r]
=\bigoplus_{\substack{\mu\vdash2r,\;\ell(\mu)\le4\\
                         \text{all parts even}}}[\mu]
 \;\oplus\!
  \bigoplus_{\substack{\mu\vdash2r,\;\ell(\mu)=4\\
                         \text{all parts odd}}}[\mu].
\tag{1}
```

Each block occurs once. Its multiplicity space in `S_mu(W)` is a
one-dimensional invariant line for `SL_2 x SL_2`, the spin cover of
`SO_4`. Thus a position-symmetric local operator acts by one scalar on
each `[mu]` block. The block dimension is the ordinary Specht dimension
`f^mu`, given by the hook formula. In particular,
`sum_mu f^mu = Catalan(r)^2`.

The decomposition is Manivel, [*A note on certain Kronecker
coefficients*](https://arxiv.org/abs/0809.3710), Theorem 1. The
`SL_2 x SL_2` / `SO_4` identification is also used in the repository's
[four-qubit invariant argument](four-qubit-orthogonal-generators.md).

## 2. Explicit spherical eigenvalues

Use the local matrix in basis `00,01,10,11`

```math
K(t,w)=\operatorname{diag}\!\left(1,
  \begin{pmatrix}t&w\\1&t\end{pmatrix},1\right),
\qquad s=t^2.
```

The invariant symmetric form on `W` is `J=epsilon tensor epsilon`, with
`epsilon=[[0,1],[-1,0]]`. A direct calculation gives

```math
\operatorname{spec}(J^{-1}K^TJK)
=\{1,1,(t+\sqrt w)^2,(t-\sqrt w)^2\},
\quad x+y=2(s+w),\quad xy=(s-w)^2,
\quad \det K=s-w.
\tag{2}
```

Write `P_nu^(2)` for the monic four-variable Jack polynomial with parameter
`2`. Its normalized zonal spherical function is

```math
\Phi_{2\nu}(s,w)=
\frac{P_\nu^{(2)}(1,1,x,y)}{P_\nu^{(2)}(1,1,1,1)},
\qquad |\nu|=r,                                      \tag{3}
```

for even blocks `mu=2 nu`. For odd blocks `mu=2 nu+(1,1,1,1)`, where
`|nu|=r-2`, the determinant twist gives

```math
\Phi_{2\nu+(1^4)}(s,w)=(s-w)\Phi_{2\nu}(s,w).       \tag{4}
```

The Jack/zonal identification is the `alpha=2` case of the classical
orthogonal spherical-polynomial theory; see [Jack's original
paper](https://doi.org/10.1017/S0080454100008517) and the account of Jack
polynomials as spherical polynomials for the orthogonal pair in [Sahi,
*Jack polynomials and the binomial formula*](https://ow3.math.rutgers.edu/~sahi/Reprints/11imrn-jack-lrbino.pdf).
Equations (3)--(4) are ordinary polynomial identities, so the square roots
in (2) disappear after writing the answer in `s,w`.

The normalized coefficient tensor `T_(A,B,C)` is the coefficient of
`s^B w^C` in the local operator family, with `A+B+C=r`. Schur-Weyl and
(1)--(4) give, over `Q`,

```math
\begin{aligned}
\operatorname{rank}_{row} T_{A,B,C}
 &=\sum_{\mu:\,[s^B w^C]\Phi_\mu\ne0} f^\mu,\\
\operatorname{rank}_{col} T_{A,B,C}
 &=\sum_{\mu:\,[s^B w^A]\Phi_\mu\ne0} f^\mu.
\end{aligned}                                         \tag{5}
```

The second line is the existing cut reciprocity
`R(T_(A,B,C))=C(T_(C,B,A))`. More generally, for any coefficient
recombination `T= sum lambda_(A,B,C) T_(A,B,C)`, replace the coefficient
test in (5) by the corresponding linear functional on `Phi_mu`; this
gives exact ordinary ranks for arbitrary pieces in the full coefficient
span. It does not assert that the resulting support optimization has a
cheap solution.

## 3. Exact controls

The standard-library part of
[`check_catalan_spherical_ranks.py`](../tools/check_catalan_spherical_ranks.py)
constructs the Jack polynomials by their triangular differential-operator
eigenproblem, evaluates (3)--(4) in rational `s,w` coordinates, and
computes all block dimensions and coefficient rank pairs through
`r=12`. Its independent word-basis control builds every eight-site full
Catalan coefficient matrix and checks both ranks modulo `101`; all 15
rank pairs agree with (5), and their modular rank cost is
`1608.377200168389`. The output is pinned in
[`catalan-spherical-ranks.json`](../evidence/catalan-spherical-ranks.json).

Run `make verify-catalan-spherical-ranks` (NumPy is used only for the
eight-site matrix control). This is a finite replay of the formula's
small-width consequences, not independent mathematical review or a
uniform factoring improvement.

## 4. Next question

The concrete remaining problem is to exploit (5) for arbitrary
recombinations: characterize the paired row/column zero sets of the
spherical-polynomial linear functionals, then find a uniform decomposition
of the endpoint evaluation with total weighted support below the target.
The natural monomial pieces are already too expensive at tiny widths, so a
useful next step must use cancellations across coefficients. The height
three modular transfer still needs a separate argument because unrestricted
path norms can acquire factors divisible by five.

This derivation and its exact checker used substantial OpenAI assistant
support. No new Lean formalization was added.
