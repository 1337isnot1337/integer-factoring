# Exact synthesis of polynomial occupation filters

## Result and scope

The identity-mode construction now admits an exact feasibility algorithm
for arbitrary polynomial occupation filters. It asks for prescribed
upper bounds on every infinity degree, prescribed lower bounds on every
zero order, and the actual crossing endpoint. It returns explicit
polynomials or proves that those bounds are impossible.

The algorithm works in every prime characteristic. Its coefficient
systems use at most `2p+1` residue patterns. If the allowed degrees are
`O(nP)`, their size is independent of the interval length `n` and
ramification `P`. Large mode scales therefore do not require long
truncated series or a large coefficient search. The output is a finite
polynomial tensor, and its ordinary degree remains charged.

This is a constructive search interface, not a new factoring exponent.
The original exponent remains `1/2`; no strict critical-pressure witness
is claimed here.

## 1. Filters, endpoint, and the exact scalar transform

Let `p` be prime, `n>=1`, and let `P` be a positive odd integer. Put

\[
x=-v,\qquad m=nP,\qquad u=x^m.
\tag{1}
\]

After the physical substitution `v=-z^2`, one has `x=z^2` and
`u=z^{2nP}`. Use the polynomial identity-block summand `Q_0` from
[the qutrit construction](truncated-symmetric-identity-filters.md).
The [prime-state diagonalization](prime-identity-filter-diagonalization.md)
applies in characteristic `p`.

Choose nonzero polynomials

\[
\alpha_j(x)\in\mathbb F_p[x],\qquad
\alpha_j(1)=1,\qquad 0\le j<p.
\tag{2}
\]

The multiplier that replaces the `p` occupation weights of `Q_0` by
these polynomials is

\[
\mathrm{Exp}_p(-Q_0)
 \sum_{j=0}^{p-1}\frac{\alpha_j(x)}{j!}Q_0^j.
\tag{3}
\]

Here `Q_0^p=0`, and every displayed factorial is a field unit. At
`x=1`, (3) equals one. Multiplying the ramified Gaussian by (3) thus
preserves the original endpoint coefficient by coefficient. This is
stronger than an endpoint identity that holds only after selecting a
particular type.

Define the finite differences

\[
\beta_k(x)=
 \sum_{j=0}^{k}(-1)^{k-j}\binom{k}{j}
 (1-x^m)^j\alpha_{p-1-k+j}(x),
 \qquad 0\le k<p.
\tag{4}
\]

At `x=1`, `beta_k(1)=(-1)^k`, so none of the `beta_k` is zero.
Write

\[
d_j=\deg_x\alpha_j,\qquad b_k=\mathrm{ord}_{x=0}\beta_k,
\qquad a=P(1-J),\qquad \kappa=P(1-2n).
\tag{5}
\]

The active zero orders and active infinity orders are exactly

\[
\{a(r-k)+2b_k:0\le r,k<p\},
\qquad
\{\kappa j-2d_j:0\le j<p\}.
\tag{6}
\]

The Gaussian spectator orders and physical energies are unchanged.
In particular, infinity uses the *full degree* of every `alpha_j`;
the number of its nonzero coefficients does not replace this cost.
The active partition product is

\[
L_p(s^{-a})
 \left(\sum_{k=0}^{p-1}s^{ak-2b_k}\right)
 \left(\sum_{j=0}^{p-1}s^{-\kappa j+2d_j}\right),
\qquad L_p(t)=\sum_{r=0}^{p-1}t^r.
\tag{7}
\]

### Laurent and rational filters

A common Laurent monomial can be cleared before using the algorithm.
Multiplying every `alpha_j` by `x^h` adds `h` to every `d_j,b_k` and
leaves (7) unchanged. The normalization at `x=1` is also unchanged.
This covers arbitrary Laurent polynomial filters after one common
clearing.

For rational filters regular at `x=1`, choose a common polynomial
denominator `D` with `D(1)!=0` and multiply by `D(x)/D(1)`. The result
is polynomial and still satisfies (2). This clearing multiplies the
two-cut partition product by

\[
s^{2(\deg D-\mathrm{ord}_0D)}.
\tag{8}
\]

Thus a nonmonomial denominator has a real cost. The synthesis algorithm
charges that cost by working with the cleared polynomials themselves.
A pole at the crossing endpoint requires a separate cancellation and
endpoint proof; (2)--(8) do not certify such a construction.

## 2. The feasibility question is linear

Supply nonnegative integer bounds `D_j,B_k`. We seek (2) with

\[
\deg\alpha_j\le D_j,\qquad
\mathrm{ord}_0\beta_k\ge B_k.
\tag{9}
\]

Expand the unknown coefficients of each `alpha_j` up to degree `D_j`.
For every `k`, set the coefficients of `beta_k` below `B_k` equal to
zero. These are homogeneous linear equations over `F_p`. The endpoint
conditions are the `p` affine equations saying that the coefficient
sum of each `alpha_j` equals one.

This immediately gives an exact but potentially large coefficient
system. The next reduction removes its dependence on the absolute
scale `m` when the degree ratios are bounded.

It also shows that extending the coefficient field cannot improve
feasibility for fixed bounds (9). The coefficient matrix and affine
right side already lie in `F_p`; their ranks do not change under a
field extension. If a solution exists over an extension, one exists
over `F_p` itself.

## 3. Residue compression

Separate every polynomial into residues modulo `m`:

\[
\alpha_j(x)=\sum_{r=0}^{m-1}x^r A_{j,r}(x^m).
\tag{10}
\]

The transform (4) depends only on `x^m`, so the homogeneous constraints
do not mix residues. For a fixed residue `r`, their bounds become

\[
D_j(r)=\left\lfloor\frac{D_j-r}{m}\right\rfloor,
\qquad
B_k(r)=\max\left(0,
 \left\lceil\frac{B_k-r}{m}\right\rceil\right).
\tag{11}
\]

A local degree bound `-1` means that the corresponding polynomial is
zero. Let `W_r` be the solution space of the local homogeneous
coefficient equations. Its endpoint image is

\[
E_r=\{(A_{0,r}(1),\ldots,A_{p-1,r}(1)): A_{\cdot,r}\in W_r\}
 \subseteq\mathbb F_p^p.
\tag{12}
\]

Since every `x^r` equals one at the endpoint, (2) is equivalent to

\[
(1,\ldots,1)\in\sum_{r=0}^{m-1}E_r.
\tag{13}
\]

The vectors in (11) change only at the integer cut positions

\[
D_j\bmod m+1,\qquad B_k\bmod m,
\tag{14}
\]

when those positions lie strictly between zero and `m`. There are at
most `2p` such cuts, hence at most `2p+1` constant patterns. Every
residue in one pattern has the same local system and the same endpoint
image. Adding a vector space to itself changes nothing. Therefore one
representative residue for each pattern suffices in (13).

This proves the compressed feasibility criterion. The implementation
keeps a coefficient lift of each endpoint basis vector, so a successful
test also constructs the polynomials in (10). A basis of the total
endpoint image has at most `p` vectors. The returned witness can thus
be supported on at most `p` representative residues.

If every `D_j<=Km`, each local coefficient vector has at most
`p(K+1)` entries. Coefficients above the largest possible transformed
degree impose no equation and are discarded. This also prevents a
very large impossible target `B_k` from creating a long list of zero
rows. For fixed `p,K`, the number of field operations is bounded
independently of `m`; the remaining integer operations use the bit
lengths of the supplied exponents.

## 4. From feasibility to a pressure certificate

The inequalities in (9) have the correct direction for constructive
search. Replacing actual `d_j` by degree upper bounds and actual `b_k`
by order lower bounds makes both sums in (7) no smaller. Therefore

\[
L_p(s^{-a})
 \left(\sum_k s^{ak-2B_k}\right)
 \left(\sum_j s^{-\kappa j+2D_j}\right)
\tag{15}
\]

is a rigorous upper bound for the active cost of any returned witness.
The full cost also includes the exact Gaussian spectators. The
physical partition uses the unchanged physical gauges. A strict
pressure inequality proved from (15) would consequently give a
concrete polynomial crossing with both cuts and its endpoint already
specified.

Feasibility is monotone: reducing degree bounds or increasing order
bounds can only remove solutions. This permits binary searches along
individual bounds, Pareto searches, or exact branch-and-bound. Floating
point pressure estimates can guide those searches; they do not replace
the final strict inequality check.

The literal tensor constructor expands (3) for small instances. Its
cost is exponential in `n`. If `D_G` is the ordinary degree of the
unramified Gaussian, a simple degree bound for the constructed tensor
is

\[
P D_G+2nP(p-1)+\max_j D_j.
\tag{16}
\]

Here `deg Q_0<=2n`. Cancellations may lower the actual degree, which the
constructor records in its coefficient table. Every specified choice
still gives a fixed finite seed; the fast profile calculation does not
make an unbounded alphabet or a growing degree free in the factoring
transfer.

## 5. Reproducible interface and controls

The implementation is
[`polynomial_state_filter_synthesis.py`](../tools/polynomial_state_filter_synthesis.py).
Its principal functions are:

- `synthesize(p,m,degrees,orders)`: return a sparse exact witness or
  `None` when the bounds are impossible.
- `profile_data(p,n,P,J,alpha)`: compute both active profiles,
  spectators, and physical mode energies, including every infinity
  degree.
- `polynomial_seed(p,n,P,alpha)`: construct the finite raw polynomial
  tensor for a small instance.

The independent checker compares the compressed oracle with full
affine systems in characteristics `2,3,5,7`, enumerates tiny polynomial
families directly, and replays the actual polynomial tensors. Run

```bash
python tools/check_polynomial_state_filter_synthesis.py \
  --output evidence/polynomial-state-filter-synthesis.json
```

The committed evidence includes 800 full affine-system comparisons,
8,964 qutrit brute-force feasibility comparisons, independent
characteristic-five degree-one enumeration, and mode scales with more
than 200 bits. Ten complete two-cut controls include `81x81` matrices;
6,885 physical endpoint coefficient equalities are checked directly.
These controls verify synthesis and profiles. They do not assert a
pressure improvement.
