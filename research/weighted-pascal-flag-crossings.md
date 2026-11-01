# Weighted Pascal flags with two different energy spectra

## Constructive extension

The [tensor-product Pascal box](tensor-pascal-box-crossings.md) extends
to arbitrary positive integer coordinate weights. The two cuts then
see different energies: the infinity spectrum uses each weight itself,
whereas the zero spectrum uses its lowest nonzero binary digit.
Both complete profiles remain explicit.

The weighted kernel also has a small, explicit Clifford matrix at
each cut. This gives an effective evaluator for arbitrary independent
physical mode heights, including heights outside the regularity cone
used for the displayed zero spectrum. The evaluator retains both
complete profiles and both independent matching-type costs.

These are exact constructor and profile results. The numerical gauge
searches require independent exact certificates before any bound could
be claimed. No seed below the half-score threshold is supplied here;
the complete factoring score remains `1/2`.

## 1. An arbitrary positive weight vector

Let `beta_1,...,beta_r` be positive integers. Choose a power of two `L`
with `sum beta_j<L`, and put

\[
 \omega_j=2^{v_2(\beta_j)},\qquad
 \kappa=L-\sum_j\beta_j>0,\qquad
 \kappa_0=L-\sum_j\omega_j>0.
\]

The mode set is `{0,1}^r`, of size `n=2^r`. Define

\[
 w(a)=\sum_j\beta_ja_j,\quad
 w_0(a)=\sum_j\omega_ja_j,\quad
 E_\infty(a)=\kappa+2w(a),\quad
 E_0(a)=\kappa_0+2w_0(a).                          \tag{1}
\]

Both energy lists have mean `L`, and complementation of all mode
bits sends either energy to `2L` minus itself. Let

\[
 P_{ab}=[b\subseteq a],\qquad c=1+u,
\]

and set

\[
 \begin{aligned}
 A&=\mathrm{diag}(c^{E_\infty(a)})P^{\mathsf T},\\
 B&=P^{\mathsf T}\mathrm{diag}(c^{2L-E_\infty(a)})P,\\
 \mathcal D&=\begin{pmatrix}A&u^{2L}P\\I&B\end{pmatrix}.
 \end{aligned}                                                    \tag{2}
\]

These polynomial matrices satisfy

\[
 AB=c^{2L}P,\qquad BPA=c^{2L}I,\qquad
 \det\mathcal D=1.                                \tag{3}
\]

Their endpoint is the same triangular crossing as in the Pascal-box
construction. All weights are positive, so every strict Pascal
inclusion raises `w(a)`. Equal physical word types remove those
strict inclusions at the endpoint. Independent matching types on
the pairs `(B,Y)` and `(X,V)` are also valid: the identity block
first fixes `(B,Y)`, and the weighted-index argument then applies
only to `(X,V)`.

The lower-right block has the useful entry formula

\[
 B_{ab}=c^\kappa
       \prod_{j:a_j=b_j=0}(c^{2\beta_j}+1).         \tag{4}
\]

This follows by multiplying the two-by-two factors
`P_2^T diag(c^(2beta_j),1)P_2`. Formula (4) is generally not a
monomial. Its exact zero valuation is

\[
 v_u(B_{ab})=2\sum_{j:a_j=b_j=0}\omega_j,            \tag{5}
\]

because `v_u((1+u)^(2beta_j)+1)=2omega_j` in characteristic two.

## 2. Explicit spectra at the natural gauges

In boundary order `(B,X,V,Y)`, use

\[
 h_0(a)=(0,E_0(a)-2L,0,E_0(a)),
\]

\[
 h_\infty(a)=(E_\infty(a)-C,2L-C,E_\infty(a)-C,-C).
                                                               \tag{6}
\]

The zero internal block `A` is regular. Equations (1), (5), and
(6) give the weighted valuation of a `B` entry as

\[
 E_0(a)+E_0(b)-2L
       +2\sum_{j:a_j=b_j=0}\omega_j
 =2\sum_j\omega_ja_jb_j\ge0.                       \tag{7}
\]

Consequently the zero row cut has two copies of `E_0(a)` as its
one-particle order list.

At infinity, let `t=u^(-2)` and write

\[
 R_j(t)=(1+t)^{\beta_j}+t^{\beta_j},\qquad
 D_{aa}=t^{w(a)},\qquad Q_{aa}=(1+t)^{w(a)}.
\]

Define the tensor-product matrix

\[
 O_{ab}=t^{\sum_j\beta_j(a_j\mathbin{|}b_j)}
           \prod_{j:a_j=b_j=0}R_j(t).               \tag{8}
\]

The two binary indices of `O` are `a,b`. Each one-bit factor is

\[
 O_j=\begin{pmatrix}R_j(t)&t^{\beta_j}\\t^{\beta_j}&t^{\beta_j}\end{pmatrix}.
\]

Direct multiplication gives

\[
 Q_jP_2^{\mathsf T}O_j=(1+t)^{\beta_j}D_jP_2.
\]

Multiplication over all coordinates therefore yields

\[
 QP^{\mathsf T}O=(1+t)^{\sum_j\beta_j}DP.            \tag{9}
\]

After removing the common factor `xi^(kappa-2C)`, where `xi=1/u`,
and the same regular row and column units as before, the infinity
matrix becomes

\[
 \begin{pmatrix}QP^{\mathsf T}D&DP\\D&(1+t)^\kappa O\end{pmatrix}.
\]

Each row `a` of `O` is divisible by `t^(w(a))`. The regular Schur
reduction, followed by (9), gives diagonal blocks

\[
 D,\qquad
 \bigl(1+(1+t)^{\kappa+\sum \beta_j}\bigr)DP=t^LDP.
\]

Thus the complete infinity one-particle spectrum is

\[
 E_\infty(a)-2C,\qquad E_\infty(a)+2L-2C.            \tag{10}
\]

At `C=L`, it is the signed list `-E_infty(a),E_infty(a)`.
The total restored physical energy and product of the two exterior
profile sums are

\[
 \begin{aligned}
 G(a)&=2(E_0(a)+E_\infty(a)-2L),\\
 \mathcal B(s)&=s^{nL}
        \prod_a(1+s^{-E_0(a)})^2(1+s^{-E_\infty(a)})^2.
 \end{aligned}                                                     \tag{11}
\]

At these gauges, the two independent matching-pair energy lists are
equal to `G/2`. Hence their paired-type pressure agrees with the
common-type pressure `B(s)<prod_a(1+s^(-G(a)/(4rho)))^(4rho)`.
Other gauges need not have equal matching-pair energies.

As an example outside the equal-spectrum box family, weights `(1,3)`
with `L=8` give

\[
 E_0=(6,8,8,10),\qquad E_\infty=(4,6,10,12).
\]

No specialization in another variable is involved in these spectra.

## 3. Small matrices for unrestricted physical mode gauges

The [Clifford Smith theorem](clifford-gaussian-smith-profiles.md)
applies even when internal Gaussian coefficients have poles after
physical gauging. The following explicit matrices make that theorem
effective for (2).

Put `F=u^(2L)`. At the zero cut, the left variables are `(B,V)`
and the right variables are `(X,Y)`. In blocks of order `n`, the
Clifford matrix is

\[
 \mathcal C_0=
 \begin{pmatrix}
 0&F^{-1}I&F^{-1}(PA)^{\mathsf T}&0\\
 P&0&0&A\\
 B&0&0&I\\
 0&F^{-1}(BP)^{\mathsf T}&F^{-1}P^{\mathsf T}&0
 \end{pmatrix}.                                      \tag{12}
\]

For verification, the cross matrix is
`K=[[0,I],[FP,0]]`; the internal alternating matrices are
`[[0,A^T],[A,0]]` and `[[0,B^T],[B,0]]`. Substitute these into the
general Clifford block formula and use (3). In particular, its
upper-left block simplifies using both identities in (3), not merely
the determinant of the full kernel.

At infinity, put all minus variables `(B,X)` on the left and all
plus variables `(V,Y)` on the right. There are no internal blocks,
and

\[
 \mathcal C_\infty=
       \mathrm{diag}(\mathcal D^{\mathsf T},\mathcal D^{-1}),
 \qquad
 \mathcal D^{-1}=\begin{pmatrix}BP&FI\\P&PA\end{pmatrix}.
                                                               \tag{13}
\]

Use `u=xi^(-1)` for its infinity valuations. Both Clifford matrices
have order `4n`; the physical Fock flattenings have order `2^(2n)`.

The two matrices also remain valid after replacing `A` by `A C^(-1)`
and `B` by `CB`, for any invertible constant flag `C`: both identities
in (3) are preserved. The endpoint is unchanged. The displayed
natural spectra in Section 2 need not survive this replacement;
the Clifford evaluator must recompute both profiles. The checker
includes such changed flags under gauges with internal poles.

For any physical mode heights `h` at a place, let `l,r` be their
restrictions to that cut's left and right variables. Multiply its
Clifford matrix by `diag(u^l,u^(-l))` on the left and by
`diag(u^r,u^(-r))` on the right, using the local valuation variable
at infinity. The complete Fock means are

\[
 \bar\nu_0=Ln+\frac12\sum h_0,\qquad
 \bar\nu_\infty=\frac12\sum h_\infty.               \tag{14}
\]

The zero contribution `Ln` is half the valuation of `det K=F^n`.
The infinity cross determinant is one. Omitting `Ln` would produce
a false gain.

If the signed Clifford orders are `beta_i`, the full cut partition
sum is obtained from

\[
 \log Z=-\bar\nu+\frac12\sum_i
              \log(e^{\beta_i/2}+e^{-\beta_i/2}).    \tag{15}
\]

For scaled heights and orders, multiply ordinary valuations by
`log(s)` in (14)--(15). The two physical matching-pair energies are
the sums of heights on `(B,Y)` and `(X,V)`, independently. Retaining
both partition functions is necessary when optimizing these gauges.

## 4. Scope and computation

The exact constructor is [weighted_pascal_flags.py](../tools/weighted_pascal_flags.py).
The finite checker reconstructs both polynomial identities, the
complete natural spectra, and the Clifford matrices. It compares
their predicted full spectra with original Fock tensors under
integer gauges, including gauges with internal poles.

Run

```sh
python3 tools/check_weighted_pascal_flags.py \
  --output evidence/weighted-pascal-flags.json
```

The optional
[unrestricted gauge search](../tools/search_weighted_pascal_flags.py)
uses the existing native matrix-profile engine and independent
matching types. It checks a local gradient against finite differences,
retains every completed run, and records any unresolved profile query.
For example:

```sh
python3 tools/search_weighted_pascal_flags.py \
  --weights 1 2 4 7 --outer 16 --restarts 6 --seconds 1200 \
  --output evidence/weighted-pascal-search-1-2-4-7.json
```

Use `--flag-shears` to generate a reproducible constant binary flag.
That option changes the actual polynomial kernel and recomputes its
full profiles; it does not reuse the natural energy-list formulas.

The arbitrary-weight family includes the original interval when
the weights are consecutive binary powers and `L=2^r`. Repeated
binary powers give the summed-coordinate boxes. Nondyadic weights
give two genuinely different cut energy lists.

The exact formulas supply a larger constructive search domain. They
do not imply that a sampled floating optimum is global, nor that
any displayed seed improves the asymptotic factoring score. The
derivation and tests used substantial OpenAI assistant support;
no new Lean formalization or external review is claimed.
