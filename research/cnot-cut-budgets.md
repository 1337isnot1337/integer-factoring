# Column-prefix budgets for polynomial crossings

## Result and scope

Consider a degree-`3t` polynomial crossing whose row Cartan degrees have
multiplicity polynomial

\[
R_t(s)=(s^3+2+s^{-3})^t.
\]
Assume its constant coefficient is nonzero. For `t=4,5,6`, its column
Smith orders at zero cannot give a strict crossing-rank improvement via

\[
\inf_{s\ge1}R_t(s)V(s)<q^{2\sqrt2},
\qquad q=2^t,\quad V(s)=\sum_{i=1}^{q^2}s^{-\nu_i}.
\]

This applies to every constant change of row bases `K(z) -> S K(z) T`
of the characteristic-two tensor-power three-CNOT path, where `S,T` are
invertible and `S K(1) T=E_q`. These matrices may mix all labels; they
need not be tensor products of smaller changes of basis. Left and right
multiplication preserve the row Cartan degrees and the polynomial degree.
The constant coefficient remains invertible in the row cut, hence nonzero
in the column cut.

The proof uses exact coefficient sums, the established exact complement
theorems, and integer polynomial certificates. It allows arbitrarily
large column orders; no numerical optimization or finite order cap enters
the proof. It excludes this zero-basepoint Smith sufficient criterion.
It does not exclude arbitrary coefficient recombinations, a different
column basepoint, or a different row profile. It does not improve the
achieved factoring exponent, which remains `1/2`.

## 1. Cutoff ranks from formal Smith factors

Let `K(z)=sum_(j=0)^L C_j z^j` be a polynomial crossing with `K(1)=E_q`.
Write `n=q^2`. Let `d_i` be its row Cartan degrees at infinity and `nu_i`
its column Smith orders at zero. For every integer `1<=k<=L`, set

\[
A_k=\sum_{j=k}^L C_j,\qquad B_k=\sum_{j=0}^{k-1} C_j,
\]

so `A_k+B_k=E_q`. Define

\[
R_k=\sum_i\max(d_i-k+1,0),\qquad
C_k=\sum_i\max(k-\nu_i,0).
\tag{1}
\]

Then

\[
r(A_k)\le R_k,\qquad c(B_k)\le C_k.
\tag{2}
\]

To prove the column statement, use a formal Smith factorization near
zero. One Smith mode has the form `z^nu a(z)b(z)^T`. Its contribution to
the first `k` coefficients, evaluated at `z=1`, is

\[
\sum_{h=0}^{k-\nu-1}
 a_h\left(\sum_{j=0}^{k-\nu-1-h}b_j\right)^T.
\]

This has rank at most `max(k-nu,0)`. Add the modes. The row statement
is the same argument applied to `w^L K(1/w)` and its first `L-k+1`
coefficients. The corresponding Smith order is `L-d_i`.

Although formal factorizations prove (2), the matrices `A_k` and `B_k`
are exact matrices over the coefficient field. Thus exact complement
bounds apply to their sum. This does not transfer an exact-only theorem
to a Laurent-series degeneration.

## 2. Exact complement input

We use the complement table from the geometric cut-rank budget theorem
in the pinned [PR 27 theorem](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md),
Section 1. It gives

\[
F_q(0)=q^2,\qquad F_q(1)=F_q(2)=q(q-1)/2.
\tag{3}
\]

For `m>=3`, a valid lower bound is the maximum of `su` over positive
integers `s,u` and `2<=A<=2m-2`, with `B=2m-A`, such that

\[
2s+Au\le2q,\qquad Bs+2u\le2q.
\tag{4}
\]

Take zero if no pair exists. These are ordinary, field-independent
complement theorems: the complement of a tensor with column rank at most
`m` has row rank at least `F_q(m)`. The cuts can be exchanged. The
geometric part follows from the integer common-subspace compression
certificates in the same PR; the one- and two-product cases use its
sharper exact small-product theorem. The replay below checks the integer
witnesses (4) used here; it does not replace those ordinary proofs.

Apply the exchanged-cut theorem to `A_k`, using any upper bound on its
row rank. Equations (1)--(2) give

\[
C_k\ge F_q\bigl(\min(n,R_k)\bigr).
\tag{5}
\]

Padding a product representation with zero terms makes the use of an
upper bound legitimate. The endpoint identity also implies that the
column polynomial matrix is generically invertible. Consequently all
`n` orders are finite, nonnegative integers, and

\[
\sum_i\nu_i\le nL.
\tag{6}
\]

A nonzero constant coefficient gives at least one zero order, so
`C_k>=k`. For `k>L`, (6) gives the further bound

\[
C_k\ge\sum_i(k-\nu_i)\ge n(k-L).
\tag{7}
\]

## 3. An exact generating identity

For `x=1/s` with `0<x<1` and every nonnegative integer `nu`,

\[
x^\nu=(1-x)^2\sum_{k=\nu+1}^{\infty}(k-\nu)x^{k-1}.
\]

Sum over all orders to obtain

\[
V(s)=(1-x)^2\sum_{k\ge1} C_k x^{k-1}.
\tag{8}
\]

All summands are nonnegative. Put

\[
g_k=\max\left(k,F_q(\min(n,R_k))\right),\qquad1\le k\le L.
\]

Using (5) on the finite prefix and (7) on the whole infinite tail yields

\[
V(s)\ge n x^L+(1-x)^2\sum_{k=1}^L g_k x^{k-1}.
\tag{9}
\]

At `s=1`, the same inequality holds by continuity: both sides are `n`.
Equation (9) therefore holds for every `s>=1` and for every admissible
integer column-order list. In particular it makes no assumption about
the largest order.

## 4. Exact certificates for the CNOT row profile

Here `L=3t` and

\[
R_t(s)=s^{-L}(1+s^3)^{2t}.
\]

Let `G_t(s)=sum_(k=1)^L g_k s^(L-k)`. Multiplying (9) by `R_t(s)` gives

\[
R_t(s)V(s)\ge
\frac{(1+s^3)^{2t}\left(ns+(s-1)^2G_t(s)\right)}{s^{2L+1}}.
\tag{10}
\]

For an integer `T`, form the integer polynomial

\[
P_{t,T}(s)=
(1+s^3)^{2t}\left(ns+(s-1)^2G_t(s)\right)-T s^{2L+1}.
\tag{11}
\]

The replay expands `P_(t,T)(1+u)` exactly, then multiplies by
`(1+u)^M`. Every resulting coefficient is strictly positive for the
following choices.

| `t` | `q` | `L` | `T` | `M` | resulting degree |
|---:|---:|---:|---:|---:|---:|
| 4 | 16 | 12 | 2600 | 11 | 48 |
| 5 | 32 | 15 | 19000 | 73 | 119 |
| 6 | 64 | 18 | 130000 | 250 | 305 |

Thus `P_(t,T)(s)>0` for all `s>=1`, proving `R_t(s)V(s)>T`.
For each row, the checker also verifies the integer inequality

\[
T^{35}>q^{99},\qquad99^2>8\cdot35^2.
\]

Therefore `T>q^(99/35)>q^(2sqrt(2))`, as required.

For audit, the nontrivial cutoff constraints at `q=64` are:

| `k` | row-tail bound `R_k` | column-prefix bound `F_q(R_k)` |
|---:|---:|---:|
| 12 | 121 | 1 |
| 13 | 42 | 24 |
| 14 | 29 | 36 |
| 15 | 16 | 70 |
| 16 | 3 | 650 |
| 17 | 2 | 2016 |
| 18 | 1 | 2016 |

The leading row mode alone already forces `C_18>=2016`. This explains
why a degree-budget target with almost all column orders near `18`
cannot work: such a target can have a very small `C_18`, even when its
total order and constant rank look favorable. The earlier cutoffs
strengthen this obstruction enough to certify the full inequality.

## 5. Replay

Run with Python's standard library:

```sh
python3 tools/check_cnot_cut_budgets.py \
  --output evidence/cnot-cut-budgets.json
```

When using this standalone checkpoint, run its adjacent checker directly.
The checker enumerates the finite integer complement witnesses, computes
the row Cartan multiplicities and the cutoff bounds, expands the exact
integer polynomials, and checks every coefficient sign. The evidence
records the source hash and a hash of each full positive coefficient
list. It contains no floating-point decision.

The proof establishes a restriction on this row profile. Constructing a
new crossing with a more favorable compatible pair of profiles remains
open. Numerical necessary-profile experiments at larger `t` are separate
from this certificate and are not asserted as constructions.
