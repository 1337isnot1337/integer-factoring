# Exact qutrit state filters on an interval Smith direction

## Result and scope

There is an explicit polynomial crossing family obtained by filtering the
three occupation states of one infinity Smith mode. Its complete profiles
reduce to nine zero orders, three infinity orders, and Gaussian spectator
modes. Three scalar valuations determine all nine zero orders.

For monomial state filters, the new rank-cost product is at least the
Gaussian product at the same gauges and scale whenever

\[
1-2n\le J\le2n.
\tag{1}
\]

The inequality is strict for every nonidentity filter. It holds at every
positive evaluation scale, for every allowed ramification, and for
unbounded integer filter exponents. It compares with the *charged*
Gaussian background; it does not assume an unproved bound for every
Gaussian gauge. At `J=0`, the existing interval frontier supplies the
absolute half-rate exclusion.

The construction and exact scalar reduction remain useful for polynomial
filters beyond monomials. The note also records why an earlier rational
positive-energy direction does not give an admissible long-interval
crossing.

## 1. The actual polynomial tensor

Work over `F3` in the truncated symmetric algebra with every physical
mode cubed equal to zero. Let `S={0,...,n-1}`, take ambient `M=2n`, and
write the all-field one-particle matrix as

\[
D_n(w)=\begin{pmatrix}A&B\\I&C\end{pmatrix}.
\]

Its row boundaries are `(V,Y)` and column boundaries `(U,X)`. The
identity block gives the polynomial rank-one summand

\[
Q_0(w)=\bigl(Y_0+(1+w)V_0\bigr)
 \left(U_0+\sum_{b=0}^{n-1}(-1)^b w^{2n-1-b}X_b\right).
\tag{2}
\]

The regular infinity reduction isolates this summand as the mode of
order `1-2n`. The derivation and both regular zero basis changes are in
[the identity-mode reduction](qutrit-identity-zero-reduction.md).

Let `P` be a positive odd integer and let `e_1,e_2` be arbitrary integers.
Set

\[
\alpha_0(v)=1,\qquad
\alpha_1(v)=(-v)^{-e_1},\qquad
\alpha_2(v)=(-v)^{-e_2},
\]

\[
a(v)=\alpha_1(v)-1,
\qquad b(v)=2\alpha_2(v)-\alpha_1(v)+2.
\tag{3}
\]

For the truncated Gaussian kernel `G_n`, define the Laurent tensor

\[
T(v)=G_n(v^P)
 \left(1+a(v)Q_0(v^P)+b(v)Q_0(v^P)^2\right).
\tag{4}
\]

The factors in (2) are linear forms whose cubes vanish, so `Q_0^3=0`.
Since `Exp_3(Q)=1+Q+2Q^2`, multiplication by (3) replaces its normalized
occupation coefficients by `1,alpha_1,alpha_2` exactly.

Every `alpha_i(-1)=1`; hence `a(-1)=b(-1)=0` and
`T(-1)=G_n(-1)` coefficient by coefficient. All matching typed tensor
powers therefore retain the original crossing. The only poles of (4)
are monomial poles at zero. Put `q=max(0,e_1,e_2)`. The tensor

\[
(-1)^qv^qT(v)
\tag{5}
\]

is polynomial and has the normalized endpoint. After `v=-z^2`, this
clearing shifts every zero order by `+2q` and every infinity order by
`-2q`, so the two partition sums' product is unchanged.

This is a fixed finite seed for every specified `n,P,e_1,e_2`. Its
ordinary polynomial degree and the full coefficient table remain part
of the construction cost. A fast profile calculation does not make a
growing alphabet or degree free in the factoring transfer.

## 2. Complete profiles with an aggregate charge

Use the natural physical gauges, multiplied by `P`. At zero, add `PJ`
to every `V` mode and subtract `PJ` from every `U` mode. This convention
is equivalent, on these balanced tensors, to adding `PJ` on `X` and
subtracting it on `Y`. The paired physical mode energies are

\[
g_i^X=P(4i-2n+2+J),\qquad
g_i^Y=P(4i-2n+2-J).
\tag{6}
\]

Put

\[
A=P(1-J),\quad B=2nP,\quad u=z^B,\quad c=1-u,
\quad \kappa=P(1-2n).
\tag{7}
\]

The zero reduction has background `Exp_3(z^A YU+z^{B-A}VX)` and

\[
Q=(z^A Y+(1-z^B)V)(U-z^{-A}X).
\tag{8}
\]

All other zero modes are spectators, of orders

\[
P(2i+1-J),\quad 1\le i<n,
\qquad
P(2i+1+J),\quad 0\le i<n-1.
\tag{9}
\]

For arbitrary Laurent functions `alpha_i`, with `alpha_0=1`, define

\[
S=\alpha_2,\qquad
D=c\alpha_2-\alpha_1,\qquad
F=1+c\alpha_1+c^2\alpha_2.
\tag{10}
\]

Let their zero valuations be `s_0,d_0,f_0`. The complete nine active
zero orders are the sorted list

\[
\boxed{
s_0+2A,\ s_0+A,\ d_0+A,\ s_0,\ d_0,\ f_0,
d_0-A,\ f_0-A,\ f_0-2A.}
\tag{11}
\]

The corresponding partition factor is particularly simple:

\[
Z_{0,\mathrm{active}}(s)
=L_3(s^{-A})
 \bigl(s^{-s_0}+s^{A-d_0}+s^{2A-f_0}\bigr),
\qquad L_3(x)=1+x+x^2.
\tag{12}
\]

The [prime-state finite-difference diagonalization](prime-identity-filter-diagonalization.md)
proves the general version.
For qutrits, the charge blocks of sizes `1,2,3,2,1` are regularly
equivalent to `(S)`, `(S,D)`, `(S,D,F)`, `(D,F)`, `(F)`, respectively.
Their scalar charge shifts give (11). The checker below verifies
independent exact polynomial factorizations of all five blocks.

At infinity, the active orders are

\[
0,\quad \kappa+2e_1,\quad 2\kappa+2e_2,
\tag{13}
\]

and the spectator modes have orders
`P(2i+1-2n)`, `1<=i<2n`. Each spectator contributes its usual factor
`L_3(s^{-order})`. Thus (9), (11), and (13) give both complete profiles.

### Exact monomial specialization

For the filters in (3), the three valuations are

\[
s_0=-2e_2,
\qquad
d_0=\begin{cases}
B-2e_1,&e_1=e_2,\\
-2\max(e_1,e_2),&e_1\ne e_2,
\end{cases}
\tag{14}
\]

\[
f_0=\begin{cases}
2B,&e_1=e_2=0,\\
-2\max(0,e_1,e_2),&\text{otherwise}.
\end{cases}
\tag{15}
\]

These formulas include the actual cancellations. In particular,
`alpha_1=alpha_2` is not a generic limit of unequal exponents. For the
identity filter, `D=-u`, `F=u^2`, and the background modes `A,B-A` are
recovered exactly.

## 3. Pointwise domination for every monomial filter in (1)

All spectator factors and `L_3(s^{-A})` cancel when comparing with the
Gaussian at the same gauges. For a real evaluation `s>1`, write

\[
Z=s^A,\quad U=s^{-B},\quad H=s^{-\kappa},\quad
R=s^{2e_1},\quad S_*=s^{2e_2}.
\tag{16}
\]

The background product remaining after cancellation is

\[
\mathcal B_G=(1+UZ+U^2Z^2)(1+H+H^2).
\tag{17}
\]

Here `H>1`, `H^{-2}<=U<H^{-1}`, and condition (1) implies

\[
H^{-1}\le Z\le U^{-1}\le H^2.
\tag{18}
\]

### Unequal exponents

If `e_1!=e_2`, (14)--(15) give the new product

\[
\mathcal B_T=
\bigl(S_*+Z\max(R,S_*)+Z^2\max(1,R,S_*)\bigr)
\left(1+\frac H R+\frac{H^2}{S_*}\right).
\tag{19}
\]

Cauchy's inequality yields

\[
\mathcal B_T\ge
(S_*+ZR+Z^2)
\left(\frac{H^2}{S_*}+\frac H R+1\right)
\ge (H+\sqrt{ZH}+Z)^2.
\tag{20}
\]

Put `y=sqrt(Z/H)`. By (18), `1/H<=y<=sqrt(H)`. The function

\[
\frac{(1+y+y^2)^2}{1+y^2+y^4}
=\frac{1+y+y^2}{1-y+y^2}
\]

is invariant under `y -> 1/y` and increases on `(0,1]`. Consequently

\[
\frac{(1+y+y^2)^2}{1+y^2+y^4}
\ge\frac{H^2+H+1}{H^2-H+1}
>\frac{H^2+H+1}{H^2}.
\tag{21}
\]

Since `U<=1/H`, equations (20)--(21) imply
`\mathcal B_T>\mathcal B_G`.

### Equal nonzero exponents

Now `e_1=e_2!=0`; put `R=S_*`, `C=1+ZU`, and `K=H+H^2`. The new
product is

\[
\mathcal B_T=
\bigl(RC+Z^2\max(1,R)\bigr)(1+K/R).
\tag{22}
\]

If `R>=1`, this is `(C+Z^2)(R+K)`, which is at least
`(C+Z^2)(1+K)>\mathcal B_G`.

If `R<1`, expand (22) as
`RC+CK+Z^2+Z^2K/R`. Its derivative in `R` is negative on `(0,1]`
because `Z^2K>C`. To verify the latter inequality, the function
`Z^2(H+H^2)-1-ZU` increases for `Z>=1/H`, and at `Z=1/H` it is
`(1-U)/H>0`. The infimum over `R<=1` is therefore its value at one,
again strictly above (17).

The remaining case `e_1=e_2=0` is exactly the Gaussian. This completes
the all-scale, unbounded-exponent comparison under (1). Because the
physical partitions in (6) are unchanged, the comparison applies to
both common physical types and independently optimized paired types.

This theorem leaves polynomial filters beyond (3), other virtual
directions, and charges outside (1) as separate questions. In particular,
it does not turn a construction-specific comparison into a factoring
lower bound.

## 4. Two checks that reject the unrepaired positive-one direction

A different infinity pivot of order `+1` gives a rational rank-one
summand entirely in `Y,X`. Its `Y_0X_0` coefficient is

\[
\frac{w^{2n-1}}{(1+w)^n}.
\tag{23}
\]

Indeed, the `(0,0)` entry after eliminating `A` is

\[
w^{2n-1}-w^{2n}\sum_{a=0}^{n-1}(1+w)^{-a-1}
=\frac{w^{2n-1}}{(1+w)^n}.
\]

The proposed coefficient `f=(1+w)/w` on its square has a genuine pole
on a retained one-copy diagonal entry. Take every physical letter to
occupy only mode zero, twice. The added coefficient is

\[
\frac{w^{4n-3}}{2(1+w)^{2n-3}}.
\tag{24}
\]

The remaining Gaussian coefficient is `V_0^2U_0^2` with coefficient
`(1+w)^2/2`; no other perturbation monomial can reach this target.
For `n>=2`, (24) has a pole of order `2n-3`. A matching-type restriction
cannot remove a pole already present on this one-copy diagonal.

There is also a separate zero-profile cost. The unweighted positive
zero order in (23) is canceled by the natural `Y_0,X_0` physical gauges.
After `w=-z^2`, the added `Y_0^2X_0^2` coefficient has weighted order

\[
2(4n-3)+2(1+1-4n)=-2.
\tag{25}
\]

Thus its multiplication map is not a regular zero unit. Ignoring those
gauges incorrectly predicts an unchanged zero profile. Equations
(24)--(25) exclude this unrepaired candidate, without excluding all
rational interior-mode constructions or possible repairs.

## 5. Exact implementation and checks

The constructor and fixed-size profile interface are in
[`truncated_symmetric_identity_filters.py`](../tools/truncated_symmetric_identity_filters.py).
The direct constructor supplies the actual polynomial tensor, including
its common monomial clearing. The profile interface has both a closed
formula and an independent multivariate-minor implementation. With
`include_spectators=False`, its active calculation uses only a fixed
number of integer operations, even for very large `n` and `P`.

Run

```sh
python tools/check_truncated_symmetric_identity_filters.py \
  --output evidence/truncated-symmetric-identity-filters.json
```

The checker verifies the five polynomial block factorizations, compares
the closed orders with every minor in the fixed blocks, reconstructs
full one- and two-mode tensors, and checks their local endpoints and both
full weighted cuts. It also checks rational instances of the domination
inequality and the two failed positive-mode calculations. The ordinary
proofs above supply the uniform statements; no floating optimization is
used as a certificate.
