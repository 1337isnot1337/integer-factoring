# Dual gauges and exact mode filters

## Results and scope

The integral Schur construction admits two explicit descriptions of its
dual freedom over a prime field.

* At a prime-power weight `h=p^r`, every dual to the **full** Schur basis
  whose parameter has weight at most `h` is its homogeneous dual plus a
  known direct sum of translation gauges. These gauges annihilate both
  modes `-h` and `+h`.
* Higher-weight parameters can filter the mode index modulo a power of
  `p`. Their crossing pieces are exactly consecutive interval sums of
  the original pieces. This construction explains how higher-weight
  duals escape the preceding creation-mode rigidity.

These are exact constructive identities. They neither improve the
known joint rank bound nor turn a finite rank improvement into an
asymptotic factoring result. In particular, the filtering operation is
an explicit realization of coefficient recombination, not a free rank
reduction.

The operators are those of the [pinned Anthropic Appendix B](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The integral lattice and prime-power homogeneous pairing are proved in
[integral-schur-crossing.md](integral-schur-crossing.md). The divided
translation identity is proved integrally in
[schur-fusion-gram.md](schur-fusion-gram.md), Section 1.

## 1. Integral divided translation

Let `Lambda_Z` be the symmetric-function lattice with Schur basis, and
let `U(p_m)=m p_(m+1)`. Write `U^(k)=U^k/k!` for its divided powers.
These preserve the integral lattice. For example, with
`H(u)=sum_(m>=0) h_m u^m`, one has

```math
e^{zU}H(u)=\frac{H(u+z)}{H(z)},
\qquad
U^{(k)}h_m=
\sum_{j=0}^k(-1)^j\binom{m+k-j}{m}h_{m+k-j}e_j.
\tag{1}
```

The first identity follows by applying `e^(zU)` to the power-sum
exponent of `H`; the second takes the coefficient of `u^m z^k`.
Multiplicativity of the first identity proves integrality on every
complete-symmetric monomial, hence on the whole lattice.

For a homogeneous parameter `c` of weight `w`, the source mode identity is

```math
(U^{(k)}c)[r]=\binom{r-w}{k}c[r].
\tag{2}
```

All binomial coefficients with negative upper argument have their
integer meaning. Equation (2) holds before reduction modulo `p`; no
factorial is inverted in the prime field.

## 2. Every bounded-weight full-Schur dual

Fix `h=p^r>0` and put `V_w=Lambda_(F_p,w)`. Define

```math
\Phi_h:\bigoplus_{w=0}^h V_w\longrightarrow V_h,
\qquad
\Phi_h(c_w)=U^{(h-w)}c_w\ (w>0),
\quad \Phi_h(1)=0.
\tag{3}
```

For `0<w<=h`, the coefficient in (2) at mode `+h` is one. At
mode `-h`, it is

```math
\binom{-h-w}{h-w}
=(-1)^{h-w}\binom{2h-1}{h-w}
\equiv1\pmod p.
\tag{4}
```

Indeed, all of the lower `r` base-`p` digits of `2h-1` are `p-1`.
For `0<=k<h`, Lucas' identity gives
`binom(2h-1,k)=(-1)^k` modulo `p`. The constant parameter has both
nonzero-index modes equal to zero. Consequently, as operators on every
state degree,

```math
(\Phi_h a)[-h]=a[-h],\qquad
(\Phi_h a)[h]=a[h].
\tag{5}
```

Since `Phi_h` is the identity on `V_h`, it is a projection. Its kernel
has the explicit direct-sum description

```math
\ker\Phi_h
=V_0\ \oplus\bigoplus_{w=1}^{h-1}
 \{c-U^{(h-w)}c:c\in V_w\}.
\tag{6}
```

Directness follows by taking the distinct lower homogeneous components.
Its dimension is `sum_(w=0)^(h-1) p(w)`.

Let `B_h` be the homogeneous vacuum Gram operator, so
`[1]a[-h]b=<B_h a,b>` for `a,b in V_h`. It is invertible when
`h=p^r`. On parameters of weight at most `h`, the full vacuum readout is
therefore exactly `B_h Phi_h`. If `a_x^hom=B_h^(-1)s_x`, then all
solutions of the full set of duality equations
`[1]a_x[-h]s_v=[x=v]` are

```math
a_x=a_x^{\mathrm{hom}}+g_x,
\qquad g_x\in\ker\Phi_h.
\tag{7}
```

In particular, every such full dual family has the same `+h` mode as
the homogeneous family. This does **not** say that its complete central
piece is fixed: the angle coefficients depend on each homogeneous
parameter weight, and the gauges can change other modes. It also does
not assert (7) for duals specified only on a proper label subspace;
their vacuum functionals have additional complementary freedom.

For `h=O(lambda^2)`, the complete parametrization, the integral
translation maps, and all vacuum equations are generated in
`2^O(lambda)` bit time. The relevant degree spaces have
`exp(O(sqrt(h)))` dimensions, and integer translation coefficients
have polynomial-in-`h` bit lengths. A generic choice can use exponential space. The companion
[prime-power involution](prime-power-gram-involution.md) now gives
polynomial-space queries for the homogeneous choice. This paragraph
does not give that space bound for every freely chosen gauge.

## 3. A finite Newton filter on all integer mode indices

Let `T=p^e`, let `a` be homogeneous of weight `A`, and let the second
parameter and the input/output states have weight `h`. Choose
`a_0 in {0,...,T-1}` with `a_0=-A-h mod T`, and define

```math
\mathcal F_Ta
=\sum_{k=a_0}^{T-1}
 (-1)^{k-a_0}\binom{k}{a_0}U^{(k)}a.
\tag{8}
```

Then

```math
(\mathcal F_Ta)[r]=[r\equiv-h\pmod T]a[r]
\quad(r\in\mathbf Z).
\tag{9}
```

To prove this, `binom(n,k) mod p` is `T`-periodic in every integer
`n` when `k<T`, because `(1+x)^T=1+x^T`. On residues
`0<=n<T`, binomial inversion gives

```math
\sum_{k=a_0}^{T-1}
 (-1)^{k-a_0}\binom{k}{a_0}\binom nk
=\binom n{a_0}(1-1)^{n-a_0}
=[n=a_0].
\tag{10}
```

The sum is empty when `n<a_0`. Apply (2) with `n=r-A`.
This proves (9), including negative modes. In particular the `-h`
mode, and hence all vacuum duality equations, is unchanged.
Applying (8) to each homogeneous part defines the same filter on an
inhomogeneous parameter.

If `T>2h`, the filter kills the `+h` mode while retaining `-h`.
Its maximum parameter weight can reach `A+T-1`; this explains why it
does not contradict the bounded-weight theorem in Section 2.

## 4. The filtered pieces are exact interval sums

Use the generalized angle from
[fermionic-piece-polynomial.md](fermionic-piece-polynomial.md), Section 1:

```math
Q_s(a,b)=\sum_i(-1)^i\binom{A+h}{i}
               a[s-h+i]b[h-s-i],
\qquad Q_a(z)=\sum_s Q_s(a,b)z^s.
\tag{11}
```

For an inhomogeneous first parameter, sum the respective homogeneous
angles. Equation (2), or the formal field translation identity, gives

```math
Q_{U^{(k)}a}(z)
=(z-1)^{A+h+k}\partial_z^{(k)}
       \bigl((z-1)^{-A-h}Q_a(z)\bigr).
\tag{12}
```

Here `partial_z^(k)` is the divided formal derivative, defined on
Laurent monomials by integer binomial coefficients before reduction.
Expanding `Q_a(z)=sum_(n>=0)q_n(z-1)^n` at `z=1`, equations
(8), (10), and (12) show that the new polynomial keeps exactly the
terms whose Taylor degree is divisible by `T`.

For every integer `s`, write `s=qT+r` with `0<=r<T`. In
characteristic `p`,

```math
z^s=(1+w)^r(1+w^T)^q,\qquad w=z-1.
```

Only the constant term of `(1+w)^r` contributes to Taylor degrees
divisible by `T`. The selected series is therefore `z^(qT)`.
The result is a finite Laurent polynomial, and

```math
Q_{\mathcal F_Ta}(z)
=\sum_s Q_s(a,b)z^{T\lfloor s/T\rfloor}.
\tag{13}
```

Thus its piece at shift `qT` is precisely
`sum_(qT<=s<(q+1)T)Q_s(a,b)`, and every other piece is zero.
The identity extends componentwise to every inhomogeneous dual.
It preserves the crossing sum because `z=1` is unchanged.

For example, a homogeneous weight-`h` family has support
`[-2h,2h]`. If `T>2h`, (13) leaves two interval pieces, the sum
of the negative shifts and the sum of the nonnegative shifts. Either
piece can have large flattening ranks. No estimate on those ranks
follows from the number of remaining pieces.

For `A,h,T=O(lambda^2)`, (8) and all its coefficients are generated
in `2^O(lambda)` bit time by the integral translation algorithm.
Alternatively, (13) constructs the new tensors directly from the old
ones without expanding the higher-weight parameters.

## 5. Verification

Run

```sh
python3 tools/check_dual_gauge_filters.py --output evidence/dual-gauge-filters.json
```

The checker compares the bounded-weight projection with independent
source vacuum and mode calculations over several prime fields. It
checks the complete kernel dimension, both exceptional modes, the
periodic Newton filters on negative and positive indices, and the
Laurent coefficient transformation. Small filtered Schur crossings
are also compared directly with interval sums of the original source
pieces. These finite checks validate the identities; the all-weight
statements are proved above.
