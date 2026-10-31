# Code quotients: preserved duality, exact degree cost, and the modular opening

## Structural result

An arbitrary nondegenerate binary code containing the all-ones word can
quotient the [folded-cube crossing](folded-bell-crossings.md) while
preserving its exact dual identity. The maximal common polynomial factor
can be removed explicitly. However, in characteristic zero the resulting
order/degree pressure bound is always at least `H_2(1/4)`. A small
covering radius by itself therefore does not improve that criterion.

Finite odd characteristic leaves a precise opening: cancellation among
minimum-weight coset leaders can raise polynomial valuations. Concrete
examples below exhibit this effect, without a factoring improvement.
Global prime-power weight divisibility only rescales the same polynomial
parameter and cannot improve optimized pressure. These are ordinary
theorems and exact controls, with no new Lean or independent proof review.
The complete factoring exponent remains `1/2`.

## 1. Nondegenerate quotients preserve the crossing transform

Let `n` be odd, and let `C<=F_2^n` contain the all-ones vector and satisfy
`C intersect C^perp={0}`. Write `r=dim C` and `V=C^perp`, of dimension
`m=n-r>=2`. All elements of `V` have even weight.
Put `D=C intersect Even(n)`. The orthogonal decompositions are

```math
C=\langle\boldsymbol1\rangle\perp D,
\qquad \operatorname{Even}(n)=V\perp D.
```

The dot product on each even summand is alternating and nondegenerate.
Hence `r` is odd and `m` is even. The quadratic form `q(v)=wt(v)/2 mod2`
on `V` gives a regular signed transform `T_V` of Bell type on
`2^(m/2)` physical labels, by the isometry and translation in the parent
note. All physical dimensions are retained.

Let `K=(n-1)/2` and let `f_w(z)` be the ambient folded polynomial. Average
over the discarded even directions:

```math
\Phi_v(z)=|D|^{-1}\sum_{d\in D}f_{v+d}(z).
```

The ambient transform factors as `T_V tensor T_D`. The uniform averaging
functional is fixed by `T_D`, because that transform is regular and
symmetric. Applying the parent dual identity therefore gives

```math
T_V\Phi(z)=z^K\Phi(1/z),\qquad \Phi(1)=\boldsymbol1.
\tag{1}
```

Thus code projection preserves the full crossing, not just a scalar
rank estimate. In any odd characteristic `|D|` is invertible.

## 2. The maximal common division and its complete degree

For `v in V`, let

```math
W_{C+v}(A,B)=\sum_{c\in C}A^{n-\operatorname{wt}(v+c)}
                                  B^{\operatorname{wt}(v+c)},
\qquad A^2+B^2=2,\quad AB=z.
```

Pairing complements gives `Phi_v=W_(C+v)/(|D|(A+B))`. Binary character
orthogonality, equivalently the direct coset Fourier expansion, yields

```math
\Phi_v(z)=2^{-K}\sum_{y\in V}(-1)^{v\cdot y}
 (1+z)^{K-\operatorname{wt}(y)/2}(1-z)^{\operatorname{wt}(y)/2}.
\tag{2}
```

This identity follows by substituting `(A+B)^2=2(1+z)` and
`(A-B)^2=2(1-z)` into the expansion of `W`. It also constructs the
path without algebraic numbers.

Let `w_max=max_(y in V) wt(y)` and `h=w_max/2`. The character matrix on
`V` is invertible over every odd-characteristic field. Consequently the
gcd of all components in (2) equals the gcd of its displayed diagonal
polynomials. The `y=0` term and a maximum-weight term show that the
**maximal** common divisor is `(1+z)^(K-h)`, up to a nonzero scalar.
There is no additional common `(1-z)` factor or another polynomial factor.

Divide by `((1+z)/2)^(K-h)`, whose value at one is one. The resulting
path has the explicit formula

```math
P_v(z)=2^{-h}\sum_{y\in V}(-1)^{v\cdot y}
 (1+z)^{h-\operatorname{wt}(y)/2}(1-z)^{\operatorname{wt}(y)/2},
\qquad T_V P(z)=z^hP(1/z),\quad P(1)=\boldsymbol1.
\tag{3}
```

Its maximum degree is exactly `h`. It has at most `h+1` coefficients.
Choosing an information set for `V`, some codeword is one on all its
`m` information coordinates, so `w_max>=m` and

```math
h\ge m/2.
\tag{4}
```

For a uniform code family with `n=O(m)`, direct enumeration, coefficient
generation and Bell factor listing cost `2^O(m)` bit operations. Degrees
and integer coefficient bit lengths are `O(n)`. A fixed finite code may
also be used as a fixed seed and tensor-powered with the parent transfer.
Larger ambient lengths require their own effective enumeration bound.

## 3. Characteristic-zero pressure cannot improve the folded limit

Let `delta_v=dist(v,C)`. In characteristic zero, the leading coefficient
of the averaged folded polynomial is a positive power of two times the
number of minimum-weight coset leaders. It cannot cancel. Division by
the common factor has nonzero constant term and preserves this order.
Thus `ord_0 P_v=delta_v`. By (3), opposite degrees are `h-delta_v`.

The [uniform order/degree bound](bell-polynomial-synthesis.md) is

```math
Q(x)=x^h\left(\sum_{v\in V}x^{-\delta_v}\right)^2,
\qquad x\ge1.
```

Choose `m` independent columns of a parity-check matrix of `C`. Every
coset has a unique representative supported on those coordinates.
Its minimum weight is at most the weight of this representative. Hence

```math
\sum_{v\in V}x^{-\delta_v}\ge(1+x^{-1})^m.
```

Combining with (4),

```math
Q(x)\ge x^{m/2}(1+x^{-1})^{2m}.
```

The minimum on `x>=1` occurs at `x=3`. Since physical label width is
`lambda=m/2`, this proves

```math
\frac{\log_2 Q(x)}{4\lambda}\ge H_2(1/4).
\tag{5}
```

This is a restriction on the displayed order/degree pressure criterion,
not a lower bound for all actual coefficient ranks, arbitrary
refactorizations, or factoring. It includes the maximal common scalar
division, so that division is not a missing source of improvement.

## 4. The finite-characteristic opening and the rescaling trap

In characteristic `p!=2`, a coset with `N_v` minimum-weight leaders has
leading coefficient equal to `N_v` times a nonzero power of two.
Its order rises above `delta_v` exactly when `p` divides `N_v`.
The maximal gcd and degree statements remain valid; positivity in
Section 3 is what can fail. A better order/degree witness within this
construction must exploit such modular cancellations, beyond merely
changing a covering radius.

One tempting global cancellation gives no gain. Let a binary generator
matrix of `V` have column multiplicities `a_g` for nonzero column types
`g in F_2^m`. For an odd prime power `p^s`,

```math
\operatorname{wt}(uG)\equiv0\pmod{p^s}\ \hbox{for every }u
\quad\Longleftrightarrow\quad
a_g\equiv0\pmod{p^s}\ \hbox{for every nonzero }g.
\tag{6}
```

Indeed, write the weight as
`sum_g a_g(1-(-1)^(u dot g))/2`. Fourier inversion over `F_2^m` is
invertible modulo `p^s`, since its determinant is a power of two.
Each nontrivial Fourier coefficient is a unit times `a_g`, proving (6).

Remove `p^s` identical copies of every nonzero column type, ignoring
zero columns and padding to odd ambient length if needed. Its dimension,
nondegenerate polar form and physical transform are unchanged, because
`p^s` is odd. In `F_p`, formula (3) and Frobenius give exactly

```math
P_v^{\rm repeated}(z)=P_v^{\rm reduced}(z^{p^s}).
```

Both the orders and degree multiply by `p^s`. Its pressure is the old
pressure with `x` replaced by `x^(p^s)`, so the optimized value is
identical. Global divisibility or full coordinate tripling alone cannot
improve the score. Nonuniform cancellation and actual rank compression
remain separate possibilities.

## 5. Exact controls and current next step

```sh
python3 tools/check_code_folded_bell.py --output evidence/code-folded-bell.json
```

The [saved controls](../evidence/code-folded-bell.json) check five
structurally specified codes in characteristics `3,5,7`, with complete
polynomial coefficients, exact duality and endpoint readback. Independent
rational coset averaging verifies (2) and the common division; information
set pressure bounds are checked at four rational scales.

The one-coordinate tripling example has 22 cosets whose valuations rise
in characteristic three. A mixed LCD example has 15. Neither supplies
a critical-pressure witness. Full tripling exactly rescales the base
orders and degree by three, providing a concrete negative control for
the theorem above. No large code search is claimed.

The conceptual target is now specific: nonuniform modular coset-leader
cancellation with an effective family and a strict complete-pressure
witness, or a rank argument outside this pressure interface. These
claims and controls used substantial OpenAI Codex assistance.
