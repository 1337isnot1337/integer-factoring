# Complete Bell-diagonal seed spaces

## Result and the remaining constructive target

There is an explicit crossing family in every odd characteristic whose two
cuts are simultaneously diagonal. Its complete local space has dimension
`d^2`, for `d=2^k`, and the change between the two eigenvalue vectors is a
regular signed Hadamard transform. It includes local tensors beyond the
three-pairing space. This is a different construction from the auxiliary
noncommuting order-three Weyl path in
[PR #31](https://github.com/1337isnot1337/integer-factoring/pull/31).

The complete symmetric tensor-power span of this larger family is still
obstructed at physical dimensions **two and four**, even with arbitrary
collective coefficient recombinations, extension fields, and any piece
count. The leading lower floors are respectively `0.774397470347...` and
`0.707133841396...`, both strictly greater than `1/sqrt(2)`.

At physical dimension eight the same argument gives only
`0.656961057491...`; it does **not** exclude the desired contraction.
Dimension eight and higher therefore remain constructive targets, rather
than implied positive results. The supported factoring exponent stays
`1/2`. There is no new Lean formalization or independent proof review.

The determinant argument adapts the local-sector argument in
[the all-field three-pairing note at 78bf3b2](https://github.com/1337isnot1337/integer-factoring/blob/78bf3b23677f74533095d581951890075e40d21f/research/all-field-pairing-barrier.md).
Here the local algebra and the two-cut identity are derived directly.
Substantial OpenAI Codex assistance was used.

## 1. Physical tensor and both eigenbases

Let `G=(F_2)^k`, `d=2^k`, and let the coefficient field `K` have
characteristic different from two. Dots below are the binary dot product;
all additions inside a dot or an index in `G` are XOR. For one vector of
indeterminates `s=(s_(a,b))_(a,b in G)`, define

```math
P_s(x,y,u,v)=\frac1d[x+y=u+v]
             \sum_{b\in G}(-1)^{b\cdot(x+v)}s_{x+y,b}.
\tag{1}
```

At `s=1` the character sum vanishes unless `x=v`. The support condition
then forces `y=u`, so `P_1=E`, the exact scale-one crossing.

Normalize the row cut to rows `(x,y)` and columns `(v,u)`. The vectors

```math
B_{(x,y),(a,b)}=[y=x+a](-1)^{b\cdot x}
\tag{2}
```

form an invertible basis: each fixed-`a` block is the character table of
`G`, and its inverse is its transpose divided by `d`. Substitution in
(1) gives `R(P_s) B=B diag(s)`.

Normalize the other cut to rows `(x,u)` and columns `(v,y)`. Its support
condition is `x+u=v+y`. Inside the block with this common value `c`,
put `r=x+v`; the matrix entry is

```math
\frac1d\sum_b(-1)^{b\cdot r}s_{c+r,b}.
```

Fourier transform in `r` diagonalizes this block. The eigenvalue at `(c,e)`
is therefore

```math
(Ts)_{c,e}=\frac1d\sum_{a,b\in G}
              (-1)^{(a+c)\cdot(b+e)}s_{a,b}.
\tag{3}
```

The column basis is again (2), now on its own pair of physical indices.
Character orthogonality proves `T^2=I` and `T1=1`. These identities hold
in every field of characteristic different from two; no square roots or
roots of unity beyond `-1` are required. In particular the construction
is valid over `F_3` without descent or a scalar decoder.

## 2. Complete tensor-power span and effective polynomial paths

Write `q=d^2` and

```math
P_s^{\otimes t}=\sum_{|\nu|=t}s^\nu A_\nu.
\tag{4}
```

In the row Bell basis, `A_nu` is the diagonal indicator of words of type
`nu`. Hence these tensors are linearly independent in every permitted
characteristic. If `F=sum_nu f_nu A_nu`, its exact row rank is

```math
\sum_{\nu:f_\nu\ne0}\binom t{\nu_1,\ldots,\nu_q}.
\tag{5}
```

In the other cut the eigenvalue on a word of type `mu` is the sum of
`f_nu` times the coefficient of `s^nu` in `(Ts)^mu`. The same integer
multinomial multiplicities give the column rank. Multiplicities count physical dimensions and must
**not** be reduced modulo the coefficient characteristic.

For any fixed polynomial vector `s(z)` with `s(1)=1`, the `O(t)` monomial
coefficients of `P_(s(z))^tensor t` sum to the crossing. Both cut ranks
can be computed from products of the local eigenpolynomials `s_j(z)` and
`(Ts(z))_j`. In a fixed finite field, all physical entries and factors can
be generated in `2^O(lambda)` bit operations at label width
`lambda=t log2(d)`. Coefficient products have degree `O(t)`. Explicit Bell
diagonal factors enumerate the eigenwords where the coefficient is
nonzero and multiply the corresponding columns of `B` and rows of `B^-1`.
This includes the cost of listing those indices and factors.

For other label widths, choose `t=ceil(lambda/log2 d)` and restrict all
legs to the first `2^lambda` words. The endpoint persists, ranks cannot
increase, and the width loss is `O(1)` for fixed `d`.

These observations supply an effective search interface, **not** a proved
uniform contraction bound. To invoke
[the factoring transfer](pieces-to-factoring.md), one still needs the
two-cut width inequality for every coefficient at every length, with
its logarithmic slack. A finite rank profile does not supply that proof.

## 3. A common determinant functional

Take `t=2q m` with `m>=1`, and define

```math
p(s)=\left(\prod_{j=1}^q s_j\right)^m
     \left(\prod_{j=1}^q (Ts)_j\right)^m.
\tag{6}
```

It is homogeneous of degree `t`, satisfies `p(1)=1`, and is invariant
under `T`. Pair its monomial coefficients with (4):

```math
L_p(F)=\sum_\nu [s^\nu]p(s)\,f_\nu.
\tag{7}
```

Since `E^tensor t=sum_nu A_nu`, `L_p(E^tensor t)=1`. In any decomposition
of a nonzero scalar multiple of the crossing into tensors in the full
span (4), at least one piece `F` therefore has `L_p(F)!=0`.

Every monomial in (6) has each exponent at least `m`. Thus such an `F`
has a nonzero row eigenvalue on a word type with every count at least
`m`. In the column coordinates `w=Ts`, (6) has exactly the same form,
because `T^2=I`. The same functional detects a column type with all its
counts at least `m`. Cancellations in a positive characteristic can
remove monomials, but cannot create a monomial violating these bounds.

Among `q` counts at least `m`, summing to `2qm`, the multinomial is
minimized at `((q+1)m,m,...,m)`: moving excess from a smaller count to a
larger count increases the factorial product. Consequently the **same
piece** has both ranks at least

```math
J_{q,m}=\frac{(2qm)!}{((q+1)m)!\,(m!)^{q-1}}.
\tag{8}
```

This proof permits every coefficient recombination in (4) and any number
of pieces. It is not merely a statement about a chosen polynomial path.
The lower bound for general `t` follows by fixing fewer than `2q` final
physical positions to `x=y=u=v=0`. This preserves `E`, gives a submatrix
of each cut, and sends the span (4) into the corresponding shorter span:
the generating tensor acquires only the scalar factor
`(d^-1 sum_b s_(0,b))^r`. No rank can increase.

## 4. The small-seed exclusion

Stirling's formula applied to (8), for fixed `q`, gives

```math
\log_2 J_{q,m}=t\left(\log_2(2q)
                  -\frac{q+1}{2q}\log_2(q+1)\right)
                  +O_q(\log(t+2)).
```

At `lambda=t log2 d`, the forced leading contraction is

```math
\gamma_q=\frac{\log_2(2q)-\frac{q+1}{2q}\log_2(q+1)}{\log_2 q}.
\tag{9}
```

For `q=4,16`, the exact integer inequalities

```math
(2q)^{2q\cdot239}>(q+1)^{(q+1)\cdot239}q^{2q\cdot169},
\qquad 2\cdot169^2>239^2
```

prove `gamma_q>169/239>1/sqrt(2)` without a floating comparison.
For `q=64` the first inequality reverses. This argument therefore stops
excluding the target at `d=8`; it does not show a construction exists there.
It also does not cover characteristic two or arbitrary nonsymmetric
position-dependent families.

## 5. Reproduction and current evidence

```sh
python3 tools/check_bell_crossing.py --output evidence/bell-crossing.json
make verify
```

The [saved replay](../evidence/bell-crossing.json) checks 297 transform
identities, 594 physical matrix eigenbasis identities, 192 independent
modular elimination ranks, and 324 literal tensor-word pencil profiles.
Physical dimensions are `2,4,8`, coefficient characteristics `3,5,7`.
The literal profiles exhaust all 81 ternary one-bit pencils through
length four. The negative controls retain integer multiplicities, including
rank 64 rather than 1 for a length-three four-symbol identity.

These finite checks support the displayed algebra. The all-length claims
are the ordinary proofs above, not an extrapolation from those checks.
