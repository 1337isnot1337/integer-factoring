# Folded Boolean cubes: an explicit all-dimension crossing and its exact limit

## Outcome

The Bell affine search leads to an explicit family on `2^k` physical
labels, with `k+1` crossing pieces, in every characteristic other than
two. Its exact row and column degree profiles are the two reversals of

```math
\left(\binom{2k+1}{0},\binom{2k+1}{1},\ldots,
                    \binom{2k+1}{k}\right).
```

The complete coefficient-span method has optimal leading contraction

```math
\rho_*=H_2(1/4)=2-\tfrac34\log_2 3
                   =0.811278124459\ldots .
\tag{1}
```

Both the upper bound and a lower bound allowing arbitrary collective
coefficient recombinations are proved below. This closes this particular
folded-cube route above the desired `1/sqrt(2)` threshold. It is a new
effective family and a scoped obstruction, **not a stronger factoring
bound**. The selected factoring exponent remains `1/2`.

The construction uses the exact [Bell two-cut interface](bell-diagonal-crossings.md).
The lower argument adapts the triangular endpoint principle in
[coefficient-basis.md](coefficient-basis.md), with different local tensors
and multiplicities. Ordinary mathematics and code used substantial OpenAI
Codex assistance. Independent proof review remains open; no Lean is added.

## 1. The folded cube and its signed transform

Put `n=2k+1`. Let `V` be the even-weight subspace of `F_2^n`; its size is
`2^(2k)`. Give it the quadratic form

```math
q(S)=|S|/2\pmod2,
\qquad \ell(S)=\min(|S|,n-|S|).
```

Since `n` is odd, complementation pairs an even set with an odd set.
Exactly `binom(n,j)` elements of `V` have level `j`, for `0<=j<=k`.
The polar form of `q` is the dot product on `V`, which is nondegenerate:
the orthogonal complement of the even subspace in the full binary space
is the all-ones line, and that line is not even.

The Gaussian sum is the integer identity

```math
\sum_{S\in V}(-1)^{q(S)}
=\operatorname{Re}(1+\mathrm i)^n
=\varepsilon_k 2^k,
\qquad \varepsilon_k=(-1)^{k(k+1)/2}.
```

Define the regular involution

```math
(T_k f)(S)=\frac{\varepsilon_k}{2^k}
       \sum_{U\in V}(-1)^{q(S+U)}f(U).
\tag{2}
```

Its rows sum to one and `T_k^2=I` by character orthogonality.

This is the physical Bell transform up to separate row and column
permutations. To see this explicitly, choose a symplectic basis for the
polar form. In its coordinates `z=(a,b)`, write
`q(z)=a dot b+B(g,z)`. Its Gaussian sign is `(-1)^(a_g dot b_g)`.
Consequently the matrix in (2) equals the Bell matrix of equation (3) in
the parent note, with its input index translated by `g`. A translation
permutes the eigenvalue vector and preserves the all-ones endpoint.
The constructor implements this basis change by binary symplectic
elimination, without an isometry oracle.

## 2. Explicit polynomials and the exact dual identity

Define integer polynomials

```math
F_0(z)=1,\quad F_1(z)=2-z,\qquad
F_m(z)=2F_{m-1}(z)-z^2 F_{m-2}(z)\quad(m\ge2).
\tag{3}
```

They have degree at most `m`, constant coefficient `2^m`, and
`F_m(1)=1`. Let `G_m(z)=z^m F_m(1/z)`; then

```math
G_0=1,\quad G_1=2z-1,\qquad
G_m=2zG_{m-1}-G_{m-2}.
```

Set

```math
f_S(z)=z^{\ell(S)}F_{k-\ell(S)}(z).
```

The key identity, over the rationals and every odd-characteristic field, is

```math
\boxed{\quad (T_k f(z))_S=G_{k-\ell(S)}(z).\quad}
\tag{4}
```

Here is an algebraic proof rather than a finite extrapolation. In a
formal algebraic extension choose `A,B,i` with
`A^2+B^2=2`, `AB=z`, and `i^2=-1`. Recurrence (3) gives

```math
F_m(z)=\frac{A^{2m+1}+B^{2m+1}}{A+B}.
```

Thus for an even set of weight `w`,

```math
f_S(z)=\frac{A^{n-w}B^w+A^wB^{n-w}}{A+B}.
```

Apply (2), using `(-1)^(|S+U|/2)=i^(|S+U|)` and selecting even
input sets by `(1+(-1)^|U|)/2`. With `U=A+iB` and `W=A-iB`, the
result, after writing `m=k-w`, is

```math
\frac{(1-i)U R^m+(1+i)W Q^m}{2(A+B)},
\qquad R=-\frac i2 U^2,\quad Q=\frac i2 W^2.
```

For completeness the scalar simplification uses
`epsilon_(m+w) i^w=epsilon_m` for even `w`, and the ratios of the
two coefficients at successive `m` are `-i/2` and `i/2`; the factors
`(UW)^w=2^w` cancel the corresponding normalization powers.
Now `R+Q=2z` and `RQ=1`. The displayed expression has initial values
`1` and `2z-1`, so it is exactly `G_m` by its recurrence. The recurrence
extends to all integer `m` with `G_(-m-1)=G_m`. This gives
`G_(k-ell(S))` also when `w>k`, proving (4).

All denominators above are legitimate in the formal extension when the
characteristic is odd; `A+B` is not identically zero. Clearing formal
denominators proves an identity over `Z[1/2]`, so specialization in any
odd characteristic is valid. Neither `i` nor `A,B` is needed by the
actual constructor, which uses (3) and the binary coordinate change.

## 3. Crossing pieces, ranks and complete cost

Use `f(z)` as the row eigenvalue vector in the Bell construction, with
the input translation from Section 1. Every component has value one at
`z=1`, so

```math
P(z)=\sum_{j=0}^k z^j A_j,\qquad
\sum_{j=0}^k A_j=E.
```

In characteristic different from two the coefficient at the beginning
of each row polynomial is `2^(k-ell)`, which is nonzero. Its valuation
is exactly `ell`. Equation (4) gives exact column degree `k-ell`, again
with nonzero leading coefficient. Therefore

```math
r_j\le\sum_{a=0}^{j}\binom na,
\qquad
c_j\le\sum_{a=0}^{k-j}\binom na.
\tag{5}
```

Over a fixed finite field, the exact ranks retain only nonzero polynomial
coefficients within these sums. The checker computes them literally.

Set label width `lambda=k`. The recurrence costs polynomial time per
polynomial; integer coefficient bit lengths are `O(k)` (for example,
the coefficient one-norm grows at most geometrically by (3)). The
physical Bell tensor and both diagonal factorizations enumerate at most
`2^O(k)` entries. Their rational denominators are powers of two with
`O(k)` bits; one common power of two clears them, with nonzero scale.
Over `F_3` the endpoint already has scale one. There are exactly `k+1`
pieces. Thus the full uniform generation and coefficient obligations
are paid; there is no compact-listing or free rank-factor oracle.

## 4. Uniform upper bound and optimality under recombination

Let `L=sum_(a=0)^k binom(n,a)3^(-a)`. Equation (5) gives

```math
r_j\le3^j L,\qquad c_j\le3^{k-j}L,
\qquad r_jc_j\le3^k L^2.
```

Since `L<=(4/3)^n`,

```math
\log_2(r_jc_j)\le
k\log_2 3+2(2k+1)\log_2(4/3)
=4H_2(1/4)k+O(1).
\tag{6}
```

This proves the effective width bound at leading contraction (1).
Ceilings add only a constant. A rational contraction strictly above (1)
therefore meets the transfer contract. Its resulting factoring bound is
weaker than the selected exponent `1/2` and is not promoted as an improvement.

The upper bound is asymptotically optimal even under arbitrary collective
recombinations of these `k+1` coefficients. The row value of `A_j` on
level `ell` vanishes for `j<ell`, and at `j=ell` it is `2^(k-ell)`.
This triangular matrix proves the `A_j` are linearly independent in
every permitted field.

For `F=sum_j b_j A_j`, let `L0,U0` be the least and greatest indices
with `b_j!=0`. Its row eigenvalue on level `U0` cannot cancel, so
`rank_row F>=binom(n,U0)`. Reversing (4), the column eigenvalue on level
`k-L0` likewise cannot cancel, so
`rank_col F>=binom(n,k-L0)`.

Every decomposition of `E` inside this coefficient span has some piece
with `b_floor(k/2)!=0`, since the sum of coefficients at that coordinate
must be one. For this same piece both ranks are at least
`binom(2k+1,floor(k/2))`: binomial coefficients increase throughout
indices `0,...,k`. Stirling's formula gives

```math
\log_2\binom{2k+1}{\lfloor k/2\rfloor}
=2k H_2(1/4)+O(\log(k+2)).
```

Together with (6), this proves exact optimal leading contraction (1)
within the displayed coefficient span, with no restriction on piece
count. In particular `3^5<2^8` implies `H_2(1/4)>4/5>1/sqrt(2)`.
The proof does not exclude other code quotients, new local tensors,
or mixtures with families outside this span.

## 5. Checks and the next structural question

```sh
python3 tools/check_folded_bell.py --output evidence/folded-bell.json
```

The [saved evidence](../evidence/folded-bell.json) includes 324 exact
rational radial identities through `k=24`, 4,776 full Bell coefficient
identities at local dimensions `2,4,8,16` in characteristics `3,5,7`,
and exact coefficient-rank extrema through `k=256`. All shown finite
profiles satisfy the proved lower and upper bounds. The low-width
normalized ranks lie below the critical threshold in some small cases;
the proved limit (1) explains why they do not imply a stronger machine.

The useful remaining question is whether a different code quotient can
replace the antipodal fold while retaining a proved polynomial dual
identity with a smaller complete degree. A small covering radius alone
does not supply such an identity or pay its physical dimensions.
