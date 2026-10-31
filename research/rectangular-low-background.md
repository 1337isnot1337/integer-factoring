# The all-low Gaussian background is a shorter interval

## Result and scope

Fix the low-only state on all four extreme physical cores of the
partial rectangular seed. The remaining Gaussian tensor has exactly
the two natural Smith profiles of an interval of length `L=M-2`.
This holds for every power of two `M>=4`, including arbitrary values
of the four bulk-charge gauges.

The scalar prefactor cancels the energy of the frozen extreme states.
Thus the conditional pressure is the ordinary Gaussian interval
pressure with `M-2` bulk modes. This identifies a preserved background
inside the partial seed. It does not prove a pressure lower bound
for the unrestricted seed, for gauges that keep other extreme states
active, or for core corrections that change this background.

## 1. Conditional Gaussian and scalar

Use the original one-particle matrix from the
[mixed repair construction](rectangular-f4-repair.md), with rows `(v,y)`
and columns `(u,x)`. The partial correction fixes output core row
zero. The top-code mixtures are inactive in a low-only state, and
the quartic requires occupied top modes on `y,v`. Consequently the
all-low restriction is independent of both code parameters and of
the quartic numerator.

In the Gaussian exterior compound, fixing the low mode as occupied
on every leg and the top mode as empty means taking a Schur
complement. The occupied-low two by two block before natural gauges is

```math
D_{00}=\begin{pmatrix}
c&w^{2M}\\1&cw^{2M-2}
\end{pmatrix},\qquad c=1+w,
\qquad\det D_{00}=w^{2M-2}.
```

Delete the two top rows and columns, then take the Schur complement
of `D_00`. Its remaining mode indices are `1,...,M-2` in each block.
The full conditional Gaussian is `det(D_00)` times the exterior
compound of this Schur complement.

After the natural physical gauges, the scalar order is `4-4M` at
infinity and zero at zero. A pair of infinity bulk charges
`lambda,mu` adds `2(lambda+mu)` to the scalar order. The zero
charges cancel on the two occupied row modes and the two occupied
column modes.

## 2. Exact infinity profile

Put `z=w^-1`, `t=z^2`, `u=1+z`, and `q=u^2`. Let `G` be the
full naturally gauged Gaussian one-particle matrix. Its known
regular diagonalization gives

```math
\kappa_j=1-2M+2j,\qquad 0\le j<2M,
\qquad\det G=1.
```

The determinant equality follows directly from the explicit regular
maps: their determinants are `u^-M` and `u^M`, while the sum of
the diagonal orders is zero.

Let `G'` delete the two top rows and columns. The explicit inverse
`G^-1=V Sigma^-1 U` has top two by two block

```math
(G^{-1})_{\mathrm{top},\mathrm{top}}
=z^{1-2M}\begin{pmatrix}u&1\\1&u\end{pmatrix}.
\tag{1}
```

For completeness, the two entries requiring a sum reduce to

```math
z+q^M z^{1-2M}=z^{1-2M},
\qquad
u z^{-1}\sum_{a=0}^{M-1}(q/t)^a=u z^{1-2M}.
```

Here `q^M=1+t^M`, since `M` is a power of two. The geometric
sum equals `t^(1-M)`, because `q+t=1`. Jacobi's complementary-minor
identity and (1) therefore give

```math
\det G'=z^{4-4M}.
\tag{2}
```

We use the following elementary DVR fact. If a matrix has sorted
invariant orders `kappa_i`, then a submatrix has its `i`th order
at least `kappa_i`. To see this, truncate a regular Smith
decomposition after its first `i-1` terms. The remainder is
divisible by `z^kappa_i`. Coordinate projections preserve both
that divisibility and the rank bound of the truncated part. A
matrix differing from rank at most `i-1` by such a remainder has
its `i`th invariant order at least `kappa_i`. Multiplying by one
global monomial first reduces Laurent matrices to integral ones.

Applying this fact to `G'` gives componentwise lower bounds by
the first `2M-2` orders of `G`. Their sum is exactly `4-4M`, which
equals (2). Every lower bound is therefore an equality.

The occupied-low block of `G'` is also the matrix in (1). Its
orders are `kappa_0,kappa_1`. Both off-diagonal rectangular blocks
adjoining it have every entry of order at least `kappa_1`.
Multiplication by the inverse occupied-low block is therefore
regular. Ordinary block elimination splits it off by regular
invertible maps.

The conditional Schur complement consequently has orders

```math
\kappa_2,\kappa_3,\ldots,\kappa_{2M-3}
=5-2M,7-2M,\ldots,2M-5.
\tag{3}
```

These are precisely `1-2L+2j` for `0<=j<2L`, where `L=M-2`.
Infinity charge gauges add `lambda+mu` to every one-particle
order. The regular transformations preserve total occupation, so
this is compatible with the complete exterior-compound profile.

## 3. Exact zero profile

Now put `z=w`, `t=z^2`, and `q=1+t`. The internal Gaussian
factors of the conditional tensor are regular units. Their
physical bulk charges cancel, since each internal edge contains
opposite charges. The two cross blocks are, up to transpose,

```math
z(D+qJ),\qquad z(DP+qJ),
\quad D=\mathrm{diag}(t^1,\ldots,t^L).
\tag{4}
```

Here `J` is the all-one matrix and `P` is the Boolean inclusion
matrix restricted to indices `1,...,M-2`. This interval of the
Boolean poset is order-convex, so `P^2=I` in characteristic two.
Every nonzero proper bit subset has an odd number of nonempty
subsets and of proper supersets. Hence `JP=J`, and right
multiplication by `P` identifies the profiles of the blocks in (4).

Write `A_L=D+qJ`. Its first diagonal entry is `t+q=1`.
Eliminating that unit pivot gives, for indices `i,j>=2`,

```math
t^i[ i=j ]+q+q^2
=t^i[ i=j ]+tq.
```

The remaining matrix is `t A_(L-1)` after relabeling. Induction
gives the exact orders `0,1,...,L-1` with respect to `t`.
Both blocks in (4) therefore have `z` orders

```math
1,3,\ldots,2L-1=2M-5.
\tag{5}
```

The zero charge sum adds `lambda+mu` to one block and its
negative to the other. Taking their exterior compound proves
the full zero profile, with the same charge dependence as an
interval of length `L`.

## 4. Pressure normalization

Relabel a retained physical mode `a` by `a'=a-1`. Its natural
total energy becomes

```math
4a-2M+2=4a'-2L+2.
```

The sum of the two frozen low-state energies is `4-4M`.
Their infinity charge contribution is `2(lambda+mu)`; their zero
charge contribution is zero. These are exactly the scalar factors
identified in Section 1. Any additional fixed extreme-core letter
weights also enter as the same scalar on the profile and physical
partition sides. They cancel from the conditional pressure equation.

Thus both profile partitions and the physical partition reduce
exactly to the Gaussian interval model at length `L=M-2`.
This can be used when a separate asymptotic argument proves that
all other extreme physical states have vanishing probability.
That additional hypothesis must be established for the family
under consideration; it does not follow from this background
identity alone.

## 5. Exact controls

The [checker](../tools/check_rectangular_low_background.py) constructs
the original polynomial Schur complement and profiles it directly.
Its [evidence](../evidence/rectangular-low-background.json) verifies
both one-particle spectra and the explicit inverse-minor identity
through `M=128`. It also restricts the original mixed four-mode
tensor to the all-low physical block and compares all sixteen
orders at both places under three charge choices.

```sh
python tools/check_rectangular_low_background.py \
  --output evidence/rectangular-low-background.json
```

The formulas above are uniform proofs. The finite checks audit
their implementation and the connection to the original tensor.
