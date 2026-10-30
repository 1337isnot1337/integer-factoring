# Complete central ranks from the small kernels of the two factors

## Result and scope

The complete homogeneous characteristic-two central piece at weight 64
has row rank **11,563** on the full 276-label regular-partition domain.
Its column rank is the same. This is within two of its 11,565 augmented
channels. The two factors separately have rank 11,564; the conclusion
concerns their complete product, including the vacuum term.

The calculation extends the previous exact complete central ranks 2, 12,
and 261 at weights 8, 16, and 32. A kernel-pairing identity avoids forming
the full `76,176` by `76,176` matrix at weight 64. The new native calculation
is checked against the existing small complete-matrix evidence and against
independent source-recursion coefficients.

These are finite ranks over `F_2`. They do not prove an asymptotic lower
bound, give a factoring lower bound, or exclude other decompositions.
The label domain is the full regular-partition space, not an asserted
binary-prefix restriction. They give no exponent above `1/2`.

The operators, crossing sum, and exchange identity are Claude (Anthropic)'s
Appendix B at [formal-math@16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The characteristic-two quotient, prime-power homogeneous duality, and
small complete-piece checks are in [modular-crossing.md](modular-crossing.md).
The new calculation and ordinary linear-algebra argument used substantial
OpenAI assistant support. No new Lean formalization is included.

## 1. A general product-rank identity

Let `U` be an `m` by `q` matrix and `V` a `q` by `n` matrix over any field.
Let `r_U=rank U`, `r_V=rank V`, and let columns of `K_U` and `K_V` be bases
of `ker U` and `ker V^T`, respectively. Then

```math
rank(UV)=r_U+r_V-q+rank(K_V^T K_U).
\tag{1}
```

Indeed, `im V=(ker V^T)^perp` for the ordinary nondegenerate coordinate
pairing. The intersection `ker U intersect im V` is the kernel of the
map `K_V^T` restricted to `ker U`. Its dimension is
`q-r_U-rank(K_V^T K_U)`. Subtract this from `r_V` to obtain (1).
No positive-definite inner product or characteristic-zero argument is used.

When both factors nearly have full channel rank, the matrix on the right
has small dimensions. Formula (1) strengthens the Sylvester lower bound;
it computes the rank exactly, including any cancellation between factors.

## 2. Apply the identity to the complete central piece

Let `h` be a power of two. Write `L_h` for all partitions of `h` into
distinct odd parts, and `M_h=L_(2h)`. Put `D=|L_h|` and `q=|M_h|`.
All coordinates below use monomial operator labels before the homogeneous
Gram inverse is applied. That inverse exists at these energies, as proved
in the modular note, and applying it to the first label is an invertible
change in each flattening. It therefore does not change either rank.

Define matrices on pair labels by

```math
U_{(a,y),nu}=[t^nu]a[h]t^y,\qquad
V_{nu,(v,u)}=[t^u]v[-h]t^nu.
\tag{2}
```

Here `a,y,u,v` belong to `L_h`, and `nu` to `M_h`. The two nonzero angle
indices are 0 and `2h`, because `binom(2h,i)` is odd only at those two
indices. The exchanged-order formula for the central piece is therefore

```math
A_0^{raw}=UV+ef^T,\quad
e_{(a,y)}=[1]a[-h]t^y,\quad f_{(v,u)}=[v=u].
\tag{3}
```

The second term is the vacuum intermediate state. In particular it must
not be omitted when computing the complete rank.
Augment the two factors as `U'=[U e]` and `V'=[V;f^T]`. Formula (1), with
`q+1` channels, computes `rank A_0` exactly.

The column flattening of the raw tensor equals the transpose of the raw
row flattening. To see this, use the original operator order instead of
the exchanged order in (3): its channel product interchanges the two
operator labels, and its vacuum term is `[a=u][1]v[-h]t^y`. In coordinates
`(a,u)` and `(v,y)` this is exactly the transpose of the row matrix.
The homogeneous dual transformation acts invertibly on the first label
in either flattening. Thus the normalized tensor also has equal row and
column ranks; it need not itself be a symmetric matrix in the chosen dual
coordinates.

## 3. Exact coefficient and kernel computation

The native helper uses the independent determinant expression for a
normal-ordered coefficient in characteristic two, already proved and
checked in the modular note. For operator parts `a`, output `pi`, and
input `sigma`, select a common unchanged subset of size
`(len(pi)+len(sigma)-len(a))/2`. The remaining parts are creations and
annihilations. The matrix entries for a part `m` and mode `k` are

```math
[(m-1)\mathbin{\&}\neg(k-1)=0]\quad\text{for creation},\qquad
[(m-1)\mathbin{\&}k=0]\quad\text{for annihilation}.
```

Take the determinant over `F_2` and sum over the unchanged subsets.
This avoids expanding the dual polynomials. All state coordinates and
operator monomials here are regular.

Stream the rows of `U'` and columns of `V'` through exact packed-bit
elimination. The resulting echelon bases give their ranks. For each free
channel coordinate, set that coordinate to one and solve the pivot equations
in ascending pivot order; these vectors form the two orthogonal complements.
Pair those small bases and perform another exact elimination for (1).

| Weight | Labels | Channels with vacuum | Factor ranks | Kernel dimensions | Kernel pairing rank | Complete central rank |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 2 | 6 | 3, 3 | 3, 3 | 2 | 2 |
| 16 | 5 | 24 | 16, 16 | 8, 8 | 4 | 12 |
| 32 | 23 | 277 | 269, 269 | 8, 8 | 0 | 261 |
| 64 | 276 | 11565 | 11564, 11564 | 1, 1 | 0 | 11563 |

The small cases show why factor ranks alone are not complete ranks.
At weight 16, for example, the kernel pairing contributes four extra
rank units beyond the Sylvester lower bound.

## 4. Verification and its boundary

```sh
python3 tools/check_modular_central_rank.py
python3 tools/check_modular_central_rank.py --full --output evidence/modular-central-rank.json
```

The default suite compiles the native helper and recomputes weights 8,
16, and 32, comparing their complete ranks with the existing independently
constructed full matrices in `modular-crossing.json`. The full suite adds
weight 64. Each case also checks 300 sampled triples and every diagonal
vacuum pairing against the source's independent normal-order recursion,
including both operator orientations and the added vacuum channel.

The result records source and input hashes. It is an author-checked native
finite computation, supported by independent coefficient tests and the
ordinary identity (1). It is not a Lean certificate or an independent
second enumeration of every large coefficient.

The next asymptotic question is whether the almost-full channel ranks and
the small kernel pairing admit a proof for growing powers of two. The
finite table alone does not answer it. A construction with a better profile
can also use different homogeneous duals, different labels, or complete
piece recombinations, as the [recombination search](modular-recombination.md)
illustrates on a different family.
