# A characteristic-two barrier for the full pairing type span

## Result

The [small-seed strengthening](pairing-small-seeds.md) now closes the three
dimensions left unresolved by the bound in this note. Together the results
exclude the full pairing-type route for every fixed seed dimension in
characteristic two, including `o(t)` separate position classes.

Let `F` be any field of characteristic two and fix a local label dimension
`d>=2`. Define the three matching tensors

```math
E=[x=v][y=u],\qquad R=[x=y][u=v],\qquad C=[x=u][y=v],
```

and put `H=E+R+C`. The crossing identity `E=R+H+C` and its base ranks are
the known three-pairing baseline; see
[PR #23, Section 3 of the seed-span note at a1d68b4](https://github.com/1337isnot1337/integer-factoring/blob/a1d68b42c287e42ce71320182f4c66af82bf1bd0/research/seed-subspace-recombination.md).
The new result concerns arbitrary collective recombinations of its tensor
power types.

For word length `t`, let `V_t` be the full coefficient span of

```math
P_t(r,h,c)=(rR+hH+cC)^{\otimes t}.
```

It contains the crossing on all `d^t` words, since
`P_t(1,1,1)=E^{tensor t}`. Write `L=binom(d,2)`. In every exact decomposition
of that crossing into elements of `V_t`, some summand has **both** matrix
ranks at least

```math
M_t=\max_{0<=a<=t}\binom ta L^{t-a}
>=\frac{(L+1)^t}{t+1}.
\tag{1}
```

This allows arbitrary recombinations, any number of summands, and any
characteristic-two extension field. It also covers every shared scalar
polynomial grading of `R,H,C`, and the eight crossing-preserving coordinate
symmetries. It is not merely an obstruction to separate pure types.

At label length `lambda=t log2 d`, (1) forces the leading contraction
parameter

```math
rho>=rho_d:=\frac{log2(L+1)}{2log2 d}.
\tag{2}
```

For every `d>=5`, this is strictly greater than `1/sqrt(2)`. The full pairing
type span therefore cannot beat the half-rate threshold at any such seed
size. The bound does not settle sizes 2, 3, or 4:

| Seed `d` | Forced `rho_d` | Half threshold excluded by this bound |
| --- | ---: | --- |
| 2 | 0.5 | no |
| 3 | 0.630929753571... | no |
| 4 | 0.701838730514... | no |
| 5 | 0.744948051202... | yes |
| 6 | 0.773705614469... | yes |

This is a restricted construction barrier, not a factoring lower bound or
an improved factoring exponent. It concerns the full word domain; arbitrary
binary label restrictions or new local tensors are outside it. The proof
and checks used substantial OpenAI assistant support, without new Lean
formalization or independent mathematical review.

## 1. Normalize the local matrices

Use the row arrangement as a matrix on `F^d tensor F^d`. Let `tau` swap the
two basis coordinates and let `w=sum_i e_i tensor e_i`. Its matrix pairings
are `E=tau`, `C=I`, and `R=ww^T`. Right multiplication by `tau` preserves
rank. The normalized pencil is

```math
M_row=rR+h(I+tau+R)+c tau.
\tag{3}
```

The column arrangement swaps `R,C` and fixes `E,H`. Thus its normalized
pencil is obtained from (3) by exchanging the scalar variables `r,c`.

For every off-diagonal unordered pair `i<j`, use the ordered basis

```math
u=e_i tensor e_j+e_j tensor e_i,\qquad v=e_i tensor e_j.
```

On that two-dimensional block, `tau=I+N`, where `Nu=0`, `Nv=u`, and `N^2=0`.
Also `R=0`. There are exactly `L` such independent blocks. Hence their two
pencils are

```math
J_row=cI+(h+c)N,\qquad J_col=rI+(h+r)N.
\tag{4}
```

These are actual constant basis changes over `F_2`, so they remain valid
over every characteristic-two extension.

## 2. The diagonal sector, including even seed sizes

On diagonal coordinates, `tau=I` and `R=ww^T`. If `d` is odd, `R` is an
idempotent. The vector `w` gives one scalar block

```math
D=r+h+c
```

in both cuts. Its trace-zero complement has `d-1` scalar blocks `c` in the
row cut and `r` in the column cut.

If `d` is even, `R` is rank-one nilpotent. On the ordered pair `(w,e_d tensor e_d)`
it is `N`, with square zero. The remaining `d-2` diagonal vectors may be
chosen in its kernel. The special two-dimensional block is

```math
K_row=cI+(r+h)N,\qquad K_col=rI+(c+h)N.
\tag{5}
```

The other diagonal blocks again have scalar pencils `c` and `r`.
The checker constructs both bases explicitly and checks every entry of
all three coefficient matrices in both arrangements. No division by `d`
is used in the even case.

## 3. A common functional on many repeated blocks

Choose `a` special blocks from Section 2 and `m=t-a` off-diagonal blocks
from Section 1, with none of the remaining scalar blocks. In the tensor
power decomposition, this gives

```math
\binom ta L^m
```

direct-sum copies of the same polynomial matrix, up to reordering its
basis factors. These are separate matrix blocks: their multiplicities are
ordinary integer counts, not scalars reduced modulo the characteristic.

The sum of the entries of the first row of one special block is `D` in
either cut. For odd `d`, the block is the scalar `D`. For even `d`, its first
rows are `(c,r+h)` and `(r,c+h)`, whose sums are both `D`.
Likewise the first-row sums in (4) are

```math
c+(h+c)=h,\qquad r+(h+r)=h.
```

Consequently, the first-row sum in each selected tensor-product block is
the same polynomial in both cuts:

```math
p_a(r,h,c)=D^a h^m.
\tag{6}
```

For a linear combination of coefficient tensors, pairing its coefficient
weights with the coefficients of (6) is thus an actual matrix-entry sum
in each arrangement. It defines one well-defined linear functional `ell_a`
on `V_t`, even if the coefficient tensors have linear relations.
Its value on the crossing is

```math
ell_a(E^{tensor t})=p_a(1,1,1)=1.
```

Some summand of every exact crossing decomposition has nonzero `ell_a`.
Every one of its selected row blocks is nonzero, and so is every selected
column block. Each has rank at least one. Direct-sum additivity forces
both ranks of that same summand to be at least `binom(t,a)L^m`.
Choose an `a` maximizing this number. The integer binomial identity
`sum_a binom(t,a)L^(t-a)=(L+1)^t` proves (1).

## 4. The width consequence and the seed-five boundary

Both factorization widths are at least the binary logarithms of their
matrix ranks. Taking logarithms in (1), their sum is at least

```math
2t log2(L+1)-2log2(t+1).
```

This gives (2). Polynomial piece counts or logarithmic width allowances
cannot change that leading restriction.

For `d=5`, exact rational logarithm bounds certify
`log2(11)/(2log2(5))>1/sqrt(2)`. For every integer `d>=6`,

```math
(d^2-d+2)^2-4d^3
=d^3(d-6)+5d^2-4d+4>0.
```

Hence `L+1>d^(3/2)`, so `rho_d>3/4>1/sqrt(2)`. This proves the claimed
exclusion for all seed dimensions at least five, without a floating-point
comparison or an unchecked monotonicity assertion.

A common polynomial substitution for `r,h,c` only takes linear combinations
of the original type coefficients, and cannot leave `V_t`. A coordinate
symmetry preserving `E` only exchanges `R,C` and fixes `H`, so it also
preserves `V_t`. Neither operation evades the theorem. Arbitrary combinations
of separately labeled nonsymmetric tensor words are a larger space and are
not covered by this type-span assertion.

## 5. Exact replay and the remaining small seeds

```sh
python3 tools/check_pairing_word_barrier.py
python3 tools/check_pairing_word_barrier.py --output evidence/pairing-word-barrier.json
make verify
```

The finite checker verifies 151,992 entries of the explicit local block
identities for seed sizes 2 through 10, including even nilpotent diagonal
sectors. It checks 88 common first-row polynomial identities, 960 integer
multiplicity estimates, exact logarithm intervals, and the polynomial
inequality used for all larger seed sizes. The all-size rank statement
follows from the direct-sum proof above.

Seeds 2, 3, and 4 remain open under this particular bound. For example,
`rho_4=log2(7)/4` is slightly below the half threshold. A favorable one-seed
rank does not settle its tensor-power coefficients: their ranks and
collective recombinations still require a uniform theorem. The separate
[characteristic-zero word barrier](word-symmetry-barrier.md) cannot be
substituted for that missing characteristic-two argument.
