# The full three-pairing type route is obstructed over every field

## Result and scope

Fix a field `F`, local dimension `d>=2`, and word length `t`. Set

```math
E=[x=v][y=u],\quad R=[x=y][u=v],\quad C=[x=u][y=v],
\qquad H=E-R-C.
```

Let `V_t` be the full coefficient span of `(rR+hH+cC)^tensor t` on
`[d]^t`. It is closed under the eight crossing symmetries. In every
crossing decomposition inside this space, some piece has both ranks
large enough to force a leading contraction strictly above `1/sqrt(2)`.
This includes arbitrary collective recombinations and any piece count.

Together with the [characteristic-two result](pairing-small-seeds.md),
the cases below exclude the entire three-pairing type route over **every
field**, including extension and function fields. A universal lower
leading contraction is

```math
gamma=\frac{log2(32)-(17/32)log2(17)}{log2(16)}
     =0.707133841396\ldots>1/sqrt(2).                    (1)
```

The seed dimension may grow: all `d>=5` estimates below are uniform;
only the three fixed small dimensions need logarithmic losses. This is
not a factoring lower bound. New local tensors, nonproduct label
relabelings, and arbitrary domain restrictions remain outside it. In
particular it does not invalidate the new relabeling seeds.

| Field/seed case | Forced leading `rho` |
| --- | ---: |
| Characteristic two, all `d` | at least 0.725919365257... |
| Other characteristics, `d=2` | 0.774397470347... |
| Other characteristics, `d=3`, except characteristic three | 0.733269522816... |
| Characteristic three, `d=3` | 0.711560037613... |
| Other characteristics, `d=4` | 0.707133841396... |
| Other characteristics, `d=5` | 0.715338279036... |
| Other characteristics, `d>=6` | greater than 3/4 |

Substantial OpenAI assistant support was used. These ordinary arguments
await independent review and add no Lean formalization. The strongest
overall factoring exponent remains `1/2`.

## 1. Local scalar and nilpotent sectors in odd characteristic

Assume `char F !=2`. Normalize the row matrix by right multiplication
with `E`. On `F^d tensor F^d`, write `tau` for swap and `w=sum_i e_ii`.
The matrices are `R=ww^T`, `C=tau`, and `H=I-R-tau`. The column cut
exchanges `r,c`.

The antisymmetric sector has dimension `A=d(d-1)/2`, with scalar pencil
`2h-c`. The symmetric sector has dimension `S=d(d+1)/2`.
If `d` is nonzero in `F`, it splits into a special line with scalar
`d(r-h)+c` and `B=S-1` scalar copies of `c`. If `d=0` in `F`, the rank-one
operator `R` is square-zero: there is one two-dimensional block
`cI+(r-h)N`, `N^2=0`, and `B=S-2` scalar copies of `c`.
These are constant basis changes. The checker verifies them directly in
several characteristics, including the nilpotent cases.

## 2. Every seed of size at least five

Ignore the special or nilpotent block. There remain `B>=A` scalar copies
of `c` and `A` scalar copies of `2h-c`. The common polynomial

```math
h^t=2^{-t}\bigl(c+(2h-c)\bigr)^t
```

is a linear combination of their tensor-type entry polynomials. The
column identity replaces `c` with `r`. Coefficient pairing with `h^t`
therefore defines the same functional in both cuts, with value one on
the crossing. A detected piece is nonzero on at least one of these
scalar types in each arrangement. Every such type has multiplicity

```math
binom(t,j) B^{t-j}A^j>=A^t.
```

This forces both ranks of the same piece to be at least `A^t`.
Vanishing binomial coefficients in the field only remove terms from the
functional expansion; they do not reduce integer block multiplicities.

For `d=5`, exact logarithm intervals show
`log2(10)/(2log2(5))>1/sqrt(2)`. For `d>=6`, `(d-1)^2>4d` implies
`A>d^(3/2)`, giving `rho>3/4`. No assumption that the characteristic
does not divide `d` was used here.

## 3. Small seeds with a scalar special sector

Suppose `d` is nonzero in `F`. Put `N=d^2` and `B=S-1`. The three local
eigenforms are independent, with multiplicities `(1,B,A)`.
For `t=2Nk`, use the homogeneous functional polynomial

```math
p=det(M_row)^k det(M_col)^k.
```

It evaluates to one on the crossing. In row eigenform coordinates every
supported monomial has exponents at least `(k,Bk,Ak)`; the same is true
in column coordinates. A detected piece has both ranks at least

```math
J_{N,k}=\frac{(2Nk)!}{((N+1)k)!(k!)^{N-1}}.              (2)
```

To see this integer bound, a type with counts `a,b,c` has multiplicity
`t!/(a!b!c!) B^b A^c`. Split `b` among its `B` physical eigenvectors
and `c` among its `A` eigenvectors, each receiving at least `k`.
The multiplicity dominates the resulting physical-word multinomial.
Among `N` counts each at least `k` summing to `2Nk`, that multinomial is
minimized by counts `((N+1)k,k,...,k)`: transferring excess to one count
maximizes the factorial product. This proves (2).

Elementary factorial asymptotics give

```math
log2 J_{N,k}=t\left(log2(2N)-\frac{N+1}{2N}log2(N+1)\right)
             +O_N(log(t+2)).
```

For general `t`, fix fewer than `2N` positions to one common label.
This preserves the crossing, cannot increase either rank, and sends the
coefficient space into that of the smaller full word family. The leading
bound is unchanged. It gives the first three small-seed entropy rates
in the table. This proof uses only the three **local** scalar sectors,
so it works in positive odd characteristic without a semisimple
symmetric-group representation assumption.

This handles `d=2,4` in every odd characteristic and `d=3` except in
characteristic three. The latter is the only remaining small case.

## 4. A four-position certificate for dimension three in characteristic three

Here the row blocks are one `K=cI+(r-h)N`, four scalars `c`, and three
scalars `beta=2h-c`. Set `alpha=r-h`. A monomial entry functional

```math
c^{t-k-j}alpha^k beta^j
```

is present in every block type with exactly `j` antisymmetric positions
and at least `k` special positions. If it is nonzero on a piece, its
row rank is at least

```math
R(t,k,j)=binom(t,j)3^j
 \sum_{a=k}^{t-j}binom(t-j,a)2^{a-k}4^{t-j-a}.            (3)
```

For a special-block count `a`, the square-zero multiplication matrix has
a nonzero monomial of degree `k`, forcing rank at least `2^(a-k)` by a
triangular minor. Summing every compatible block gives (3).

The common polynomial `p=r c h^2` has degree four, is invariant under
`r,c` exchange, and equals one at `r=h=c=1`. In the row variables,

```math
p=2c^4+2c beta^3+c^3 alpha+2c^2 alpha beta+c alpha beta^2.
```

The five values of (3) are respectively

```text
1296, 648, 520, 912, 540.
```

Thus a piece detected by this one functional has both ranks at least
520. The same certificate repeats: compatible block choices on disjoint
position sets inject into the global choices, giving

```math
R(t_1+t_2,k_1+k_2,j_1+j_2)
 >=R(t_1,k_1,j_1)R(t_2,k_2,j_2).
```

Every monomial in `p^m` therefore has floor at least `520^m`. Its value
on the crossing remains one. Fixing at most three leftover positions
proves both ranks at least `520^floor(t/4)`. The leading contraction is
`log2(520)/(8log2(3))=0.711560037613...>1/sqrt(2)`.
All steps hold over every extension field of characteristic three.

## 5. The combined conclusion and replay

The smallest table entry is the odd-characteristic dimension-four value
in (1); exact rational bounds verify its strict half-threshold gap.
The large-seed estimates are uniform, so growing the seed does not
evade this full-type obstruction. The usual binary prefix for a **fixed**
seed also contains a full `[d]^(t-1)` fiber; fixing its remaining position
preserves the leading lower bound. Arbitrary restrictions need separate
analysis.

```sh
python3 tools/check_all_field_pairing.py
python3 tools/check_all_field_pairing.py --output evidence/all-field-pairing-barrier.json
make verify
```

The checker replays local scalar/nilpotent bases, the explicit
characteristic-three polynomial and block floors, their repetition
inequalities, 64 independent global characteristic-three rank cases,
integer determinant-sector bounds, and exact logarithm
intervals. The all-size conclusion is the written argument above.
