# Mixed-field partial repair of the rectangular residual core

## Result and scope

The four-by-three residual core admits an explicit height-triangular repair
over every characteristic-two field whenever

```math
\delta=1+\alpha\beta\ne0.
```

The repair preserves the all-low Gaussian background. Its rational local
tensor has denominator at most `1+w`, uniformly over all finite mode sets
allowed by the construction. Independent matching types on the two
crossing pairs cancel this possible pole and give polynomial crossing
tensors of every admissible word length.

In the four-element field, `alpha=1` and `beta=omega` activate both mixed
directions while satisfying the condition. The accompanying exact
constructor supplies four finite examples, including alphabets `(16,12)`.
This establishes a construction and its endpoint obligations. It does not
establish a better rank-pressure bound or a factoring exponent above `1/2`.

The predecessor is [rectangular residual repair](rectangular-residual-repair.md),
especially Sections 5 and 6. Its coefficient matrices, coordinate maps,
and ramified Gaussian proof are retained. The contribution here is the
explicit two-parameter partial repair, its uniform simple-pole bound,
and the exact mixed-field constructor. The initial source is PR #36 at
`dde77eccab94e558aa134621d8f6d7d256bec8bb`.

## 1. The sparse repair

Index the twelve core states by `v+4*y`. The core heights are

```math
h_X=(0,0,1,1),\qquad h_Y=(0,0,1),\qquad
h(i)=h_X(i\bmod4)+h_Y(\lfloor i/4\rfloor).
```

Use the residual matrix `F(alpha,beta)` of the predecessor. The third
input state is `top + alpha*both`; the third output covector is
`top + beta*both`. Both summands have the same index height.

For `delta != 0`, define the rows of `C=C(alpha,beta)` as follows.
Here `e_j^T` is coordinate row `j`; every unspecified coefficient is zero.

```math
\begin{aligned}
C_{i,*}&=e_i^T &&(0\le i\le7),\\
C_{8,*}&=\delta^{-1}e_8^T,\qquad
C_{9,*}=\delta^{-1}e_9^T,\\
C_{10,*}&=\delta(e_0^T+e_7^T+e_{10}^T)
             +\beta\delta e_3^T+\alpha e_9^T,\\
C_{11,*}&=\alpha e_8^T+\delta e_{11}^T.
\end{aligned}
\tag{1}
```

Direct multiplication gives

```math
\begin{aligned}
CF=U={}&I+E_{0,7}+E_{0,10}+E_{1,6}+E_{4,9}+E_{7,10}\\
       &+\alpha(E_{0,9}+E_{3,10})
        +(\beta/\delta)(E_{8,11}+E_{9,10}).
\end{aligned}
\tag{2}
```

Each off-diagonal entry in (2) has `h(row)<h(column)`. Thus every extra
transition strictly lowers height. The determinant of `C` is one: in
height order it is lower triangular with diagonal entries consisting of
eight ones, two `delta^-1` entries, and two `delta` entries. The repair
is therefore invertible.

Equation (2) is a rational-function identity, not a conclusion inferred
from testing finitely many field elements. The checker independently
verifies `delta*C*F=delta*U` entry by entry in the polynomial ring
`F2[alpha,beta]`. At `alpha=beta=0`, (1) reduces to the predecessor's
two-entry partial repair.

For `alpha=1`, `beta=omega`, and `omega^2=omega+1`, one has
`delta=omega^2`, `delta^-1=omega`, and `beta*delta=1`. In particular,

```math
\begin{aligned}
C_{10,*}&=\omega^2(e_0^T+e_7^T+e_{10}^T)+e_3^T+e_9^T,\\
C_{11,*}&=e_8^T+\omega^2e_{11}^T.
\end{aligned}
```

This point is unavailable over `F2`: its only doubly active point has
`alpha*beta=1` and a singular middle principal height block.

## 2. Dressed correction and the all-length endpoint

Let `S` be any finite mode set containing `0` and maximum `T>0`, and let
the power of two `M` be larger than `T`. Retain the one-particle matrix
and the quartic term of [coordinate cores](coordinate-core-crossings.md).
Write `c=1+w`. Before correction the restricted tensor is

```math
K=\mathrm{restrict}_{\alpha,\beta}\bigl(G(1+Q/c^2)\bigr),
\qquad Q=x_0y_Tu_0v_T.
```

On core output pairs put `H_ii=c^h(i)`. Apply

```math
K'=H^{-1}CH K.
\tag{3}
```

The predecessor proves that `w=1+t^T`, together with the actual occupied
index gauge on each coefficient, makes the original tensor regular in
`t` with endpoint `F` on the core and crossing on the bulk modes.
Because the mixed coordinate vectors preserve index height, that proof
holds for the present constants `alpha,beta`. After (3), the regular
endpoint is `U` on the core and the same bulk crossing.

Fix independent word types for the `x,v` pair and the `u,y` pair. The
total index height cancels on each pair. Therefore the restricted tensor
power is regular at `w=1`. Its only possible finite pole was at that
point, so every restricted entry is a polynomial in `w`.

Every term in the endpoint tensor power using an extra entry of (2)
strictly decreases total height. Matching paired types have equal input
and output total heights. Such a term cannot survive. The remaining
endpoint is exactly

```math
[x=v][u=y].
\tag{4}
```

For a fixed local seed, the degree of these polynomial entries grows at
most linearly with word length. Thus this is a genuine all-length
polynomial crossing family. The ordinary effective coefficient and
factorization procedures for a fixed finite seed apply; no seed rank
bound is supplied by the endpoint argument.

The output row for the all-low core pair is unchanged by (1). The
original quartic is absent from that row. Restricting all four core
labels to their low state consequently gives exactly the undeformed
Gaussian background, for every value of `w`, including all bulk entries.

## 3. Uniform denominator at most `c`

The safe height-dressed formula initially gives denominator `c^4`.
In fact, at least three of these factors cancel for every allowed mode
set and all constants with `delta != 0`.

For a fixed pair of input letters and fixed bulk output letters, write
`G_j` for the Gaussian contribution in output core row `j`, and `Q_j`
for its quartic contribution before division by `c^2`. The latter
vanishes except in rows `10` and `11`, since the quartic requires both
top output modes. Formula (3) therefore reads

```math
\begin{aligned}
K'_j&=G_j &&(0\le j\le7),\\
K'_8&=\delta^{-1}G_8,\qquad
K'_9=\delta^{-1}G_9,\\
K'_{10}&=\delta G_{10}
  +\delta c^{-2}(Q_{10}+G_0)
  +c^{-1}(\beta\delta G_3+\delta G_7+\alpha G_9),\\
K'_{11}&=\delta G_{11}+\delta c^{-2}Q_{11}
  +\alpha c^{-1}G_8.
\end{aligned}
\tag{5}
```

At `w=1`, the one-particle matrix is

```math
\begin{pmatrix}0&P\\I&0\end{pmatrix},
\qquad P_{ab}=[b\subseteq a],
```

where the rows are `(v,y)` and columns are `(u,x)`. The low row of `P`
can use only its low column. A row indexed below `T` cannot use column
`T`. The `u`-to-`y` block is the identity. These facts give the exact
residue identities

```math
Q_{10}(1)=G_0(1),\qquad Q_{11}(1)=0.
\tag{6}
```

For the first identity, a surviving `Q_10` minor has the four quartic
variables deleted. Its remaining `u` and `y` occupations must agree.
The mixed both-state contributions would leave an unmatched low or top
`u/y` mode, and hence vanish. On the `x/v` side, a remaining top `x`
column cannot be matched below `T`. The only surviving low-core case
is therefore exactly the bulk minor occurring in `G_0`, whose low
`x/v` and `u/y` matches have coefficient one. This also covers zero
minors for all other input letters.

For the second identity, deleting the quartic from output row `11`
leaves the low `v` row but removes the low `x` column. The low row of
`P` then vanishes. These arguments use only the minimum and maximum
indices, so gaps in `S` do not affect them.

Both numerators in the `c^-2` terms of (5) are divisible by `c` by (6).
Every other term has at most a simple pole. Thus `c*K'` is polynomial,
as claimed. The finite bulk examples below have actual simple poles,
so the local tensor should not be treated as a polynomial before the
paired-type restriction.

For profile calculations, use the numerator after exact cancellation
and its recorded denominator degree. Dividing by `c` adds one to each
infinity order and leaves the zero orders unchanged. Restore that
factor once.

## 4. Exact examples and reproduction

The field is encoded as `0,1,omega,omega+1 = 0,1,2,3`. Each polynomial
is represented by two binary coefficient planes. Tensor entries use
`x+A*(y+B*(u+B*v))`, with alphabets `(A,B,B,A)`.

The following table uses `alpha=1,beta=omega`.

| Modes | `M` | `(A,B)` | Numerator degree | Remaining denominator | Local endpoint checks | Two-copy paired entries |
|---|---:|---|---:|---|---:|---:|
| `0,1` | 2 | `(4,3)` | 8 | `1` | 144 | 420 |
| `0,1,2` | 4 | `(8,6)` | 25 | `c` | 2,304 | 7,920 |
| `0,1,2,3` | 4 | `(16,12)` | 33 | `c` | 36,864 | 136,896 |
| `0,2,5` | 8 | `(8,6)` | 49 | `c` | 2,304 | 7,920 |

Run:

```sh
python tools/check_rectangular_f4_repair.py --full \
  --output evidence/rectangular-f4-repair.json
```

The `--seed-directory` option exports complete coefficient arrays,
logical heights, denominator metadata, and coefficient hashes. The
[recorded evidence](../evidence/rectangular-f4-repair.json) includes:

- The symbolic polynomial identity proving (2) over every
  characteristic-two field.
- Literal exterior endpoint checks for all sixteen F4 parameter pairs;
  all thirteen nonsingular repairs and all three singular rejections.
- 4,096 independent scalar-convolution checks of polynomial
  multiplication in the two-plane representation.
- Every local ramified coefficient and every indicated two-copy paired
  type entry, including exact denominator division.
- Independent literal exterior expansions of all numerator entries for
  the three examples with at most three modes.

The finite checks support the formulas and implementation. Sections 1–3
give the all-length argument. A scoring improvement still requires an
exact pressure or coefficient-rank certificate for these or other seeds,
followed by the existing effective transfer. Fixed-field descent can
be applied after the complete typed tensor power; there is no need to
tensor a conjugate copy merely to remove F4 coefficients.
