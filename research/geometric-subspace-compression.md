# Universal subspace compression from integer Jacobian certificates

## Result

Let `P_1,...,P_m` be arbitrary endomorphisms of a `d`-dimensional
vector space. Work over an algebraic closure of the coefficient field.
The certificates below guarantee subspaces

```math
\dim U=s,\qquad \operatorname{codim}W=q,
\qquad P_iU\subseteq W\quad(1\le i\le m).             (1)
```

There are two explicit families.

1. **Split-map certificates.** Choose positive integers `a,b` with
   `a+b=m`. Condition (1) holds whenever

   ```math
   s+a q\le d,\qquad b s+q\le d.                    (2)
   ```

2. **Three-map certificates and extensions.** For `m=3+a+b`, with
   `a,b` nonnegative integers, condition (1) holds whenever

   ```math
   s+(a+1)q+\lceil q/2\rceil\le d,\qquad
   (b+1)s+q+\lceil s/2\rceil\le d.                  (3)
   ```

For three maps, this reduces to

```math
s+q+\max(\lceil s/2\rceil,\lceil q/2\rceil)\le d.    (4)
```

All these statements hold in every characteristic. Their certificates
are signed permutation minors or four explicitly displayed integer
minors with determinant `1` or `-1`. They require no assumption about
the generic classification of matrix tuples.

For a product sum `B=sum_(i=1)^m P_i tensor Q_i`, equation (1) implies

```math
\operatorname{rank}(E_d-B)\ge s q.                  (5)
```

The same geometric bounds survive a formal degeneration of the
crossing, even when the individual summands have poles. This last
statement does not import the separate exact-only bound
`d(d-1)/2` for one or two product terms.

These are lower bounds on candidate crossing decompositions. They do
not improve the factoring exponent, which remains `1/2`.

## 1. An explicit algebraic certificate for universality

Put `ell=d-q` and `p=d-s`. A standard pair of subspaces has
`U=K^s` in the input and `W=K^ell` in the output. A map carrying
`U` into `W` has blocks

```math
P_i=\begin{pmatrix}A_i&B_i\\0&D_i\end{pmatrix},
\qquad A_i\in\operatorname{Mat}_{\ell\times s},
\quad D_i\in\operatorname{Mat}_{q\times p}.           (6)
```

Define the integer-linear Jacobian map

```math
\mathcal J_{A,D}:\operatorname{Mat}_{p\times s}
  \oplus\operatorname{Mat}_{q\times\ell}
 \longrightarrow\bigoplus_{i=1}^m\operatorname{Mat}_{q\times s},
\qquad (X,Y)\longmapsto(D_iX-YA_i)_{i=1}^m.           (7)
```

**Certificate lemma.** If integer matrices `A_i,D_i` make (7) have
a full-row square minor with determinant `1` or `-1`, then every
`m`-tuple of `d` by `d` matrices over every algebraically closed
field admits (1).

Here is a chart proof, including the passage from a generic tuple to
every tuple. Write moving subspaces as graphs of `X` and `Y`. For a
general block matrix `P_i=[[A_i,B_i],[C_i,D_i]]`, containment is the
equation

```math
C_i+D_iX-YA_i-YB_iX=0.                              (8)
```

Thus a chart of the incidence variety is affine space with free
coordinates `A_i,B_i,D_i,X,Y`, and

```math
C_i=YA_i+YB_iX-D_iX.                                (9)
```

At `X=Y=0`, the differential of the projection to all matrix entries
is the identity on the `A,B,D` coordinates; its remaining part is
the negative of (7). A unimodular minor therefore gives a nonzero
full-row Jacobian minor of this polynomial projection in every
characteristic.

This implies dominance. To see this directly, suppose there were a
nonzero polynomial relation of minimum degree among the output
coordinates. Differentiating and using full Jacobian rank over the
function field forces every composed partial derivative to vanish.
Any nonzero partial derivative would therefore be a relation of
smaller degree. If all partial derivatives were the zero polynomial
in positive characteristic, perfection of an algebraically closed
field would make the relation a power of a lower-degree relation.
Both cases contradict minimality.

The full incidence variety is the closed subset

```math
\{(P,U,W):P_iU\subseteq W\text{ for all }i\}
 \subseteq (\operatorname{Mat}_d)^m
        \times\operatorname{Gr}(s,d)\times\operatorname{Gr}(\ell,d).
```

The product of Grassmannians is projective. Consequently its
projection to the affine matrix-tuple space is closed. The image is
both dense, by the chart calculation, and closed; it is the entire
matrix-tuple space. This proves the certificate lemma. The only
geometric facts used are the usual Grassmannian charts and the closed
projection property of a projective variety.

For the background facts, the [Stacks Project's Grassmannian
construction](https://stacks.math.columbia.edu/tag/089R) gives the affine
charts used here. The Pluecker embedding realizes each Grassmannian as
a closed projective variety. A projective morphism is proper by
[Stacks Project, Section 29.44, Lemma 29.44.5](https://stacks.math.columbia.edu/tag/01W7),
and proper morphisms are universally closed by
[Section 29.42, Definition 29.42.1](https://stacks.math.columbia.edu/tag/01W0).
This is the closed-projection implication in the preceding paragraph.

Subspaces need not descend to the original field. This causes no
problem in (5), because matrix rank is unchanged by extending a field.

## 2. Split-map certificates

Assume (2). The `p=d-s` rows of `X` contain `a` disjoint blocks of
`q` rows. For the first `a` maps, set `A_i=0` and let `D_i` select
one such block. Their outputs in (7) are precisely those `a` blocks
of `X`.

Similarly, the `ell=d-q` columns of `Y` contain `b` disjoint blocks
of `s` columns. For the remaining `b` maps, set `D_i=0` and let
`A_i` select one such block. Their outputs are the negatives of
those `b` blocks of `Y`.

Choose the corresponding `m q s` input coordinates. The resulting
Jacobian minor is a signed permutation matrix, hence unimodular.
The certificate lemma proves (2).

The case `a=1,b=m-1` gives

```math
\operatorname{rank}(E_d-B)
 \ge\max_{s\ge0}s\bigl(d-(m-1)s\bigr),              (10)
```

whenever `m>=2`. Thus the split family includes the earlier
two-map-flag quotient. Other allocations can be stronger.

For instance, with `a,b>=2`, set

```math
d_0=ab-1,\qquad s_0=a-1,\qquad q_0=b-1.
```

Both inequalities in (2) are equalities. Scaling these dimensions
gives the asymptotic complement coefficient

```math
\frac{(a-1)(b-1)}{(ab-1)^2}.                         (11)
```

No scaling or divisibility assumption is needed to use the exact
integer inequalities (2).

## 3. Four compatible three-map minors

There are two types of `A` blocks. Type `A_1` has input dimension
one and output dimension two; type `A_2` has input dimension two
and output dimension three. Their three matrices are

```math
A_1:
\quad\begin{pmatrix}0\\1\end{pmatrix},
\quad\begin{pmatrix}1\\1\end{pmatrix},
\quad\begin{pmatrix}1\\1\end{pmatrix};
```

```math
A_2:
\quad\begin{pmatrix}0&0\\1&1\\0&1\end{pmatrix},
\quad\begin{pmatrix}0&0\\1&0\\1&1\end{pmatrix},
\quad\begin{pmatrix}0&1\\0&0\\0&1\end{pmatrix}.       (12)
```

There are two corresponding `D` types. Type `D_1` has input
dimension two and output dimension one; type `D_2` has input
dimension three and output dimension two:

```math
D_1:
\quad\begin{pmatrix}1&0\end{pmatrix},
\quad\begin{pmatrix}1&1\end{pmatrix},
\quad\begin{pmatrix}1&0\end{pmatrix};
```

```math
D_2:
\quad\begin{pmatrix}0&1&0\\1&0&0\end{pmatrix},
\quad\begin{pmatrix}0&1&0\\0&0&1\end{pmatrix},
\quad\begin{pmatrix}0&0&1\\1&1&1\end{pmatrix}.         (13)
```

Order the input coordinates of (7) by all entries of `X` in row
order, followed by all entries of `Y` in row order. Order outputs
by map index, then row, then column. For each pair `(A_i,D_j)`,
the minor formed by its first `3ij` columns has the following
determinant:

| `A` type | `D` type | Full Jacobian size | Selected minor | Determinant |
|---|---|---|---|---|
| `A_1` | `D_1` | `3` by `4` | first `3` columns | `-1` |
| `A_2` | `D_1` | `6` by `7` | first `6` columns | `1` |
| `A_1` | `D_2` | `6` by `7` | first `6` columns | `1` |
| `A_2` | `D_2` | `12` by `12` | all columns | `1` |

These are finite integer certificates. The replay reconstructs each
matrix from (12)--(13), checks the determinant by exact elimination,
and checks an integer inverse of each selected square minor.

To obtain arbitrary input dimension `s`, take `floor(s/2)` blocks
of type `A_2` and, if needed, one `A_1` block. Their total output
dimension is `s+ceil(s/2)`. Independently, take enough `D_2,D_1`
blocks for output dimension `q`. Their total input dimension is
`q+ceil(q/2)`.

Use these block sums for all three maps. The variables `X` and `Y`
split into rectangular blocks indexed by one `A` block and one
`D` block. Equation (7) acts separately on each such pair. Selecting
the certified minor in every block gives a block-diagonal unimodular
minor. Extra dimensions are filled with zeros. This is possible
precisely under (4), proving the three-map claim.

In particular, every triple of five by five matrices admits
`dim U=2,codim W=2`, giving complement rank at least four. At
dimension nine, `s=3,q=4` gives complement rank at least twelve.

## 4. Adding maps to the three-map system

Starting from the block system in Section 3, add `a` maps of the
`D`-only type and `b` maps of the `A`-only type. The former use
`a q` extra rows of `X`; the latter use `b s` extra columns of `Y`.
Their selected variables are disjoint from the three-map core.

The required input-complement and output dimensions are

```math
d-s\ge q+\lceil q/2\rceil+a q,
\qquad
d-q\ge s+\lceil s/2\rceil+b s.
```

These are exactly (3). The combined Jacobian minor is the direct
sum of the three-map core minors and signed identity blocks. It
is unimodular, so the certificate lemma applies.

Ignoring only integer roundings, the extra family supplies split
allocations with coefficients `a+3/2` and `b+3/2`, whose sum is
`m`. The exact inequalities (3), rather than a continuous
approximation, are used in all finite bounds.

## 5. The complement-rank consequence

Let `B=sum_i P_i tensor Q_i`, and choose `U,W` from (1). Write
`pi:V -> V/W` for the quotient and restrict the input to
`U tensor V`. Then

```math
(\pi\otimes I)B\big|_{U\otimes V}=0.
```

The crossing swaps its input factors. Hence

```math
(\pi\otimes I)E_d(u\otimes v)=\pi(v)\otimes u,
```

whose rank is exactly `q s`. Restriction and projection cannot
increase rank, proving (5).

Any matrix of opposite-cut rank at most `m` has a representation
as a sum of at most `m` products. Padding with zero terms allows
these bounds to be applied at any larger product count as well.

Define `F_geo(d,m)` as the maximum `s q` over (2) and (3), or zero
if there is no positive feasible pair. One may also take the
monotone envelope `max_(k>=m) F_geo(d,k)`. The latter remains valid
by zero padding; it is not an assumption about cancellation in a
chosen product representation.

Equivalently, for `m>=2` this maximum ranges over integers
`2<=A<=2m-2`, `B=2m-A`, and positive `s,q` satisfying

```math
2s+Aq\le2d,\qquad Bs+2q\le2d.                       (14)
```

Even `A,B` give (2); odd `A,B` give (3), including its roundings.
Swapping `A,B` and `s,q` permits the smaller range `2<=A<=m`.
The checker exposes this exact optimized function as
`geometric_complement(d,m)`. The exceptional exact bounds for
`m=0,1,2` are deliberately left to the calling theorem.

The finite ranges differ from the older flag-only formula. The
split family can stay positive through `m=2d-2`, while both
families vanish for `m>=2d-1`. Indeed, a positive pair in (2)
forces `a,b<=d-1`; a positive pair in (3) forces
`a,b<=d-3`. Thus a weighted minimum using only these geometric
complement bounds can stop at `max(m_min,2d-1)`. The earlier
cutoff `d+1` must not be reused.

For exact crossings, the separate one- and two-product bound
`d(d-1)/2` can be combined with `F_geo`. That stronger small-count
bound has a different proof and is not asserted for formal limits.

## 6. Formal degeneration with moving subspaces

Suppose `T(t)` has entries in `F[[t]]`, with `T(0)=E_d`, and

```math
T(t)=A(t)+B(t),\qquad
B(t)=\sum_{i=1}^m P_i(t)\otimes Q_i(t),              (15)
```

where all individual factors and pieces may have Laurent-series
poles. Apply a geometric compression certificate over an algebraic
closure of `F((t))`. Its finite collection of subspace coordinates
lies in some finite algebraic extension `L` of `F((t))`.

Choose the extension of the `t`-adic valuation to `L`, and let
`O` be its discrete valuation ring. The lattices

```math
U_O=U_L\cap O^d,\qquad W_O=W_L\cap O^d
```

are saturated submodules of a finite free module over a DVR.
They are direct summands, with ranks `s` and `d-q`. Choose an
integral inclusion `iota:O^s -> O^d` with image `U_O` and an
integral quotient `pi:O^d -> O^q` with kernel `W_O`. Their
reductions have full ranks `s` and `q`.

All of `B`, including its poles, is killed exactly:

```math
(\pi\otimes I)B(\iota\otimes I)=0.
```

The same projection of `A` equals that of the integral tensor
`T`. Reducing this integral matrix modulo the maximal ideal gives

```math
(\pi_0\otimes I)E_d(\iota_0\otimes I),
```

of rank `s q`. Therefore the generic rank of `A` is at least
`s q`. Movement of both subspaces is harmless because their
saturated lattices keep their dimensions after reduction.

This proves the formal version of (5) for the geometric
certificates. It does not use a specialization of an unprojected
Laurent matrix, where poles could invalidate the argument.

## 7. Replay

Run

```sh
python tools/check_geometric_subspace_compression.py \
  --output evidence/geometric-subspace-compression.json
```

The standard-library replay reconstructs the four basic integer
minors and their inverses, builds larger compatible block systems,
checks split and extended certificates, verifies the polynomial
chart identity, and compares the optimized integer bounds with
brute-force feasible-pair enumeration. Finite-field complement
tests are included as controls. The chart and closed-projection
proof establishes universality; random tests do not replace it.

The [flag quotient note](flagged-cut-rank-budget.md) gives the earlier
subset-budget application. The geometric construction and exact replay
used substantial OpenAI assistant support.
