# Odd full-conservation pieces cost at least q

## Result and scope

Let `q >= 3` be odd and work with labels in `Z/qZ`. For a scalar kernel
`K:(Z/qZ)^2 -> F`, define

```math
T_K(x,y,u,v)=[x+y=u+v]\,K(x+y,x+u).
\tag{1}
```

The two standard cuts are `(x,y)|(u,v)` and `(x,u)|(y,v)`. Their
ordinary ranks are the row and column convolution-rank sums

```math
r(T_K)=\sum_{S\in G}\operatorname{rank}C_{K(S,\cdot)},\qquad
c(T_K)=\sum_{w\in G}\operatorname{rank}C_{K(\cdot,w)},
\tag{2}
```

where `C_f(a,b)=f(a+b)`. The identity is `T_I` with `I(S,w)=[S=w]`.

Every nonzero piece has

```math
\sqrt{r(T_K)c(T_K)}\ \ge\ q.
\tag{3}
```

Consequently an exact decomposition `I=sum_i K_i` with `k` nonzero
pieces has rank-sum cost at least `kq`. At the finite-seed threshold
`q^sqrt(2)` this leaves:

| q | target | forced number of nonzero pieces |
|---:|---:|---|
| 5 | 9.7385 | 1, but the only one-piece identity costs `q^2` |
| 7 | 15.6729 | exactly 2 |
| 11 | 29.6997 | exactly 2 |
| 13 | 37.6144 | exactly 2 |
| 17 | 54.9690 | 2 or 3 |
| 19 | 64.3325 | 2 or 3 |

In particular, **the q=13 general free-conservation search reduces to
two pieces**. This is a reduction of the search space, not an
odd-alphabet rank-sum obstruction. At q=17 and q=19 three pieces remain
possible; for larger q the floor allows up to about `q^(sqrt(2)-1)`
pieces.

**Scope correction.** The two-piece reduction is not the decisive
statement for the fixed conservation family. The identity `T_I` is the
crossing tensor `E_q`, and the
[minimum-cut rank budget](minimum-cut-rank-budget.md) proves, for every
field and every exact decomposition of `E_d`,
`sum sqrt(r_i c_i) > d^sqrt(2)` for all `2 <= d <= 21`. Therefore the
whole conservation-kernel family is already closed for `q <= 21`, with
arbitrary piece count, independently of the floor below. The first
unexcluded conservation alphabet is `q=23`; the q=13 two-piece search is
closed. See [§3](#3-piece-count-reduction) for the exact table and the
remaining q=23 scope.

The characteristic-two theorem in
[the full-conservation barrier](full-conservation-rank-sum-barrier.md)
uses the Boolean degree filtration. That argument does not transfer to
odd q: the character basis is a basis of convolution rank one. The
present note therefore proves only the per-piece floor and records a
finite two-piece control. It makes no factoring-exponent or complete-bound
claim.

## 1. The support lemma

For a nonzero function `f:Z/qZ -> F` with `t` nonzero entries, let
`C_f(a,b)=f(a+b)`.

**Support lemma.** `rank C_f >= q/t`.

Every row of `C_f` is a translate of `f`, and the row space is closed
under translation. Choose `r=rank C_f` rows spanning it. If a column is
zero in all chosen rows, it is zero in every translate of `f`, hence in
every column of `C_f`. But every column is a nonzero translate of `f`.
Thus the supports of the chosen rows cover all `q` columns. Each has
size `t`, so `r t >= q`. The argument is valid over every coefficient
field and does not need a `q`-th root of unity.

The same lemma was used in the
[incidence-kernel barrier](incidence-kernel-rank-sum-barrier.md).

## 2. A nonzero piece has cost at least q

For a nonzero kernel `K`, let:

- `u` be its number of nonzero rows,
- `v` its number of nonzero columns,
- `m` its total number of nonzero entries,
- `t_S` the support size of row `K(S,\cdot)`,
- `s_w` the support size of column `K(\cdot,w)`.

Rows of zero support contribute nothing to `r(T_K)`. By (2) and the
support lemma,

```math
r(T_K)\ge q\sum_{S:t_S>0}\frac1{t_S}
       \ge q\frac{u^2}{\sum_S t_S}
       =q\frac{u^2}{m}.
```

The second inequality is Cauchy--Schwarz. The column calculation gives

```math
c(T_K)\ge q\frac{v^2}{m}.
```

Since `m <= u v`, because the nonzero entries lie in the `u` by `v`
incidence rectangle,

```math
r(T_K)c(T_K)
\ge q^2\frac{u^2v^2}{m^2}
\ge q^2.
```

Taking square roots proves (3).

This is a floor for every individual nonzero piece. It does not bound a
sum of pieces below by anything stronger than `kq`; cancellation between
pieces is not used.

## 3. Piece-count reduction

Suppose `I=sum_i K_i`, discard zero pieces, and let `k` be the number of
remaining pieces. By (3),

```math
\sum_i\sqrt{r(T_{K_i})c(T_{K_i})}\ge kq.
\tag{4}
```

For `q=7,11,13` the target `q^sqrt(2)` lies strictly between `2q` and
`3q`. A one-piece decomposition is `K_1=I`, whose cost is `q^2`, above
the target. Hence a crossing decomposition, if one exists, must have
exactly two pieces.

For `q=5` the target is below `2q`, so only the one-piece identity is
possible, and it fails. For `q=17,19` the target lies between `3q` and
`4q`, so at most three pieces are possible. At `q=23` the target
`84.2896` lies between `3q=69` and `4q=92`, so at most three pieces are
possible.

The per-piece floor is only a reduction. For every `q <= 21`, the
stronger [minimum-cut rank budget](minimum-cut-rank-budget.md) closes the
family outright: it excludes every exact decomposition of `E_q` for all
`2 <= q <= 21`, with no piece-count restriction. Thus the q=7,11,13
two-piece reductions and the q=17,19 three-piece windows are superseded
as crossing questions. The two-piece catalogue below remains a finite
control, not an open search target.

The script
[`check_full_conservation_piece_count.py`](../tools/check_full_conservation_piece_count.py)
recomputes the table, checks the cut-rank identities (2) against dense
tensors at q=3,5,7, and checks (3) on all nonzero q=3 kernels over `F_2`
and `F_5`, plus random kernels at q=5,7,13. It also verifies the exact
piece-count arithmetic through q=29.

## 4. Finite two-piece catalogue

If a q=13 crossing has exactly two pieces, write it as

```math
I=K+(I-K).
```

The remaining question is a genuine two-piece classification, not a
one-piece search. The checker records a finite catalogue of natural
kernels:

- one nonzero entry;
- one row or one column supported on a character;
- a product of two characters;
- a chirp `K(S,w)=omega^(aSw)`;
- a cyclic permutation;
- diagonal masks, including intervals and random subsets.

For every catalogue element at q=7,11,13,17,19, the total
`cost(K)+cost(I-K)` is at least `q^2`. At each q the minimum in this
catalogue is attained by a single entry and equals exactly `q^2`:
`49,121,169,289,361`. These are finite controls, not a two-piece
theorem. In particular, the catalogue does not cover arbitrary dense
kernels, mixed supports with cancellation, or nonlinear recombinations.

The recorded q=13 minimum is `169`, compared with target `37.6144`.
The trivial entry-plus-complement pair already shows that the coarse
floor `2q` is far from the observed catalogue scale, but it does not
prove a two-piece lower bound. The minimum-cut budget now supplies the
required closure through `q=21`; a structural two-piece classification
is no longer on the critical path.

## 5. Verification

```sh
python3 tools/check_full_conservation_piece_count.py \
  --output evidence/full-conservation-piece-count.json
```

The evidence records 20 dense cut-rank comparisons, 1,954,315
piece-floor checks, the piece-count table, and the five finite
two-piece catalogues. A pass checks the implementation and the stated
finite controls. It does not establish an odd-alphabet obstruction or
a stronger complete factoring exponent.

The closure through `q=21` is the separate, already-checked
[minimum-cut certificate](minimum-cut-rank-budget.md):

```sh
python3 tools/check_minimum_cut_rank_budget.py \
  --output evidence/minimum-cut-rank-budget.json
```

That checker verifies 71,354 rank profiles, 577,536 feasible dual
vertices, 498,130 unrounded endpoint inequalities, and the exact
comparison `L^408 > d^577` for every `2 <= d <= 21`; the checked
run reports PASS and a first uncertified profile at `d=22`. The identity
`T_I=E_q` is immediate from (1): `[x+y=u+v][y=u]=[x=v][y=u]`. No q=22
or q=23 conservation closure is claimed here.
