# The second coefficient as a compressed transposition class sum

For the [coherent bounded-spin crossing](bounded-spin-crossing.md),
the second row-cut coefficient has a uniform representation-theoretic
formula. It also identifies the primary blocks of the existing
[eight-site row algebra](bounded-spin-rank-algebra.md). It does not
give a favorable all-length cut-rank bound or a new factoring score.
No new Lean formalization is used.

## 1. The operator identity

Work over `Q` or a field of characteristic greater than three, with one
common selected spin-zero path frame whose paired norms are units.
Let `n=2k`. Use the original normalized row cut `xy|uv`, so that
`C_0=I` for

```math
G_n(z)=\sum_{j=0}^n C_jz^j,
\qquad K(z)=I+z(S-I).
```

For each pair of sites `i<j`, let `P_ij` swap those two sites in one
binary tensor copy. Define its compressed matrix by

```math
(T_{ij})_{ab}=\langle v_a,P_{ij}v_b\rangle/\nu_b,
\qquad \mathcal T=\sum_{i<j}T_{ij}\otimes T_{ij}.
```

Then the exact uniform identity is

```math
\boxed{3C_2=\mathcal T+k(k+1)I.}                   \tag{1}
```

The matrix `T_ij` is a compression. The selected path set need not be
stable under each transposition for (1) to hold. If the quotient action
is identified, `mathcal T` is its usual diagonal transposition class sum.
The proof of (1) does not assume that identification.

## 2. The two-site transition calculation

For two spin-zero vectors, contract every site except `i,j` to form
a four-by-four transition matrix `R(v,w)`. Algebraic `SL_2` invariance
makes it commute with the two-site action. In the stated
characteristics the centralizer is spanned by the identity and the
two-site swap `P`. Write `R=aI+bP`. If

```math
g=\operatorname{tr}R=\langle v,w\rangle,\qquad
h=\operatorname{tr}(PR)=\langle v,P_{ij}w\rangle,
```

then `4a+2b=g`, `2a+4b=h`, and hence

```math
a=(2g-h)/6,\qquad b=(2h-g)/6.
```

For another transition matrix with traces `g',h'`, this gives

```math
\operatorname{tr}(RR')
=(gg'+hh')/3-(gh'+hg')/6.                          \tag{2}
```

Swapping the two copies at exactly sites `i,j` contracts these two
transition matrices by (2). After the prescribed input-norm division,
the resulting operator is

```math
\frac13(I\otimes I+T_{ij}\otimes T_{ij})
-\frac16(T_{ij}\otimes I+I\otimes T_{ij}).          \tag{3}
```

The one-site scalar-transition identity from the
[three-jet proof](bounded-spin-three-jets.md) says that a single
copy swap compresses to `I/2`. In the coefficient of `z^2`, the two
single-swap subtractions cancel the added identity for each site pair.
Thus `C_2` is precisely the sum of (3) over `i<j`.

## 3. Evaluate the single-copy class sum

Let `E_ab` be the collective matrix-unit operators on all `n` binary
sites. Direct expansion on the tensor power gives

```math
\sum_{a,b\in\{0,1\}}E_{ab}E_{ba}
=2nI+2\sum_{i<j}P_{ij}.
```

On spin-zero vectors the two diagonal operators are both `kI`, and
the off-diagonal collective raising and lowering operators annihilate
the vector. Consequently

```math
\sum_{i<j}P_{ij}=k(k-2)I,
\qquad \sum_{i<j}T_{ij}=k(k-2)I.                  \tag{4}
```

Compressing preserves this identity for any selected frame. Summing
(3), its identity coefficient is
`(binom(n,2)-k(k-2))/3=k(k+1)/3`. This proves (1).

## 4. Interpret the existing eight-site row certificate

At `n=8` over `F_5`, the scalar term in (1) vanishes. Let `R` be the
global swap of the two path labels. The crossing-normalized row
matrices in the existing certificate are `A_j=C_j R`; therefore

```math
R=A_0,\qquad \mathcal T=3A_2A_0.                  \tag{5}
```

Evaluating (5) in the published polynomial model gives the following
ordinary primary labels. The dimensions and original nilpotent rank
are inherited from that source certificate; the new calculation
identifies their operators.

| Replica swap `R` | Class sum `mathcal T` | Dimension | Nilpotent rank |
| --- | --- | --- | --- |
| `+1` | `4` | 20 | 0 |
| `+1` | `3` | 1 | 0 |
| `-1` | `3` | 43 | 0 |
| `+1` | `2` | 70 | 13 |
| `-1` | `4` | 35 | 0 |

These five joint residue characters are distinct. On the double-root
70-dimensional component, `mathcal T` has derivative four with respect
to the certified algebra generator, so it has a nonzero nilpotent part.
The replica swap has zero derivative there, as follows from `R^2=I`.
The Chinese remainder decomposition therefore shows that `R,mathcal T`
generate the entire six-dimensional certified row algebra, including its
square-zero component. No inference of simultaneous diagonalizability
is used.

This is a finite interpretation of the eight-site model. Formula (1)
does not assert that these two operators generate every larger-length
coefficient algebra, identify the other cut, or improve the existing
finite-block obstruction. Uniform primary multiplicities and same-piece
two-cut ranks remain open.

## 5. Reproduce the exact controls

Run `make verify-bounded-spin-second`. The
[checker](../tools/check_bounded_spin_second.py) constructs every
single-copy transposition from literal integral path vectors and
compares 198 class-sum matrix entries with (4). It compares 104 entries
of (1) at lengths four, six and eight against the separate magnetic-state
coefficient constructor. It then checks the polynomial swap identity,
joint primary labels and nonzero derivative in the published row model.
It does not rebuild the dense rank certificate in that last step.

The [receipt](../evidence/bounded-spin-second.json) binds those controls
to hashes of the existing constructor and rank model. The construction,
norm conventions and finite rank model are credited to the bounded-spin
lane; this lemma and its independent entry controls continue the
[Discussion #51 exchange](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19324321).
