# Polynomial crossing profiles with bounded bit weights

## Result and scope

There is a new effective family of characteristic-two crossing polynomials
with degree linear in the binary label width. Consequently, **the entire
symmetry span can be searched exactly within `2^O(lambda)` construction
cost**, and an optimal decomposition in that span needs only `O(lambda)`
pieces. This goes beyond bounded-size recombinations without violating the
factor-generation budget.

The family is a deformation of the folded Pfaffian entries. It is not a
regrading of the original source operators, and it asserts no replacement
for their exchange or intermediate-energy identities. Its crossing sum is
proved directly at `z=1`; its factors come from exact matrix elimination.

No uniform improved width bound is proved. Finite full-span optima at label
widths 1, 2, and 3 have maximum rank products 4, 36, and 256. These do not
establish a factoring exponent above `1/2`. The all-size construction and
cost statement are ordinary mathematical arguments, not Lean proofs.

The inspiration is Claude (Anthropic)'s Appendix B at
[formal-math@16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth)
and the [folded Pfaffian formula in PR #22](https://github.com/1337isnot1337/integer-factoring/blob/735b3ddb37d2e63dc469cf833687355e42097dd4/research/modular-pfaffian.md).
The deformation, complete-span optimization argument, and finite searches
below used substantial OpenAI assistant support and retain Apache-2.0.

## 1. Labels and explicit polynomial entries

For `D=2^lambda`, put `H=4D`. Label `a=0,...,D-1` by the distinct odd pair

```math
mu(a)=(2a+1,H-2a-1).
```

These pairs are disjoint from one another and have the same ordinary sum
`H`. Their integers have only `lambda+2` bits.
Fix an absolute integer `B>=1`, choose `w_0=1` and
`1<=w_j<=B` for `1<=j<=lambda+1`, and define

```math
l_w(n)=sum_{j:bit_j(n)=1} w_j,\qquad
f_w(n,z)=prod_{j:bit_j(n)=1}(1+z^{w_j}).
```

The Hamming specialization has all weights one. Form an eight-by-eight
alternating polynomial matrix with two vertices in each block `U,A,B,Y`.
For parts `t in u`, `n in x`, `m in v`, and `t in y`, use:

| Blocks | Entry |
| --- | --- |
| `U_t,A_n` | `z^(l_w(n OR t))(1+z)` |
| `A_n,Y_t` | `z^(l_w(n AND t)-1)(1+z)` |
| `A_n,B_m` | `[bits(n-1) subset bits(m-1)] z^(l_w(n))` |
| `U_t,B_m` | `[bits(m-1) subset bits(t-1)] f_w(m,z)` |
| `B_m,Y_t` | `[m AND t=1] z^(l_w(t)) f_w(m,z)` |
| `U_t,Y_s` | `[t=s] z^(l_w(s))` |

Within-block entries vanish. Since all parts are odd and `w_0=1`, the
exponent in the second row is nonnegative. Define `P_lambda(x,y,u,v;z)`
to be its Pfaffian and take its polynomial coefficients as pieces.

## 2. Exact crossing at z=1

At `z=1`, all blocks except `AB` and `UY` vanish. The `UY` determinant
is one exactly when the two disjoint pair labels agree. To compute the
`AB` determinant, replace each mode `n` by `(n-1)/2`. A label becomes
`{a,bar(a)}`, where the complement uses `lambda+1` bits and the leading
bit of `a` is zero. Between labels `a,b`, the subset matrix is

```math
\begin{pmatrix}
[a subset b]&[a subset bar(b)]\\
0&[b subset a]
\end{pmatrix}.
```

The lower-left entry is zero because `bar(a)` has leading bit one and
`b` has leading bit zero. Its determinant is `[a=b]`.
The Pfaffian at `z=1` is the product of these two determinants, hence

```math
sum_i [z^i] P_lambda(x,y,u,v;z)=[x=v][y=u].
```

This proof is independent of any homogeneity of the modified entries.
The exact scale is one, so the pieces extend to every commutative
`F_2`-algebra.

## 3. Degree and the full construction cost

Every pair label has weighted size `g=2+sum_(j>=1) w_j` because its two
parts have complementary non-low bits and both have low bit one.
Each displayed entry has degree at most the sum of its endpoint sizes.
For example, `l_w(n OR t)+1<=l_w(n)+l_w(t)` because `n AND t` contains
the low bit. Therefore every perfect matching has degree at most

```math
4g<=4(2+B(lambda+1))=O_B(lambda).
```

There are `O_B(lambda)` coefficient pieces, and an entry of the complete
profile takes polynomial time in `lambda` to compute: the matrix has just
eight vertices and polynomial-length binary coefficient strings.
All `D^4` entries of every piece cost `2^O(lambda)` time.

Take all eight crossing-preserving coordinate symmetries of every piece
and let their span over `F_2` be `S_lambda`. Its dimension `m` is at most
`8(4g+1)=O_B(lambda)`. Enumerate all `2^m` tensors in that span, calculate
both exact matrix ranks, and sort by rank product. At the first cutoff
whose candidate span contains the crossing, solve for a decomposition.
It uses at most `m` tensors. A rank-`r` matrix admits an explicitly computed
rank factorization via independent rows and columns, padded to binary
width `ceil(log2 r)`.

All enumerations, ranks, sorting, span tests, and factor generation have
fixed polynomial cost per tensor in `D` and `m`. The total is therefore
`2^O_B(lambda)`, including factor values and metadata. This is a single
uniform algorithm, not advice choosing an uncomputable low-rank subset.
Increasing fixed constants absorbs finitely many small widths.

For any fixed `B`, one can additionally enumerate all `B^(lambda+1)`
weight assignments and select the best resulting decomposition. This
still costs `2^O_B(lambda)`, with only `O_B(lambda)` output pieces. Allowing
unbounded weights such as `w_j=j` changes the degree bound; the full-span
cost argument must not be applied to those experiments unchanged.

## 4. Exact optimization target and the factoring connection

Let `Q_B(lambda)` be the selected minimum maximum rank product. If one proves

```math
log2 Q_B(lambda)<=4 rho lambda+O(log lambda)
```

for fixed `1/2<=rho<1`, then rounding adds at most two width bits.
All effective modular-pieces conditions are already satisfied, so the
[scale-one transfer](scale-one-modular-transfer.md) gives factoring
exponent `log2(1/rho)`. This is a conditional implication; the required
uniform rank bound is not supplied by the finite data.

For comparison, every polynomial-piece decomposition with maximum product
`Q` satisfies `D^2<=M Q`: the crossing has row rank `D^2`, each nonzero
piece has row rank at most its rank product, and rank is subadditive.
Thus a leading width coefficient below `2lambda` is impossible in this
model. This is a restriction on the decomposition, not on factoring.

## 5. Finite evidence and remaining question

For all weights one, the symmetry spans at widths 1, 2, and 3 have dimensions
3, 13, and 21. Full enumeration gives cutoffs 4, 36, and 256. The stored
witnesses are replayed directly; the latter two improve substantially on
using the original coefficient pieces, whose maxima are 108 and 2209.
The cutoff 256 at width three does not prove a growing saving.

The uncombined all-one profiles through width five attain maximum binary
width sums `2,8,12,16,20`. Bounded-size recombinations alone do not settle
the optimized full-span behavior. A uniform rank theorem, or a scalable
explicit decomposition drawn from this span, remains the next target.

## 6. Replay and verification scope

```sh
python3 tools/check_hamming_pfaffian.py
python3 tools/check_hamming_pfaffian.py --full --output evidence/hamming-pfaffian.json
make verify
```

The full suite checks 1,118,480 crossing entries in the all-one families
through width five, plus 4,352 entries for two nonconstant bounded weight
assignments. It records the complete coefficient rank profiles. Native
full-span enumeration checks all 7, 8,191, and 2,097,151 nonzero span tensors
at widths 1, 2, and 3, respectively. Every selected witness rank is checked
with independent dense elimination and its complete crossing sum is replayed.

The native classifier is additionally compared with Python on 131 span
indices per small domain, including zero, generators and deterministic
random samples. At width two, a separate complete Python enumeration of
all 8,191 span tensors confirms the optimum. The default suite replays the
stored width-three witness and oracle checks; `--full` repeats its full
native search and the larger raw profiles. Source and input hashes bind
the recorded evidence to the exact construction.

The native helper supports the recorded domains up to eight labels. The
uniform all-length algorithm is given explicitly in Section 3 and in the
Python span generator; its asymptotic cost does not imply that larger
finite enumerations are practical. Neither threshold clipping in the
small searches nor a positive result at a few widths proves the missing
uniform width bound.
