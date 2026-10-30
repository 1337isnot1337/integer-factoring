# Joint label relabelings and a product-sixteen seed

## Finite result

The [subsequent complete lower-bound check](hamming-relabeling-optimum.md)
proves exact cutoff four in this 102-dimensional span, with a 20-piece
witness. The six-piece cutoff-16 seed below remains useful as a different
piece-count tradeoff and for its smaller algebraic block model.

The earlier [joint one/two-pair span](hamming-joint-optimum.md) has exact
cutoff 25 over `F_2`, at binary width two. Applying the same label
permutation to all four coordinates preserves the crossing and both
ranks. Combining different such relabelings enlarges the source span,
so its optimized cutoff need not stay the same.

For all 24 permutations of the four labels, the joint orbit of the old
33-dimensional space has dimension 102. A six-piece decomposition in
this larger space has ranks

```text
(1,16), (16,1), (4,4), (4,4), (4,4), (4,4).
```

Its maximum rank product is **16**, and maximum rounded width sum is
**four**, improving the old values 25 and six. This is an upper bound in
the larger space; no optimum claim is made for that space. The old
optimality theorem remains valid in its precisely smaller span.

The complete entry tables and source relabelings are retained in
[the witness](../evidence/hamming-relabeling.json). The checker regenerates
the original one/two-pair sources, their canonical 33-element basis, every
label image, and each witness's source representation. It independently
checks dense ranks and all 256 crossing entries. The witness also works
over every characteristic-two extension field.

This finite gain is not an improved asymptotic factoring theorem. The
strongest complete exponent remains `1/2`. New proofs used substantial
OpenAI assistant support and await independent review; no Lean work was
added.

## 1. Separate the new seed

Write the first two pieces as the usual matching tensors

```math
R=[x=y][u=v],\qquad C=[x=u][y=v],
```

and the last piece as

```math
X=[x\mathbin\oplus y\mathbin\oplus u\mathbin\oplus v=0],
```

where the labels are the two-bit vectors `0,1,2,3`. The three remaining
pieces `B_1,B_2,B_3`, each of ranks `(4,4)`, sum to

```math
A=[x,y,u,v\text{ are all distinct}].
```

Indeed `A=E+R+C+X`: on this four-element group a zero xor is either a
paired tuple or all four distinct labels. Thus

```math
E=R+C+B_1+B_2+B_3+X=R+C+A+X.                            (1)
```

The six-piece version has cutoff 16. The four-piece aggregate has
component ranks `(1,16),(16,1),(6,6),(4,4)` and is useful for an exact
block analysis. Aggregating is a change of analysis; it does not preserve
the six-piece finite cutoff.

## 2. A constant block model

Consider the aggregate pencil `rR+cC+aA+hX`. Right multiplication by the
crossing permutation preserves rank. A constant `F_2` change of basis
splits its row matrix into two scalar blocks `c`, one block

```math
K_{row}=cI_2+(r+h)N,\qquad N^2=0,
```

and three identical four-dimensional blocks

```math
J_{row}=cI_4+(c+a)U+(h+a)UV,
\qquad U^2=V^2=0,\ UV=VU.                              (2)
```

The latter is the regular Boolean algebra representation on
`F_2[U,V]/(U^2,V^2)`, in the order `(UV,U,V,1)`.
Explicit matrices are

```text
U  = rows (4,8,0,0)
UV = rows (8,0,0,0)
```

with binary row encodings. The column arrangement simply exchanges
`r,c`. The same constant basis works. The checker constructs this basis
from the actual seed and verifies all 2,048 entries of all four
coefficient matrices in both cuts; no variable-dependent inverse or
unproved simultaneous diagonalization is used.

For an algebraic construction of the basis, put `N_0=C E+I`,
`R_0=R E`, and `X_0=X E`. Their pairwise products vanish. Their joint
image has dimension seven, and joint kernel dimension nine.
`L=N_0+R_0` identifies a seven-dimensional domain complement with the
image. On that complement let `P=L^(-1)R_0`, `Q=L^(-1)X_0`. Exact
calculation gives

```math
P^2=P,\qquad PQ=QP=P,\qquad Q^2=P,
\qquad rank P=1,\quad rank(Q+P)=3.
```

Split off the image of `P`. On its six-dimensional kernel, `Q` is the
direct sum of three square-zero chains. Lifting these chains gives (2);
two remaining vectors give the scalar blocks. These finite constant
identities imply the block statement over every characteristic-two field.

## 3. Exact tensor-power ranks and effective new families

For word length `t`, each block type consists of `i` special blocks,
`j` four-dimensional blocks, and `k` scalar blocks, with `i+j+k=t`.
Its matrix is

```math
c^k K^{\otimes i}\otimes J^{\otimes j},
```

repeated `t!/(i!j!k!) * 3^j * 2^k` times. Polynomial coefficient
extraction and any coefficient recombination act inside every repeated
block. The exact rank is the sum of these multiplicities times the
corresponding block ranks. The column formula swaps `r,c`. At `t=2`,
104 independent direct global-matrix computations replay this formula.
The largest block has dimension `4^t`, versus global matrix dimension
`16^t`.

Either version of (1) gives effective polynomial tensor-power families.
Assign its constant number of pieces fixed scalar polynomials with value
one at `z=1`, and take the `t`-fold tensor product of their sum. The
crossing scale remains one, degree is `O(t)`, and the full coefficient
symmetry span has dimension `O(t)`. At label width `lambda=2t`, its entire
span can therefore be enumerated, ranked, optimized and factored in
`2^O(lambda)` time, with `O(lambda)` retained pieces. An odd width is
handled by restricting one binary bit in the next even-width family.

These families are outside the universal three-pairing type obstruction:
they have additional local tensors. A uniform rank-product bound below
the half threshold remains unproved. Tensoring a finite witness alone
leaves exponentially many pieces; collecting words into coefficients can
increase rank, so the finite cutoff 16 must not be extrapolated.
The [reciprocal choice](sl2-coefficient-ranks.md)
`r=1,c=z^2,a=z,h=1+z+z^2` has since been solved: full-span optimization
has exact leading contraction `0.884638857886...`, so that particular
choice does not supply an improvement. The other gradings remain open.

```sh
python3 tools/check_hamming_relabeling.py
python3 tools/check_hamming_relabeling.py --output evidence/hamming-relabeling.json
make verify
```
