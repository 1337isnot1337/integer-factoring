# Pencil algebras and low-rank correction modules

## Structural target

The [bounded-support factorization](bounded-spin-support-ranks.md) gives
a coefficient tensor `T=T_(n/2,0,0)` with rank at most

```math
a_{n/2}=O(g^n),\qquad g=\sqrt{2+\sqrt2}<\varphi^{4/3}.
```

This suggests using `T` as a generator of a small correction module,
rather than using corner coefficient tensors alone. The general lemma
below makes the needed hypotheses precise. Exact eight- and ten-site
algebras support that target, but the simplest linear-span version
already fails at twelve sites. No uniform favorable decomposition is
proved here.

## 1. A correction-module lemma

Let `V` be a finite-dimensional vector space, let `B` be a unital
matrix algebra on `V` of dimension `d`, and let `T` have rank at most
`a`. Put

```math
M=B(\operatorname{im}T)
 =\sum_{j=1}^d B_j(\operatorname{im}T),
```

where the `B_j` form any basis of `B`. Then:

1. `dim M<=d a`.
2. `M` is invariant under the algebra `A=Alg(B,T)`.
3. For every `X in A`, some `Y in B` has `rank(X-Y)<=d a`.
4. If `B` is commutative, the induced algebra on `V/M` is commutative.

The first assertion follows from the displayed sum. Multiplication by
`B` preserves it because `B` is an algebra. Multiplication by `T` sends
it into `im T`, which is contained in `M` since `I in B`. This proves
invariance. In a word in `B` and `T`, every term containing `T` has
image in `M`; the terms containing no `T` belong to `B`. Subtracting
their sum proves the third assertion and also the quotient assertion.

Apply this to the normalized row operators of the coherent spin family.
If the original pencil generates an algebra of polynomial dimension,
and all six-vertex coefficient operators belong to `Alg(B,T)`, then
the correction module has dimension `poly(n)*g^n`, below the desired
rank-cost growth. These are substantive **uniform hypotheses**, not
consequences of the scalar invariant-ring dimension.

The [cut involution](bounded-spin-projector.md) identifies the two
operator spaces and fixes the extreme tensor. The same argument can
therefore be made in each cut. The resulting approximations need not
come from the same tensor in `B`. Simultaneous two-cut lifting and a
crossing decomposition remain additional requirements; a small invariant
module alone does not establish a factoring improvement.

## 2. A verified ten-site algebra

At length ten there are 34 path labels, and each normalized cut matrix
has size 1156. The original coefficient tensors have dimension eight.
Use the basis

```math
(C_0,C_2,C_3,C_4,C_5,C_6,C_7,C_{10}).
```

Their row operators form a commutative algebra isomorphic to

```math
F_5^6\times F_5[\epsilon]/(\epsilon^2).              \tag{1}
```

The [model certificate](../evidence/bounded-spin-ten-algebra-model.json)
gives a unit, multiplication operators, seven orthogonal idempotents,
a radical generator, and two generating linear combinations. All full
matrix closure identities for those two generators are checked. Their
actions on the small coefficient algebra generate the entire
eight-dimensional space from its unit; hence these identities prove
the full algebra assertion, rather than only behavior on a sampled
vector.

This algebra has seven distinct maps to `F_5`. It **cannot** be generated
by one element over `F_5`: a character of a one-generator algebra is
determined by the generator's value, giving at most five possibilities.
The monogenic eight-site method consequently needs a real structural
change at ten sites. The two-generator certificate avoids adjoining
a field merely for the verification.

The primary dimensions on the 1156-dimensional module are:

| Component | Module dimension | Radical rank |
| --- | ---: | ---: |
| Scalar | 34 | 0 |
| Scalar | 1 | 0 |
| Scalar | 225 | 0 |
| Scalar | 300 | 0 |
| Scalar | 35 | 0 |
| Dual numbers | 406 | 56 |
| Scalar | 155 | 0 |

Their sum is 1156. On the dual-number component, a nonzero scalar part
gives rank 406; a zero scalar part and nonzero radical coefficient give
rank 56; both zero give rank zero. The other components contribute
their whole dimension precisely when their character is nonzero.
This is an exact rank formula over every extension of `F_5`.

## 3. Full-span finite optimization remains above the target

The eight character/radical forms and their cut-involuted counterparts
give sixteen linear forms in eight coordinates. Enumerating their
9,805 distinct common kernels gives the sharp minmax rank product

```math
\boxed{50625=225^2>34^3.}                          \tag{2}
```

Every direction below this product lies in
`3*w_3+3*w_5=0`, while the crossing's coordinates
`(1,2,0,3,1,4,2,1)` have value one. This proves the lower bound over
all characteristic-five extensions using the same common-kernel
argument as the [eight-site certificate](bounded-spin-rank-algebra.md).
An explicit base-field spanning list, scaled to sum to the crossing,
attains the bound. The replay checks both actual cut ranks of those
pieces directly, using the column tensor rather than only the involution.

Thus the complete ten-site original span still fails even the weaker
finite `L^(3/2)` comparison. This is not an all-length obstruction:
polynomial slack in a uniform rank theorem cannot be decided from one
finite threshold.

## 4. The vertex module and the twelve-site distinction

At ten sites the full six-vertex coefficient space has dimension nine
and is exactly `B+F_5 T`. The extreme tensor has rank 116. Direct
calculation of `sum B_j im T` gives dimension **170**, and the lemma
proves its invariance under every coefficient operator. Thus the
quotient has dimension 986. This is different from the smaller
commutator-generated module observed in the earlier prototypes.

At twelve sites, twenty explicit entry evaluations certify:

```text
original pencil span:             10
original pencil plus the vertex:  11
full six-vertex span:             14
```

The entry matrix supplies the lower bounds. The upper bound ten is the
uniform three-jet bound, and the full bound fourteen follows by counting
degree-twelve monomials in the degree `2,4,4,6` generators of the
[odd-characteristic invariant ring](four-qubit-orthogonal-generators.md).
The old-plus-vertex space is at most eleven-dimensional and the same
entry matrix attains eleven.

Consequently, `F=B+F_5 T` is not the uniform formula to prove. The
weaker algebra-containment target `F subset Alg(B,T)` allows nonlinear
operator words and remains open. So do a polynomial bound on `dim B`
and a compatible two-cut decomposition using the quotient and module.
These are the next structural questions; a larger version of the same
fixed-block optimization is not supplied as their answer.

## Reproduction and arithmetic

Run `make verify-bounded-spin-growth`. This optional target needs NumPy,
a C compiler, `pkg-config`, and the FLINT development library. The
entry constructors used by the factoring-family interface retain their
standard-library implementations.

The replay batches signature blocks with the same difference weight.
After deleting fixed coordinates, each block has the same matrix indexed
by constant-weight subsets. Concatenate their encoder columns, apply
that repeated local matrix, then perform one final matrix product per
monomial. [`bounded_spin_batched_blocks.py`](../tools/bounded_spin_batched_blocks.py)
implements this grouping with exact `F_5` arithmetic through the
[small FLINT wrapper](../tools/native_gf5_flint.py).
The FLINT functions used are documented in its
[modular matrix interface](https://flintlib.org/doc/nmod_mat.html).

The batched eight-site tensor is compared in full with the independent
integer block construction. Sixteen ten-site entries are compared with
the magnetic constructor; full basis, endpoint, algebra, idempotent,
radical, module and separator identities are then checked. Native small
matrix operations are compared with independent NumPy/elementary
field computations. The
[receipt](../evidence/bounded-spin-growth.json) pins the sources and
records the FLINT version. No large tensor cache is needed in the repo.

These proofs and computations used substantial Codex AI assistance.
The complete factoring score remains `1/2`; no new Lean work or
independent acceptance is claimed.
