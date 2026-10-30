# Local tensor powers and weight symmetry of the crossing deformation

## Result and limits

The [bounded-bit-weight crossing polynomial](hamming-pfaffian.md) has an
exact constant-size local description at every label width. Expanding its
Pfaffian gives a sum of 23 products of four-leg polynomial tensors of
local dimension two. This description is effective for arbitrary positive
bit weights with low weight one.

There is also an exact symmetry under permuting **all** non-low bit
positions, including the leading position used to choose a label
representative. Consequently, for fixed `B`, optimization over weights in
`{1,...,B}` needs only

```math
\binom{\lambda+B}{B-1}
```

weight multisets. Testing every ordered weight assignment is unnecessary.
The full-span construction budget was already `2^O_B(lambda)`; this lemma
reduces its outer search to polynomial size.

Neither statement supplies the missing uniform rank bound. In particular,
the generic ranks of the local tensors cannot be substituted directly for
the ranks of their polynomial coefficients. These are ordinary arguments,
with no new Lean formalization. The underlying Pfaffian construction and
its source attribution remain as stated in the linked note. The local
formula and symmetry argument used substantial OpenAI assistant support.

## 1. Separate one bit of a Pfaffian matching

Write the two parts of label `a` as `2a+1` and `H-2a-1`, where
`H=2^(lambda+2)` and `0<=a<2^lambda`. The low bit of both parts is one.
Their remaining `lambda+1` bits are complements: the first part has the
`lambda` bits of `a` followed by a leading zero, and the second has their
complement followed by a leading one.

Number the eight vertices `U0,U1,A0,A1,B0,B1,Y0,Y1`. A matching may only
join different blocks. There are 60 such perfect matchings. For a free
bit with input values `(x,y,u,v)`, vertex `(block,sheet)` has bit
`input(block) XOR sheet`, where the block inputs are `(u,x,v,y)`.
For an edge from a block on the left to a block on the right, call its
two bits `p,q`. Its local polynomial in `Z` is:

| Edge | Local factor `e(p,q;Z)` | Low-bit factor in `z` |
| --- | --- | --- |
| `UA` | `Z^(p OR q)` | `z(1+z)` |
| `UB` | `indicator(q<=p)(1+Z)^q` | `1+z` |
| `UY` | `[p=q] Z^q` | `z` |
| `AB` | `[p<=q] Z^p` | `z` |
| `AY` | `Z^(p AND q)` | `1+z` |
| `BY` | `[pq=0] Z^q(1+Z)^p` | `z(1+z)` |

This table follows entry by entry from the explicit folded matrix. For
example, the `AY` low-bit contribution to `l_w(n AND t)` is canceled by
its minus one, leaving the displayed low factor `1+z`. The `BY` condition
`m AND t=1` imposes `pq=0` at every non-low bit and imposes no additional
condition at the low bit. The `UY` equality is likewise a product of
bitwise equalities.

For a matching `M`, define its low-bit scalar and local tensor by

```math
L_M(z)=\prod_{e\in M}\operatorname{low}_e(z),\qquad
A_M(X,Y,U,V;Z)=\prod_{e\in M}e(p_e,q_e;Z).
```

In characteristic two the Pfaffian is a sum over matchings without sign
changes. If bit `j` has weight `w_j`, its complete contribution is obtained
by substituting `Z=z^(w_j)`. The leading bit always has four input bits
zero. Thus the exact formula is

```math
P_{\lambda,w}(x,y,u,v;z)=
\sum_M L_M(z) A_M(0,0,0,0;z^{w_{\lambda+1}})
\prod_{j=0}^{\lambda-1}
 A_M(x_j,y_j,u_j,v_j;z^{w_{j+1}}).
\tag{1}
```

For all weights one, this is a scalar polynomial times `A_M(z)^tensorlambda`
for each matching. Identical local tensors are combined by adding their
low-bit scalars. Terms with identically zero leading entry vanish for
every positive leading weight. Enumerating the 60 fixed matchings leaves
23 nonzero local terms. The source code generates these constants; the
evidence file records every local matrix and scalar, so the formula does
not depend on an omitted table or an unspecified choice.

## 2. Permuting bit positions preserves complete profiles

Let `pi` permute the `lambda+1` non-low positions and put
`w'_j=w_{pi(j)}`. Apply that same bit permutation to every mode in every
label. An unordered antipodal pair stays an unordered antipodal pair.
If its chosen first mode now has leading bit one, complement its non-low
bits to choose the other member as representative. This defines a
bijection `phi_pi` on the `2^lambda` labels.

Every displayed folded entry is preserved: the bitwise subset, equality,
and disjointness conditions are invariant, and the permuted weights keep
all weighted degrees and products unchanged. Reordering the two vertices
inside a block changes no Pfaffian sign in characteristic two. Therefore

```math
P_{\lambda,w}(x,y,u,v;z)=
P_{\lambda,w'}(\phi_\pi(x),\phi_\pi(y),
                       \phi_\pi(u),\phi_\pi(v);z).
\tag{2}
```

This includes permutations exchanging a free bit and the leading bit.
Restricting to permutations of the free bits alone would give a weaker
reduction and would miss equivalent choices of leading weight.

Taking polynomial coefficients in (2) simultaneously permutes the rows
and columns in either matrix arrangement. It preserves both ranks,
commutes with each of the eight crossing-preserving coordinate
symmetries, and preserves the crossing tensor. Hence the entire
symmetry-span optimization problem, including its optimum cutoff and
witness piece count, is identical for weight assignments with the same
multiset.

There are `lambda+1` non-low positions. Their multisets on `B` weight
values are counted by weak compositions of `lambda+1` into `B` parts,
giving `binom(lambda+B,B-1)`. Enumerating these compositions and choosing
one sorted representative per composition is a uniform polynomial-time
outer loop for fixed `B`.

## 3. What the local ranks do and do not imply

Over the rational function field `F2(z)`, the 23 local tensors have these
pairs of matrix ranks:

| Row rank, column rank | Number of local terms |
| --- | --- |
| `(4,1)` | 2 |
| `(1,4)` | 2 |
| `(3,3)` | 18 |
| `(4,4)` | 1 |

The ranks are calculated by exact polynomial minors. In the all-one
specialization, tensor powers multiply these generic ranks. But the
crossing pieces are polynomial coefficients of their sum. Coefficient
extraction can increase matrix rank relative to a rational-function
factorization, and cancellation between different terms can decrease it.
Neither operation is accounted for by simply multiplying the table's
rank pairs.

In particular, the single `(4,4)` local term is `z^2 E_bit`, with its
all-one scalar `z^6`. Its tensor power contributes
`z^(2lambda+6) E_lambda`. The rest of the profile must be analyzed together
with this full-rank contribution. Dropping it would destroy the crossing
sum. The formula makes this obstruction explicit and supplies a fixed
collection of small matrices for further analysis; it does not claim to
resolve it.

## 4. Exact replay

```sh
python3 tools/check_hamming_tensor_power.py
python3 tools/check_hamming_tensor_power.py --output evidence/hamming-tensor-power.json
make verify
```

The checker compares 2,944 complete polynomial profiles with the original
eight-vertex Pfaffian evaluator through width six, including nonconstant
weights. It also checks 5,024 complete-profile permutation identities
through width four, using every permutation of the non-low positions at
these widths and exercising exchanges with the leading position. These
finite checks supplement the entrywise all-size proofs above. Exact
polynomial minors verify all local generic ranks independently of the
large-domain rank routines. The complete finite-span optimum claims stay
in [hamming-pfaffian.md](hamming-pfaffian.md).
