# Joint physical conjugations of the four-letter reflection seed

A constant invertible map acting jointly on two physical legs can change
the opposite-cut Smith profile while preserving the first cut exactly.
For the four-letter reflection seed, all binary maps that also preserve
total physical height and commute with swap form a finite family of
1,536 matrices. Every member has pressure at least that of the original
seed at every scale. Thus this complete family gives no improvement,
including after optimization of the common physical type.

The result is specific to these binary maps and the inherited gauges.
It does not concern larger alphabets, maps with other coefficient fields,
or conjugations that change the height charges.

## 1. A new local crossing with one unchanged cut

Use the degree-twelve scalar seed `T(z)` and the heights

$$h=(0,3,1,4)$$

from [the reflection row-fiber calculation](triangular-auxiliary-row-fiber.md).
View its row-cut matrix as an operator on the pair of physical labels.
Let `S` be swap, and let `D(z)` be diagonal with entries
`z^(h(x)+h(y))`. Choose a constant invertible binary matrix `U` with

$$US=SU,\qquad UD(z)=D(z)U.$$

Define

$$T'(z)=U T(z)U^{-1}. \tag{1}$$

Since `T(1)=S`, equation (1) has the same literal crossing endpoint.
Its entries remain polynomials of degree at most twelve.

The restored zero matrix is `Z=z^(-4)D(z)T(z)`. Commutation with
`D` gives `Z'=UZU^{-1}`, so every zero Smith order is unchanged.
The regular infinity matrix before reshuffling is
`W(w)=w^4 T(w^(-1))D(w)`, and likewise `W'=UWU^{-1}`.
These identities preserve polynomial regularity at both endpoints.

The opposite cut is the reshuffle of `W'`. Joint similarity on the
physical pair need not be an invertible row/column operation on this
reshuffle. Its Smith profile can therefore change, and the checker
computes it directly.

As in the original seed, the zero scalar and infinity regularization
cancel in the restored product. The type energies remain `G_j=4h(j)`.
There is no new normalization charge in (1).

## 2. The binary family is complete

The pair heights range from zero through eight. Four fibers have dimension
one. Four more fibers have dimension two, each with swap acting as the
two-coordinate exchange. Its invertible binary centralizer consists of
the identity and that exchange.

The height-four fiber has basis

$$((0,3),(3,0),(1,2),(2,1)).$$

Here swap is the direct sum of two exchange matrices. Every commuting
matrix consists of four blocks of the form

$$\begin{pmatrix}a&b\\b&a\end{pmatrix},\qquad a,b\in\mathbb F_2.$$

There are 256 such matrices, of which 96 are invertible. This can also
be counted as `|GL_2(F_2)|*2^4=6*16`: in characteristic two the exchange
is `I+epsilon`, with `epsilon^2=0`.

The full family therefore has `96*2^4=1536` members. Enumeration of
these blocks exhausts all binary `U` satisfying the two commutation
conditions, rather than sampling a few elementary conjugations.

## 3. Exact prefix inequalities give an all-scale conclusion

The original restored profile at both endpoints is

$$a=(0,1,3,4,5,6,7,8,8,9,10,11,12,13,15,16). \tag{2}$$

Across the 1,536 conjugations there are 22 distinct opposite-cut
profiles `b`, each in increasing order. The exact computation proves

$$\sum_{i=1}^k b_i\le\sum_{i=1}^k a_i
       \qquad(1\le k\le16). \tag{3}$$

Here is an elementary pressure consequence of (3). For every real
`c>=0`,

$$\sum_i(c-b_i)_+
 =\max_{0\le k\le16}\left(kc-\sum_{i=1}^k b_i\right)
 \ge\sum_i(c-a_i)_+.$$

For `t>0`, multiply by `t^2 exp(-tc)` and integrate over `c>=0`.
All orders are nonnegative, and

$$t^2\int_0^\infty e^{-tc}(c-u)_+\,dc=e^{-tu}.$$

Consequently

$$\sum_i s^{-b_i}\ge\sum_i s^{-a_i}\qquad(s>1). \tag{4}$$

The zero profile stays (2), and the energies do not change. Equation
(4) proves that neither the pressure product nor its common-type
criterion can improve at any scale. Equality of the profiles recovers
the original seed.

The computation gives a sharper diagnostic. Exactly 32 conjugations
retain opposite-cut rank one at the endpoint. Every one of those
commutes with the entire polynomial `T(z)` and makes no tensor change.
They are also exactly the matrices commuting with `W(0)`. All other
conjugations raise the opposite-cut endpoint rank to at least three.
These are finite classification statements, not asserted for larger
reflection seeds.

## 4. Replay

[`joint_height_conjugations.py`](../tools/joint_height_conjugations.py)
constructs the complete family, both restored endpoints and their
profiles. The checker resolves all 1,536 zero profiles independently
as well as all opposite profiles. It checks the endpoint identities,
all inequalities in (3), and 66 independent block-Toeplitz ranks for
the 22 distinct profiles.

```sh
python3 tools/check_joint_height_conjugations.py \
  --output evidence/joint-height-conjugations.json
```

The [evidence](../evidence/joint-height-conjugations.json) records each
distinct profile, its multiplicity, a representative joint map, all
prefix slacks, and source hashes. The replay uses only exact binary
and integer arithmetic. The conclusion also holds for these same
binary matrices after coefficient-field extension in characteristic
two, because their Smith profiles do not change.
