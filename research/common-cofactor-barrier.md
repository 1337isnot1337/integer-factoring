# A shared projected cofactor forces large actual ranks

## Result and scope

A projected cofactor supplies a stronger common rank certificate than the
determinant quotients for the asymmetric quadratic family. For every
normalized ternary binary-alphabet quadratic whose two-cut rank product
is at most two at `0` and `2`, the exact catalogue contains one polynomial
that is a projected `3` by `3` minor in **both** cuts and is nonzero at the
crossing endpoint.

An ordinary representation-theoretic argument amplifies each certificate
to every tensor length. It forces

```math
\rho\ge\frac{\log3}{2\log2}=0.792481250\ldots
   >\frac34>\frac1{\sqrt2}.                         \tag{1}
```

Thus this entire ternary quadratic family cannot improve the existing
half-rate factoring result through its full word coefficient spans.
The complete finite classification covers **32,265 normalized seeds**,
representing **12,852,225 actual coefficient seeds**. The earlier
rank-product-three cutoff is equivalent here, so it is covered as well.

The amplification theorem holds over every field whenever the displayed
common-minor certificate exists. The existence classification here is
for coefficients in `F_3`, and remains valid over their extensions.
It does not assert existence for arbitrary coefficients in other fields.
It does not cover new physical code restrictions or recombinations of
spaces from different seeds. No factoring exponent above one half is
established, and no Lean formalization is added.

## 1. A coefficient of a large repeated representation

Let `V` have dimension `n`, put `k=n-1`, and write `C_k(M)` for the
matrix of `k` by `k` minors of an `n` by `n` matrix `M`. For fixed
nonzero vectors `a,b`, a projected cofactor is

```math
p=a^T C_k(M)b.                                      \tag{2}
```

Every `k`-vector in an `n`-dimensional space is decomposable when
`k=n-1`. Therefore (2) is also an actual projected determinant
`det(L M R)` for fixed full-rank matrices `L,R` of sizes `k` by `n`
and `n` by `k`.

The induced polynomial module of highest weight `(t^k,0)` is

```math
\operatorname{Sym}^t(\mathop{\bigwedge}^{n-1}V)
 =\operatorname{Sym}^t(V^*)\otimes\det^t.            \tag{3}
```

This is the projective-space case of the Schur construction: the
Grassmannian of hyperplanes is projective space, so its Pluecker
coordinates have no further relations. The polynomial `p^t` is a matrix
coefficient of (3), by applying a pure symmetric power of `b` and
evaluating at a pure power of `a`. This remains an identity in positive
characteristic; no factorial is inverted.

By the dual Weyl filtration in Doty's
[*Schur-Weyl duality in positive characteristic*, arXiv:math/0610591v3](https://arxiv.org/pdf/math/0610591v3),
Lemma 13, and its characteristic-independent formal characters, the
module (3) occurs in `V^(tensor kt)` with multiplicity

```math
h_k(t)=\frac{(kt)!\prod_{j=0}^{k-1}j!}
                   {\prod_{j=0}^{k-1}(t+j)!}
       =\frac{k^{kt}}{\operatorname{poly}(t)}.       \tag{4}
```

The [joint determinant proof](joint-determinant-barrier.md) gives the
elementary polynomial bound and explains why these integer multiplicities
are not reduced as field scalars. As there, a finite base field can be
extended to its algebraic closure without changing ranks.

In a constant basis respecting that filtration, every polynomial matrix
`M(z)` acts block triangularly. If a coefficient functional `L` has
`L(p^t)!=0`, the induced operator on each of the `h_k(t)` selected
diagonal factors is nonzero. A minor selecting one nonzero entry from
each triangular block yields

```math
\operatorname{rank}L(M^{\otimes kt})\ge h_k(t).      \tag{5}
```

This uses a filtration, with no semisimplicity assumption. Tensoring on
`r` extra positions preserves the conclusion whenever
`L(p^t ell^r)!=0` for an entry polynomial `ell` of `M`.

## 2. One shared polynomial bounds both cuts of one piece

Let `K(z)` be a fixed polynomial crossing seed on `q` physical letters,
so each cut matrix has order `n=q^2`. Suppose its cut matrices admit
projected cofactors with exactly the same polynomial `p`, nonzero at
the crossing endpoint. Take word length

```math
m=(n-1)t+r,\qquad 0\le r<n-1,
```

and choose an entry `ell` that is nonzero at the endpoint.
The polynomial `p^t ell^r` is a linear combination of entries of
`K^(tensor m)` and is nonzero at the endpoint. For pieces

```math
A_i=L_i(K^{\otimes m}),\qquad \sum_i A_i=\kappa E_m,
\qquad\kappa\ne0,
```

some `L_i` must have `L_i(p^t ell^r)!=0`: apply the entry covector
representing this polynomial to the displayed tensor identity.
Because the polynomial is a projected cofactor in both cuts, (5) applies
to both ranks of that same piece. Hence

```math
\boxed{\max_i r(A_i)c(A_i)\ge
 h_{n-1}(\lfloor m/(n-1)\rfloor)^2.}                \tag{6}
```

No bound on the number of pieces is needed for (6). Its leading
contraction floor is `log(n-1)/log n`. For binary letters, `n=4`,
and (1) follows. Both strict comparisons in (1) can be checked without
floating logarithms: `3^2>2^3` gives the first, and `18>16` the second.

The theorem applies to arbitrary coefficient recombinations and therefore
to Taylor jets and homogeneous infinity coefficients. It concerns the
full word alphabets; unrelated physical restrictions are outside its proof.

## 3. Complete finite certificate coverage

The [low-rank quadratic catalogue](multipoint-low-rank-seeds.md) gives
nine first-value representatives and every one of 3,585 second values.
Each pair determines a unique quadratic with strict crossing at one.
Constant invertible leg changes and the allowed physical permutations
preserve existence of projected-cofactor certificates: exterior powers
transform the projection vectors invertibly, while a cut swap exchanges
the two certificates. The independently checked orbit coverage therefore
transfers the 32,265 certificates to every pair of the 3,585 allowed
finite values, giving `3585^2=12,852,225` actual coefficient seeds.

For each normalized seed the saved record specifies:

- Four exact ternary projection vectors, one left and one right in each cut.
- The common degree-at-most-six polynomial, normalized to value one at one.
- The exact first-value and second-value indices in the complete catalogue.

The fast constructor fixes the left vector to a coordinate vector and
solves a linear intersection of cofactor-row spaces. This covers **32,263**
cases. Two cases, indexed `(7,29)` and `(7,56)`, need general left
projections. A complete projective search supplies them over `F_3`.
There are no remaining cases.

The verifier uses a different construction from the search. It lifts each
four-vector to a literal `3` by `4` matrix with those Pluecker coordinates,
then expands `det(L M R)` directly. It checks **64,530 projected determinant
identities** and all **32,265 nonzero endpoint identities**. It does not
infer polynomial equality merely from agreement at the three field points.
Independent literal tensor-power controls test the representation rank
bound in characteristics two, three and five, including nonsemisimple cases.

```sh
python3 tools/check_common_cofactor_barrier.py \
  --output /tmp/common-cofactor-replay.json
python3 tools/check_common_cofactor_barrier.py \
  --generate /tmp/common-cofactor-catalogue.json \
  --output /tmp/common-cofactor-fresh-replay.json
cmp evidence/common-cofactor-quadratics.json /tmp/common-cofactor-catalogue.json
```

The first command independently verifies every saved certificate. The
second reconstructs the complete deterministic search, verifies it, and
permits byte-for-byte comparison with the saved catalogue. Source and
catalogue hashes are retained in the validation record.

## 4. Actual finite jets and the next constructive requirement

The accompanying actual-jet search checks one witness for each of the
195 unequal-determinant signatures in the earlier quadratic catalogue.
It tries every minimal degree-budget multiplicity allocation and its
exact Hermite coefficients, omitting zero interpolation coefficients.
It retains both actual ranks of every local jet and reconstructs the
complete crossing with the winning schedule.

| Tensor copies | Smallest maximum rank product over these witnesses |
| --- | ---: |
| 1 | 4 |
| 2 | 64 |
| 3 | 816 |

These are finite schedule optima over the retained witnesses, not all
coefficient tensors within a signature. The saved computations reconstruct
**851,760 crossing entries**. A separate polynomial-Kronecker and column-rank
oracle checks selected witnesses in all three lengths. The common-cofactor
theorem supplies the uniform obstruction that the finite comparisons alone
could not establish.

```sh
python3 tools/search_asymmetric_quadratic_jets.py --copies 1 2 3 \
  --output /tmp/asymmetric-quadratic-jets-replay.json
```

The previous determinant quotient floor `2/3` left room numerically.
The projected cofactor closes that opening for this full ternary quadratic
family. A new construction must escape the shared projected-minor
polynomials, in addition to supplying asymmetric pieces and complete costs.
Generic increases in numerical rank are not an upper-bound mechanism;
the next useful target is a structured family in which the two projected-
minor images separate while actual coefficient cancellation persists.

This proof and its computations used substantial OpenAI Codex assistance.
