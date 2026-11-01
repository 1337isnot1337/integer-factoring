# Five-label twisted-pairing preflight

A bounded exact search tested scalar entry spaces outside the previously
excluded standard three-pairing and restricted-subalphabet families.
Every space in the specified catalogue was excluded before auxiliary
matrix enumeration. This is a finite catalogue result, not a theorem
about arbitrary twisted pairings.

## 1. The narrow five-label target

The [universal observable detector](local-alphabet-and-pairing-spectral-barriers.md)
forces spectral base at least `125/13` on five labels. An initially
suggested target `19/2` was therefore impossible and was discarded
before the preflight. A consistent sufficient target is `97/10`:

```math
\frac{125}{13}<\frac{97}{10}<5^{\sqrt2}.
```

The comparisons are exact: `97*13>125*10`, `24^2<2*17^2`, and
`97^17<10^17*5^24`. Thus there is only a narrow interval for an
improving five-label individual-word spectral construction.

## 2. A finite catalogue of new scalar spaces

Work over F2. Define row and column rank-one tensors

```math
R(A,B)_{xyuv}=A_{xy}B_{uv},\qquad
C(A,B)_{xyuv}=A_{xu}B_{yv}.
```

Use five-by-five matrices from the catalogue `I,P,I+P,D,A`, where
`P_(i,j)=[j=i+1 mod 5]`, `D=diag(1,1,0,0,0)`, and
`A=e_(0,1)+e_(1,0)+e_(2,3)+e_(3,2)`.
There are 25 ordered row tensors and 25 ordered column tensors.
Choose two distinct row tensors `U_1,U_2` and one column tensor `V`,
giving `binom(25,2)*25=7500` catalogue choices for spaces

```math
S=\operatorname{span}_{\mathbf F_2}(E,U_1,U_2,V).
```

The intended observable template was

```math
T_0=\begin{pmatrix}
0&0&U_1&0\\0&0&U_2&0\\0&0&E&V\\0&0&0&0
\end{pmatrix},\qquad c=d=e_3,
```

where `e_3` is the third standard basis vector. Every middle-to-middle
word uses the middle self-loop at each step, so its all-word endpoint
is exactly the crossing tensor power. Conjugation by `W` preserves
the endpoint with boundaries `c^T=e_3^T W`, `d=W^{-1}e_3`.
There are 20,160 matrices in `GL_4(F2)`.

## 3. Exact preflight result

Different catalogue choices can span the same scalar space. The count
7,500 refers to distinct index triples, not distinct subspaces.

All 7,500 displayed tensor bases were independent. The checker evaluated
all fifteen nonzero coefficient vectors in each space: **112,500 tensor evaluations
and 225,000 exact cut-rank evaluations**. Every tensor with crossing coefficient
one obeyed

```math
r(t)c(t)\ge121.
```

The smallest observed product was exactly 121. One attaining space uses
row indices `0,1` and column index `0`, with coefficient vector
`(1,1,0,1)`. Its tensor is the standard `E+R+C`, whose two ranks are
`5*4/2+1=11`; this witness also follows from the general pairing formula.
Hence the scalar linear
map extracting the crossing coefficient is a rank detector with
multiplier eleven. The [detector lemma](linear-rank-detectors-for-observable-tables.md)
gives spectral base at least **11** for every finite observable table
whose entries lie in any of these spaces. Since
`11^7>5^10` and `sqrt(2)<10/7`, this exceeds the improvement threshold.

No linear program or auxiliary conjugation search was needed. The
exclusion applies uniformly in auxiliary dimension over F2, but only
to the specified finite catalogue and the individual-word rank spectral
method. It does not lower-bound actual ranks after grouping words.

## 4. Reproduction

[Catalogue checker](../tools/check_twisted_pairing_catalogue.py) and
[deterministic catalogue evidence](../evidence/twisted-pairing-catalogue-certificate.json);
[independent dense controls](../tools/check_twisted_witnesses_dense.py) and
[their evidence](../evidence/twisted-witnesses-dense.json).

`check_twisted_pairing_catalogue.py` recomputes the complete catalogue
using exact bitwise elimination and no dependencies outside Python's
standard library. Its deterministic evidence records all 7,500 distinct
index triples, their minimum products, the attaining witness, and the
SHA256 of the complete ordered rank-table stream. The full scan completed
in about ten seconds in the session environment.

An independent dense check also verifies the attaining `(11,11)` witness.
For comparison, the tensor `E+C(I,I)+R(J_01,J_01)+R(J_23,J_23)`,
where `J_ij=e_ij+e_ji`, has cuts `(8,18)` over F2. These row matrices
lie outside the catalogue. Thus the product bound must not be replaced
by a claim that each cut is at least eleven.
`check_twisted_witnesses_dense.py` constructs both literal flattenings
and verifies these controls by ordinary dense elimination.

No linear programming or auxiliary conjugation code is needed to verify
this result. More general catalogue searches remain open. A failure to
find a detector in a future search would leave a space unresolved; it
would not establish an improving construction.
