# A rank barrier for position-invariant word families

## The all-size result

Fix an alphabet size `d>=2` and a field of characteristic zero. Labels are
all words of length `t` in `[d]`, so the label count is `D=d^t`. Suppose
an exact crossing decomposition

```math
E_D(x,y,u,v)=[x=v][y=u]=sum_i F_i(x,y,u,v)
```

has every summand invariant under simultaneously permuting the `t` word
positions in all four labels. The number of summands is arbitrary; they
may be arbitrary linear recombinations of invariant tensors.

Put `q=d^2`. If `t=2qk`, some summand has **both** matrix ranks at least

```math
B_{q,k}=
\frac{(2qk)!}{((q+1)k)!\,(k!)^{q-1}\,(2qk+q)^{q(q-1)/2}}.
\tag{1}
```

For general `t>=2q`, the same assertion holds with
`k=floor(t/(2q))`, by restricting fewer than `2q` positions. Consequently,
writing `lambda=t log2 d`, some summand has total binary channel width at
least

```math
4rho_d lambda-O_d(log lambda),\qquad
rho_d=\frac{log2(2q)-\frac{q+1}{2q}log2(q+1)}{log2 q}.
\tag{2}
```

In particular:

| Seed alphabet `d` | Forced contraction parameter `rho_d` | Rate ceiling under the existing transfer |
| --- | ---: | ---: |
| 2 | 0.774397470347... | 0.368853854985... |
| 3 | 0.733269522816... | 0.447584517591... |
| 4 | 0.707133841396... | 0.499944790681... |

All three forced parameters are strictly larger than `1/sqrt(2)`.
For seed four the positive margin is about `0.00002706021`; the strict
comparison is certified with rational logarithm bounds, not rounding.
Thus these characteristic-zero invariant word families cannot reach the
half rate, even with collective recombinations beyond ordinary type pieces.

This is a restriction on the stated families and the usual width-to-rate
transfer, not on factoring algorithms. It does not cover characteristic
two, arbitrary label maps, binary subsets of a nonbinary word domain, or
families breaking the position symmetry. It supplies no improved factoring
exponent. The proof and checks used substantial OpenAI assistant support,
without new Lean formalization or independent mathematical review.

## 1. Invariant matrices and their large repeated blocks

Regroup the local four-leg tensor space in the row arrangement as
`End(U)`, where `U` has dimension `q=d^2`. A position-invariant tensor of
word length `t` is a matrix commuting with the place-permutation action
of `S_t` on `U^t` (tensor power notation is used below).

We use three standard characteristic-zero results: Schur-Weyl decomposition,
the determinant twist of Schur modules, and the hook dimension formula.
They are given in Etingof et al.,
[*Introduction to Representation Theory*](https://math.mit.edu/~etingof/reprbook.pdf),
Sections 5.17--5.19 and Proposition 5.22.2. The argument below combines them
for the two crossing cuts; that combination is the new claim here. The
[consulted source manifest](../evidence/word-symmetry-sources.json) pins the
PDF hash and exact standard results used.

Schur-Weyl gives

```math
U^{\otimes t}=\bigoplus_{mu\vdash t,\,len(mu)<=q}
 [mu]\otimes S_mu(U).
```

Here `[mu]` has dimension `f^mu`, and the commuting matrix of an invariant
tensor acts as

```math
T=\bigoplus_mu I_{[mu]}\otimes T_mu,\qquad
rank T=\sum_mu f^mu\,rank T_mu.
\tag{3}
```

Hence any nonzero block `T_mu` forces rank at least `f^mu`.
The same statement applies in the column arrangement.

## 2. A common determinant functional for both cuts

For a local tensor `A` on `[d]^4`, let `R(A)` and `C(A)` be its two
`q`-by-`q` matrix arrangements. At word length `t=2qk`, define the
homogeneous polynomial

```math
P(A)=det R(A)^k\,det C(A)^k.
\tag{4}
```

The functional can be defined explicitly from polynomial coefficients. Each
orbit of local-entry tuples is specified by a degree-`t` monomial `m` in
the `q^2` local matrix entries; an invariant tensor has one value `F_m` on
that orbit. Set `Phi_P(F)=sum_m [m]P F_m`. Then
`Phi_P(A^tensor t)=P(A)`. Pure tensor powers span the invariant space in
characteristic zero, so this normalization uniquely determines the functional.
The definition is the same in either regrouping of the four local legs.

In the row regrouping, write `P(A)=det R(A)^k Q(A)`, where `Q` has degree
`qk`. Decompose the functional of `Q` using (3) at degree `qk`. Its values
on a pure power are sums of matrix coefficients of the Schur representations
`S_nu(R(A))`, with `nu` a partition of `qk` into at most `q` parts.

Multiplication of such a coefficient by `det R(A)^k` is a coefficient of

```math
(det)^k tensor S_nu(U)=S_(nu+(k,...,k))(U).
```

Thus, at degree `t`, the row-block support of `Phi_P` consists only of
partitions `mu` with `mu_q>=k`. Equality on pure powers proves the support
assertion for the entire invariant space. This is also a direct way to see
the needed part of the Cauchy decomposition of matrix polynomials.

The column factor in (4) proves the identical support assertion for the
column arrangement. Therefore

```math
Phi_P(F)!=0  ==>  rank_row F>=min_(mu_q>=k) f^mu
                    and rank_col F>=min_(mu_q>=k) f^mu.
\tag{5}
```

This is one common linear functional, not two unrelated tests whose
nonzero summands could be different.

## 3. The crossing has nonzero functional value

Let `e(x,y,u,v)=[x=v][y=u]` on the local alphabet. The global crossing on
all words is exactly `e^tensor t`. Both `R(e)` and `C(e)` are the same
permutation matrix exchanging the two local alphabet coordinates. Their
determinants have equal sign, so their product is one. Equation (4) gives

```math
Phi_P(E_D)=P(e)=1.
```

Linearity now forces `Phi_P(F_i)!=0` for some summand of every exact
invariant decomposition. That particular summand satisfies both inequalities
in (5). The argument places no bound on the number of summands and survives
all linear recombinations within the invariant space.

## 4. An explicit lower bound for the Specht dimensions

Pad `mu` to `q` parts. The hook formula in Vandermonde form is

```math
f^mu=\frac{t!\prod_{i<j}(mu_i-mu_j+j-i)}
              {\prod_{i=1}^q(mu_i+q-i)!}.
```

The numerator product is a positive integer, hence at least one. Also

```math
\prod_i(mu_i+q-i)!
<=\left(\prod_i mu_i!\right)(t+q)^{q(q-1)/2}.
```

Under `mu_i>=k` and `sum mu_i=2qk`, the product `prod_i mu_i!` is maximized
by `((q+1)k,k,...,k)`. Indeed, transferring one unit from a nonfirst part
above `k` to the largest part multiplies that product by
`(mu_1+1)/mu_j>=1`; repeat until all other parts equal `k`.
These inequalities give exactly (1).

Stirling's formula, with `q` fixed, yields

```math
log2 B_(q,k)=t H_q-O_q(log t),\qquad
H_q=log2(2q)-\frac{q+1}{2q}log2(q+1).
```

Apply this to both ranks in (5), and use `lambda=t log2 d`, to obtain (2).
Ceilings and any logarithmic width allowance cannot remove this leading
restriction.

If `t` is not a multiple of `2q`, fix the last fewer than `2q` word
positions to the same alphabet symbol in all four labels. The crossing
restricts to the shorter full-word crossing, every summand remains invariant
under permutations of the remaining positions, and neither matrix rank
increases under restriction. This proves the general-length statement with
the same leading coefficient.

## 5. Permutations within several position classes

The same leading bound holds if the positions are partitioned into `b`
classes of lengths `t_1,...,t_b` and every summand is invariant under
`G=product_j S_(t_j)`. The class partition may depend on `t`. It is enough
that `b=o(t)`; full symmetry across different classes is unnecessary.

In class `j`, keep `m_j=2q floor(t_j/(2q))` positions and fix the rest to
one common symbol. Fewer than `2qb` positions are lost. On the remaining
classes, take the tensor product of the determinant functionals above.
It has value one on the crossing. The representation of the direct product
`G` decomposes into tensor products of the Schur-Weyl blocks, so any block
seen by this one common functional has repeated dimension at least

```math
prod_(j:m_j>0) B_(q,m_j/(2q)).
```

A summand with nonzero functional therefore has both ranks at least that
product. This argument uses the product functional on the whole invariant
tensor; it does not choose potentially different summands class by class.

Stirling's estimate in Section 4 gives the lower logarithmic rank

```math
H_q t-O_q(b+sum_j log2(t_j+1))
>=H_q t-O_q(b log2(2+t/b)).
```

The last step is concavity of the logarithm. If `b=o(t)`, then
`(b/t) log(2+t/b)` tends to zero, so the error is `o(t)` and the same
leading restriction `rho>=rho_d` follows. Fixed bounded weight classes
are one instance of this product symmetry, over characteristic zero.
A number of classes comparable to `t` is outside this conclusion.

## 6. Scope relative to other current barriers

Type coefficients of a tensor power of any fixed seed are invariant under
position permutations. The theorem therefore applies to arbitrary collective
recombinations of those type coefficients on the full word domain. It is
stronger in this particular characteristic-zero setting than checking just
one pure repeated type or the crude sum-of-square-roots estimate.

It is distinct from the sparse-defect invariant barrier in PR #23: here
all words are labels, with no restriction to few defect positions. It does
not automatically apply to the characteristic-two Pfaffian constructions.
Equation (3) uses semisimplicity and its repeated irreducible blocks; that
step must not be imported into positive characteristic without a new proof.
The fixed-characteristic rank phenomena already show that field extensions
can change individual separating-rank assertions.

For seed `d=5`, the present bound has `rho_d<1/sqrt(2)`. This argument alone
therefore gives no half-threshold exclusion from seed five onward. It also
makes no claim about restrictions to an arbitrary binary subset of `[d]^t`:
such a restriction can destroy the required position symmetry.

## 7. Exact arithmetic replay

```sh
python3 tools/check_word_symmetry_barrier.py
python3 tools/check_word_symmetry_barrier.py --output evidence/word-symmetry-barrier.json
make verify
```

The checker compares two exact hook-dimension formulas and verifies the
lower estimate on 11,084 partitions in finite shifted regions. Rational
logarithm series with an explicit geometric remainder certify the three
strict half-threshold comparisons, and a seed-five negative control prevents
extending the numeric conclusion past its scope.

A small independent representation check uses local matrix size two and
word length four: an explicit invariant Johnson-scheme matrix has rank two
and nonzero value for the `det(A)^2` functional, attaining the corresponding
Specht multiplicity. These checks validate the arithmetic and a finite block
example. The all-size support argument and Schur-Weyl step are supplied by
the proof above, not inferred from the finite examples.
