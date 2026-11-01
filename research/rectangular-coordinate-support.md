# Coordinate support with different alphabets and any piece count

**Ordinary all-field result, 2026-11-01.** The existing
[coordinate-support proof](coordinate-support-rank.md) extends to different
physical alphabets and gives a rank-sum bound without a piece-count loss.
This concerns support-preserving decompositions. It is not a bound on
general crossings, all supporting angles or factoring algorithms. The
complete factoring score remains **1/2**. No new Lean formalization is used.

Let `SX,SY` be nonempty sets of distinct binary words of the same length,
with sizes `mX,mY`, and put `m=min(mX,mY)`. Use `x,v in SX`, `y,u in SY`.
At each bit, permit the crossing pairing `x=v,y=u`, and one of
`x=y,v=u` or `x=u,v=y`, selected by a fixed mask `J` for that tensor.
Let `r(B),c(B)` be the two complete matrix ranks and let

```
d(B) = |{(x,y) in SX*SY : B(x,y,y,x) != 0}|.
```

Then, over **every field**,

```math
r(B)c(B)\ge d(B)^2/m.                                      \tag{1}
```

If finitely many such tensors sum to a nonzero scalar multiple of the
crossing, each may have a different fixed mask. For any number of pieces,

```math
S:=\sum_i\sqrt{r(B_i)c(B_i)}
\ge \frac{m_Xm_Y}{\sqrt{\min(m_X,m_Y)}}.                    \tag{2}
```

Thus the balanced finite-block rank-sum parameter satisfies

```math
\rho_*:=\frac{\log S}{\log m_X+\log m_Y}
\ge 1-\frac{\min(\log m_X,\log m_Y)}
              {2(\log m_X+\log m_Y)}\ge\frac34,            \tag{3}
```

when `mX*mY>1`. There is no assumption of common physical types or equal
alphabet sizes. At equal sizes this recovers the predecessor's exponent;
unequal sizes strengthen this normalized rank-sum restriction.

## Proof

Every supported entry satisfies `x xor y=v xor u`. For each XOR difference
`a`, put `Xa={x in SX : x xor a in SY}`. Both matrices split into blocks
with row and column indices `Xa`. In particular `|Xa|<=m`.

Let `Da` be the first coordinates of nonzero crossing-diagonal entries
in that block and `Ta=support(a xor J)`. Exactly as in the credited
square proof, nonzero row entries preserve projection onto `Ta` and
nonzero column entries preserve projection onto its complement. Each
occupied projection class has a nonzero diagonal entry, so

```
r_a >= |projection_Ta(Da)|,
c_a >= |projection_complement_Ta(Da)|,
r_a*c_a >= |Da|.
```

The two projections determine the word. Rank additivity across XOR
blocks and Cauchy-Schwarz give

```
sqrt(r(B)*c(B)) >= sum_a sqrt(r_a*c_a)
                >= sum_a sqrt(|Da|)
                >= sum_a |Da|/sqrt(m) = d(B)/sqrt(m).
```

Squaring proves (1). Every one of the `mX*mY` crossing-diagonal entries
of a nonzero crossing sum occurs in at least one piece. Therefore
`sum_i d(B_i)>=mX*mY`. Summing the unsquared inequality proves (2), and
taking logarithms proves (3). All inequalities compare ordinary integer
ranks and counts; no ordered-field structure is assumed for coefficients.

## Source audit: diagonal regularized transpose words

For an invertible diagonal `D` over a characteristic-two field define

```
S_D(M)=D M^T D^(-T),       F_D(z)=I+(1+z) S_D.
```

Each factor has `F_D(1)=I`. Its four-index tensor is supported on just
the crossing `x=v,y=u` and the transpose pairing `x=u,y=v`. Products
continue to preserve each unordered coordinate pair, so the same two
global pairings support every polynomial coefficient. Tensor powers,
coordinate label restrictions, monomial physical gradings and coefficient
linear recombinations retain the corresponding fixed coordinate mask.
Consequently (2) excludes an improved **balanced rank-sum decomposition
whose pieces remain in that coefficient support**, including different
coordinate alphabets on the two pairs and any piece count.

This audit does not automatically constrain a procedure producing pieces
outside the support. Non-diagonal transpose channels, arbitrary physical
retractions, new polynomial generators, mixed masks inside one piece,
and an all-angle recursion analysis require separate arguments. In
particular no general weighted-Smith pressure obstruction is inferred
merely from zeros in the original source tensor.

The explicit restriction to diagonal factors matters. For any invertible
`A`, `S_A(M)=A M^T A^(-T)` is an involution and `I+S_A` is square zero.
Its image has dimension `q(q-1)/2`. Products of regularized **non-diagonal**
channels need not stay in the support above. A low-rank product at zero
remains low rank under further left/right multiplication, while an odd
number of full transpose channels has rank-one opposite realignment at
the leading coefficient. These are construction facts to test, not a
favorable profile or escape certificate.

## Finite controls and credit

[The checker](../tools/check_rectangular_coordinate_support.py) independently
checks the one-bit support truth table, 600 arbitrary unequal-alphabet
tensors over F2, F3 and F5, and complete coefficients of four regularized
diagonal-transpose products over F4. The latter use binary restriction of
scalars for their ranks, literal operator products and full endpoint
checks. Their coefficient diagonal counts cover the whole crossing.
The mixed-mask negative control violates the claimed inequality when
the support hypothesis is removed, as required.

The [receipt](../evidence/rectangular-coordinate-support-controls.json)
records the alphabets, masks, ranks, diagonal counts and source hashes.
The all-field/all-size conclusion follows from the ordinary proof, not
from the finite cases. The original support theorem and its XOR/projection
argument are credited directly; the contribution here is the unequal-pair
and all-piece-count formulation plus the source audit. Directed by Cinna
Davis / 1337isnot1337 with substantial OpenAI assistant support. No
independent mathematical acceptance or stronger factoring bound is claimed.
