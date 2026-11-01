# Smooth represented norms give short relation words directly

## Result and scope

The [bounded smooth-norm interface](bounded-smooth-norm.md) characterizes
membership in a quotient of two sets of small reduced smooth forms. Its
balanced-divisor search is needed for that exact characterization. It is
not needed by the [positive-kernel relation construction](positive-relation-kernel.md).

A supplied primitive representation of a smooth integer already gives a
short word for the target class. The relation algorithm can use that word
without producing two reduced forms. More generally, its useful budget is
the number of prime factors counted with multiplicity, rather than an
upper bound on the represented integer or on two balanced divisors.

This note proves the direct converter and a relaxed, sufficient search
contract. It removes the previous `2^O(b/log b)` divisor-enumeration step
from this contract and accepts witnesses that the balanced converter
rejects. It does not supply a fast algorithm for finding a witness, and
does not improve the complete factoring exponent `1/2`.

The base repository revision is
`a5e05feb756b4ad35ae1877d54a5f2e4eb387c19`. The ideal/form correspondence
and signed prime-form factorization are classical; see Lenstra--Pomerance,
*A rigorous time bound for factoring integers*, JAMS 5 (1992), Sections
2.7--2.8, in the [primary paper](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf).
The new work here is the weaker extraction contract and its integration
with the existing diagonally dominant relation construction.

## 1. Direct conversion of a represented norm

Let `Delta<0` be a discriminant. Let `Q=(A,B,C)` be a primitive positive
definite form of that discriminant. Fix one orientation of every split
prime form of norm `p<=Y` with `p` coprime to `Delta`; these are the
factor-base generators. An inverse orientation is represented by a
negative exponent of the chosen generator.

Suppose a supplied primitive vector `(x,y)` satisfies

```math
n=Q(x,y)>0,\qquad \gcd(n,\Delta)=1,
\qquad n=\prod_{p\le Y}p^{e_p}.
```

Then one can construct a signed exponent vector `v` such that

```math
[Q]=\prod_{p\le Y}[P_p]^{v_p},
\qquad |v_p|=e_p,
\qquad \|v\|_1=\Omega(n)=\sum_p e_p.
```

If the factor base keeps both orientations as distinct coordinates,
use nonnegative exponents instead. Aliases between distinct coordinates
do not affect the norm bound or the relation argument.

### Proof and algorithm

Choose integers `u,v` with `x*v-y*u=1`, using the extended Euclidean
algorithm. The proper change of variables with first column `(x,y)`
has middle coefficient

```math
B'=2Axu+B(xv+yu)+2Cyv.
```

It transforms `Q` into `(n,B',C')`, with
`B'^2-4nC'=Delta`. For each `p` dividing `n`, take the prime form

```math
P'_p=\left(p,b_p,\frac{b_p^2-\Delta}{4p}\right),
\qquad b_p\equiv B'\pmod{2p}.
```

The displayed third coefficient is integral. This also covers `p=2`:
if `2` divides `n` and is coprime to `Delta`, then `Delta` is `1` modulo
`8`, and the residue of `B'` modulo `4` specifies the orientation. No
modular square-root search is required for these supplied witnesses.

The primitive ideal corresponding to `(n,B',C')` has norm `n`. Since
`n` is coprime to `Delta`, its prime-ideal factors are unramified and
coprime to the conductor. Primitiveness prevents both conjugate
orientations above one rational prime from occurring: their product is
the rational ideal `(p)`, which would give integer content. Consequently
this ideal is the product of the oriented prime ideals above `p`, each
to exponent `e_p`. Passing to ideal classes proves
`[Q]=product_p [P'_p]^e_p`. Choosing the sign relative to the stored
orientation yields the stated word.

A nonprimitive input vector may first be divided by `gcd(x,y)`. All
conditions and the multiplicity budget must be checked on the resulting
primitive norm. The zero vector is rejected.

### Cost

Let `b=bits(|Delta|)` and let the normalized represented integer have at
most `H` bits. For a reduced input form, positive definiteness bounds
the coordinates by `2^O(b+H)`, since

```math
(2Ax+By)^2+|\Delta|y^2=4An.
```

Extended Euclid, change of variables, orientation selection and form
reduction therefore have polynomial bit cost in `b+H`. Trial division
through `Y` verifies smoothness and produces the exponents in

```math
Y\,\mathrm{poly}(b+H)
```

bit operations. If a valid prime factorization is supplied, its
verification and the conversion instead have polynomial bit cost in
`b+H+Y` (or in the factorization length when primality certificates are
also supplied). Efficient standard form arithmetic gives the same cost
for an optional independent check of the resulting word.

There is no divisor enumeration and no balanced-partition subroutine.
The finite control implementation uses the repository's elementary
form-composition routine; it is not a benchmark of the efficient
asymptotic composition algorithm used in this cost statement.

## 2. A sufficient relation-search contract

Let `g_1,...,g_m` generate the class group. The list includes the
factor-base prime forms needed by the direct converter. Suppose a
finite sample `U` is supplied with known words `w_u` satisfying
`||w_u||_1<=W`.

Fix an integer `L>=0` and set

```math
K=W+L+1.
```

For each `i`, the relaxed extraction task is to find an element `u` of
the sample and a primitive represented norm of the reduced class

```math
Q_i=[g_i]^K u
```

such that the norm is coprime to `Delta`, is `Y`-smooth, and has
`Omega(n)<=L`. The direct converter supplies a word `v_i` for `Q_i`
with `||v_i||_1<=L`. Thus

```math
r_i=K e_i+w_u-v_i
```

is a relation. Its correction to `K e_i` has norm at most `W+L=K-1`.
The square matrix with rows `r_i` is strictly diagonally dominant:

```math
|r_{ii}|-\sum_{j\ne i}|r_{ij}|\ge1.
```

It is nonsingular, so its determinant annihilates the whole group by
the [existing adjugate argument](arithmetic-torsion-extraction.md).
Each row has norm at most `2K-1`; Hadamard's inequality gives

```math
0<|\det R|\le(2K-1)^m.
```

Hence the annihilator has at most
`1+m*ceil(log2(2K-1))` bits. This explicit size bound does not require
a saturated relation lattice or an exact class number.

If every sample is one prime-form generator or its inverse, take `W=1`
and `K=L+2`. This statement is conditional only on the stated sampling
and extraction hypotheses; it does not assume an unproved sampling
theorem.

## 3. Why the old density guarantee still supplies a witness

Let `S_R` be the set from [bounded-smooth-norm.md](bounded-smooth-norm.md):
reduced forms of leading coefficient at most `R`, coprime to `Delta`
and `Y`-smooth, with `4R^2<|Delta|`. Suppose its density is `delta>0`
and the sample has absolute nonprincipal character bias at most
`epsilon`, where `epsilon<delta/(1-delta)` when `delta<1`.

The positive-kernel lemma gives, for every target `t`, an element `u`
in the sample with `t*u=a*b^-1` for `a,b` in `S_R`. The ideal-content
argument in the bounded-norm note gives a primitive represented norm

```math
n\le R^2,
\qquad \Omega(n)\le\lfloor\log_2(R^2)\rfloor.
```

Thus choosing `L=floor(log2(R^2))` proves existence for the relaxed
contract. The algorithm can compute this integer using bit lengths;
no real logarithm is needed. The case `R=1` has `L=0` and is harmless
when the stated density and sampling hypotheses hold.

The search may accept many more witnesses than this proof guarantees.
In particular, `Omega(n)<=L` implies only `n<=Y^L`, so the candidate
region can include represented integers much larger than `R^2` and
ones with no balanced divisor pair. Their bit length is bounded by
`1+L*ceil(log2 Y)`, which is polynomial on the proposed arithmetic
parameter scales. This enlargement requires no new density assumption.

## 4. Two controls that matter

For `Delta=-999`, the reduced form `(2,-1,125)` represents `125` at
`(0,1)`. With `Y=7` and `R=12`, there is no divisor `d` of `125` for
which both `d` and `125/d` are at most `12`. The old balanced converter
therefore rejects this witness. The direct converter accepts it: `125`
is coprime to `999`, and its multiplicity is `3`. The budget supplied
by the old existence argument is `floor(log2(144))=7`.

The norm or word-length budget cannot simply be dropped. Given a target
`g_i^K u`, its already-known input word has length at most `K+W`.
Returning that same word makes `K e_i+w_u-v_i` zero. Smoothness alone
does not supply a full-rank relation matrix. The strict inequality
`||v_i||_1<K-W` is the step that proves progress.

## 5. Remaining task and reproduction

The arithmetic search must still locate one of the guaranteed primitive
smooth represented norms in the sampled target family. Enumerating all
integer vectors in the larger region, all smooth norms, all reduced
forms, or all short words can take exponential time in `b`. None of
those enumerations is claimed to meet the desired bound.

The direct conversion, strict diagonal dominance, determinant bound and
relaxed existence guarantee are ordinary proved statements. Finite
controls are in [check_smooth_norm_words.py](../tools/check_smooth_norm_words.py),
with a retained [record](../evidence/smooth-norm-words.json). Run:

```sh
python3 tools/check_smooth_norm_words.py
```

The checker validates supplied witnesses and finite relation controls;
it does not implement the missing asymptotic search. This work used
substantial OpenAI Codex assistance. No new Lean formalization is added.

The retained replay checks 2,366 primitive supplied representations and
4,732 normalized scaled controls across 17 discriminants. Among the
accepted representations, 592 exceed the old bound `R^2`, and another
244 satisfy that bound but have no old balanced divisor split. Both
orientations above `2` occur. Fourteen nonempty factor bases also have
explicit dominant relation matrices whose determinants annihilate every
listed generator; three discriminants have no factor-base prime at this
small test bound and are reported as such. These finite enumerations are
validation of the converter and relation algebra, not evidence of fast
witness search.
