# The global divided square of quartic charge creation

The quartic creation operator in [PR #97](https://github.com/1337isnot1337/integer-factoring/pull/97)
has an exact global simplification at the eight-color alphabet. Its full
product is linear in the quartic coefficient, rather than merely linear
through a fixed jet order. At larger alphabets there can be a quadratic
term, with a closed binary relation-code formula. These are construction
identities, not cut-rank estimates or a stronger complete factoring bound.
The proven score remains **1/2**.

The charge-fusion construction and its first-layer mechanism belong to the
[primary worker's proposal](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19357910)
and PR #97, which retain the earlier homogeneous Gaussian and worker-4
color-class sources. Sol worker 2 supplied the complementary global
identity and controls below, first reported in
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19358338).
The original scope's degree-five polynomial was superseded by PR #97's
degree-four-per-occupied-vertex construction. The identity here works for
either polynomial and, more generally, any common coefficient `H`.

## 1. The operator expansion

Let the eight spatial vertices be `V={0,...,7}`. Physical colors are
`F_2^s`, with zero as vacuum and `q=2^s`. Coefficients lie in an arbitrary
commutative ring `R` of characteristic two. For every spatial four-set `S`,
the operator `C_S` replaces vacuum on `S` by four distinct nonzero colors
whose XOR is zero, summing each assignment to the four fixed vertices
once. It annihilates inputs that are nonvacuum anywhere in `S` and leaves
other vertices unchanged. This is precisely the creation convention of
PR #97; it uses the additive color group, without requiring an embedding
of `F_q` in `R`.

Creations on overlapping spatial supports compose to zero; those on
disjoint supports commute. There cannot be three disjoint four-sets in
eight vertices. Consequently, for `H in R`,

```math
\prod_{|S|=4}(I+H C_S)=I+HU+H^2D,\qquad
U=\sum_{|S|=4}C_S,\qquad
D=\sum_{\{S,V\setminus S\}} C_S C_{V\setminus S}.       \tag{1}
```

The last sum is over **35 unordered** complementary pairs. Its operator
`D` annihilates every input except the all-vacuum word. At an output word
`c=(c_0,...,c_7)`, its coefficient is the parity of the number of partitions
of the eight positions into two unordered four-blocks, each carrying four
distinct nonzero colors of XOR zero. A zero color makes this number zero.

Although `U^2=0` in characteristic two, this does **not** imply `D=0`:
the ordinary square counts every disjoint pair twice. We call `D` the
divided-square term to emphasize the unordered-pair convention, not to
identify it with the ordinary square.

## 2. Repeated-color cancellation and the exact q8 identity

Suppose two positions `i,j` carry the same nonzero color. Every valid
partition puts them in different blocks, since each block has distinct
colors. Fix these two positions independently of the partition and
exchange them between the blocks. Both color multisets remain unchanged,
so the result is another valid partition.

This is an involution with no fixed unordered partition. A block changes
by the two-element symmetric difference `{i,j}`, so it does not stay the
same. Nor can it become its complement: a four-set and its complement
have symmetric difference of size eight. Thus valid partitions pair off
and the coefficient of `D` is zero. If a color occurs three or more times,
there are no valid partitions in the first place.

There are only seven nonzero colors when `q=8`. Every output with eight
nonzero colors therefore repeats a color, proving the **operator identity**

```math
\boxed{\prod_{|S|=4}(I+H C_S)=I+H\sum_{|S|=4}C_S
       \quad(q=8).}                                  \tag{2}
```

It holds on all input words, over every characteristic-two coefficient
ring, without truncating `H` or the spatial variables. The same vanishing
holds whenever there are fewer than eight nonzero colors. Multiplying
(2) by the Gaussian pair-creation operator gives the exact q8 charge-fused
tensor, not merely its initial jet.

The common coefficient `H` is essential. If a single gate on `S_0` has
weight `H+Delta` and the others still have weight `H`, its q8 product is

```math
I+HU+\Delta C_{S_0}+H\Delta C_{S_0}C_{V\setminus S_0}.
```

The last term has coefficient `H*Delta` on the output with colors
`1,2,4,7` in each block of that partition. Thus (2) must not be used for
arbitrary spatial weights. For PR #97's nonzero polynomial `H`, the
endpoint-preserving variation
`Delta=(xyz)^5(x+y)` in the later
[spatial-quartic family in PR #100](https://github.com/1337isnot1337/integer-factoring/pull/100) gives a nonzero
term. Its order is at least nine in each cut, so this scope distinction
does not affect that family's order-five obstruction. This was reported
in the [focused inspection](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19358929).

## 3. Eight distinct colors and the binary code formula

Now suppose the eight output colors are distinct and nonzero. If their
total XOR is nonzero, no two XOR-zero blocks can partition them. Otherwise
form the binary matrix with columns `(1,c_i)`, and put

```math
k=8-\operatorname{rank}_{\mathbb F_2}[(1,c_i)]_{i=0}^7.
```

Its kernel is the code of even-cardinality zero-XOR subsets of positions.
The all-one word belongs to this code because all eight colors XOR to
zero. Distinctness excludes weight two. Complementation by the all-one
word then excludes weight six. The only possible weights are zero, four
and eight. Exactly one codeword has weight zero and one has weight eight,
so there are `2^k-2` weight-four codewords. Complementation pairs them
without a fixed point. The integer number of valid unordered partitions is

```math
\frac{2^k-2}{2}=2^{k-1}-1.                            \tag{3}
```

Here `k>=1`. Thus `D[c]=1` precisely when the colors are distinct and
nonzero, their total XOR is zero, and `k>=2`, equivalently the augmented
rank is at most six. These conditions determine every coefficient of `D`
uniformly in `s`. A coefficient needs only 35 literal quartet checks or
binary elimination on eight columns of `s+1` bits; its bit cost is `O(s)`.
Enumerating the entire local vector still costs `O(q^8 s)` by this method.
This entry formula does not give a faster global coefficient contractor.

Two guards distinguish this result from simpler invalid formulas:

- At `q=16`, the affine three-plane with colors `8,...,15` has augmented
  rank four and seven valid partitions, so `D` is nonzero.
- The colors `1,2,4,8,16,32,64,127` have total XOR zero but augmented rank
  seven. Their only nonempty even zero-XOR subset is the full set. They
  have no valid quartet partition, so total XOR alone is insufficient.

## 4. Independent finite controls and the next obligation

```sh
make verify-charge-creation-square
```

The standard-library [checker](../tools/check_charge_creation_square.py)
directly expands every q8 disjoint gate pair. It enumerates 168 ordered
quartets and 35 spatial partitions, accumulating **987,840** integer terms.
The 229,320 supported output words all have even multiplicity: 211,680 have
multiplicity four and 17,640 have multiplicity eight. This calculation uses
neither the repeated-color involution nor the code formula.

A separate literal partition count tests the relation-code formula on all
6,435 eight-element subsets of the fifteen q16 nonzero colors. Of these,
6,000 have nonzero total XOR, 420 have three valid partitions, and fifteen
have seven. The checker also tests the fixed-pair involution on repeated
colors and the two guards above. The [receipt](../evidence/charge-creation-square.json)
pins the checker and deterministic record hashes. These finite controls
support the implementation; the all-alphabet statements follow from the
proofs, rather than extrapolation of the panel.

The [combined full default receipt](../evidence/charge-reciprocity-full-verification.json)
and [compressed log](../evidence/charge-reciprocity-full-verification.stdout.txt.xz)
record exit-zero `make verify` at clean frozen source `a21c955` in 888.40
seconds, inside a marker-checked private `/tmp` namespace. It includes
PR96/97/98/99; the later common-H caution and PR100 are outside that exact
pin. The caution has a literal eight-partition check and documentation
inspection, while this checker is byte-identical to the full source.

This checkpoint simplifies the full q8 representation and controls the
extra term for growing alphabets. It does not determine higher Smith
orders, the normalized pressure with the restored degree cost, or a
complete factoring transfer. Continued work must address one of those
structural obligations; a larger field or generic-profile sweep is not
justified merely by (2). No new Lean work or independent mathematical
review is claimed.
