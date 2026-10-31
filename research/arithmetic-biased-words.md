# Biased words without a class-group decomposition

## Result and input model

Suppose g_1,...,g_m generate a finite abelian group G, m>=1, and an integer
B>=1 satisfies |G|<=2^B. Given 0<epsilon<1, an explicit multiset of exponent
words with absolute nonprincipal character bias at most epsilon can be
constructed in polynomial time in m, B and 1/epsilon, using deterministic
group operations. Its size is polynomial in those parameters, and its word
l1 norm W satisfies

```math
\log_2(1+W)=O(B+\log(mB)+\log(1/\epsilon)).
```

No group order, invariant factors, discrete logarithms or coordinate oracle
are needed. The known *upper bound* B is sufficient. Repeated elements are
retained with their multiplicity. This closes the biased-word obligation in
the [positive relation kernel](positive-relation-kernel.md), once actual
generators are supplied. It does not construct those generators or extract
smooth pairs.

The argument combines an elementary powering seed below with the published
derandomized-squaring construction of
[Rozenman--Vadhan](https://people.seas.harvard.edu/~salil/research/derand_squaring.pdf)
and generic bias amplification, for example
[Chen--Moore--Russell, Theorem 4](https://arxiv.org/abs/1304.5010).
It is an application of these standard methods, with no novelty claim for
explicit Cayley expanders.

## A seed with a known gap

Begin with the following list of words: 2m(B+1) copies of the zero word and
one copy of each +2^j e_i and -2^j e_i, for 1<=i<=m and 0<=j<=B. Its length
is

```math
D=4m(B+1).
```

Every character average is real and nonnegative because the distribution
is symmetric and has identity mass 1/2.

For a nonprincipal character chi, some chi(g_i)=exp(2*pi*i*a/q) is not 1,
with 2<=q<=|G|<=2^B. Repeated doubling of the angle produces a j<=B for
which its distance from the nearest integer is between 1/4 and 1/2: before
that interval is reached, doubling its distance has no wrap-around. The
real part of chi(g_i^(2^j)) is therefore at most zero. The corresponding
positive/negative pair makes the character average at most 1-2/D.

Consequently the Cayley averaging operator has absolute nonprincipal norm
at most 1-2/D. The constructor [arithmetic_bias_seed.py](../tools/arithmetic_bias_seed.py)
uses only m and B. Evaluating the seed needs O(mB) repeated squarings and
inversions; constructing its explicit exponent vectors has polynomial
output cost. Its maximum word norm is 2^B, with only B+1 bits per exponent.

## Amplification, with all costs charged

Rozenman--Vadhan Observation 4.3 preserves Cayley words under derandomized
squaring. Their Theorem 4.4 gives the spectral bound

```math
\lambda_{new}\le\mu+(1-\mu)\lambda^2.
```

Their Lemma 5.1 supplies uniform explicit constant-degree index expanders
on Q^k vertices, of degree Q and norm mu<=1/100, for one fixed constant Q.
Pad the seed with identities to the next power Q^k. This increases D by
less than the fixed factor Q, preserves nonnegative eigenvalues, and leaves
gap at least 2/(QD). Every subsequent list size is again a power of Q, so
the index-graph sizes match exactly. No unproved all-sizes expander is used.
For gap gamma=1-lambda<1/4, the new gap is at least 3gamma/2. Hence
O(log D) iterations give lambda<=3/4. The list size and original-word
length each grow by a constant factor per iteration, so both are polynomial
in D. A fixed number of ordinary squarings then gives bias <=1/10, still
with polynomial size and word length.

Chen--Moore--Russell Theorem 4 constructs an epsilon-biased multiset from
this constant-bias multiset in polynomial time and size, using group
multiplication alone; its stated size overhead is O(epsilon^-11). The
construction consists of product words in the previous lists. Its word
length is bounded by the same polynomial construction cost. Thus the final
l1 norm is at most 2^B times a polynomial in mB and 1/epsilon, proving the
displayed logarithmic bound. All explicit word-vector arithmetic also has
polynomial bit cost. Equal group elements need not be identified or
deduplicated.

The index expanders are constructed on the *word labels*, not by enumerating
G. Their sizes depend on m, B and epsilon, rather than |G|. This is what
makes the input-model claim meaningful for an unknown class group.

## Effect on the arithmetic contract

Choose epsilon<=delta/2 for a proved smooth-set density lower bound delta.
The positive-kernel sample has polynomial size in m, B and 1/delta. Its
possibly large numerical word norm does not cause a large relation
encoding: K=W+2L+1 has only O(log(W+L)) bits. A full m-by-m relation matrix
therefore still has polynomial encoded size in the stated parameters and
the available smooth-word bit lengths.

For class groups, elementary reduced-form counting supplies an explicit
B=O(log|Delta|) bound without knowing the class number. With m and 1/delta
at exp(O(sqrt(b log b))) scale, this sampler has that same scale with a
different constant. It would fit any fixed logarithmic-saving score. The
unresolved arithmetic tasks are actual generator construction and efficient
smooth-pair witness extraction. No improved deterministic factoring bound
follows yet.

## Source audit and finite evidence

Jalan--Moshkovitz [Theorem 1.1](https://arxiv.org/abs/2105.01149) states a
stronger near-optimal support bound. We do not need that sharper bound here.
In the consulted v1 proof, Proposition 3.9 invokes CMR Theorem 4 for an
initial weakly biased set, whereas the consulted CMR theorem states
amplification from a constant-bias set. The powering-seed argument above
supplies an explicit starting gap and avoids inferring an input model from
that citation. This observation is about the dependency used here, not a
claim that the near-optimal expander theorem is false.

Run `python3 tools/check_arithmetic_bias_seed.py`. Exact integer controls
check all 130,816 nontrivial cyclic character angles with denominator at
most 512. Exact rational controls check 502 characters in elementary
two-groups, the seed gap, and the scalar gap-growth inequality. The sparse
expander construction itself is supplied by the cited theorems and is not
implemented or independently reverified by this finite checker. The
[record](../evidence/arithmetic-bias-seed.json) preserves this scope.

This work used substantial OpenAI Codex assistance and adds no Lean code.
