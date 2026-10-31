# A positive Fourier kernel for arithmetic relation search

## A proved translate-hitting lemma

Let G be a finite abelian group and let S be a subset of density delta>0.
With averages normalized by |G|, define

```math
C(g)=\mathbb E_x[\mathbf1_S(x)\mathbf1_S(xg)].
```

If mu is a probability distribution with absolute nontrivial character bias
at most epsilon, then for every target t in G,

```math
\mathbb E_{u\sim\mu} C(tu)
\ge\delta^2-\epsilon(\delta-\delta^2).
```

For 0<delta<1, epsilon<delta/(1-delta) guarantees a support element u for
which tu=a b^-1 with a,b in S. Taking epsilon<=delta/2 is sufficient
uniformly for 0<delta<=1; the case delta=1 is immediate.

This is a positive existence/hitting statement for **quotients of two
elements of S**. It does not say that the sample hits S itself, or that the
pair a,b can be extracted quickly. Those distinctions are essential for a
factoring application.

## Fourier proof

The normalized Fourier coefficient of C at a character chi is
|hat(1_S)(chi)|^2, so all its coefficients are nonnegative. The principal
coefficient is delta^2. Their total sum is C(1)=delta, and hence the
nonprincipal sum is delta-delta^2. Multiplying each coefficient by the
sampled character value and using its absolute bias gives the inequality.
Translation by t only changes unit-modulus character phases.

The key gain is the Fourier l1 bound delta, independent of |G|. The
unsmoothed indicator 1_S need not have a comparable Fourier l1 norm.
Standard Fourier facts prove this lemma; no novelty claim is made for
positive-definite autocorrelation. Its proposed arithmetic use, and the
necessary extraction interface below, are the research direction.

## Deterministic relation rows without a rank-growth lottery

Suppose g_1,...,g_m generate G. Every element a in S has an available
integer word v_a in these generators with l1 norm at most L. Every sample
element u has a known word w_u of norm at most W. Assume a deterministic
procedure can extract u,a,b from the guaranteed intersection for any target.

Put K=W+2L+1. For each i apply the procedure to t=g_i^K and form

```math
r_i=K e_i+w_u-v_a+v_b.
```

It is a relation. Its correction to K e_i has l1 norm at most W+2L<K.
Thus the m-by-m relation matrix is strictly diagonally dominant and has
nonzero determinant. A maximum-coordinate argument proves nonsingularity:
if Av=0 and |v_i| is maximal and nonzero, its diagonal term strictly
dominates the sum of all other terms in that row.

The [two-part extraction theorem](arithmetic-torsion-extraction.md) then
uses the determinant as an annihilator. The relation matrix need not be a
saturated presentation and its relation lifts need not be uniformly random
modulo twice the hidden kernel. The remaining generator and pair-extraction
requirements are explicit assumptions, not consequences of smoothness
density.

## Class-group interpretation

For the negative-discriminant class group, let S consist of reduced forms
(a,b,c) with gcd(a,Delta)=1 and B-smooth first coefficient a. Its prime
factorization gives a short word in the factor-base prime forms; its l1
length is O(log|Delta|). The recognizable-smooth subset in
Lenstra--Pomerance supplies a lower-density result for the discriminants
satisfying their hypotheses.

The pair-extraction target is to find two such small smooth forms whose
quotient is t times one of the sampled classes. This is a problem about
smooth ideal/norm representations. The positive Fourier bound proves that
a pair exists for a suitable small-bias sample. It does not implement a
fast arithmetic witness extractor. Enumerating the entire class group or
all smooth forms would generally have exponential bit cost and would not
improve the factoring score.

Likewise, a small-bias set for a group supplied in known coordinates is
not automatically an efficient construction for an unknown class group.
The input assumptions of [Jalan--Moshkovitz's abelian Cayley expander work](https://arxiv.org/abs/2105.01149)
and related amplification constructions must be audited before using them
as this algorithmic primitive.

## Exact finite pipeline

The [checker](../tools/check_positive_relation_kernel.py) uses G=(Z/2Z)^4
and a six-element quadratic small-bias set with epsilon=1/3. A five-element
test set S has density 5/16. Exact rational Fourier coefficients verify
the autocorrelation identity, positive coefficients, and the lower bound
5/192 for every translate; the actual minimum is 1/24.

The sample deliberately misses S itself. It nevertheless hits every
translate of S-S and yields all sixteen pair witnesses. These witnesses
are found by exhaustive finite search; this is not the missing fast
arithmetic extractor.

The checker constructs a strictly diagonally dominant relation matrix
with determinant 5,120 for the group of order 16. Halving those relation
rows reaches only two torsion elements. Determinant projection and the
generic two-group basis routine recover the full rank-four two-torsion.
This validates the presentation-to-torsion replacement even with a very
unsaturated relation lattice.

The [record](../evidence/positive-relation-kernel.json) retains the set,
sample, exact biases, kernel bound, relation matrix and extraction result.
Run `python3 tools/check_positive_relation_kernel.py` to reproduce it.

These proved components and finite controls do not establish a new
deterministic factoring bound. The generator, bias-construction and smooth
pair-extraction obligations remain open. This work used substantial
OpenAI Codex assistance and adds no Lean formalization.
