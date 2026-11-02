# High two-power roots as relation certificates

## Question and result

Can supplied high two-power roots of row residues replace known odd
annihilators while retaining **every** two-primary direction of the supplied
subgroup? Yes, with a size bound for the entire ambient group and enough
root depth relative to the determinant. The construction below verifies
the roots and projects the generators in polynomial charged time. It
does not acquire roots or improve the complete factoring exponent `1/2`.

This is an alternative certificate for the
[odd-residual relation interface](odd-quotient-relations.md), not a larger
set of admissible residues. At sufficient depth, a residue has such a root
if and only if it has odd order. The possible algorithmic difference is
acquiring a root without acquiring an odd annihilator.

The next research question is consequently an acquisition question. The
guaranteed-positive specialization below removes row selection if a suitable
deep-root engine can be acquired; more generally, roots and independent
rows could suffice. Every cost must fit the desired complete bound.
Further fixed-group root enumeration would not resolve this question.
Sol worker 2's
[deep-root investigation](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19362022)
is coordinated separately; no factoring or class-group oracle is available
to this construction.

## Projection theorem

Let `G` be a finite abelian group with `|G|<=2^b`, and let
`H=<g_1,...,g_m>`, with `m>=1`. Let `R` be a nonsingular integer `m` by `m`
matrix. Its row residues are

```
h_i = product_j g_j^(R_ij).
```

Write `d=det R=2^s u`, where `u` is odd and may be negative. Suppose
explicit witnesses `z_i` in **G** satisfy

```
z_i^(2^B) = h_i,       B>=b+s.                         (1)
```

The ambient size bound must cover the witnesses, even when they are
outside `H`. Let `A=adj R` and define

```
Z_j = product_i z_i^(A_ji),
w_j = g_j^u * Z_j^(-2^(B-s)).                          (2)
```

Then the `w_j` generate exactly the Sylow two-subgroup `P_H` of `H`, and
every `w_j` has order dividing `2^s`.

**Proof.** The adjugate identity `A R=d I` and commutativity give
`Z_j^(2^B)=g_j^d`. Therefore `(w_j)^(2^s)=1`.
Also `B-s>=b`, so raising any element of `G` to `2^(B-s)` kills its
two-primary component. Thus the correction in (2) has odd order.
Writing `pi_2` for the canonical two-primary projection in the finite
abelian group, we obtain

```
w_j = pi_2(w_j) = pi_2(g_j)^u.
```

The `pi_2(g_j)` generate `P_H`, and raising to the odd power `u` is an
automorphism of that group. Hence (2) generates all of `P_H`, regardless
of whether the individual `z_i` lie in `H`. This also proves that the
output lies in `H`. No discrete logarithms or odd-order factorization
enter the construction. QED.

For `s=0`, the conclusion is `P_H={1}`. Root existence then forces the
supplied subgroup to have odd order; a factor is not promised from such
a subgroup.

## Precision and ambient scope

The depth condition is sufficient and cannot simply be omitted from this
general guarantee. In additive `C4`, take `g=1`, `R=[2]`, `B=1`, `z=1`.
The root equation `2z=2g` holds, but formula (2) gives `w=g-z=0`, losing
the entire two-part. The verifier rejects this depth.

A subgroup size bound is also insufficient when roots can be outside it.
In additive `C16`, take `H=<4>`, `R=[2]`, `B=3`, `z=1`. The root equation
`8z=2g=8` holds. Falsely using `b=2` from `|H|=4` passes the numerical
inequality, yet (2) gives `w=4-4z=0`. The ambient bound requires `b>=4`.
The black-box verifier cannot certify the size of an unknown group from
these identities: membership in the common ambient abelian group and
its size bound remain explicit promises.

The tautological choice `R=2^B I` supplies roots `z_i=g_i`, but has
`s=mB`. It fails `B>=b+s` on any nontrivial ambient group and supplies no
factoring algorithm. If the matrix depends on the chosen depth, its
**actual** determinant valuation must be checked; selecting `B` first
does not establish the inequality for a later matrix.

## Guaranteed-positive deep roots remove row selection

Sol worker 2 contributed the following
[specialization](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19362022).
For any supplied `g` in `G`, put `h=g^(2^b)`. This element has odd order,
so there **exists** a root `z` satisfying `z^(2^(2b))=h`. Apply the theorem
separately to each generator with the one-row matrix `R=[2^b]`, root
depth `B=2b` and `s=b`. For **any** returned root,

```
w = g / z^(2^b) = pi_2(g).                             (3)
```

All displayed powers cost `O(b)` squarings per generator. Thus a charged
algorithm for these guaranteed-positive deep-root queries would retain
all factoring directions without selecting relation rows at all. This
is not the rejected equal-depth choice `R=[2^B]`: the root depth in (3)
is twice the row depth.

The known root `g` at the shallower depth `b` is not a root at depth `2b`
in general. Nor does an arbitrary chain of square-root choices guarantee
continued divisibility: choices can introduce a two-primary component.
No procedure for acquiring the deeper roots is supplied here. The
specialization changes the precise open acquisition target, rather than
solving it.
The [complete all-input transfer](deep-class-root-projection.md) gives the
precise engine promise and polynomial overhead using the existing
generation import and early-stop factor extraction.

## Relation to odd-order certificates

For `B>=b`, the image of the power map `x -> x^(2^B)` is exactly the
odd-order subgroup of `G`. This follows by decomposing `G` into its
two-primary and odd-order factors: the map kills the first and is an
automorphism on the second.

If a positive odd integer `M_i` is supplied and verified to annihilate
`h_i`, choose `e_i` with `2^B e_i=1 mod M_i` and set `z_i=h_i^e_i`.
For `M_i=1`, take `z_i=1`. These are certificates (1).
Conversely, (1) implies that `h_i` has odd order, but it does not supply
an efficiently computed integer annihilator of that order. No faster
root acquisition is inferred from existence. The current
[smooth-square decision](smooth-square-decision-reduction.md) and
[odd-pivot](odd-pivot-two-primary-elimination.md) notes are separate conditional
interfaces, not an acquisition algorithm for (1).

## Algorithm and full invoice

[The implementation](../tools/high_root_projection.py) reuses the
[odd-residual word evaluator](../tools/odd_quotient_relations.py) and the
existing exact fraction-free determinant. It validates the dimensions,
nonsingularity, precision inequality and every root equation, computes
the adjugate by signed minors, and evaluates (2). It returns generators,
not a discrete-log table or a complete two-group basis.

Let `L>=1` bound the bit length of matrix entries and put
`K=O(m(L+log(m+1)))`. Hadamard's bound gives `O(K)` bits for the
determinant, its odd part and each adjugate entry. The signed-minor
implementation uses `O(m^5)` exact integer arithmetic operations on
polynomial-size integers; no fast matrix arithmetic is required.
Evaluating all row residues and adjugate words, taking signed powers,
and checking/projecting the roots costs

```
O(m^2 K + m B) group operations,
```

in addition to the polynomial integer work. In particular, **B squarings
are charged**, not just `log B` bits of its encoding. The implementation
does not materialize the enormous integer `2^B`. A polynomial-size
certificate algorithm must keep the numerical depth polynomial in its
stated parameters. For a fixed matrix, choosing `B=b+s` suffices.

For the quadratic class-group factoring application at discriminant
`-4N`, the existing [early ambiguous extraction](early-stop-class-factoring.md)
can be applied to the `w_j` with annihilator `2^s`. Its polynomial cost
in the input, supplied generator count and exponent length remains
charged, and no factorization of an odd annihilator is needed. The
uniform generator/analytic premises and all preprocessing costs of the
existing reduction remain necessary. A useful two-primary subgroup
suffices for a factor; the theorem does not assert that an arbitrary
supplied `H` contains one. No claim of polynomial generic full-basis
extraction is made.

## Checked evidence and credit

```
make verify-high-root-projection
```

The [checker](../tools/check_high_root_projection.py) and
[receipt](../evidence/high-root-projection.json) contain three fixed
abstract product-group controls. They check a negative determinant and
non-diagonal adjugate, an outside-subgroup root, and the trivial two-part
case. Expected projections are explicit, and the existing independent
two-group routine checks their order logs. Guards retain the low-depth
and false-ambient examples, reject a false root and singular/tautological
rows, and demonstrate linear charged work at depths 64 and 2048.
There is no root-enumeration or class-group acquisition campaign.

The [full default receipt](../evidence/high-root-full-verification.json)
and [complete compressed stdout](../evidence/high-root-full-verification.stdout.txt.xz)
record exit-zero `make verify` at frozen
`be4aae16f0d1e1432b5a142e0ed39b109220992f`, tree
`c4dcf4c9b82dd9f6f5c01080688555747dc8ae9e`, in 1,300.61 seconds.
A fresh marker probe verified private `/tmp`; the host root was read-only,
networking disabled, and only the frozen checkout writable. Tracked state
was clean before and after. The receipt pins raw and compressed stdout
hashes, and the complete log contains the new checker. PR107's merge tree
is identical to that source. The later guaranteed-positive corollary and
this archive follow-up have focused/documentation checks; no executable
source changed after the full replay.

Sol worker 1 contributed the derivation in
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360821);
Sol worker 2 checked the adjugate projection and necessary scope in
[the response](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19361354).
The integration reuses Sol worker 2's odd-residual implementation and
the predecessor determinant, group arithmetic and extraction work, with
their attribution retained. This is ordinary author-checked mathematics
with substantial OpenAI Codex assistance, not a novelty claim for
localization at two, independent mathematical review or new Lean work.
The complete proven factoring bound remains `1/2`.
