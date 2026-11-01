# A hereditary adjugate cap excludes the retained five-letter target

Every scalar strict crossing over every coefficient field, with at least
two zero-energy labels on each normalized paired alphabet, obeys

```
4 S3(r) + 3 S2(s) <= 12 d,
4 S3(s) + 3 S2(r) <= 12 d.
```

Here `r,s` are its complete normalized Smith profiles, `S_k` is a prefix
sum, and `d` is the normalized minimum-order deficit. This follows by
applying the **existing graded adjugate theorem to a two-letter physical
restriction**, then restoring the core minima and comparing minors.
It does not require parity conservation or fixed local frames.

The retained target at paired energies `(0,0,12,16,18)` and deficit six
has `S2=6`, `S3=15`. The new cap demands `78<=72`, a contradiction.
Thus the previously retained target and its entire dominating average-
prefix sector are not realizable over any field. This supersedes the
parity condition in the [elementary core obstruction](zero-energy-parity-core.md).
It does **not** exclude every target at these energies, every asymmetric
criterion, or the entire [rational output family](joint-six-to-five-retractions.md).
The factoring score remains `1/2`.

## 1. Restrict before applying the old theorem

The input theorem is equation (1) of the
[graded adjugate note](https://github.com/1337isnot1337/integer-factoring/blob/3670b7df8db2bca84bc1dc1ea492c7044fa93e60/research/graded-adjugate-caps.md).
Its proof, graded cofactor identity and rank-one adjugation lemma are
retained as inputs and credited to that contribution. The new step is
closure under physical restrictions and the transfer back to the full
profile. Checking the old inequalities only on the full alphabet misses
this consequence.

Choose two zero-energy symbols in the X,V alphabet and two in Y,U.
Restrict all four indices to these sets. The result is a scalar strict
two-letter tensor, so both cost matrices are generically invertible
four-by-four matrices. Its paired energy lists are `(0,0)`.

Compute its profiles relative to the original full-tensor normalization,
calling them `u,v`. Every core minor is a full-tensor minor, hence

```
S_k(u) >= S_k(r),       S_k(v) >= S_k(s),       k<=4.
```

Let `a=min u>=0`, `b=min v>=0`. Normalizing the restricted tensor gives
profiles `u-a,v-b` and deficit `d_core=d-a-b`. A nonzero matched endpoint
entry and its support bound give `a+b<=d`, so this deficit is nonnegative.

## 2. The core adjugate inequality transfers with a favorable sign

Apply the old theorem with `q=2`, `n=4`, `k=4`. All core cell energies
vanish, and its adjugate shift is `C=3 d_core`. The inequality becomes

```
4 S3(u-a) + 3 (v_2-b) <= 12 (d-a-b).
```

Since `v_2-b=S2(v)-2b`, expansion cancels the unknown zero-core minimum:

```
4 S3(u) + 3 S2(v) + 6 b <= 12 d.
```

The extra term is nonnegative. Dropping it and using the minor-prefix
comparisons proves the first displayed full-profile cap. Exchanging the
places proves the second. This argument permits arbitrary splitting of
the paired energies among the individual leg charges, exactly as the
input theorem does. Integer and rational heights are covered by the
polynomial and ramified formulations.

For `gamma=(r+s)/2`, adding the two caps gives the necessary linear
constraint

```
4 S3(gamma) + 3 S2(gamma) <= 12 d.
```

It therefore belongs in the averaged-profile relaxation as well. A target
with `d=6`, `S2(gamma)>=6`, `S3(gamma)>=15` violates the cap by at least
six. No assumption that the actual two profiles coincide is needed.

## 3. What is excluded and what remains

The target retained in the input note passes its full-alphabet cell,
shared and graded adjugate caps. Its first orders are `0,6,9,9`, so it
fails the inherited two-letter cap. Any row-lattice target guaranteeing
both of those prefix sectors is also excluded, regardless of its frame
choice or polynomial degree. Parity mixing and same-parity zero labels
avoid the elementary parity lemma but cannot avoid this stronger cap.

The exact [six-to-five rational module](joint-six-to-five-retractions.md)
remains a valid construction interface; its output profiles must satisfy
the new inequality. Other prefix allocations, other energy shapes and
restrictions to larger or nonzero-energy cores are separate questions.
This result does not turn a failed target into a stronger factoring bound.

The [checker](../tools/check_zero_energy_hereditary_adjugate.py) verifies
the normalization algebra and the exact `78>72` rejection. It separately
checks 32 dense F4 two-letter tensors with no parity constraint, including
256 complete orders. Eight literal five-letter tensors are restricted
and independently profiled at both places: 464 full and core orders.
Their core prefixes dominate the full prefixes and satisfy the stronger
inequality retaining `6b` in each orientation. The
[receipt](../evidence/zero-energy-hereditary-adjugate.json) pins the inputs
and controls. These examples check the implementation; the restriction
argument and input theorem prove the all-field result.

Run `make verify-zero-energy-hereditary`. This focused result postdates
the full-suite source `f2a16c0`. No Lean work or independent review is
claimed.
