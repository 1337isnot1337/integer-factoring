# Correct the entry checker's lower-tail sign

## Error and corrected scope

At source `89aeec9`, both `lower_tail_record` and `required_constant` in
`tools/check_entry_moment_cycles.py` formed `good*W-h*count`, then selected
negative deviations. This computes the **upper** tail while labeling it
as the lower tail. The correction changes both deviation definitions to
`h*count-good*W`, so the existing negative-tail filter computes exactly
`E(mu-X)_+^(2q)`.

The minimal asymmetric control is the flag `(1,0,0)` with period three,
`W=q=1`. Its true lower moment is `2/27`; the old helper returned the
upper moment `4/27`. The corresponding normalized required constants
are `2/9` and `4/9`. The new regression failed on the old code and passes
after the correction. A separate existing implementation,
`check_cofactor_orbit_moments.moment_record`, already used the correct
convention and is unchanged.

## Evidence repaired without new walks

The [canonical tables](entry-moment-cycles.md) and
[receipt](../evidence/entry-moment-cycles.json) are recomputed from the
identical two archived flag files. No input, cutoff, phase, order, window,
period, hit/entry count or composition budget is changed. The new
`--saved-flags` replay checks the existing flag hashes and reprofiles them;
the historical 219306-walk count is checked as metadata, not executed.

Exactly 83 scalar fields change: 42 required-constant decimals,
10 exact tail numerators, 16 `C=1` flags and 15 `C=2` flags. All ten
nonwrapping sufficient-window inequalities still pass. All 32 corrected
matched-mean tests pass `C=2`; 30 pass `C=1`. The two `C=1` failures are
the raw mean-eight tests at the smaller order on each cycle.
The [correction receipt](../evidence/entry-tail-sign-correction.json)
pins the old/new sources, saved flags, field counts and separate lower-tail
cross-checks.

The former upper-tail numbers remain accessible in the
[exact pre-repair receipt](https://github.com/1337isnot1337/integer-factoring/blob/89aeec9ffc8a8b34af19ef45650b8ff520b33c40/evidence/entry-moment-cycles.json).
They are not discarded or promoted to lower-tail evidence.

The former claim that raw matched-mean lower tails generally fail `C=2`
is withdrawn. The older tables measured upper deviations. Also,

```math
\mathbb E(\mu-X)_+^{2q}\le\mu^{2q},\qquad
C_{\rm req}\le\mu/q,
```

so small lower-tail constants at low mean/high order are not by themselves
evidence of independent hits or a uniform empty-window estimate. The
unchanged empty fractions and gaps still distinguish raw hits and entries.
In fact this universal bound already implies all 32 matched-mean `C=2`
passes and 28 of the 30 `C=1` passes. The two remaining `C=1` passes are
the entry mean-eight tests at the smaller order. Those are finite findings,
not an asymptotic moment theorem.

The abstract lower-tail moment theorem is unaffected. The independent
clustered-cycle lower-tail checker and saved arithmetic words are
unaffected. Sol1's cumulant and word-boundary work imports only cyclic
counts/event flags from these modules, and its own exact calculations do
not use the two erroneous tail helpers. The complete factoring score
remains **1/2**.

## Verification and historical receipts

```sh
python3 -m unittest discover -s tests -p test_entry_moment_tails.py -v
make verify-entry-moment-saved
make verify-docs
```

The saved-flag replay invokes no group-composition routine. The original
`make verify-entry-moment-cycles` remains available for a complete arithmetic
replay with the corrected statistics; it is not rerun to regenerate this
correction. Historical full receipts remain intact at their source pins:
they establish replay of the old code, not correctness of its tail labels.
The regression and independent-helper checks are meaningful semantic
checks in addition to matching a stored receipt.

This is a correctness repair, not a new factoring result or independent
mathematical review. The original selector, registered experiment and
Grok/independent-cycle contributions retain their source credit.
Developed with substantial OpenAI Codex assistance. No Lean work.

The [clean frozen default replay](../evidence/entry-tail-sign-full-verification.json)
passes at exact source `de0c57c7d803fc17a32611a8116e011a607bff7f` in
2494.28 seconds. Its [complete compressed log](../evidence/entry-tail-sign-full-verification.stdout.txt.xz)
is archived. The namespace has private `/tmp`, read-only host root and
only the frozen checkout writable. The later automatic-bound clarification
and current-main integration retain both corrected source hashes in the
correction receipt, and pass saved-flag/boundary integration and docs checks
separately. This is code verification, not independent mathematical review.
