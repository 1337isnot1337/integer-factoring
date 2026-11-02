# Path C: a uniform density bridge for changes of the relation row

**Result.** Count an event when a smooth hit changes the integer relation
row from the preceding exponent, or follows a miss. These events retain
at least a `1/(2K_p-1)` fraction of the original hit density, uniformly.
This removes the separate raw-hit-cluster-length assumption by changing
the event process, without restricting the original accepted smooth set.
Uniform subgroup density and moment estimates remain open. The complete
factoring score remains **1/2**.

This combines the [fixed-row multiplicity lemma](clustered-hits-and-pooled-rank.md)
with the existing [entry/moment framework](cofactor-orbit-moment-criterion.md).
The [registered scope](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19342034)
uses only previously saved rows and an abstract cyclic control, with no
new arithmetic orbit or cutoff search.

## 1. Define a periodic event, despite the row's drift

Fix a target prime `p` and its cyclic class orbit of period `h`. Let
`H` be the existing reduced-form height, `K=K_p=min{k:p^k>H}` and
`L=2K-1`. At each smooth phase choose the same deterministic, oriented
norm decoder as the selector. Its integer word `v(k)` is periodic:
`v(k+h)=v(k)`, and

\[
 \prod_q q^{|v_q(k)|}\le H,\qquad |v_p(k)|\le K-1.
\]

For all integer exponents with smooth hits, write

\[
 r(k)=k e_p-v(k),\qquad r(k+h)=r(k)+h e_p.
\]

The integer row labels themselves are not periodic. Define an event at
`k` if it is a smooth hit and the preceding exponent is either not a
smooth hit or has a different row. The equality test is

\[
 r(k)=r(k-1)\quad\Longleftrightarrow\quad
                 v(k)-v(k-1)=e_p.
 \tag{1}
\]

The right side uses only periodic data, so the event flag is periodic.
In particular, it handles the boundary between phases `h-1` and zero
without incorrectly comparing `r(h-1)` to `r(0)`.

If a phase produces a proper gcd, count it as an event immediately,
without requiring a decoded word. Treat such a phase separately from a
smooth phase when comparing consecutive words. The original
full-smooth-or-proper-gcd predicate is preserved.

## 2. The density loss is at most polynomial

For a fixed integer row `r`, every exponent at which it occurs satisfies

\[
 r_p-(K-1)\le k\le r_p+(K-1).
\]

Consequently there are at most `L` occurrences of that row over the
entire integer line. Every consecutive segment of smooth hits with
identical row therefore has length at most `L`. Each such segment starts
with exactly one event. A proper-gcd phase is an event of its own.

There is at least one event per cycle whenever the original predicate
has a hit. Otherwise, repeatedly following a smooth hit backwards would
give infinitely many occurrences of its same row, contradicting the
finite multiplicity bound. A proper-gcd hit would already be an event.
This argument also handles an all-good cycle: its rows must change,
because they drift by `h e_p` after a period.

Partition the hits on the cycle into the segments starting at events.
Each has at most `L` points, including segments crossing the period
boundary. If `delta` is the original admissible density and `kappa` the
new event density, this proves the exact uniform inequality

\[
 \boxed{\kappa\ge\frac{\delta}{2K_p-1}.}
 \tag{2}
\]

No assumption on raw hit-run lengths, independent samples, or rational
rank is used. Different events may still have dependent rows. Events
separated by a miss can even repeat the same row; (2) does not require
global uniqueness. The event process is locally computable from two
consecutive norm tests and their decoded words, and requires no known
group order.

The bound is sharp for the abstract periodic-word contract. For example,
in `Z/3` take the target generator one, prime `p=2`, height `H=2`, and
canonical words `v=(0,1,-1)` at the three phases. All phases are smooth,
`K=2`, and the change-event pattern is `(0,0,1)`. Thus
`delta=1`, `kappa=1/3`, and `L=3`. Pure hit-run entries would instead
have density zero. This is a cyclic model, not a new factoring input.

## 3. The maximum zero gap increases by only an additive polynomial

Let `G_raw` be the longest cyclic run of inadmissible phases and `G_event`
the longest cyclic run with no change event. Set the longest zero run to
zero for an all-good flag. Then

\[
 \boxed{G_{\rm raw}\le G_{\rm event}
                   \le G_{\rm raw}+2K_p-2.}
 \tag{3}
\]

The first inequality follows because every event is admissible. For the
second, consider two consecutive events. After the first, its identical-row
hit segment has length `ell<=L` (a proper-gcd segment has length one).
Until the next event, this is followed only by a miss run of length
`G<=G_raw`: a later hit after a miss would itself be an event. Hence
the intervening event-zero run has length `ell-1+G<=L-1+G_raw`.

Thus if every raw window of length `W` hits, every event window of
length `W+2K_p-2` hits. This is a statement about the worst case over
starting phases, not the difference between first-hit times at an
individual start. The all-good cyclic control above attains the additive
gap overhead exactly.

## 4. What remains in the moment route

Every event is an admissible query. The existing sliding-window moment
proof applies to any periodic zero-one flag; it needs no inversion
symmetry. If a uniform theorem supplies original density
`delta>=1/D`, then (2) supplies event density

\[
 \kappa\ge1/((2K_p-1)D).
\]

Use this denominator in the existing moment-window formula, and prove
the lower-tail estimate for **this changed event process**. It then
gives an admissible power for the same selector and its existing complete
factoring transfer. The raw selector can stop at any admissible hit;
it need not evaluate event status. Since `K_p=O(b)`, the density bridge contributes only
`O(log b)` to the logarithm of the cost, preserving every fixed saving
exponent supported by the remaining uniform estimates.

Neither the original subgroup-density bound nor the event moment bound
is proved here. A moment estimate for raw hits, or for pure hit-run
entries on other inputs, is not automatically an estimate for these
events. The new lemma supplies one previously separate density-loss
step; it does not establish a stronger complete runtime.

## 5. Exact saved-data check

The [checker](../tools/check_relation_change_events.py) reconstructs the
canonical words from all 670 saved rows on the existing larger cycle,
and fills the small initial powers using the proved reduced-ideal
criterion `4 p^(2k)<|Delta|`. It compares (1) across every phase,
including the periodic boundary. The original group/order/row/factor
checker is its prerequisite. This postprocessing performs no new class
composition or arithmetic cycle.

| Target | Original hits | Change events | Maximum segment / proved cap | Match previously measured entry flags |
| --- | ---: | ---: | ---: | --- |
| 5 | 343 | 61 | 17 / 17 | Yes, every phase |
| 7 | 343 | 73 | 13 / 13 | Yes, every phase |

Thus the existing finite entry-moment measurements also apply to the
changed process on these two saved directions. This equality is a
finite observation, not a general identity between the two processes.
The maximum segments include the identity arc across phase zero, outside
the earlier `k>=K` nonwrapping multiplicity table.
The existing exact maximum-zero-run records give raw/event gaps
`17939/17939` for target five and `4506/4508` for target seven,
consistent with (3). No new phase enumeration is required for these values.
The checker includes the sharp all-good cyclic control and a
factor-producing-phase control. Its
[receipt](../evidence/relation-change-events.json) pins semantic upstream
data while omitting runtime-only fields.

```sh
make verify-relation-change-events
```

The density proof is an ordinary argument, not an inference from the
finite panel. Existing selector, Grok data/entry analysis and repository
moment credits are retained. No new Lean or independent mathematical
review is claimed.

The [clean frozen replay](../evidence/relation-change-events-full-verification.json)
passes the default and focused gates at exact source
`fd5970d76735c4d4e0f041bc908f3589c5933ee2`, with its complete compressed
log archived. The later gap corollary changes prose only; the subsequent
main integration retains identical event-checker source and has separate
focused and documentation checks recorded in the receipt.
