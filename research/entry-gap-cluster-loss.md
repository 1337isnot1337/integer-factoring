# Entry gaps and the cluster-length loss in the Path C moment criterion

## Result and scope

Let `f` be the cyclic admissible-class indicator of one generator orbit,
with period `h` and `s` good classes, and put

\[
 e(k)=f(k)\bigl(1-f(k-1)\bigr). \tag{1}
\]

The one entries of `e` are the starts of the hit runs.  This note records
an exact structural fact about them:

> **Entry-gap lemma.**  Let `L_r` be the length of hit run `r` and `G_r`
> the number of zeros immediately after it.  The cyclic gaps between
> consecutive entries of `e` are `d_r=L_r+G_r`.  They satisfy
> `d_r>=2` and `sum_r d_r=h`.  For every integer `W>=1`, a cyclic window
> of length `W` with no entry exists if and only if
> `max_r d_r >= W+1`.  Consequently, if `W>=max_r d_r`, every cyclic
> window contains an entry and therefore a hit.

This is elementary, but it changes the interpretation of the lower-tail
moment criterion in [the cofactor window note](cofactor-orbit-moment-criterion.md).
That criterion applied to `e` is a sufficient certificate for
`W>=max d_r`; it is not a way to obtain a window smaller than the maximum
entry gap.  In particular the exact entry-gap window is never smaller
than the raw maximum-zero-run window, and it exceeds it by at most the
maximum cluster length.

The second point is the **cluster-length loss**.  Write
`D_r=ceil(h/s)` and `D_e=ceil(h/runs)`, where `runs` is the number of hit
runs.  Since `mean_cluster=s/runs`,

\[
 D_e=\left\lceil \frac{h}{runs}\right\rceil
    \le D_r\,\left\lceil\frac{s}{runs}\right\rceil
    =D_r\lceil mean\_cluster\rceil. \tag{2}
\]

Thus an entry moment bound gives a window at most a factor
`ceil(mean_cluster)` larger than the raw-density window.  If
`mean_cluster` is polynomial in the input length, this is a lower-order
`O(log b)` term in the log-cost (15) of the existing criterion.  It does
not by itself improve the factoring exponent: the uniform density and
moment hypotheses are still unproved.

The finite checks below are independent of the two cycles reserved by the
owner and worker 5.  They support the structural premise, not a uniform
theorem.  No factoring score, Lean result or independent mathematical
acceptance is claimed.  This ordinary source audit used OpenAI assistant
support; the owner's Path C scope and the coordination record in
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51)
are retained.

## 1. The exact entry-gap equivalence

Work on the cyclic word `f` of length `h`.  An entry is a position whose
value is one and whose predecessor is zero.  Number the entries cyclically
and let `d_r` be the forward distance from entry `r` to the next entry.
Between those two entries there are exactly `d_r-1` positions, all
non-entries: the `L_r` positions of the run beginning at entry `r`
including the entry itself, followed by `G_r` zero positions, with
`L_r+G_r=d_r`.  Distinct runs are separated by at least one zero, so
`L_r>=1` and `G_r>=1`; hence `d_r>=2`.  Summing the distances around the
cycle gives `sum_r d_r=h`.

If `W<=d_max-1`, choose the gap with `d_r=d_max`.  The `d_max-1` positions
strictly between its two entries contain no entry, so a window of length
`W` inside that block is empty.  Conversely, an entry-free window of
length `W` lies between two consecutive entries; the distance between
them is at least `W+1`, so `d_max>=W+1`.  This proves the equivalence,
including windows that wrap around the period.

The raw indicator `f` has zero runs of lengths `G_r`.  Therefore

\[
 \max_r G_r \le \max_r(L_r+G_r) \le \max_r G_r+\max_r L_r. \tag{3}
\]

The left inequality says the entry-gap window is never smaller than the
raw zero-run window; the right inequality is the exact cluster-length
loss.  The loss is additive in the maximum cluster, not multiplicative,
when the exact gaps are used.  The moment route replaces this exact gap
by a larger sufficient certificate; the finite table below shows how
large that certificate can be.

## 2. The conditional cluster-loss criterion

Apply the lower-tail criterion of the [cofactor window
note](cofactor-orbit-moment-criterion.md) to the entry word `e`.  Its
density is `delta_e=runs/h`, its density denominator is `D_e=ceil(h/runs)`,
and its period is still `h`.  For positive integers `q,C` and a certified
period bound `Hbar>=h`, define

\[
 T=\left\lceil\bigl((2q+1)(Cq)^q\overline H\bigr)^{1/(q+1)}\right\rceil,
 \qquad W_e=D_e(T+1). \tag{4}
\]

If the entry lower-tail moment satisfies

\[
 \frac1h\sum_x\bigl(\delta_e W_e-E_x(W_e)\bigr)_+^{2q}
   \le (CqW_e\delta_e)^q, \tag{5}
\]

then the proof of (8)--(10) in the cited note shows that no length-`W_e`
window has zero entries.  Hence every length-`W_e` window has a hit.
Combining (2) with (4),

\[
 W_e \le \bigl\lceil mean\_cluster\bigr\rceil
       \,D_r(T+1), \tag{6}
\]

so the entry route pays at most the polynomial cluster factor relative to
the raw criterion.  This is the precise content of the proposed
"polynomial cluster-length loss plus entry moment" mechanism.  The
remaining obligations are uniform bounds on `mean_cluster` and on the
entry lower-tail moment (5), not a new endpoint identity.

The exact-gap lemma gives a sharper alternative.  If one can prove

\[
 \max_r(L_r+G_r)\le G_{entry}(b) \tag{7}
\]

uniformly, then `W=G_entry(b)` already works, with no moment hypothesis.
The finite data below is consistent with a polynomial right side, but it
does not prove (7).  The raw zero-gap bound `max G_r` remains the existing
obstruction; the entry formulation can only add the maximum cluster to it.

## 3. Finite independent controls

The [checker](../tools/check_entry_gap_cluster_loss.py) recomputes the
nine complete orbits already recorded in
[the cofactor moment evidence](../evidence/cofactor-orbit-moments.json),
excluding the owner-named cycle `N=1000036000099, Y=19, g=5` reserved by
worker 5.  The retained cases are `N=1022117` (two generators),
`N=100160063` (four generators), and `N=1000036000099` with `g=7` and two
reached residual ideals.  All nine periods are certified by return to the
identity; there are no incomplete periods.

The exact gap statistics are:

| `N` | bits | generator | period `h` | `delta` | runs | `D_e` | mean cluster | max cluster | raw max zero | max entry gap |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1022117 | 20 | prime 3 | 58 | 0.3966 | 3 | 20 | 7.667 | 13 | 19 | 24 |
| 1022117 | 20 | residual 49 | 58 | 0.4138 | 14 | 5 | 1.714 | 3 | 5 | 8 |
| 100160063 | 27 | prime 2 | 4230 | 0.7565 | 366 | 12 | 8.743 | 44 | 13 | 50 |
| 100160063 | 27 | prime 3 | 423 | 0.7400 | 54 | 8 | 5.796 | 21 | 4 | 23 |
| 100160063 | 27 | residual 59 | 8460 | 0.7537 | 1518 | 6 | 4.200 | 29 | 5 | 30 |
| 100160063 | 27 | residual 677 | 4230 | 0.7565 | 804 | 6 | 3.980 | 22 | 4 | 24 |
| 1000036000099 | 40 | prime 7 | 198850 | 0.2115 | 16949 | 12 | 2.481 | 16 | 64 | 71 |
| 1000036000099 | 40 | residual 9169 | 198850 | 0.2115 | 32823 | 7 | 1.281 | 7 | 57 | 58 |
| 1000036000099 | 40 | residual 167 | 198850 | 0.2115 | 34607 | 6 | 1.215 | 6 | 39 | 40 |

Every row passes the exact checks: entries are isolated, the entry gaps
sum to the period, every gap is at least two, the raw maximum-zero-run
matches, a window of length `d_max-1` is empty, and a window of length
`d_max` is nonempty.  The table is consistent with
`max entry gap <= 2b` and `mean cluster <= b` on these orbits, but those
inequalities are finite observations, not uniform bounds.  The largest
cluster is 44 at 27 bits; the largest exact entry gap is 71 at 40 bits.
The exact entry-gap window is often much smaller than the moment window
(4).  For example, at `N=1000036000099, g=7`, the exact entry gap is 71,
while (4) with `q=bit_length(b+1)=6`, `C=1` and `Hbar=2N` gives
`W_e=4620`.  The moment criterion is a loose certificate for a gap bound,
not a numerical optimization.

For the conditional moment check, use the owner's predeclared parameters
`q=bit_length(b+1)` and `q^2`, `C=1,2`, and `Hbar=2N`.  In the
non-wrapping cases at `q=bit_length(b+1)`, `C=1`, the entry lower-tail
bound (5) holds in all 6 cases; the raw bound holds in all 7 non-wrapping
raw cases.  The entry windows are larger by the cluster factor, and three
entry windows (both `N=1022117` generators and `N=100160063, g=3`) wrap
and are marked trivial.  The reserved `g=5` case is deliberately not used here; the
owner's existing note already records that its raw `C=1` bound fails at
its window, while worker 5 owns the exact entry comparison.

The synthetic controls exhaust all short cyclic words of period at most
12 (8,166 nonconstant words) and 44 clustered words of periods 17--64.
They verify the equivalence for every `W`, the lower-tail Lipschitz
controls, the cluster-loss inequality (2), and the absence of adjacent
entries.  These are exact finite checks, not a proof of a uniform
cluster bound.

## 4. Raw hit-prefix replay

The four `tools/window_lemma/data/scan_*.jsonl` files contain 70,465
distinct `(N,cutoff,p)` records on 328 integers and 1,103 primes.  The
checker verifies the exact SHA256 of each file, distinct keys, sorted and
unique hit offsets inside the scanned prefix, and the absence of
conflicting prefixes.  There are no validation errors; 13 records are
censored at the scan cap.  The median number of recorded hits is 2, the
median first-hit offset is 201, the largest first-hit offset is 3,394,411,
and the median starting exponent `K_p` is 5 (maximum 81).

The prefix is not a complete period and is stopped after a fixed number
of hits, so its density is a biased diagnostic.  With the prefix density
`len(hits)/scanned`, the empty-window fractions for non-wrapping windows
are:

| mean count `c` | `W=ceil(c/prefix density)` | records used | windows | empty fraction | records with any empty window |
|---:|---:|---:|---:|---:|---:|
| 1 | median 250 | 70,462 | 256,151,017 | 0.43613 | 69,639 |
| 2 | median 500 | 70,460 | 133,400,520 | 0.14343 | 1,305 |
| 4 | median 1,301 | 2,803 | 2,803 | 0.0 | 0 |
| 8 | -- | 0 | 0 | -- | 0 |

At the actual starting phase `K_p` (the beginning of the prefix), the
length-`W` window is empty in 35,072 of 70,462 cases for `c=1`
(0.49774), 334 of 70,460 for `c=2` (0.00474), and none for `c=4` or
`c=8`.  The `c=4` and `c=8` windows are mostly longer than the prefix
(67,659 and 70,462 skipped), so those rows are not informative.  The
terminal position also biases the prefix: the scan usually stops at a
hit, and the final zero gap can be long in the censored cases (maximum
1,861,137).  These numbers therefore show only that the raw hit prefixes
are consistent and that the fixed-hit stopping rule leaves substantial
empty-window mass at mean count 1; they are not a uniform density or gap
estimate.

## 5. What changed and what remains

The exact entry-gap lemma separates two questions that the moment
criterion had merged:

1. **Exact window.**  The cofactor closure needs a bound on the maximum
   cyclic gap between entries.  This is at least the raw maximum zero run
   and at most the raw maximum zero run plus the maximum cluster length.
   A uniform polynomial bound on either quantity gives a much smaller
   window than the moment certificate.
2. **Moment certificate.**  If the maximum gap is not directly bounded,
   the entry lower-tail moment (5) can certify a larger window.  The
   entry process is isolated, so its lower tail is more regular than the
   raw hit process, but the price is the factor `ceil(mean_cluster)` in
   (6).  This is a proof technique, not an asymptotic improvement over a
   successful raw max-gap bound.

The finite data is consistent with a polynomial cluster length and a
polynomial entry gap on the independent complete orbits, and with the
entry `C=1` moment bound at the predeclared parameters.  The next
theoretical obligation is a uniform bound on the maximum entry gap
(equivalently, on the cluster and zero-gap lengths) or a uniform proof of
(5).  The raw zero-gap bound remains the existing obstruction; the entry
route cannot remove it.  The complete factoring score remains `1/2` until
one of these uniform estimates, together with the density bound, is
proved.

```sh
python3 tools/check_entry_gap_cluster_loss.py \
  --output evidence/entry-gap-cluster-loss.json
```
