# Exact entry lower-tail moments on two Path C cycles

**2026-11-02 correction:** the earlier checker selected the upper tail while
labeling it as the lower tail. Both sign sites and the tables below are
corrected from the identical archived flags, with zero new group compositions.
The former raw-versus-entry lower-tail contrast is withdrawn. See the
[exact correction and regression](entry-moment-tail-sign-correction.md).

- **Base commit:** `443bd08` (current `main` when this checkpoint was written).
- **Registered scope:** the bounded experiment at
  [Discussion #51, 2026-11-01T20:41:39Z](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338636).
- **Status:** exact finite diagnostic. It proves no uniform density, moment or
  window theorem, and it does not improve the complete factoring score.
- **Tool:** [`tools/check_entry_moment_cycles.py`](../tools/check_entry_moment_cycles.py).
- **Evidence:** [`evidence/entry-moment-cycles.json`](../evidence/entry-moment-cycles.json).

## 1. What was tested

The owner's registered mechanism separates smooth-hit *clusters* from cluster
*entries*. For the predicate `f` of
[the cofactor orbit criterion](cofactor-orbit-moment-criterion.md), put

```text
e(k) = f(k) (1 - f(k-1)).
```

Every entry is a smooth hit (or a proper-gcd hit), but an entry is counted only
once per hit run. If `delta` is the density of `f` and `kappa` that of `e`, then
exactly

```text
kappa = delta / mean_hit_run.
```

A window containing an entry contains a hit. Thus the same lower-tail moment
criterion can be applied to `e`, with density `kappa`, and a window theorem for
entries implies a window theorem for hits. The diagnostic asks whether removing clusters reduces lower deviations; the
old upper-tail calculation did not test that assertion. The cost is the factor `1/kappa = mean_hit_run / delta`.

The predeclared diagnostic used:

- the full-smooth-or-proper-gcd predicate with cofactor cap `1`;
- two complete base cycles, `N=450296711, Y=23, g=2` and
  `N=1000036000099, Y=19, g=5`;
- the exact discrete index of `g=7` in the second cycle;
- raw and entry lower-tail moments at `q=bit_length(b+1)` and `q^2`;
- constants `C=1,2` and the certified input-only period bound `2N`;
- empty-window fractions at matched mean counts `1,2,4,8`;
- the actual `K_p` starting phase.

The current replay enforces an aggregate cap of 250000 base-walk
compositions. The two periods use exactly 219306 steps; the generator-seven
index is recorded during the existing generator-five walk and costs no
extra walk. The earlier separate index helper added 66509 steps, making
its aggregate 285815; its `composition_cap` field did not certify that
aggregate. That bookkeeping repair preserved the base flags and computations. The later
sign correction changes only moment statistics, not either run profile.

For a window `W`, with cyclic window counts `A_x(W)` and `mu=delta W`, the
tested inequality is

```text
(1/h) sum_x (mu - A_x(W))_+^{2q} <= (C q W delta)^q.
```

The sufficient window is the integer from
[equation (11)](cofactor-orbit-moment-criterion.md):

```text
W = D ( ceil( ((2q+1)(Cq)^q * 2N)^(1/(q+1)) ) + 1 ),   D = ceil(1/delta).
```

When `W >= h` the window wraps a full period and the criterion is trivial; the
tool marks that case separately.

## 2. Exact cycle controls

Both cycles closed inside the `250000` composition cap. The period bound `h <= 2N`
was checked for each.

| `N` | bits | `Y` | `g` | `H` | `K_p` | period `h` | good | `delta` | entries | `kappa` | mean run |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 450296711 | 29 | 23 | 2 | 12251 | 14 | 20456 | 2821 | 0.1379057489 | 453 | 0.0221450919 | 6.2274 |
| 1000036000099 | 40 | 19 | 5 | 577360 | 9 | 198850 | 343 | 0.0017249183 | 61 | 0.0003067639 | 5.6230 |

Controls checked by the tool:

- the identity is reached, so each period is complete;
- `h <= 2N`;
- `f(0)=1`, and the inverse symmetry `f(k)=f(-k)` holds on both cycles;
- every entry is a hit, and the entry count equals the number of hit runs;
- `kappa = delta / mean_hit_run`;
- the gcd branch contributed zero hits in both cycles;
- `g_7 = g_5^66509` in the second cycle. Since `K_7=7`, the actual `g_7` start
  phase in the `g_5` cycle is `7*66509 mod 198850 = 67863`.

## 3. Corrected sufficient-window lower tails

The decimals below use the true lower tail. Exact integer inequalities are
in the corrected receipt; all previously specified windows are retained.

### Cycle A: `N=450296711`

| process | `q` | `C` | `W` | `C_req` | empty windows |
|---|---|---|---|---|---|
| raw | 5 | 1 | 1432 | 0.69338 | 0/20456 |
| raw | 5 | 2 | 2536 | 0.357367 | 0/20456 |
| raw | 25 | 1 | 464 | 0.473949 | 0/20456 |
| raw | 25 | 2 | 896 | 0.24697 | 0/20456 |
| entry | 5 | 1 | 8234 | 0.669454 | 0/20456 |
| entry | 5 | 2 | 14582 | 0.295602 | 0/20456 |
| entry | 25 | 1 | 2668 | 0.156032 | 0/20456 |
| entry | 25 | 2 | 5152 | 0.220861 | 0/20456 |

### Cycle B: `N=1000036000099`

| process | `q` | `C` | `W` | `C_req` | empty windows |
|---|---|---|---|---|---|
| raw | 6 | 1 | 223300 | wraps | not used |
| raw | 6 | 2 | 403680 | wraps | not used |
| raw | 36 | 1 | 46400 | 0.0635417 | 0/198850 |
| raw | 36 | 2 | 90480 | 0.0138492 | 0/198850 |
| entry | 6 | 1 | 1255100 | wraps | not used |
| entry | 6 | 2 | 2268960 | wraps | not used |
| entry | 36 | 1 | 260800 | wraps | not used |
| entry | 36 | 2 | 508560 | wraps | not used |

All ten nonwrapping sufficient-window inequalities still pass. The six
wrapped sufficient windows retain their previous coverage status. On the
sparse 40-bit cycle, every entry sufficient window still wraps the period;
that is not a nontrivial entry-window certificate.

## 4. Corrected matched-mean lower tails

Keep `W=ceil(m/delta)` and the previous empty-phase fractions. The changed
numbers are the lower-tail constants and their pass/fail flags.

### Cycle A

| process | `m` | `W` | empty fraction | `C_req(q=5)` | `C_req(q=25)` |
|---|---|---|---|---|---|
| raw | 1 | 8 | 0.7193 | 0.206578 | 0.0435521 |
| raw | 2 | 15 | 0.5922 | 0.372581 | 0.0810277 |
| raw | 4 | 30 | 0.4022 | 0.691646 | 0.159567 |
| raw | 8 | 59 | 0.2274 | 1.23693 | 0.306745 |
| entry | 1 | 46 | 0.3093 | 0.161122 | 0.0388788 |
| entry | 2 | 91 | 0.1229 | 0.26513 | 0.0741243 |
| entry | 4 | 181 | 0.0088 | 0.334459 | 0.132677 |
| entry | 8 | 362 | 0.0000 | 0.218159 | 0.0709386 |

### Cycle B

| process | `m` | `W` | empty fraction | `C_req(q=6)` | `C_req(q=36)` |
|---|---|---|---|---|---|
| raw | 1 | 580 | 0.8207 | 0.161339 | 0.0276382 |
| raw | 2 | 1160 | 0.7362 | 0.316889 | 0.0551098 |
| raw | 4 | 2319 | 0.6559 | 0.621634 | 0.10982 |
| raw | 8 | 4638 | 0.5510 | 1.2102 | 0.218578 |
| entry | 1 | 3260 | 0.6134 | 0.153639 | 0.0274047 |
| entry | 2 | 6520 | 0.4659 | 0.293507 | 0.054392 |
| entry | 4 | 13040 | 0.1708 | 0.497941 | 0.105794 |
| entry | 8 | 26079 | 0.0000 | 0.13043 | 0.0294155 |

All 32 matched-mean lower-tail tests pass `C=2`. Thirty pass `C=1`;
the two exceptions are the raw `m=8` tests at the smaller order on each
cycle. The previous claim that raw tests generally fail `C=2` was about
the upper tail and is withdrawn.

A universal bound explains much of the small constants:
`E(mu-A)_+^(2q)<=mu^(2q)`, so `C_req<=mu/q` in this normalization.
Low-mean, high-order passes therefore do not establish independence or
short empty-window control. Entry constants are sometimes smaller, but
the old claimed orders-of-magnitude lower-tail contrast is unsupported.
The empty fractions and their genuine difference between raw and entry
flags are unchanged.

## 5. Actual starting phases

The actual selector starts at `K_p`, not at the identity. Both starts lie in
initial gaps:

- cycle A: `K_2=14`, `f(14)=0`, first hit at offset `7`. The exact one-mean
  window `W=floor(1/delta)=7` is empty; the ceilings `W=8` (raw) and `W=46`
  (entry) contain one hit and two entries respectively.
- cycle B: `K_5=9`, `f(9)=0`, first hit at offset `664`. The raw `m=1`
  window `W=580` is empty; the raw `m=2` window `W=1160` and the entry
  `W=3260` window contain hits (nine and three respectively).
- The **actual `g_7` power walk** has flags
  `f_7(j)=f_5(66509*j mod 198850)` and starts at exponent `K_7=7`.
  Its first hit is at offset **3823**. It has 343 raw hits and 73 entries;
  maximum raw and entry zero gaps are 4506 and 4508. Its raw windows
  `W=580,1160` and its one-mean entry window `W=2724` are empty at this
  start. The two-mean entry window `W=5448` contains two entries.
- A separate **shifted generator-five walk** begins at base phase `67863`
  but advances by one base step, not by `66509`. Its first hit offset is
  `1265` and its entry window `W=3260` contains one entry. Those old
  measurements are retained under explicit shifted-base labels. The
  actual generator-seven entry window of length `3260` is empty.

Thus both raw and entry windows can miss at an actual starting phase.
The derived control records its own stride and density; the moment tables
above continue to concern the two original base directions. The broader
[saved-cycle replay](clustered-hits-and-pooled-rank.md) independently checks
the actual generator-seven direction and its moments.

## 6. Verdict and remaining obligation

The corrected finite findings are:

1. Both raw and entry matched-mean lower tails satisfy `C=2`. Much of the
   decrease with `q` follows from the universal `mu/q` bound; it is not an
   independent-hit certificate.
2. The cluster-length loss is small on these cycles (`mean_hit_run` about
   `5.6` and `6.2`), so the loss of density is a small constant factor here.
3. At every sufficient window that fits inside the period, the predeclared
   `C=1,2` bounds hold exactly.

They do **not** establish a uniform theorem. The missing uniform inputs are:

- a uniform bound `mean_hit_run <= poly(b)` on the reached cycles;
- a uniform entry lower-tail moment bound with `C=O(1)` at the required
  high order `q` and window `W`;
- the usual density and generation inputs of the arithmetic selector.

The [relation-change event follow-up](relation-change-event-density.md)
replaces pure entries by changes of the integer row. It proves a polynomial
density-loss bound for that changed process, without assuming a raw-cluster
length bound. Its own uniform event moments and original subgroup density
remain open. The finite flags of the two processes agree on the saved
generator-five and generator-seven directions; this agreement is not a
general identity.

On the sparse 40-bit cycle, the entry density loss makes the `C=1,2` sufficient
window exceed the period. That is a small-scale obstruction to the finite
certificate, not evidence that the entry constants fail: the constants at
smaller windows are below the predeclared thresholds. A uniform proof must
control the period-root/density tradeoff for the actual parameter scales, not
only the matched-mean constants.

No complete factoring exponent changes. The experiment isolates the next
analytic target: a uniform moment estimate for the chosen event process,
together with its stated original density input.

## 7. Replay

```sh
python3 tools/check_entry_moment_cycles.py --replay evidence/entry-moment-cycles.json --saved-flags
make verify-entry-moment-cycles
make verify-docs
```

The checker uses only integer arithmetic for the cycle predicates, periods,
entries, window counts, empty-window counts and exact moment inequalities. The
displayed `C_req` decimals are diagnostic; the exact inequalities are the
certificates.

Credit: the registered experiment and the cluster-entry mechanism are the
owner's. This checkpoint is an independent exact implementation and replay.
The arithmetic orbit code is the repository's existing prime-form machinery.
Substantial AI assistance was used; no human review, Lean formalization or
stronger score is claimed.

The [clean full replay](../evidence/entry-control-correction-full-verification.json)
passes the default and focused entry/event gates at exact source
`059a215cab000d216bcaff4838b0c4e6a6a913f4` in 748.87 seconds, with
clean tracked state before and after. Its complete compressed log is
archived. The receipt also records the instrumented 219306-call replay
and the forced-cap regression; later main integration checks are scoped
separately.

The historical full replay above predates the sign repair. It checked the
then-current code and receipt, whose tail labels were incorrect; it is not
evidence for the corrected statistics. The correction note records the
new saved-flag replay and asymmetric regressions separately.
