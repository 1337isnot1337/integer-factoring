# Exact entry lower-tail moments on two Path C cycles

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
entries implies a window theorem for hits. The candidate gain is that the
entry process removes the heavy lower tail created by clustered hits. The cost
is the factor `1/kappa = mean_hit_run / delta`.

The predeclared diagnostic used:

- the full-smooth-or-proper-gcd predicate with cofactor cap `1`;
- two complete base cycles, `N=450296711, Y=23, g=2` and
  `N=1000036000099, Y=19, g=5`;
- the exact discrete index of `g=7` in the second cycle;
- raw and entry lower-tail moments at `q=bit_length(b+1)` and `q^2`;
- constants `C=1,2` and the certified input-only period bound `2N`;
- empty-window fractions at matched mean counts `1,2,4,8`;
- the actual `K_p` starting phase.

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

## 3. Sufficient windows

Exact required constants `C_req`, where the bound holds at the displayed `W`
with `C_req^q` equal to the exact left/right ratio.

### Cycle A: `N=450296711`

| process | `q` | `C` | `W` | `C_req` | empty windows | start empty? |
|---|---|---|---|---|---|---|
| raw | 5 | 1 | 1432 | 0.6358 | 0/20456 | no |
| raw | 5 | 2 | 2536 | 0.5116 | 0/20456 | no |
| raw | 25 | 1 | 464 | 0.4479 | 0/20456 | no |
| raw | 25 | 2 | 896 | 0.2187 | 0/20456 | no |
| entry | 5 | 1 | 8234 | 0.6533 | 0/20456 | no |
| entry | 5 | 2 | 14582 | 0.2920 | 0/20456 | no |
| entry | 25 | 1 | 2668 | 0.1920 | 0/20456 | no |
| entry | 25 | 2 | 5152 | 0.2254 | 0/20456 | no |

All eight windows fit inside the period. The predeclared `C=1` and `C=2`
bounds hold, and the required constant is below `1` in every case.

### Cycle B: `N=1000036000099`

| process | `q` | `C` | `W` | result |
|---|---|---|---|---|
| raw | 6 | 1 | 223300 | wraps |
| raw | 6 | 2 | 403680 | wraps |
| raw | 36 | 1 | 46400 | `C_req=0.1319`, 0 empty, start nonempty |
| raw | 36 | 2 | 90480 | `C_req=0.1432`, 0 empty, start nonempty |
| entry | 6 | 1 | 1255100 | wraps |
| entry | 6 | 2 | 2268960 | wraps |
| entry | 36 | 1 | 260800 | wraps |
| entry | 36 | 2 | 508560 | wraps |

Here the entry density loss matters: `D=3260` for the entry process, so every
entry sufficient window at the predeclared `C,q` exceeds the period. Those
criteria are trivial, not positive evidence. The raw high-moment window at
`q=36, C=1` fits and is nonempty at the actual start.

## 4. Matched mean counts

`W = round(m/delta)` for mean `m`. The table gives the empty-window fraction and
the exact required constant at `q=bit_length(b+1)` and `q^2`.

### Cycle A

| process | `m` | `W` | empty fraction | `C_req(q=5)` | `C_req(q=25)` |
|---|---|---|---|---|---|
| raw | 1 | 7 | 0.7385 | 4.167 | 1.335 |
| raw | 2 | 15 | 0.5922 | 6.774 | 2.660 |
| raw | 4 | 29 | 0.4128 | 7.027 | 3.729 |
| raw | 8 | 58 | 0.2318 | 7.462 | 3.890 |
| entry | 1 | 45 | 0.3157 | 0.508 | 0.274 |
| entry | 2 | 90 | 0.1257 | 0.516 | 0.249 |
| entry | 4 | 181 | 0.0088 | 0.594 | 0.282 |
| entry | 8 | 361 | 0.0000 | 0.338 | 0.142 |

### Cycle B

| process | `m` | `W` | empty fraction | `C_req(q=6)` | `C_req(q=36)` |
|---|---|---|---|---|---|
| raw | 1 | 580 | 0.8207 | 16.666 | 6.039 |
| raw | 2 | 1159 | 0.7363 | 21.503 | 6.906 |
| raw | 4 | 2319 | 0.6559 | 19.954 | 6.446 |
| raw | 8 | 4638 | 0.5510 | 12.845 | 3.505 |
| entry | 1 | 3260 | 0.6134 | 1.137 | 0.377 |
| entry | 2 | 6520 | 0.4659 | 1.427 | 0.434 |
| entry | 4 | 13039 | 0.1708 | 0.591 | 0.158 |
| entry | 8 | 26079 | 0.0000 | 0.257 | 0.0767 |

The qualitative separation is clear. Raw windows of one or a few mean counts
have a heavy lower tail and fail `C=2` (the only exception is cycle A,
`m=1,q=25`, at `C_req=1.335`). Entry windows satisfy `C=2` everywhere in the
table and satisfy `C=1` from `m=4` onward on cycle B and everywhere on cycle A.
Increasing `q` lowers the entry constant by a factor of about two to four.
So the cluster-entry decomposition is doing real work on these cycles.

## 5. Actual starting phases

The actual selector starts at `K_p`, not at the identity. Both starts lie in
initial gaps:

- cycle A: `K_2=14`, `f(14)=0`, first hit at offset `7`. The raw `m=1`
  window `W=7` is empty; the entry `W=45` window contains two entries.
- cycle B: `K_5=9`, `f(9)=0`, first hit at offset `664`. The raw
  `W=580` and `W=1159` windows are empty; the entry `W=3260` window contains
  three entries.
- `g_7` (cycle B): index `66509`, start phase `67863`, `f=0`, first hit at
  offset `1265`. The raw `W=580` and `W=1159` windows are empty; the entry
  `W=3260` window contains one entry.

This is the concrete phenomenon the registered experiment was aimed at: small
raw windows can be empty at the actual phase, while the entry process at a
matched or slightly larger window is nonempty.

## 6. Verdict and remaining obligation

The finite data support the owner's mechanism in the following precise sense:

1. The entry lower-tail moment constant is `O(1)` at matched mean counts on
   both cycles, and it improves with `q`.
2. The cluster-length loss is small on these cycles (`mean_hit_run` about
   `5.6` and `6.2`), so the loss of density is a small constant factor here.
3. At every sufficient window that fits inside the period, the predeclared
   `C=1,2` bounds hold exactly.

They do **not** establish a uniform theorem. The missing uniform inputs are:

- a uniform bound `mean_hit_run <= poly(b)` on the reached cycles;
- a uniform entry lower-tail moment bound with `C=O(1)` at the required
  high order `q` and window `W`;
- the usual density and generation inputs of the arithmetic selector.

On the sparse 40-bit cycle, the entry density loss makes the `C=1,2` sufficient
window exceed the period. That is a small-scale obstruction to the finite
certificate, not evidence that the entry constants fail: the constants at
smaller windows are below the predeclared thresholds. A uniform proof must
control the period-root/density tradeoff for the actual parameter scales, not
only the matched-mean constants.

No complete factoring exponent changes. The experiment isolates the next
analytic target: an entry-process analogue of the cofactor orbit moment
criterion, with a polynomial cluster-length bound and a `q`-uniform constant.

## 7. Replay

```sh
python3 tools/check_entry_moment_cycles.py --output /tmp/entry-moment-cycles.json
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
