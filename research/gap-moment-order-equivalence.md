# Moment order, maximum gaps, and the cofactor window

## Result and scope

The [cofactor moment criterion](cofactor-orbit-moment-criterion.md) uses a
lower-tail moment at a sufficient window (11). This note records two
elementary facts about that criterion and one finite diagnostic.

1. At the period-scale order `q+1 >= log_2(Hbar)`, the lower-tail criterion
   is equivalent, up to an absolute constant in `C`, to a maximum
   inadmissible-gap bound. The moment language removes a small constant,
   not the gap obligation.
2. At smaller orders the moment criterion is strictly stronger than the
   gap bound: it must control windows with few hits, not only empty ones.
   The previously checked order `q = bit_length(input_bits+1)` is in this
   smaller-order regime.
3. On the saved generous-predicate cycles the maximum gap is below
   `2 D log_2 h` on every complete orbit of the recorded panel. A clustered
   abstract orbit of the same density violates the bound by a factor 85.
   Raising `q` to the input-bit scale shrinks the sufficient window by a
   factor 4.6 and lowers the least lower-tail constant from 2 to 1 on the
   saved generator-five cycle.

The new statements are ordinary inequalities. The finite data are not a
uniform theorem. No factoring-score improvement is claimed; the complete
bound remains **1/2**.

The concrete target that emerges is a uniform gap law
`G = O(D log h)` for the original smooth set along target cycles, with
`D = ceil(1/delta)`. It is a target-orbit pseudorandomness statement
stronger than the [all-coset density theorem](small-index-smooth-density.md),
and it is not implied by density or coset density alone. It is weaker than
the ambient [positive-minorant](near-identity-positive-minorants.md)
condition that the [uniform obstruction](uniform-ambient-minorant-obstruction.md)
excludes.

## 1. Notation

Let `f` be a nonzero periodic zero-one flag on a cycle of length `h`, with

```math
s=\sum_x f(x),\qquad \delta=s/h,\qquad D=\lceil h/s\rceil .
```

Thus `1 <= delta D < 2`: the lower bound is `D >= 1/delta`, and the upper
bound follows from `D < h/s + 1`. Write `G` for the longest cyclic run of
zeros, and for a positive integer `W` put

```math
A_x(W)=\sum_{j=0}^{W-1}f(x+j),\qquad \mu=\delta W,
\qquad M_q(W)=\frac1h\sum_x\bigl(\mu-A_x(W)\bigr)_+^{2q}.
```

The criterion is

```math
M_q(W)\le (CqW\delta)^q=(Cq\mu)^q. \tag{1}
```

Fix `Hbar >= h` and

```math
T=\left\lceil\bigl((2q+1)(Cq)^q\overline H\bigr)^{1/(q+1)}\right\rceil,
\qquad W=D(T+1). \tag{2}
```

This is equation (11) of the moment note. Since `delta D >= 1`,
`mu >= T+1`.

## 2. The two implications

### 2.1 The criterion forces the gap bound

If `A_x(W)=0`, the sliding-window lower bound (8) of the moment note gives

```math
M_q(W)\ge\frac{\mu^{2q+1}}{(2q+1)h}.
```

Together with (1) this would force
`mu^{q+1} <= (2q+1)(Cq)^q h`. But (2) gives
`mu >= T+1 > ((2q+1)(Cq)^q Hbar)^{1/(q+1)}`, hence
`mu^{q+1} > (2q+1)(Cq)^q Hbar >= (2q+1)(Cq)^q h`, a contradiction.
Therefore:

```math
\boxed{\text{(1) at the window (2)}\quad\Longrightarrow\quad G<W.} \tag{3}
```

This is the contrapositive of the empty-window contradiction already used
to build (11). It is stated here because it is the exact gap consequence.

### 2.2 A gap bound gives the moment bound when `mu <= Cq`

If `G<W`, every window has at least one good phase, so `A_x(W)>=1` and

```math
\bigl(\mu-A_x(W)\bigr)_+ < \mu .
```

Hence `M_q(W) < mu^{2q}`. The universal lower-tail bound
`M_q(W)<=mu^{2q}` and its consequence `C_req<=mu/q` are already recorded
in the [tail-sign correction](entry-moment-tail-sign-correction.md); the
gap condition supplies the additional `A_x(W)>=1` needed to compare the
same window with the same `C`. If also `mu <= Cq`, then

```math
mu^{2q}\le (Cq\mu)^q,
```

and (1) follows **at the same `W` and with the same `C`**. Thus

```math
\boxed{G<W\ \text{and}\ \mu\le Cq\quad\Longrightarrow\quad
\text{(1) at }(C,q,W).} \tag{4}
```

The condition `mu <= Cq` is the only place where the order enters. It is
not a density assumption; it compares the mean window count with the
moment order.

## 3. The period-scale threshold

Use `mu < 2(T+1)` and bound the ceiling in (2). If
`q+1 >= log_2(Hbar)`, then `Hbar^{1/(q+1)} <= 2`. For `q>=1` also
`(2q+1)^{1/(q+1)} <= 2`. Therefore

```math
T \le 4Cq+1,\qquad \mu<2(T+1)\le 8Cq+4\le 12Cq.
```

So (4) applies with `C'=16C` (the constant is deliberately loose):

```math
\boxed{G<W(C)\quad\Longrightarrow\quad
\text{(1) at }(16C,q,W(C))\quad
(q+1\ge\log_2\overline H).} \tag{5}
```

Combining (3) and (5):

**Order-gap equivalence.** For `q+1 >= log_2(Hbar)`, the lower-tail
criterion (1) at the window (2) is equivalent, up to replacing `C` by
`16C`, to the maximum-gap bound `G < W(C)`. At `Hbar=2N` and an `n`-bit
input, the threshold is `q >= n`. The window at that order is

```math
W=D(T+1)=O(CDn).
```

At the smaller order `q=(log n)^alpha` used in the conditional instance
(18) of the moment note, `Hbar^{1/(q+1)}` is `exp(Theta(n/(log n)^alpha))`,
which is much larger than one. There the gap bound `G<W` does **not**
imply (1) with a comparable constant. The moment criterion is strictly
stronger: it also controls windows with one or a few hits. The crude
sufficient constant from the gap would be
`O(C Hbar^{1/(q+1)})`, not `O(C)`.

## 4. A conditional gap theorem and its cost

Suppose the following uniform arithmetic input.

> **Gap law.** There are an absolute `K` and a density bound `D_unif` such
> that every generated target `g` satisfies
> `G_g <= K D_g log_2 h_g`, where `D_g=ceil(1/delta_g)<=D_unif`.

Since `h_g<=2N` and `N` has `n` bits, `log_2 h_g<=n+1`. Scanning
`W_g=G_g+1<=K D_unif(n+1)+1` consecutive powers from the selector's start
therefore produces an admissible hit. The closure uses at most
`(m+K_B)W_max` norm queries, with `m<=Y+1` initial generators and
`K_B=B bitlength(B)` residual keys. The moment window and the
period-root factor are not used. Under the inherited generation,
preprocessing and extraction contracts, the cost is

```math
\log T_{\rm factor}
=O\bigl(\log n+\log Y+\log B+\log D_{\rm unif}+\log K\bigr). \tag{6}
```

If `log Y+log B+log D_unif=O(n/log n)`, this is
`exp(O(n/log n))`, i.e. saving exponent `alpha=1`. This is the same
ceiling as the current density/factor-base scale; the gap law alone does
not exceed it. It does remove the need for a lower-tail moment at every
phase. Neither the gap law nor (6) is proved here.

## 5. Density and coset density do not give the gap law

The bound (3) is not reversible from density. Let `h=2m`, let `S` be an
interval of `m` consecutive phases, and let `f=1_S`. Then

```math
\delta=1/2,\qquad D=2,\qquad G=m=h/2,
```

while `2D log_2 h = 4 log_2(2m)`. The ratio `G/(D log_2 h)` tends to
infinity. So a set of constant density can cluster in one interval.

The same example is compatible with the all-coset density conclusion of
the small-index theorem. A coset of a subgroup of index `q` is a
`q`-spaced arithmetic progression. An interval of length `delta h` meets
it in about `delta h/q` points, so its density in the coset is about
`delta` as long as the interval is long compared with `q`. Thus even
density in every coset of every small-index subgroup does not bound the
cyclic gaps. The gap law is a genuinely stronger interval statement about
the arithmetic smooth set, not a corollary of the density proof.

The [positive-minorant criterion](near-identity-positive-minorants.md) is
a different sufficient condition with a window `(2K_p-1)D`. Its uniform
ambient version is excluded by the
[arithmetic norm-gap obstruction](uniform-ambient-minorant-obstruction.md).
The gap law is weaker: it asks only for the longest bad run on each target
cycle, not for a nonnegative minorant of the whole smooth set. That is why
it remains a candidate after the ambient obstruction.

## 6. Finite panel

The tool reads the two existing complete-orbit receipts and computes
`D`, `G`, and `G/(D log_2 h)` from the recorded integers. No new arithmetic
walk is used for this table. All ten generous-predicate complete orbits
satisfy the integer bound `G < 2 D bit_length(h)`.

| Input | Generator | `h` | `s` | `D` | `G` | `G/(D log2 h)` |
|---|---:|---:|---:|---:|---:|---:|
| 1022117 | prime 3 | 58 | 23 | 3 | 19 | 1.081 |
| 1022117 | residual 49 | 58 | 24 | 3 | 5 | 0.285 |
| 100160063 | prime 2 | 4230 | 3200 | 2 | 13 | 0.540 |
| 100160063 | prime 3 | 423 | 313 | 2 | 4 | 0.229 |
| 100160063 | residual 59 | 8460 | 6376 | 2 | 5 | 0.192 |
| 100160063 | residual 677 | 4230 | 3200 | 2 | 4 | 0.166 |
| 1000036000099 | prime 5 | 198850 | 42053 | 5 | 100 | 1.136 |
| 1000036000099 | prime 7 | 198850 | 42053 | 5 | 64 | 0.727 |
| 1000036000099 | residual 167 | 198850 | 42053 | 5 | 39 | 0.443 |
| 1000036000099 | residual 9169 | 198850 | 42053 | 5 | 57 | 0.648 |

The strict generator-seven orbit gives a second density scale:
`h=198850`, `s=343`, `D=580`, `G=4506`, ratio `0.441`. The strict
generator-five orbit is excluded from this panel because its entry
analysis was reserved by another worker; the existing
[entry-moment note](entry-moment-cycles.md) records its strict raw
window at `q=36`.

The abstract clustered control (`h=4096`, `s=2048`, `D=2`, `G=2048`) has
ratio `85.33` and violates the integer bound. It is the density-compatible
counterexample to the gap law.

## 7. Order scan on the saved cycles

The generous predicate is the one in the saved larger input:
`N=1000036000099`, `Y=19`, residual cap `10000`, `h=198850`,
`s=42053`, `D=5`. The exact lower-tail check uses `Hbar=2N`.

| Generator | `G` | `q` | `C` | `W` | min hits | least `C` | holds |
|---|---:|---:|---:|---:|---:|---:|---|
| 5 | 100 | 6 | 1 | 1925 | 318 | 2 | no |
| 5 | 100 | 20 | 1 | 405 | 42 | 1 | yes |
| 5 | 100 | 40 | 1 | 415 | 43 | 1 | yes |
| 5 | 100 | 60 | 1 | 490 | 53 | 1 | yes |
| 7 | 64 | 6 | 1 | 1925 | 343 | 1 | yes |
| 7 | 64 | 40 | 1 | 415 | 45 | 1 | yes |

At `q=40` the window is 4.6 times shorter than at `q=6`, and the
generator-five lower-tail constant drops from 2 to 1. The mean count at
`q=40` is `mu=87.8`, which is larger than `q`; the observed `C=1` is
therefore stronger than the crude gap implication (5), which would allow
`C'=16`. The actual window distribution has far fewer one-hit windows than
the worst-case bound assumes.

The strict predicate (`Y=19`, residual cap `1`, `h=198850`, `s=343`,
`D=580`) behaves differently. For generator seven, `G=4506` and the raw
`q=40`, `C=1` window is `48140<198850`; the exact lower-tail bound holds
with least constant 1 and minimum count 30. The entry process has 73
entries, `D_entry=2724`, maximum entry gap 4508, and `q=40`, `C=1`
window `226092>198850`: it wraps the period. Raising `q` therefore does
not rescue the entry process on this cycle. The entry density loss
dominates the order gain, extending the `q=36` observation of the
[entry-moment note](entry-moment-cycles.md) to `q=40`.

The empty-window direction has a useful qualification. On the strict
generator-seven cycle, a window of length `G=4506` is empty, but the
lower-tail bound still holds at `q=40`, `C=1` because `mu=7.77<=q`. The
contradiction only appears at the larger prescribed window (2), where
`mu>=T+1`. The abstract clustered control, whose gap exceeds its own
prescribed window (`102<2048`), fails the bound with an empty window.
Thus (3) is specific to the prescribed window, not to every window longer
than the gap.

## 8. Remaining obligation and next mechanism

The missing arithmetic input is now explicit:

> Prove, or refute, a uniform bound `G_g=O(D_g log h_g)` for the
> original smooth set along the generated target cycles.

This is stronger than the all-coset density theorem and weaker than the
ambient positive minorant. It is a concrete target-orbit pseudorandomness
statement; the finite ratios above are consistent with a constant below
two, but ten saved orbits are not a uniform theorem. A bounded next check
is to compute the same ratio on already saved complete orbits with new
sizes or densities, looking for a ratio that grows with `h` or exceeds two.
That would falsify the law without a new parameter search. A proof would
need interval discrepancy for the even-product construction, not only
subgroup or coset density; the clustered control shows why the existing
Fourier density argument does not supply it.

## 9. Replay and scope

```sh
python3 tools/check_gap_moment_order.py \
  --output evidence/gap-moment-order.json
make verify-gap-moment-order
```

The tool recomputes only the saved generous generator-five and
generator-seven orbits and the strict generator-seven orbit from the
existing prime forms. It checks the integer gap bound, the abstract
clustered violation, the exact lower-tail scans at `q=6,20,40,60`, the
strict raw/entry comparison at `q=40`, and the two empty-window controls.
The panel table is read from the existing receipts.

Credits: the moment criterion, the window (11), the universal `mu/q`
lower-tail bound, the density theorem and the entry/event notes are prior
repository work. The equivalence
inequalities (3)-(5), the gap-law formulation and the order scan are the
new checkpoint. Substantial AI assistance was used. No Lean, independent
mathematical review or stronger complete bound is claimed.
