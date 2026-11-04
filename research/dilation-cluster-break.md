# Exponent dilation and the strict cluster-entry wrap

## Result and scope

The [gap-moment note](gap-moment-order-equivalence.md) reduced the
period-scale moment criterion to a maximum-gap bound on the original
target cycle. The [entry-moment cycle note](entry-moment-cycles.md)
recorded one finite exception: on the strict generator-seven cycle the
raw window fits the period, but the cluster-entry window wraps.

This note records a finite structural repair and a new exact identity.

1. A coprime exponent dilation `k -> a k` changes the cyclic order while
   preserving the good-class set. On the saved strict generator-seven
   cycle, step `a=1` gives 73 entry clusters and an entry window
   `226092 > 198850`; step `a=3` gives 185 clusters and window `89225`;
   step `a=13` gives 343 isolated good classes and window `48140`. The
   entry lower-tail moment holds with `C=1` at each tabulated step.
2. For **any** zero-one cycle of length `h` with `s` ones, the average
   number of entry clusters over all `h` cyclic dilations is exactly
   `s(1-s/h)`. The average over the unit group is the same expression
   minus the excess of unit-difference pairs over the uniform baseline.
   On the strict cycle the exact unit average is `342.409` out of `343`
   good classes.
3. The finite data therefore identify a new target: unit-difference
   equidistribution of the smooth set. If the number of ordered pairs
   with unit difference is close to `(phi(h)/h)s^2`, a random unit
   dilation has near-maximal cluster count. This is a weaker and
   different arithmetic statement than the maximum-gap law.

The exact identity is an ordinary finite-group calculation. The dilation
repair is a finite diagnostic. Neither proves a uniform cluster-breaking
lemma, an entry-gap theorem, or a complete bound better than **1/2**.

## 1. Dilation of the exponent order

Let `f` be a zero-one flag on `Z/hZ` with support `S` and `s=|S|`. For a
unit `a` modulo `h`, define the dilated flag

```math
f_a(k)=f(ak),\qquad S_a=a^{-1}S.
```

The relation row at exponent `ak` is the same row as before, with the
power `g^{ak}` in place of `g^k`:

```math
(ak)e_p-v=\sum_j u_j c_j .
```

If `a` is fixed, the scan generates `g^{ak}` by multiplying by `g^a` at
each step, so the number of power compositions and smoothness tests is
unchanged. Choosing one `a` costs one extra initial power. Trying a set
`A` costs `|A|` scans. The density `s/h` and the period `h` are unchanged;
the gap and the number of clusters generally change.

For the saved period `h=198850`, the smallest odd prime not dividing `h`
is `3`. This is the deterministic step used in the finite check. For a
general period, use the smallest odd prime not dividing it.

## 2. Clusters and the entry window

The entry flag is

```math
e(k)=f(k)\bigl(1-f(k-1)\bigr).
```

Its support has one point at the start of every cyclic run of ones in
`f`. Write `R=sum_k e(k)` and

```math
D_{\rm entry}=\left\lceil\frac{h}{R}\right\rceil .
```

The entry moment at order `q`, constant `C`, and period bound `Hbar` uses
the window (2) of the gap-moment note with `D_entry` in place of `D`:

```math
W_{\rm entry}=D_{\rm entry}(T+1),
\qquad
T=\left\lceil\bigl((2q+1)(Cq)^q\overline H\bigr)^{1/(q+1)}\right\rceil .
```

The entry process fits the period exactly when `W_entry < h`. This is the
condition that failed on the strict cycle at step one.

On the strict cycle (`N=1000036000099`, cutoff `19`, cap `1`, generator
`7`, `h=198850`, `s=343`, `D=580`), the exact data are:

| step `a` | raw gap | runs `R` | `D_entry` | entry gap | `W_entry` | fits |
| ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| 1 | 4506 | 73 | 2724 | 4508 | 226092 | no |
| 3 | 3949 | 185 | 1075 | 3949 | 89225 | yes |
| 7 | 8013 | 301 | 661 | 8013 | 54863 | yes |
| 13 | 5535 | 343 | 580 | 5535 | 48140 | yes |

At `a=13` every good class is isolated, so the entry flag equals the raw
flag and `D_entry=D`. The raw gap increases from `4506` to `5535`, but the
entry window still fits because the entry density is maximal. At `a=3`
the raw gap decreases and the cluster count more than doubles; the entry
window is `89225`, less than half the period. The entry lower-tail moment
at `q=40`, `C=1` passes in all four rows; the failure at `a=1` is the
scan-window fit, not the moment inequality.

The same check on the two dense large-period cycles
(`N=1000036000099`, cutoff `19`, cap `10000`, generators `5` and `7`)
gives:

| generator | step | raw gap | runs | `D_entry` | `W_entry` | fits |
| ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| 5 | 1 | 100 | 14329 | 14 | 1162 | yes |
| 5 | 3 | 40 | 28113 | 8 | 664 | yes |
| 5 | 7 | 37 | 33063 | 7 | 581 | yes |
| 7 | 1 | 64 | 16949 | 12 | 996 | yes |
| 7 | 3 | 50 | 29921 | 7 | 581 | yes |
| 7 | 7 | 37 | 33107 | 7 | 581 | yes |

Thus the strict cycle is the only saved large-period case in which the
entry window wraps at step one; the smallest-prime dilation repairs it,
and the dense cycles were already inside the period.

## 3. The exact dilation average

For a shift `a` modulo `h`, let

```math
R_a=\sum_k f(k)\bigl(1-f(k-a)\bigr)
    =s-\bigl|S\cap(S+a)\bigr| .
```

Summing over all shifts and counting ordered pairs,

```math
\sum_{a\in\mathbb Z/h\mathbb Z}R_a
=hs-\sum_a\bigl|S\cap(S+a)\bigr|
=hs-s^2
=hs\Bigl(1-\frac{s}{h}\Bigr).
```

Therefore the average over all `h` shifts is exactly

```math
\boxed{\ \frac1h\sum_a R_a=s\Bigl(1-\frac{s}{h}\Bigr).\ } \tag{1}
```

For a sparse set this is close to `s`; the average number of clusters is
large even when one particular shift is badly clustered. The interval
control of the gap-moment note has `s=h/2` and all-shift average `h/4`,
but its identity step has one run. The average is not the obstruction;
the unit restriction and the entry-gap are.

For the unit group `U=(Z/hZ)^*`, the same calculation gives

```math
\frac1{\varphi(h)}\sum_{a\in U}R_a
=s-\frac{N_U}{\varphi(h)},
\qquad
N_U=\#\{(x,y)\in S^2:x-y\in U\}. \tag{2}
```

The uniform baseline for `N_U` is `(varphi(h)/h)s^2`. If the difference
set is equidistributed over units and non-units, the unit average is
close to (1). On the strict cycle the exact values are

```math
h=198850,\quad s=343,\quad \varphi(h)=76800,
\quad N_U=45380,
```

so

```math
\frac1{\varphi(h)}\sum_{a\in U}R_a
=343-\frac{45380}{76800}
=342.409114583\ldots
```

and the all-shift average is `342.408353...`. The unit average is within
`0.0008` of the all-shift average, and the unit-difference pair fraction
is `45380/117649=0.38572`, against the baseline
`76800/198850=0.38622`. The strict set is not concentrated on non-unit
differences.

For the dense generator-five cycle the exact all-shift average is
`33159.589`; a deterministic sample of 32 units gives `33129.813` with
minimum `32635` and maximum `34321`. For generator seven the sample
average is `33112.938`. The sampled unit averages are within `0.15%` of
the all-shift baseline. These are samples, not exact unit averages.

## 4. What a uniform cluster-breaking lemma would need

The finite data suggest the following target.

> **Unit-difference equidistribution.** For every generated target cycle
> of length `h` and support size `s`, there is a small error `epsilon`
> such that
> `|N_U-(varphi(h)/h)s^2| <= epsilon s^2`,
> with `epsilon` small compared with the density margin.

If `epsilon` is small, (2) gives `R_avg >= s(1-s/h)-epsilon s`. A positive
fraction of unit dilations then has `R=Omega(s)`, so `D_entry=O(D)`. To
make the entry window fit, one also needs the entry gap below
`W_entry`. The entry gap is at most the raw gap plus the longest run of
ones, so a sufficient target is:

> **Dilated gap plus cluster bound.** There is a fixed small set `A` of
> units such that for every target cycle some `a in A` satisfies
> `R_a=Omega(s)`, `G_a=O(D log h)`, and the longest run of ones is
> `O(D log h)`.

The first condition is new and combinatorial-arithmetic. The second and
third are gap/cluster bounds for the dilated raw flag. Together they give
`D_entry=O(D)` and an entry gap `O(D log h)`, hence an entry window
`O(D log h)` that fits inside the period up to a constant factor. The
finite strict data satisfy the direct entry-gap bound without separating
these two terms.

The all-shift identity (1) is general. The unit restriction is not. A set
invariant under every unit dilation has `R_a=R_1` for every unit `a`, so
the unit average is just the generator-order cluster count and gives no
improvement. The interval control is not unit-invariant, but its identity
step has one run; dilation breaks it. The arithmetic input is therefore a
statement that the smooth set is far from both a union of multiplicative
cosets and an interval. The
all-coset density theorem controls cosets of small-index subgroups;
the non-unit indicator has Fourier support on the small-index subgroups
(and on the large prime factors of `h`, whose cosets have few residues). This suggests a route from the existing density
character estimates to unit-difference equidistribution. It is not
carried out here.

## 5. Scope and next check

The exact identity (1) and the finite tables are proved or computed
exactly by `tools/check_dilation_cluster_break.py`; the unit average for
the strict cycle is exact, and the dense-cycle unit averages are samples.
The finite repair does not prove that a fixed small dilation works for
every target, and it does not bound the number of trials needed to find
one. It also does not improve the complete exponent: the density
index barrier and the uniform dilated-gap obligation remain.

The next check is a bounded arithmetic test of unit-difference
equidistribution on the saved cycles: compare `N_U` with
`(phi(h)/h)s^2` for the exact sparse cycles, and test whether the
small-index Fourier terms of the non-unit indicator are controlled by
the all-coset density receipt. If that test passes, the next note should
state a precise cluster-breaking lemma and its cost; if it fails, the
dilation repair is only a finite accident and should not be pursued as a
uniform route.

The finite data are a continuation of worker 2's hit-gap lane and use
the saved complete orbits of the gap-moment and entry-moment receipts.
No new arithmetic walk or factoring score is claimed.
