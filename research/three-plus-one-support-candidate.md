# A three-plus-one cut support diagnostic and its two-alternative obstruction

## Result and scope

The [asymmetric-scale floor](asymmetric-degenerate-scaling-floor.md)
identifies a non-rectangular `3+1` or `1+3` leg partition as a possible
escape from the rectangular zero-scale identity.  This note constructs the
smallest explicit q4 instance from the Gaussian endpoint tensor and applies
the existing asymmetric support expression to its four Smith channels.

Cut the four old legs as

```math
(x,y,u)\mid(v),\qquad (x)\mid(y,u,v).                 \tag{1}
```

Both cut matrices are `64 x 4` and have rank four.  With zero physical
letter heights, the infinity/zero Smith orders give the following
diagnostic values.  The ordinary q4 `2+2` cuts are replayed as a
calibration.

| intrinsic scale | calibration `2+2` (balanced) | diagnostic `3+1` (balanced) | diagnostic `3+1` (`a=0.2`) |
| ---: | ---: | ---: | ---: |
| `1e-4` | `1.4147906914` | `0.7078281542` | `0.5906193862` |
| `1e-3` | `1.4199894432` | `0.7143228067` | `0.5971204403` |
| `1e-2` | `1.4724312944` | `0.7794964984` | `0.6629313224` |

The calibration is above one at every tested scale.  The `3+1` expression
is below one and increases with the intrinsic scale.  Using the same
`2 log q` normalization as the four-channel calibration, its formal
zero-scale limit is

```math
R_{\text{rect}}=a+b,\qquad
R_{3+1}=\frac{a+b}{2}.                                \tag{2}
```

This is only a diagnostic of the Smith orders under the unusual cut.  It
is not a crossing criterion.

## The cut-interface lemma

Let all four alphabets have size `q` and let

```math
E(x,y,u,v)=[x=v][y=u].                                \tag{3}
```

For any bipartition of the four legs, the flattening rank of `E` is `q^2`
exactly when both matched pairs `(x,v)` and `(y,u)` are split across the
two sides; otherwise it is `q`.  Indeed, a split pair contributes one free
diagonal label to the image, while an unsplit pair contributes a free
label on each side but no additional cross-relation.

Consequently a **one-piece** endpoint seed has rank-sum cost `q` for a
pair in which neither cut splits both pairs, `q^{3/2}` when exactly one
cut does, and `q^2` when both do.  Since `q<q^{sqrt(2)}<q^{3/2}<q^2` for
`q>1`, only the first case is trivially certified by the endpoint.
This is a necessary interface condition, not a lower bound.

## The two-alternative obstruction

The pairing of the four legs matters.  Under the **two alternative cuts**
`(x,y,u)|(v)` and `(x)|(y,u,v)`, `E` has rank `q` on both sides: the
`(x,y,u)|(v)` matrix has one nonzero entry in each `v` column, and the
`(x)|(y,u,v)` matrix has one nonzero entry in each `x` row.  Therefore the
single trivial piece `A_1=E` has rank-sum cost `q`.

For every `q>1`, `q < q^{sqrt(2)}`.  A rank-sum criterion that admitted
this **pair** of cuts would certify the trivial one-piece endpoint as a
strict crossing, independently of any cancellation or auxiliary
construction.  The low diagnostic value is an artifact of changing both
flattenings used by the factoring interface, not evidence of a new seed.
The interpretation of the table as a support value below one is
withdrawn.

The ordinary cuts `(x,y)|(u,v)` and `(x,u)|(y,v)` both have rank `q^2`
for `E`, so the one-piece seed has cost `q^2` and is correctly excluded.

## Mixed-cut endpoint control

A mixed pair with one standard cut and one `1+3` cut is the remaining
support-only possibility.  At q4, with `i` the standard `(x,u)|(y,v)`
cut at infinity and `j` the `(x)|(y,u,v)` cut at zero, the Gaussian
diagnostic and the **one-piece endpoint** have the following values at
scale `1e-4`:

| `a` | Gaussian `i` + `1+3` | endpoint `i` + `1+3` | Gaussian `3+1` + `j` | endpoint `3+1` + `j` |
| ---: | ---: | ---: | ---: | ---: |
| `0.05` | `0.5926434308` | `0.5493746089` | `1.0993104311` | `1.0993104311` |
| `0.20` | `0.7030219537` | `0.6898979486` | `1.1528292710` | `1.1528292710` |
| `0.50` | `0.9383361895` | `0.9330127019` | `1.1885206326` | `1.1885206326` |
| `0.95` | `1.1089322899` | `1.1061249500` | `0.8595744515` | `0.8595744515` |

The Gaussian profile is **above** the endpoint profile in the first
orientation and **equal** to it in the second.  At q8 the same
comparison approaches equality as `t` decreases: at `a=0.05`, `t=1e-6`,
the first-orientation values are `0.549382364` (Gaussian) and
`0.549374609` (endpoint).  Thus the mixed support value is the endpoint
value up to a positive finite-scale correction, not a cancellation gap.
The endpoint is one piece with ranks `q` and `q^2`, whose rank-sum cost
`q^{3/2}` is above `q^{sqrt(2)}`; it is a valid negative control, not a
seed.  This closes the **support-only mixed-cut** lane.  A mixed leg cut
still needs a different transfer that charges the unequal widths; the
support value alone does not supply one.

## What remains open

The obstruction is specific to the **both-alternative** pair.  A mixed
pair with one standard cut and one `3+1`/`1+3` cut does **not** admit the
one-piece endpoint: the two endpoint ranks are `q` and `q^2`, so the
trivial cost is `q^{3/2} > q^{sqrt(2)}`.  Such a mixed pair is not tested
here, and its channel/letter invoice and zero-scale normalization are
open.

The following questions are also separate and are not advanced by the
table:

- whether any non-rectangular cut can be used with a new transfer that
  charges unequal widths; the support-only mixed lane is now closed;
- the correct channel/letter invoice for a grouped side;
- effective rank factors, degree, coefficient bit growth, and the
  finite-alphabet recursion.

The finite q4 value is not a complete-bound improvement.  No score
changes; the complete exponent remains `1/2`.

## Reproduction

```sh
make verify-three-plus-one-support
```

The target runs the checker and the documentation link check.  The
[checker](../tools/check_three_plus_one_support.py) builds the q4 tensor,
the ordinary and alternative cuts, exact F4 Smith profiles, and the
support values; the [receipt](../evidence/three-plus-one-support.json)
records the six scale/angle cases and the formal zero-scale identities.
The Gaussian tensor, matrix oracle and asymmetric support expression are
existing repository machinery; the alternative cut construction and the
two-alternative obstruction are the new content.
