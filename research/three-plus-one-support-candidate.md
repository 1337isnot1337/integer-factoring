# A three-plus-one cut support diagnostic and its one-piece obstruction

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

## Why it is not a crossing criterion

The endpoint tensor itself is

```math
E(x,y,u,v)=[x=v][y=u].                                \tag{3}
```

Under either `3+1` or `1+3` cut, `E` has rank `q`, not `q^2`: the
`(x,y,u)|(v)` matrix has one nonzero entry in each `v` column, and the
`(x)|(y,u,v)` matrix has one nonzero entry in each `x` row.  Therefore
the **single trivial piece** `A_1=E` has both cut ranks `q` and rank-sum
cost `q`.

For every `q>1`, `q < q^{sqrt(2)}`.  A criterion that admitted the
`3+1`/`1+3` cuts would therefore certify the trivial one-piece
decomposition of the endpoint as a strict crossing, independently of any
cancellation or auxiliary construction.  The low diagnostic value is an
artifact of changing the two flattenings used by the factoring interface,
not evidence of a new seed.

The ordinary two cuts `(x,y)|(u,v)` and `(x,u)|(y,v)` are essential: both
have rank `q^2` for `E`, so the one-piece seed has cost `q^2` and is
correctly excluded.  Any proposed non-rectangular escape must either
preserve the two-cut interface and change the channel algebra inside a
piece, or prove a new lower-bound criterion that does not admit the
trivial endpoint.  The `3+1` profile above supplies neither.  In
particular, the earlier interpretation of this table as a support value
below one is withdrawn.

## What remains open

The table is still a reproducible diagnostic of the rank-four Smith
profile under an alternative leg grouping.  The following questions are
separate and are **not** advanced by it:

- which non-rectangular grouping, if any, survives inside a piece while
  retaining both standard flattenings;
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
one-piece obstruction are the new content.
