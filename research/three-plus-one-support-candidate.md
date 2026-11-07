# A three-plus-one cut support candidate

## Result and scope

The [asymmetric-scale floor](asymmetric-degenerate-scaling-floor.md)
identifies a non-rectangular `3+1` or `1+3` cut, with a rectangular Smith
matrix, as a genuine escape from the rectangular zero-scale identity.  This
note constructs the smallest explicit instance from the q4 Gaussian
endpoint tensor and evaluates the existing asymmetric support equation.

Cut the four old legs as

```math
(x,y,u)\mid(v),\qquad (x)\mid(y,u,v).                 \tag{1}
```

Both cut matrices are `64 x 4` and have rank four.  With zero physical
letter heights, the Smith orders at the infinity and zero ends give the
following support values.  The ordinary q4 `2+2` cuts are replayed as a
calibration.

| intrinsic scale | calibration `2+2` (balanced) | candidate `3+1` (balanced) | candidate `3+1` (`a=0.2`) |
| ---: | ---: | ---: | ---: |
| `1e-4` | `1.4147906914` | `0.7078281542` | `0.5906193862` |
| `1e-3` | `1.4199894432` | `0.7143228067` | `0.5971204403` |
| `1e-2` | `1.4724312944` | `0.7794964984` | `0.6629313224` |

The calibration is above one at every tested scale, as the rectangular
zero-scale identity requires.  The `3+1` candidate is below one throughout
this range and is increasing with the intrinsic scale; its infimum is the
zero-scale value.  For four `X` and four `Y` letters the zero-scale limits
are

```math
R_{\text{rect}}=a+b,\qquad
R_{3+1}=\frac{a+b}{2}.                                \tag{2}
```

Thus the candidate has a strict gap at the leading scale, not a small
finite-scale correction.

## What this does and does not establish

This is a **numerical support diagnostic**, not a transfer theorem.  The
existing support equation is applied to the four Smith channels of the
rank-four rectangular cut.  The ordinary rectangular formula has
`n=AB` channels and assumes a full-rank profile; the `3+1` cut is
rank-deficient in the `64 x 4` ambient matrix, and the grouped side has
external alphabet `64` rather than `4`.  The following obligations are
therefore open:

- prove that the two-cut transfer accepts a rank-four rectangular cut with
  unequal external alphabets;
- derive the correct channel/letter invoice for the grouped side and check
  that the endpoint scale is unchanged;
- charge the effective rank factors, degree, coefficient bit growth and the
  resulting finite-alphabet recursion;
- reproduce the candidate at q8, where the `3+1` matrix is `4096 x 16` and
  the current matrix oracle's `1024`-order limit does not apply directly.

The finite q4 value is a reason to pursue the construction, not a
complete-bound improvement.  No score changes; the complete exponent
remains `1/2`.

## Reproduction

```sh
make verify-three-plus-one-support
```

The target runs the checker and the documentation link check.  The
[checker](../tools/check_three_plus_one_support.py) builds the q4 tensor,
the ordinary and candidate cuts, exact F4 Smith profiles, and the support
values; the [receipt](../evidence/three-plus-one-support.json) records the
six scale/angle cases and the zero-scale identities.  The Gaussian tensor,
matrix oracle and asymmetric support equation are existing repository
machinery; the `3+1` cut construction is the new diagnostic.
