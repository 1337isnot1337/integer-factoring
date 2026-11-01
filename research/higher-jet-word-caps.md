# Higher adjugate jets pay for ordered Smith-channel words

The scalar second-jet bound extends to every jet below the physical
alphabet size. The crossed adjugate profile is bounded below by the
ordered-word costs in the nonleading source Smith channels, truncated at
`q*s2`. This holds over **every field**, including singular leading channel
factors. Graded cell/shared constraints inherit the bound after physical
restriction and restoration of its minima.

These necessary conditions exclude the twelve-letter point that survives
the [preceding linear and second-jet checks](inherited-second-jet-caps.md#5-a-successor-point-survives-the-declared-nonlinear-checks).
These word conditions alone do not exclude every twelve-letter crossing
or improve the factoring score **1/2**. The separate earlier
[box-budget entropy certificates](bare-seven-to-twelve-letters.md) already
give uniform bare coverage through twelve letters. The later
[historical thirteen-letter point](joint-bare-profile-feasibility.md) passes
both those original budgets and the word conditions, but the subsequent
[box-adjugate composition](box-adjugate-word-caps.md) excludes it.

## 1. The integral low-jet identity

Work over a DVR `O` with uniformizer `t`. A regular crossed matrix with
normalized Smith profile `s1=0,s2,...,sn`, `n=q^2`, has a unit pivot.
Its exact pivot rank-one lift gives an ordinary tensor operator `R` such
that the crossed remainder `E=T-R` has zero pivot row and column and
nonzero Smith orders exactly `s2,...,sn`. All entries of `E` have order
at least `k=s2`. When `k=0` the conclusion below is immediate.

Write `R(M)=A M^T B`, put `D=det(A)det(B)`, and define the regular
transpose channel `S(M)=adj(B)^T M^T adj(A)^T`. On the generic invertible
locus,

```
det R = epsilon D^q,   R^(-1)=D^(-1)S,
epsilon=(-1)^(q(q-1)/2).
```

Let `e_a(W)` be the sum of all principal `a`-minors of `W`, with `e_0=1`.
Multiplying the determinant and inverse expansions in an independent
variable `x` gives the homogeneous coefficient

```
[x^ell] adj(R+xE)
 = epsilon D^(q-1-ell) sum_(a+b=ell) e_a(SE)(-SE)^b S,
                                      0<=ell<=q-1.     (1)
```

Every exponent of `D` in this range is nonnegative. Thus (1) is a
polynomial identity over the integers, initially proved on the generic
characteristic-zero locus and valid after specialization over every field.
No inverse of a singular channel is used. Terms of higher `E` degree
contain at least `q` entries of `E` and are divisible by `t^(qk)`.

Group (1) by `b`. This expresses the low-jet part as
`sum_(b=0)^(q-1) f_b(E,D)(SE)^b S`, with each `f_b` a **regular scalar**.
The elementary-symmetric factors therefore do not create further channel
rank or negative valuations.

## 2. Ordered words bound the complete lower profile

An integral Smith expansion of the crossed remainder gives

```
E = sum_(j=2)^n t^(s_j) Q_j,
```

where each `Q_j` is a transpose channel. Every summand of `(SE)^b S`
is a composition of `2b+1` transpose channels, hence another transpose
channel, whose crossed matrix has rank at most one over `O`. Its cost is
`s_(i1)+...+s_(ib)` for an **ordered** word `(i1,...,ib)`. Distinct words
are retained with multiplicity, even when their costs coincide.

Form the multiset consisting of the empty-word cost zero and all word
costs of lengths `1,...,q-1` on the alphabet `s2,...,sn`. Add arbitrarily
many copies of `qk`, sort, and let `h_i(s)` be the first `n` costs.
Costs at or above `qk` can be discarded before adding those copies.
Then the crossed adjugate Smith orders satisfy

```
eta_i >= h_i(s),               1<=i<=n.               (2)
```

Indeed, modulo `t^a` for `a<=qk`, all channels with word cost at least
`a` and the remainder vanish. The remaining rank-one lifts require at
most one cyclic generator per word. The number of adjugate Smith orders
below `a` is at most that number, proving (2). No independence among the
word channels is assumed.

For `q=2`, (2) is exactly the credited
[second-jet bound](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/adjugate-seven-letter-pressure-floor.md#8-a-stronger-universal-second-jet-constraint)
`eta_i>=min(s_i,2*s2)`. For larger `q`, the only words below `2k` are
the empty word and single letters, so (2) dominates that earlier bound.
For example, at `q=3` and source `(0,1,3,3,3,3,3,3,3)`, it gives
`(0,1,2,3,3,3,3,3,3)`, which a literal control attains.

## 3. Graded restriction and local minima

Put `H_k(s)=sum_(i<=k) h_i(s)`. Substitute this prefix into the existing
[graded adjugate cell/shared proof](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/graded-adjugate-caps.md).
On a square `m`-letter restriction, use `n=m^2` and
`C=(n-1)d+m*(sum X+sum Y)`. The necessary inequalities are

```
k S_(n-1)(r) + H_k(s) <= k C - TopCell_k,         1<=k<=n,
(l-1)S_n(r) + S_(n-l)(r) + H_l(s)
  <= l C - max(TopX_l+l*Ymax,TopY_l+l*Xmax),      1<=l<=m.  (3)
```

For the raw restricted second profile `v` with minimum `b>=0`, define the
same word-cost function using `v2,...,vn` and the empty word zero. Shifting
to `v-b` changes each letter by `b`, each permitted word by at most
`(m-1)b`, and the cap by `mb`. Order statistics consequently change by at
most `mb`; the first order remains zero. Thus

```
H_k(v-b) >= H_k(v) - (k-1)m b.                        (4)
```

The first-profile coefficients in (3) sum to `k(n-1)` for its row size
`k`. Restoring the restricted minima leaves an extra term at least
`[k(n-1)-(k-1)m]b>=0`. Dropping it, and using individual Smith
interlacing and coordinatewise monotonicity of word costs, proves (3) for
the first `n` orders of the full tensor. Smallest paired-energy subsets
give the strongest right sides as before. Both cut orientations apply.
Rational heights follow by ramification.

These functions are nonlinear. The [compiler](../tools/higher_jet_word_caps.py)
checks the two proposed profiles separately and makes no substitution of
their average. A heap computes only the first `m^2` words below the cap;
the complete exponential word set is not enumerated.

## 4. Literal checks and the twelve-letter rejection

Run `make verify-higher-jet-words`. The
[checker](../tools/check_higher_jet_word_caps.py) expands determinant
coefficients independently by a subset dynamic program over truncated
polynomial rings. It checks 1,650 identities from (1) over F2 and F10007,
including singular left and right channel factors. It independently
interpolates 259 complete cofactors and uses ordinary truncated Toeplitz
ranks to recover their lower profiles; three controls attain (2), and a
singular-channel control also satisfies it.

Six exhaustive word enumerations compare against the heap. Another 348
rational local-shift prefixes check (4) and its favorable restoration
sign. Nine independently graded literal tensors give 36 complete profiles
and 744 inherited inequalities. The
[compressed packet](../evidence/higher-jet-word-controls.json.xz) retains
matrices, coefficients, Toeplitz ranks, profiles and source hashes.

The prior twelve-letter point still passes its declared linear and
second-jet relaxation, with exact pressure below `-1/100`. Nevertheless it
violates **258 oriented word caps**. The whole-alphabet determinant case
demands `96841739<=96597072`, a contradiction of **244667** in its integer
normalization. The point and its historical weaker-relaxation receipt are
retained; it is no longer a construction target under (3).

The written argument supplies the all-field claim; the finite controls
test the implementation. A [clean full replay](../evidence/higher-words-full-verification.json)
at exact `269efb8` passed defaults and fifteen named joint/cofactor targets:
356 Python command blocks, 66 unit tests, five crossing tests and 988
Markdown files. Later main integration and integer targets are outside that
source pin. No uniform twelve-letter obstruction, realizing tensor, new Lean
proof or independent mathematical review is claimed. Substantial OpenAI
assistant support was used.
