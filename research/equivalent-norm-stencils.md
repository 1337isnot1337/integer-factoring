# Exact equivalent-norm stencils and a charged finite comparison

## Result

The [new constructor](../tools/equivalent_norm_stencils.py) obtains several
valid norm candidates from each supplied prime-form power. It uses proper
integer changes of variables and passes the resulting unreduced forms to
the same oriented residual presentation machinery. Every extra norm,
including a duplicate or unusable norm, is counted.

The [complete finite comparison](../evidence/equivalent-norm-stencils.json)
finds no improvement at equal norm budgets. In all eight tested budgets
across the retained 64- and 80-bit inputs, the leading-only method has
at least the relation rank of the two- and four-norm alternatives. The
experiment does not change the active factoring policy or score exponent.

Proper quadratic-form equivalence and represented-norm sieving are
classical. The repository already has a larger
[coordinate norm sieve](arithmetic-norm-sieve.md) and a
[complete bounded ellipse sieve](complete-smooth-norm-sieve.md).
This checkpoint tests a small explicit stencil in the newer unfactored
residual-ideal presentation, with complete costs and exact replay.

## 1. Four norms in the same ideal class

Let `Q=(A,B,C)` be a reduced primitive positive form. The four proper
changes, written by their columns, are

$$
U_0=\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
U_1=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
U_2=\begin{pmatrix}1&0\\1&1\end{pmatrix},\quad
U_3=\begin{pmatrix}1&0\\-1&1\end{pmatrix}.
$$

Every determinant is one. The corresponding leading coefficients are

$$A,\qquad C,\qquad A+B+C,\qquad A-B+C.$$

For a general matrix with columns `(x,y)` and `(u,v)`, direct substitution
gives the transformed coefficients

$$
\bigl(Ax^2+Bxy+Cy^2,
2Axu+B(xv+yu)+2Cyv,
Au^2+Buv+Cv^2\bigr).
$$

The discriminant and primitivity are unchanged, and proper equivalence
preserves the ideal class. Thus a power `g_p^e=[Q]` supplies an exact
class relation from any of these leading norms without computing a new
class-group power. The extra norm arithmetic, smooth-part tests, known
prime divisions, atom refinement, and presentation work remain real costs.

The transformed norm need not satisfy the reduced leading-norm bound.
In particular `C=(B^2-Delta)/(4A)` can be much larger than `A`.
Consequently this experiment makes no diagonal-dominance claim using
the original height threshold. It uses exact rational selection of
original relations and the existing torsion-image extractor instead.

## 2. Decoding an unreduced primitive form

The decoder first checks the discriminant, primitivity, positivity, and
coprimality of the leading norm with the factoring input. Known prime
factors are stripped completely. Their orientations come from the middle
coefficient. An unfactored residual is represented by its odd unramified
norm and a canonical root. It is not subjected to primality testing or a
factoring oracle.

These orientation facts do not require reduction. At an odd unramified
prime, a root branch lifts uniquely through its whole prime-power
component. For `Delta=-4N` with `N=1 mod 4`, an even leading coefficient
in a primitive form has exactly one factor of two: writing `B=2b`,
primitivity makes `b` and `C` odd, and `AC=b^2+N=2 mod 4`.
The existing ramified-two class therefore remains the correct factor.
For odd discriminant, an even leading norm requires `Delta=1 mod 8`;
two is split and its orientation is determined by `B mod 4`.

A proper gcd between any tested norm and `N` is recorded directly. A
gcd equal to `N` is recorded as unusable. The latter occurs, for example,
when a proper transform exposes the discriminant-sized coefficient of
a principal form. Both cases remain counted in the norm budget.

## 3. Independent verification

The constructor is separate from the replay routine. Replay computes
each supplied target power again, evaluates the transformed quadratic
form at three independent coordinate vectors, reconstructs all known
prime words and root signs, and checks every resulting class identity.
It then independently verifies every atom's literal ideal word, all four
integer graph maps, and every compact unit substitution. Exact rational
selection checks the same retained original rows and ranks.

The [checker](../tools/check_equivalent_norm_stencils.py) includes 550
small reduced forms, 2,200 proper transforms, 8,800 independent coordinate
evaluations, and 1,882 complete class identities. Its 481 valid leading
cases reproduce the entire old decoder payload. Four invalid matrices
and eight corrupted result records are rejected.

The 24 fixture cases independently replay 33,600 supplied powers,
57,600 declared norm tests, and 57,514 class identities. Every case is
also checked after a JSON round trip. The receipt binds 25 source files
and the two original complete fixture archives by SHA-256.

## 4. Equal-budget observations

Each budget is tested with one, two, and four norms per supplied power.
The declared budgets are 256, 512, 1,024, and 2,048 for the 64-bit input,
and 1,024, 2,048, 4,096, and 8,192 for the 80-bit input. No budget was
selected by looking for a successful stencil result.

| Fixture | Norm budget | Stencil | Supplied powers | Exact-key rank | Refined rank | Factor found |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| 64 bits | 2,048 | Leading | 2,048 | 6 | 17 | 3,332,910,829 |
| 64 bits | 2,048 | Leading and trailing | 1,024 | 3 | 11 | No |
| 64 bits | 2,048 | Four norms | 512 | 2 | 6 | No |
| 80 bits | 8,192 | Leading | 8,192 | 1 | 3 | No |
| 80 bits | 8,192 | Leading and trailing | 4,096 | 2 | 2 | No |
| 80 bits | 8,192 | Four norms | 2,048 | 2 | 2 | No |

The larger stencils can improve rank at a fixed number of supplied
powers. That comparison spends more norm tests. It does not establish
an advantage under the equal norm budgets above, or an end-to-end
factoring speedup.

The measured analysis processes every declared norm and builds one final
presentation. The original power generation, archive loading, and
independent replay have separate costs. In particular, these timings
must not be compared with a fresh search that repeatedly checkpoints.
The replay timing is recorded separately for every case.

## 5. Reproduction and remaining scope

```sh
python3 tools/check_equivalent_norm_stencils.py \
  --output /tmp/equivalent-norm-stencils-replay.json
```

The returned data include the complete in-memory presentation and its
verification. The retained receipt stores the replay counts, exact
budgets, outcomes, hashes, and timings; it does not duplicate all 24
large intermediate presentations. The source archives and deterministic
checker reconstruct them.

No larger stencil or interval sieve is excluded by these finite results.
Such an extension would need its own charged candidate construction and
comparison. The tested four-norm policy has not earned a change to the
active factoring search.

This implementation and analysis used substantial OpenAI assistant support.
