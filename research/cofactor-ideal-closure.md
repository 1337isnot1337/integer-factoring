# Composite ideal closure without residual factorization

## Constructive result

The [cofactor collision selector](../tools/cofactor_collision_relations.py)
waits for repeated residual keys. A different construction promotes
each residual ideal to a generator and asks for one admissible power
of that generator. Every selected row is strictly diagonally dominant
under fixed norm weights on the **entire augmented generator set**.
The resulting rows are independent without a rational-rank search.

No residual norm is factored or tested for primality. A residual norm
may contain repeated prime factors. Its middle-coefficient root already
specifies the required invertible ideal.

At residual cap `B`, at most `B bitlength(B)` residual keys can be
introduced. If each reached generator has an admissible power in a
window of length `W`, the closure requires at most one successful
selection per generator. Its cost is polynomial in the number of
introduced generators, `W`, and the input bit length. A uniform
improved bound on those parameters remains open.

The implementation is
[cofactor_ideal_closure.py](../tools/cofactor_ideal_closure.py). Its
finite controls include a 40-bit factorization in 301 norm queries
with 70 new residual generators. These are finite results, with no
new factoring-score claim.

## 1. Fixed prime base and composite ideal generators

Let `N` be odd and use the existing discriminant choice

\[
 \Delta=\begin{cases}-N,&N\equiv3\pmod4,\\-4N,&N\equiv1\pmod4.
 \end{cases}
 \qquad
 H=\left\lfloor\sqrt{|\Delta|/3}\right\rfloor.
 \tag{1}
\]

Start with the known prime-form base, including the ramified-two
generator when appropriate. A base prime `p` has recorded norm `d=p`.
The norm is a label attached to the original ideal; reduction of its
form need not retain that leading coefficient.

A residual key has the form

\[
 \kappa=(R,C),\qquad
 1<R\le B,\quad \gcd(R,2\Delta)=1,\quad
 C^2\equiv\Delta\pmod{4R},
 \tag{2}
\]

where `C` is the canonical root up to sign modulo `2R`. Introduce the
class of

\[
 (R,C,(C^2-\Delta)/(4R))
 \tag{3}
\]

as a generator `h_kappa`, and record its norm label `d=R`. Its reduced
form is computed by the existing exact reduction algorithm.

The ideal specified by (3) is primitive and invertible because of the
coprimality in (2). This statement does not require `R` to be prime or
squarefree. Conjugating the root replaces the ideal class by its
inverse. Distinct keys can represent the same class; the construction
retains them as separate generators with their own fixed norm labels.
It never changes an existing generator's weight.

## 2. One admissible query gives one row

For a generator `g_i` with recorded norm `d_i>1`, start at

\[
 K_i=\min\{k:d_i^k>H\}.
 \tag{4}
\]

Inspect reduced powers `g_i^k` for `k=K_i,...,K_i+W-1`. Proper gcds
between their leading norms and `N`, and useful target involutions,
give factors immediately. Otherwise a leading norm `A<=H` is coprime
to `N` and has the exact split

\[
 A=S R,\qquad \gcd(S,R)=1,
 \tag{5}
\]

where `S` contains all factors from the fixed known-prime base. The
modular primorial-power test computes `S`; division is only by those
known primes. An admissible query has `R<=B`.

The partial ideal identity is

\[
 g_i^k=\prod_{p\text{ in base}}g_p^{v_p}\,h_\kappa^\epsilon,
 \qquad \epsilon\in\{-1,+1\},
 \tag{6}
\]

with the residual term omitted when `R=1`. The smooth prime word obeys

\[
 \prod_p p^{|v_p|}=S.
 \tag{7}
\]

If the residual key is new, append its generator to the queue before
recording the relation. The selected relation row is

\[
 r_i=k e_i-\sum_p v_p e_p-\epsilon e_\kappa.
 \tag{8}
\]

This uses the whole residual norm as one ideal generator. Repeated
factors inside `R` do not produce a hidden prime decomposition.
After selecting one row for `g_i`, proceed to the next queued
generator. A generator is never processed twice.

## 3. Weighted dominance includes self-reference and conjugation

Give every generator the fixed positive column weight

\[
 w_j=\log d_j.
 \tag{9}
\]

Let `c` be the full exponent vector of the right side of (6), after
collecting contributions in the augmented columns. Equations (5)--(7)
give

\[
 \sum_j|c_j|w_j
 \le\log S+\log R
 =\log A\le\log H<k\log d_i=k w_i.
 \tag{10}
\]

For `R=1`, use `log R=0`. If several contributions enter the same
column, collecting them can only decrease the left side.

Since `r_i=k e_i-c`, inequality (10) implies `k>|c_i|` and hence a
positive diagonal entry. More precisely,

\[
 \boxed{\displaystyle
 (r_i)_i w_i-\sum_{j\ne i}|(r_i)_j|w_j
 \ge k w_i-\sum_j|c_j|w_j>0.
 }
 \tag{11}
\]

This proves strict dominance even when the target generator appears
in its own right-side decomposition. A same-orientation residual
subtracts one from the diagonal; a conjugate residual adds one. Both
are covered by (10)--(11). Negative smooth-word exponents also cause
no change to the proof.

When all queued generators have a selected row, the full augmented
matrix is square. Scaling column `j` by `w_j` makes it strictly row
diagonally dominant, so it is nonsingular. The integer relation rows
are independent over the rationals.

At an intermediate stage, restrict the columns to the already
processed generators. Removing other columns only strengthens (11).
Thus every selected prefix is independent as well. The generic
two-primary relation extractor can safely try partial presentations
before closure, as the implementation does.

The proof uses logarithms only to establish independence. The program
does not approximate them or construct the potentially enormous
integer `d_i^k`; `k>=K_i` certifies the required strict power inequality.

### Integral elimination keeps the extractor in the original base

When a residual generator first appears, its coefficient in the new
row is `-epsilon`, a unit in the ordinary integers. Eliminate that
generator immediately. If `W_i` is the exponent word for its parent
in the original prime base, its definition is

\[
 W_{\mathrm{new}}=\epsilon(k W_i-v).
\]

Initial generators have the standard unit words. Every new definition
has one previously existing parent and a known-prime correction word,
so these definitions form a forest. No residual factorization is
introduced by this class-word identity.

A query whose residual key already exists emits the original-base row
`k W_i-v-epsilon W_residual`. A smooth query emits `k W_i-v`.
Queries defining a new key emit no original-base row: they supply the
unit relation used to eliminate that new column.

These are integer presentation transformations. After `s` augmented
rows and `C` definitions, exactly `s-C` rows remain on the original
base. They are independent: reversing the unit-column eliminations
would otherwise give a dependence among the augmented rows, which
(11) excludes. Thus the implementation needs no rational-rank test,
and its extractor has the original dimension `m`. At full closure
there are exactly `m` independent original-base rows.

The coefficient bits must be charged. Put
`M=max(2, all selected exponents)` and `V=ceil(log_2 H)`. Each smooth
correction coordinate has absolute value at most `V`. Along a
definition path of depth `D`,

\[
 \|W_{\mathrm{child}}\|_\infty
 \le M\|W_{\mathrm{parent}}\|_\infty+V,
 \qquad
 \|W\|_\infty\le(V+1)M^D.
\]

A resulting relation therefore has coefficient bit length

\[
 O\bigl(\log(n+1)+(D+1)\log M\bigr),
 \qquad D\le C,\quad M\le O(n)+W.
\]

The implementation stores the forest and expands only endpoints of
emitted relations. There are at most `m` such relations and at most
two endpoint expansions per relation. Caching endpoint words bounds
the number of traversed definition nodes by `2mC`; each traversal
updates at most `m` integer coordinates with the bit bound above.
The augmented matrix is retained as an audit certificate, not used as
the extractor's linear-algebra dimension.

When literal partial-ideal verification is enabled, every selected
identity is checked before it enters the forest. Exact integer
substitution then establishes all emitted original-base relations.
The extractor may use its already-established-relations contract,
avoiding repeated evaluation of every long original-base word at each
extraction. The independent checker replays the integer substitutions.

## 4. Finite closure and cost

Replace `B` by `min(B,H)` for counting purposes. By the
[residual-key bound](cofactor-collision-quota.md), the total number of
new keys is at most

\[
 K_B=B\,\mathrm{bitlength}(B).
 \tag{12}
\]

Let `m` be the initial base size and `C` the number of keys actually
introduced. Then

\[
 C\le K_B,\qquad r=m+C,
 \qquad \text{norm queries}\le r W.
 \tag{13}
\]

Assume each generator reached by this deterministic queue has an
admissible query in its allotted window, unless a factor is found
earlier. The queue must then close after at most `r` selections. The
same condition could be imposed on all allowed residual keys, which
is a stronger sufficient hypothesis.

Each group operation, split, and form reduction has polynomial bit
cost in the input size and the exponent's binary length. Augmented
relation entries have `O(n+log W)` bits. Eliminated original-base
entries obey the explicit depth bound above. Independence is supplied
by (11) and integer elimination, and the existing relation extractor
has polynomial cost in `m` and those coefficient bit lengths.
Thus the overall cost has the form

\[
 r W\,\mathrm{poly}(n,\log W)
   +\mathrm{poly}(r,n,\log W).
 \tag{14}
\]

All rejected queries are included in `W`. Use a deterministic ordered
key dictionary for a worst-case comparison-cost implementation.
Finite Python controls use ordinary dictionaries.

At `B=H`, every coprime leading norm is admissible, so `W=1` suffices.
This gives an explicit finite closure with `r<=m+H bitlength(H)`.
That bound is exponential in the binary input length and does not
improve the current score. Smaller caps trade a smaller potential
generator set for an unresolved admissible-hit problem.

## 5. Factoring and score contract

The initial prime forms must generate the intended class group for
the existing full-presentation factoring theorem. Adding residual
classes preserves that generation property. Once closure is complete,
the independent genuine relations meet the generic extractor's
contract, without a supplied annihilator or an unknown-prime
factorization routine.

The usual primality and perfect-power preprocessing, discriminant
choice, and useful-involution result remain part of the upstream
factoring reduction. A small experimental cutoff need not satisfy
the generating-base hypothesis. Consequently failure of a finite
small-base run is not a primality certificate.

For a polynomial-size generating base, bounds

\[
 B(n),W(n)\le
 \exp\!\left(O\!\left(\frac n{(\log n)^\alpha}\right)\right)
 \tag{15}
\]

would imply the same exponent parameter `alpha` through (12)--(14).
A bound on the actual reached-key count `C(n)` can replace the cap
bound. Polynomial powers of `r` change the constant in the exponent,
not `alpha`.

The missing theorem is a uniform admissible-hit estimate along the
powers of the generated ideal classes, together with a suitable cap
or reached-key bound. Global class-group density alone does not
provide this estimate. The construction removes residual
factorization and rational-rank selection from this closure mechanism;
it does not assert the remaining analytic bound.

## 6. Exact finite controls

[check_cofactor_ideal_closure.py](../tools/check_cofactor_ideal_closure.py)
replays selected powers, partial ideal identities, and every complete
relation word. It verifies (11) by integer products and independently
checks augmented and original-base row ranks by exact rational
elimination. All admitted candidates were also checked during the runs.

| Input | Bits | Prime cutoff | Norm queries | Residual generators added | Selected rows / active generators |
|---|---:|---:|---:|---:|---:|
| `1009 * 1013` | 20 | 5 | 6 | 4 | 6 / 6 |
| `1019 * 1021` | 20 | 5 | 8 | 6 | 8 / 9 |
| `1031 * 1033` | 21 | 5 | 5 | 3 | 5 / 5 |
| `10007 * 10009` | 27 | 7 | 11 | 10 | 11 / 13 |
| `10037 * 10039` | 27 | 7 | 16 | 15 | 16 / 17 |
| `1000003 * 1000033` | 40 | 19 | 301 | 70 | 71 / 73 |

All six runs use residual cap 10,000 and window 32. Each returns a
proper factor. Partial rows are sufficient in several cases, as their
torsion image already contains a useful involution. These traversals
introduce new generators, so the table is not a matched-traversal
comparison with the fixed-base collision selector.

The controls retain unfactored norms `49=7^2`, `121=11^2`,
`413=7*59`, `91=7*13`, and `77=7*11`. These displayed factorizations
are verification data and are never supplied to the algorithm. The
`N=1013` control contains a conjugate self-reference in a residual
generator's selected relation and completes without returning a factor.

Three budget controls exhaust a smooth-only window, a query cap, and a
generator cap. Another control, `N=31*37` with cutoff 7, closes its
small subgroup presentation without finding a factor. It records the
need for the generating-base hypothesis rather than treating closure
as an unconditional prime/composite decision.

Full records are in
[cofactor-ideal-closure.json](../evidence/cofactor-ideal-closure.json).

## 7. Larger controls and word-bit accounting

The integral-elimination implementation also factors the existing
92-bit, 124-bit, and 144-bit control inputs:

| Input bits | Prime cutoff | Residual cap | Norm queries | Original base size | Original-base relation rank | Active augmented generators | Maximum definition depth | Maximum relation-entry bits |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 92 | 1,021 | `1021^2` | 81,430 | 85 | 81 | 1,062 | 26 | 156 |
| 124 | 4,093 | `4093^2` | 1,391,826 | 283 | 280 | 5,321 | 57 | 435 |
| 144 | 8,191 | `8191^2` | 3,395,784 | 513 | 510 | 10,245 | 71 | 540 |

The returned factors are `2147483647`, `2305843009213693951`, and
`4108225714090758818561`, respectively. The runs select 1,058, 5,318,
and 10,242 augmented rows. Their extractor dimensions are only 85,
283, and 513.
The recorded coefficient sizes explicitly charge the growth caused by
the definition chains.

Every selected partial ideal identity was verified during each run.
An independent replay recomputes every selected power and partial
identity, checks each weighted dominance inequality using integer
power thresholds, and repeats all exact forest eliminations. This
certifies the original-base relations and their independence without
performing a dense rank computation on the augmented matrices.

The [large-control manifest](../evidence/cofactor-ideal-closure-large.json)
records source hashes, full-record hashes, parameters, and replay
counts. The complete witnesses are
[cofactor-ideal-closure-92.json](../evidence/cofactor-ideal-closure-92.json)
[cofactor-ideal-closure-124.json](../evidence/cofactor-ideal-closure-124.json),
and [cofactor-ideal-closure-144.json](../evidence/cofactor-ideal-closure-144.json).
For example, replay the last with

```text
python3 tools/check_cofactor_ideal_closure.py \
  --replay evidence/cofactor-ideal-closure-144.json
```

These finite experiments establish the implemented reduction and its
bit accounting on those inputs. The uniform estimate in Section 5 is
still the open step.
