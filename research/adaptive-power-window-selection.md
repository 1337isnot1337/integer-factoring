# Adaptive power windows and a fresh 144-bit factorization

The fixed-base power selector can retain its successful rows and extend
only unfinished targets. This gives an exact deterministic search whose
window need not be chosen in advance. A fresh 144-bit control factors
after 8,232,730 candidate norms. The missing uniform short-window theorem
remains open; the complete factoring score is still one half.

This note builds on [smooth leading-power selection](smooth-leading-power-selection.md)
and [determinant-free relation halving](determinant-free-ambiguous-extraction.md).
The new implementation is
[adaptive_leading_power_windows.py](../tools/adaptive_leading_power_windows.py).

## 1. Continue the exact candidate sequence

Use the fixed prime forms, discriminant, height `H`, and starting powers
`K_p` from the fixed-base selector. An exhausted phase of absolute window
`W` has attempted every target. Each target either has an accepted smooth
leading norm or has tested all powers

$$K_p,K_p+1,\ldots,K_p+W-1.$$

To extend to `W'>W`, retain every accepted row unchanged. For each remaining
prime, construct the first untested form by binary exponentiation:

$$g_p^{K_p+W}.$$

Then test successive powers through exponent `K_p+W'-1`, stopping that
target at its first accepted norm. Coprime CRT composition, the HNF
fallback, and the primorial smoothness filter are unchanged.

This schedule has no repeated candidate within a target and no skipped
exponent in its tested prefix. Accepted targets are not searched again.
The global order differs from one wide scan, so early extraction can
return after a different set of targets. Every accepted row remains the
same first successful row for its target at the covered exponent.

For a norm with prime word `v`, the retained row is `k e_p-v` and satisfies

$$k\log p>\log H\ge\sum_q |v_q|\log q.$$

Target primes are distinct. Thus the designated minor stays strictly
dominant after logarithmic column weighting. Adding later rows preserves
their independence. Partial relation halving has exactly its earlier
torsion-image scope; an unsuccessful phase is not a primality certificate.

## 2. Validate retained state without assuming failed candidates

Before resuming, the implementation rebuilds the prime base and checks
its roots and forms against the stored base. It independently reconstructs
each accepted target power, verifies its written factorization over the
known base, decodes its word, and checks the exact relation row. It also
checks that every target has an attempt record and that each unfinished
target reached the old boundary.

The factor exponents are bounded by `bits(H)` before any power product is
formed. Validation therefore does not create huge integers from an
arbitrary binary-encoded exponent in a malformed norm record.

The total candidate count equals the sum of the per-target attempt
counts. The maximum tested exponent follows from those counts. For a
successful nonunit norm, the number of known-prime divisions is its
prime multiplicity plus the one-based index of its largest prime factor;
the unit norm uses one division. These identities also validate the
division count of a completed phase using the primorial filter.

This validation proves the stored accepted relations. It does not prove
that every earlier rejected candidate was actually tested. That history
comes from executing the preceding phase or repeating its complete
search. The implementation and checker state this distinction explicitly.

## 3. Doubling removes the need to know a successful window

For each target prime `p`, index its candidate sequence by positive
integers `j`, with exponent exactly `K_p+j-1`. Define its first-hit index

$$w_p=\min\{j\ge1:\text{the leading norm of }g_p^{K_p+j-1}
       \text{ passes the fixed-base norm checks}\}.$$

The checks include coprimality to the input and complete smoothness over
the known base. Suppose a positive integer `W_*` satisfies `w_p<=W_*`
for every target. Equivalently, each target has an acceptable norm among
the `W_*` exponents `K_p,...,K_p+W_*-1`.

In an actual execution, let `t_p` count all candidates tested for target
`p` across every phase. They form one prefix, with no repetitions. The
algorithm retains the first successful row and never resumes that target,
so `t_p<=w_p` at every time. A factor found by a gcd or partial extraction
can stop the execution sooner; an unvisited target has `t_p=0`. Hence
the cumulative number of candidate norms, including all rejected values,
satisfies the sharper bound

$$\sum_p t_p\le\sum_p w_p\le mW_*.$$

This argument applies to any increasing absolute window schedule that
eventually covers the first hits. It does not charge the unused part of
a phase's nominal window: testing a target stops at its retained hit.

To avoid knowing `W_*` in advance, run absolute windows `1,2,4,8,...`.
The boundary `2^ceil(log_2 W_*)` covers every first hit. Thus at most
`1+ceil(log_2 W_*)` phases and `ceil(log_2 W_*)` resumptions are needed;
when `W_*=1`, exactly one phase suffices. The candidate bound remains
`mW_*`, despite rounding the final phase boundary upward.

All setup work must still be counted. Each resumed phase rebuilds the
base and validates the retained state. If `C_base` is the base setup
cost and `C_test` bounds one norm test and its successor update, a bound
for the full computation is

$$O(mW_*C_{\rm test})
 +O(\log(W_*+2)\,C_{\rm base})
 +\mathrm{poly}(m,b,\log Y,\log(W_*+2)).$$

The last term includes retained-state validation, binary exponentiations
at phase boundaries, and relation extraction. The prime product is
reconstructed when a phase starts; its written size is charged in the
setup and modular-reduction cost. No rejected prefix is replayed merely
to reach the next boundary.

Under the [proved generator theorem](polylog-prime-form-generators.md),
one can choose a polynomial-size base. Consequently a uniform bound

$$W_*(b)=\exp\!\left(O\!\left(\frac{b}{(\log b)^\alpha}\right)\right)$$

would give score `alpha`, with the usual primality, perfect-power, and
recursive-splitting steps. The algorithm would not need to know the
constant or the successful window beforehand. This is a conditional
consequence. No such uniform bound is proved here.

## 4. A fresh 144-bit control

The input primes come from a specified SHA-256 seed and have complete
recursive Lucas certificates. The search receives only their product
and its public search parameters:

$$\begin{aligned}
N&=17355670471043633865640786569418441084204301\\
 &=4108225714090758818561\cdot4224614633883335075341.
\end{aligned}$$

Both factors have 72 bits. The prime cutoff is 8191 and the base has
513 generators. The small experimental cutoff is not asserted to
satisfy the much larger uniform analytic cutoff.

| Phase | Absolute window | New norm tests | Retained rows | Outcome |
| --- | ---: | ---: | ---: | --- |
| Initial search | 50,000 | 8,077,846 | 494 | Window exhausted |
| Continue unfinished targets | 200,000 | 154,884 | 499 | Proper factor found |

The total is **8,232,730 norm tests and 182,539 known-prime divisions**.
Every leading norm has at most 73 bits. The second phase retains the
494 earlier rows and adds five more. The final presentation has
deficiency fourteen; this is a successful partial extraction, not a
claim that every such partial presentation must succeed.

The prime certificates are used to certify the test input. They are not
passed to the norm-search algorithm. Neither phase calls a routine to
factor an unknown leading norm. The accepted norms are divided only by
the known base primes after passing the exact primorial test.

## 5. Reproduction

The complete stored state is in
[adaptive-leading-power-windows.json](../evidence/adaptive-leading-power-windows.json).
Run:

```sh
python3 tools/check_adaptive_leading_power_windows.py
```

The default verifies the input-prime certificates, all accepted prefix
relations, and a fresh continuation of the 154,884 later candidates. It
also runs ten small searches, including six resumed phases, checks the
accepted prefixes with independent HNF composition and literal division,
and rejects malformed schedules and work counters.

To repeat the complete 8,232,730-candidate search:

```sh
python3 tools/check_adaptive_leading_power_windows.py --full
```

Use `--regenerate` to rebuild the input certificates and stored evidence.
The original two-phase search was executed in full during this session.
Passing these controls establishes this finite result and the stated
implementation contract. It does not supply the missing worst-case
window bound.
