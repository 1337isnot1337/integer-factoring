# Constructive bounded-width class-group navigation

## Status

This is a deterministic finite search that can cross the
[strict height-descent traps](prime-form-height-traps.md). It retains a
bounded frontier, records exact generator words, and checks every returned
answer. A pruned search can miss an existing word. Such a miss is reported
as **inconclusive**, not as the negative answer required by
[the square-decision reduction](smooth-square-decision-reduction.md).

The implementation navigates the explicit height trap with width 16 and
constructs complete dominant relation matrices for the same 51-bit and
60-bit inputs used by [the exact square-box search](factored-square-box-search.md).
Its costs are polynomial in its explicit width and depth parameters.
No polynomial width or uniform success theorem is proved, so these
constructions do not improve the complete factoring score.

## 1. Two bounded frontiers

Let \(g_1,\ldots,g_m\) be supplied classes, and fix a width \(W\),
left depth \(L\), and right depth \(R\). All classes are stored as
canonical reduced forms. Build a frontier from the identity using the
\(2m\) moves \(g_i^{\pm1}\), retaining at most \(W\) new class values
at each depth. Store a generator word and its path depth for every
retained class. The resulting identity-side index has at most \(1+WL\)
entries and can be reused for several targets.

For a target \(h\), build the same kind of frontier starting at \(h\),
up to depth \(R\). Before pruning a newly generated state, look for its
class in the identity-side index. If the meeting class has words
\(a,c\), then

\[
\prod_i g_i^{a_i}=h\prod_i g_i^{c_i},
\qquad h=\prod_i g_i^{a_i-c_i},
\qquad \sum_i|a_i-c_i|\le L+R.
\tag{1}
\]

Both meeting identities and the final word are checked. A class retained
at an earlier depth is not retained again at a later depth; its earlier
word gives at least as much remaining path budget. Within one layer,
equal canonical classes are deduplicated while retaining an explicit
word.

The implementation uses sorted arrays, sorting, merge operations, and
binary search. Its deterministic accounting does not depend on expected
hash-table behavior. It also uses no translation-invariance assumption
for the ordering on reduced forms: ordering serves only exact lookup and
frontier selection.

## 2. Retaining different height ranges

The baseline `height` policy retains the \(W\) smallest leading
coefficients, with the remaining form coefficients breaking ties. This
allows increases between layers, but can still concentrate the search
near a restricted set of small-height classes.

The `stratified` policy divides the new states by the bit length of their
leading coefficient. Within each bucket it uses canonical form order.
It then takes the first state from each nonempty bucket, the second state
from each bucket, and so forth, stopping after \(W\) states. This keeps
some larger-height candidates available for later compositions.

The policy is a concrete search rule, not a theorem about the geometry of
class-group paths. Increasing width does not supply a general monotonic
success guarantee for a pruned heuristic. The finite controls compare
the actual returned paths and misses at stated parameters.

## 3. What a failed search proves

If neither frontier ever pruned a class, the left index is the complete
ball of radius \(L\), and the right search visits the complete target
ball of radius \(R\). Any word of length at most \(L+R\) can be split
between these two balls. Therefore failure in this unpruned case is an
exact negative answer for that word-length bound.

If either side pruned a class, the implementation returns
`incomplete_search` and sets `negative_answer_is_exact` to false. It does
not promote the frontier's exhaustion to an arithmetic nonrepresentation
claim. In particular, this implementation cannot simply be substituted
for the exact variable-cap decision oracle.

The full-word bound here is \(\ell_1\), whereas the exact square-box
search uses separate coordinate caps. Both are useful concrete searches,
but their negative-answer contracts differ.

## 4. Explicit cost

The identity-side expansion uses at most \(2mWL\) class compositions.
Each target expansion uses at most \(2mWR\). These bounds include the
partial layer in which a meeting is found. Stored words have \(m\)
integer coordinates, each of magnitude at most the corresponding depth.
Candidate sorting, exact lookups, and merging the retained arrays add
polynomial factors in \(b,m,W,L+R\), where
\(b=\lceil\log_2|\Delta|\rceil\).

Thus the entire search is deterministic and polynomial in these explicit
parameters. For a collection of \(k\) targets, the expansion count is at
most

\[
2mW(L+kR).
\]

The evidence's `group_compositions` fields count only these frontier
expansions. Setup inversions, the polynomial final-word checks,
determinant computation, and factor extraction are separate. The data
are operation counts, not complete wall-clock benchmarks.

The retained index can have \(O(W(L+R))\) classes across both sides, and
one layer has \(O(mW)\) temporary candidates. These are explicit memory
costs. The depth is charged by value, not merely by its binary encoding.
A sufficiently large width can make the search exhaustive, but that
width can be exponential in \(b\). Establishing success with a useful
uniform width remains open.

## 5. Dominant relations from the paths

Fix \(T=L+R\), then choose \(K=T+1\) before defining the targets. For
each factor-base generator, search for a word for \(g_i^K\) of length at
most \(T\). A returned word \(w^{(i)}\) gives

\[
R_i=K\mathbf e_i-w^{(i)},\qquad
\prod_jg_j^{(R_i)_j}=1.
\]

The diagonal is positive and exceeds the sum of the absolute
off-diagonal entries by at least \(K-T=1\). Hence the matrix is
nonsingular, and its absolute determinant annihilates every supplied
generator. The driver verifies every row and every annihilation
identity, then applies [early-stop extraction](early-stop-class-factoring.md).

The known length-\(K\) word for the target cannot itself satisfy the
search budget: \(K>T\). Thus successful rows are not obtained by merely
reprinting the prescribed target powers.

The driver receives only \(N,m,W,T\) and the frontier policy. It
deterministically constructs the first \(m\) split unramified prime
forms of discriminant \(-4N\). No factorization, class number,
annihilator, or complete class enumeration is supplied. The finite bases
are not asserted to generate the full class group. A returned proper
divisor is checked directly; a missing relation or an unhelpful subgroup
does not certify primality.

## 6. Completed finite constructions

| Quantity | 51-bit input | 60-bit input |
| --- | ---: | ---: |
| \(N\) | 1621982935400141 | 1000000016000000063 |
| Prime forms \(m\) | 9 | 12 |
| Width \(W\) | 1024 | 2048 |
| \(L=R\) | 32 | 40 |
| Diagonal \(K\) | 65 | 81 |
| Retained identity-side classes | 30855 | 78137 |
| Identity-side expansion compositions | 536958 | 1826136 |
| Total target-side expansion compositions | 290230 | 2394099 |
| Maximum returned word length | 23 | 28 |
| Minimum row dominance slack | 42 | 53 |
| Derived annihilator | 19727896832033528 | 94161325139492316876528 |
| Annihilator bits | 55 | 77 |
| Returned factorization | \(32452843\cdot49979687\) | \(1000000007\cdot1000000009\) |

The identity-side frontiers prune 360196 and 1362181 class occurrences,
respectively. Their success is therefore certified by the returned
words, not by an exhaustive search guarantee. The target-side meetings
occur by depths 9 and 23, respectively.

A fixed-parameter calibration uses the 51-bit input, width 256, and
total depth 64. The smallest-height policy finds one of nine targets;
the stratified policy finds four. All unsuccessful searches are correctly
marked inconclusive. The wider successful run provides actual bounded
words for every target, demonstrating that those earlier misses cannot
be treated as exact negative answers.

For the explicit height trap
\(N=1022117\), \(P=\{3,7,11\}\), and
\(h=(17,10,60126)\), width 16 and \(L=R=6\) find

\[
h=g_7^5g_{11}^{-3}.
\]

This gives an eight-move path to the identity. Its first move necessarily
increases the leading coefficient, consistent with the proved trap.
The same depth budget at width one misses the word and reports an
inconclusive result.

## 7. Verification and reproduction

The small unpruned controls compare 975 decisions against independent
complete bounded-word expansion for discriminants of absolute value at
most 128. They include 798 positive and 177 exact negative decisions.
The trap control checks the complete returned path to the identity and
the narrow-frontier miss.

The independent factoring replay checks both sides of every meeting,
every word-length bound, all 21 dominant relations, and the annihilation
identities. It computes the determinants by rational elimination rather
than the construction's Bareiss routine and replays the exact factor
extraction. The returned factors are trial-prime certified only after
extraction.

The evidence consists of
[the 51-bit construction](../evidence/bounded-class-factoring-51bit.json),
[the 60-bit construction](../evidence/bounded-class-factoring-60bit.json),
and [the controls and independent replay](../evidence/bounded-class-word-search.json).

```sh
python3 tools/run_bounded_class_factoring.py --n 1621982935400141 --primes 9 --width 1024 --depth 64
python3 tools/run_bounded_class_factoring.py --n 1000000016000000063 --primes 12 --width 2048 --depth 80
python3 tools/check_bounded_class_word_search.py
```

The implementation is
[`bounded_class_word_search.py`](../tools/bounded_class_word_search.py).
This supplies an executable global path search beyond strict height
descent. The unresolved issue is a uniform success guarantee at a width
small enough to improve the factoring bound.
