# Linear space for the regular modular factors

## Result and scope

Fix a prime `p`. The original regular-partition crossing family has
deterministic factor, piece, width, and index queries with bounds

```math
\text{time}=2^{O_p(\lambda)},\qquad
\text{working space}=O_p(\lambda+1)
\quad\text{bits}.                                      \tag{1}
```

The factor values lie in the fixed prime field. The space includes the
current compressed labels, intermediate partition, channel index, and arithmetic
workspace. It does not include the calling evaluator's stack or its
coefficient-ring elements. Thus (1) is a local interface; a global
factoring space bound requires separate recursive accounting.

This improves the unspecified polynomial-space bound in
[the characteristic-two factor routine](modular-factor-space.md) and
[its fixed-characteristic extension](modular-factor-space-odd.md).
It retains their `2^O(lambda)` time and the exact half-rate width
estimate from [the regular modular construction](modular-crossing.md)
and [the exact endpoint analysis](rough-input-transfer.md).
No smaller contraction coefficient is asserted.

The source operators and crossing identities are Anthropic's Appendix B,
pinned to
[`16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The improvements here are ordinary exact arithmetic and streaming
arguments; no new Lean formalization is used.

## 1. Extract a coefficient without storing its polynomial

Let `P(S)` be a polynomial over `F_p` of degree at most `L`. Choose

```math
q=p^r>L+1.
```

It suffices to take the least such power, so `q<=p(L+1)` unless `L=0`,
which is handled directly. Construct `F_q` as a polynomial quotient of
degree `r` over `F_p`. For `0<=i<=L`, the exact identity is

```math
[S^i]P(S)
 =-\sum_{a\in\mathbf F_q^*}P(a)a^{-i}.                 \tag{2}
```

Indeed, `sum_(a!=0) a^k` is zero unless `q-1` divides `k`, when it is
`q-1=-1` in the field. For the exponents `k=j-i`, with `0<=j,i<=L`,
only `k=0` is divisible by `q-1`. The strict inequality `q-1>L`
is essential. In characteristic two the minus sign in (2) disappears.
The final value belongs to the embedded prime field because it is the
original coefficient.

Enumerate field elements by their `r` base-`p` coordinates. The routine
keeps one element, one polynomial value, and one running sum.
Exponentiation computes inverses and powers in `poly(log q)` field
operations and `O_p(log q)` bits of workspace. No table of field
elements, powers, evaluations, or coefficients is needed.

For completeness, the extension field can be constructed in polynomial
time in `q` and `O_p(log q)` space. Enumerate monic degree-`r`
polynomials. Test irreducibility using successive modular `p`th powers
and polynomial gcds: require `X^(p^r)=X` modulo the candidate and
`gcd(X^(p^j)-X,f)=1` for `1<=j<=floor(r/2)`. A reducible degree-`r`
polynomial has an irreducible factor of degree at most `r/2`, so these
tests are sufficient. An irreducible polynomial of degree `r` exists.
Each polynomial and gcd workspace has `O_p(r)` bits, and there are
`p^r=q` candidates. The same construction is reused throughout a query.

## 2. Evaluate the folded factor at one field element

The existing characteristic-two proof gives the complete angle profile
for an operator label `mu` of weight `h` and particle count `ell`:

```math
P(S)=(1+S)^{h+\ell}
       \sum_\eta\operatorname{per}M_\eta(S),           \tag{3}
```

where `eta` ranges over the allowed unaffected subsets of the input
and output states. The matrix has `ell` rows and `ell` signed columns.
Its entries, for operator part `n` and state mode `k`, are

```math
S^{(n\mathbin{\mathrm{AND}}k)-1}\quad\text{for a creator},
\qquad
S^{n\mathbin{\mathrm{AND}}(\mathrm{NOT}\ k)}
   \quad\text{for an annihilator}.                    \tag{4}
```

All modes are odd, so the creator exponent is nonnegative. The degree
of (3) is at most `2h`. These are exactly the folded source entries;
the expanded dual polynomial is never constructed.

At a particular `a in F_q`, enumerate `eta` and then use Ryser's formula

```math
\operatorname{per}M_\eta(a)
 =\sum_{J\subseteq[\ell]}
       \prod_{n\in\mu}\sum_{j\in J}M_\eta(a)_{n,j}.
                                                               \tag{5}
```

For each `J`, accumulate one row sum and one product. Compute each
entry in (4) from its two mode indices by field exponentiation.
There is no need to store even the `ell` by `ell` matrix. The signed
column list occupies `O(ell log(h+2))` bits, and the two subset counters
occupy `O(sqrt h)` bits.

The input and output weights are at most `4h`. Distinct odd parts with
particle count `s` have total weight at least `s^2`. Consequently all
labels and states in one factor query have `O(sqrt h)` parts. There
are `2^O(sqrt h)` choices of `eta` and `J`. Equations (3)--(5) therefore
give one value `P(a)` in time
`2^O(sqrt h) poly(h)` and space `O(sqrt h log(h+2))`.

Apply (2) with `L=2h`. Its `O(h)` evaluations multiply the time by only
a polynomial. The angle factor of any specified index is obtained in

```math
2^{O(\sqrt h)}\operatorname{poly}(h)\ \text{time},
\qquad O(\sqrt h\log(h+2))\ \text{bits}.               \tag{6}
```

The bottom factor is a scalar permanent of the usual binary creation
and annihilation tests. Apply the same subset loops directly over `F_2`;
no interpolation is needed. Zero-particle cases and angle indices
outside `[0,2h]` are handled directly.

### Every fixed prime

For fixed odd `p`, use the unaffected multiplicities, nonzero factorial
weights, and folded one-particle polynomials from
[the fixed-characteristic proof](modular-factor-space-odd.md).
Evaluate each folded polynomial directly at `a`; it is a sum of at most
`h` terms whose binomial coefficients can be computed modulo `p`.
Keep the usual signs in Ryser's formula and the minus sign in (2).

Every multiplicity is below `p`. A retained partition has
`O_p(sqrt h)` particles, so the unaffected-multiplicity iterator and
Ryser subsets still number `2^O_p(sqrt h)`. Their lists and counters
have `O_p(sqrt h log(h+2))` bits. The degree remains at most `2h`.
Thus (6) holds for every fixed `p`, with constants depending on `p`.
No denominator divisible by `p` is introduced by this improvement.

## 3. Compress the partitions and remove the recursion stack

A list of `O(sqrt h)` separate binary mode indices would still cost an
extra factor `log h`. The following elementary code removes it. For
strictly increasing positive parts `j_1<...<j_ell`, of total weight at
most `H`, put `g_i=j_i-j_(i-1)>=1`, with `j_0=0`. Then

```math
\sum_{i=1}^{\ell}(\ell-i+1)g_i\le H.
```

Weighted AM--GM and `ell! >= (ell/e)^ell` give

```math
\sum_i\ln g_i
\le \ell\ln(H/\ell)-\ln(\ell!)
\le \ell\ln(eH/\ell^2)
\le2\sqrt{H/e}.                                      \tag{7}
```

The last inequality follows by maximizing the preceding expression over
positive real `ell`. Also `ell<=sqrt(2H)`. Encode each positive gap by
its Elias-gamma code, of length `2 floor(log_2 g_i)+1`, and place one
initial sentinel bit before the concatenation. Its total length is at most

```math
1+\sqrt{2H}+\frac4{\ln2}\sqrt{H/e}<1+5\sqrt H.           \tag{8}
```

Scanning the code decodes any requested part, or computes its particle
count, weight, and last part. This needs only `O(log(H+2))` temporary
bits. Appending a gap is a shift and an insertion. Popping the last gap
scans to find its code length and shifts it away. Repeated scans cost
only polynomial time in `H`.

At a fixed odd prime, encode the increasing *distinct* mode list this
way and attach each multiplicity using a fixed `O_p(1)` number of bits.
There are `O_p(sqrt H)` such modes. Thus a regular partition at every
fixed prime has an explicit `O_p(sqrt H+1)`-bit code. No count table is
used to encode it.

Enumerate distinct-odd partitions by an iterative depth-first search.
Keep the current gap code, its total weight, last part, and one next-part
candidate. Appending an allowed candidate descends to a child. When no
candidate fits, decode and pop the last part; its successor is the next
candidate at the parent. At a leaf of the target weight, emit the code
and perform the same pop. No per-level frame stack is kept. The stored
code always describes a partition of weight at most `H`, so (8) applies
throughout the search, including partial states.

For fixed odd `p`, use the same rule on distinct modes and their bounded
multiplicities. Upon returning from a child, the popped mode and
multiplicity determine the next candidate pair at the parent. Candidate
scans introduce at most a polynomial factor in `H`.

Each partial list is itself one regular partition of weight at most the
target. Hence the number of search nodes is at most a polynomial factor
times the cumulative partition count through that target. The ordinary
partition bound gives `2^O_p(sqrt h)` nodes for target at most `4h`.
Enumerating candidates in a fixed order therefore takes
`2^O_p(sqrt h) poly(h)` time and `O_p(sqrt h+1)` space.

The permanent routine can also keep its partitions compressed. At `p=2`,
represent `eta` by a bit mask on the output's positions. Reject masks
whose selected parts do not occur in the input or have the wrong particle
count. Generate the signed columns by scanning the output and input
codes and skipping the selected unchanged modes. For each row and
Ryser subset, scan the codes again to recover the needed mode index.
There is no stored signed-column list or matrix. The masks have
`O(sqrt h)` bits, and the field/scanning variables have `O(log h)` bits.

For fixed odd `p`, use a bounded-multiplicity mask in a fixed base, and
generate repeated signed columns by scanning the coded modes and their
multiplicities. Its length remains `O_p(sqrt h)` bits. All inversions
and factorial factors are the same prime-field units as in the existing
fixed-characteristic formula. This reduces the factor space in (6) to

```math
\text{space}=O_p(\sqrt h+1)\ \text{bits},               \tag{9}
```

with unchanged exponential-time order. Only a fixed number of encoded
partitions and masks is live at once.

## 4. Label, channel, and width queries

Count by incrementing an `O_p(sqrt h)`-bit counter. Unrank by scanning
until the desired counter value is reached. Repeating these scans in a
polynomial number of weight or angle loops leaves the time bound
unchanged. This deliberately trades time already allowed by the local
interface for smaller space.

The least label weight with at least `2^lambda` regular partitions is
found in the same way. The explicit equal-weight block construction in
[the modular note](modular-crossing.md), Section 4, bounds the search by
`h<=8lambda^2` at `p=2`, and `h<=12lambda^2+lambda` for fixed odd `p`.
One may stop each count at `2^lambda`. This gives `h=O_p(lambda^2)`
without a count table or a monotonicity assumption on partition numbers.

For the source channel sets, enumerate the angle and intermediate state:

```math
(i,\pi),\quad 0\le i\le2h,\quad
|\pi|=s+i\quad\text{or}\quad |\pi|=2h-s-i.            \tag{10}
```

All nonnegative weights in (10) are at most `4h`. The number of channels
has `O_p(sqrt h+log(h+2))` bits. Its exact count, binary width, and
unranking therefore fit (9). The current channel contains only one
compressed intermediate partition.

An entry of a complete piece is the sum over these channels of two
scalar factor entries. Enumerate the channels sequentially. The product
of two `2^O_p(sqrt h)` time bounds is still `2^O_p(sqrt h)`, and the
space is the maximum workspace plus a fixed number of compressed labels.
Thus complete-piece queries satisfy the same time bound and (9). When an individual channel
count is capped by the trivial identity factorization on the binary
label prefix, the newly required piece-entry query has this bound.

Substituting `h=O_p(lambda^2)` in (9) proves (1), including the metadata,
binary prefix, and individual-width cap interfaces.

## 5. Exact checks

Run

```sh
python3 tools/check_near_linear_factor_space.py \
  --output evidence/near-linear-factor-space.json
```

The production module
[`linear_modular_factors.py`](../tools/linear_modular_factors.py)
evaluates one requested angle coefficient by field summation and streamed
scalar Ryser loops on gap-coded partitions. It has no polynomial
coefficient array, matrix table, expanded mode list, or partition-recursion
stack. The checker also retains the uncompressed intermediate routine
as a comparison implementation. Its verification code compares both with
the independent source recursion and the prior polynomial-profile
routine, checks field-sum extraction independently, and checks compressed
partition/channel indexing against enumerated reference sets. It verifies
the gap-code length inequality on its enumerated finite partitions. Finite
tests check the formulas; the uniform bit-space claim is the algorithmic
argument above.
