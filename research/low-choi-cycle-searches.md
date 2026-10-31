# Two finite searches for regular binary q8 cycle sources

## Outcome and scope

These searches test sufficient Smith profiles for asymptotic coefficient tensor powers. They are distinct from the already known exclusion of completed native eight-label finite seeds. A small source can only improve the score if its larger tensor-power constructions satisfy the full transfer requirements.

Two explicit permutation families were screened over `F_2` at physical size
`q=8`. Every retained permutation is one cycle of length 64. Thus `S=I+U`
is a regular nilpotent operator and gives the complete nilpotent-resolvent
crossing with ideal column valuations `(0,64,...,64)`.

Neither finite family has a source passing the first row-Smith gate. These
are complete statements about the two stated finite classes. They do not
exclude other nonlinear label gauges, other feedback functions, nonpermutation
operators, or larger physical sizes. The factoring exponent is unchanged.

For `r=rank Phi(U)`, the exact number of row-Smith orders equal to one is

\[
k=\operatorname{rank}\begin{pmatrix}A&0\\A+C&A\end{pmatrix}-2r,
\qquad A=\Phi(U),\quad C=\Phi(U^2).
\]

The sufficient Smith cost is at least `64r+63k`. Success requires this to be
less than `64^sqrt(2)`, approximately 358.36. Leading rank at least six
already fails. The native programs therefore compute exact first jets only
for the remaining small-rank cases and use conservative rank bounds to skip
cases that cannot improve a recorded minimum.

## 1. Nonlinear de Bruijn feedback

On six binary bits, consider the invertible feedback register

\[
U(x)=(x\mathbin{\mathrm{div}}2)+32\bigl((x\bmod2)+g(x\mathbin{\mathrm{div}}2)\bigr),
\]

where addition inside the last parentheses is in `F_2`. The five-variable
Boolean function `g` is the sum of a polynomial of degree at most two and
the indicator of the zero input. There are 16 multilinear monomials of degree
at most two, so the program checks all `2^16=65,536` feedback functions.

It retains the 1,952 full 64-cycles. For each, it tests ten balanced partitions
of the six bit positions between the two physical labels. Complementary
partitions are identified because they transpose both physical labels and
preserve the relevant ranks. This gives 19,520 full-cycle/partition cases.

| Leading Choi rank | Cases |
|---:|---:|
| 3 | 80 |
| 4 | 752 |
| 5 | 1,406 |
| At least 6 | 17,282 |

The best first-jet floor among cases of leading rank at most five is 761,
with `(r,rank Phi(U^2),k)=(5,8,7)`. Its feedback truth table is the 32-bit mask
`4244438269`, with the row label using bit positions `{0,2,3}` and the column
label using `{1,4,5}`. The fast checker reconstructs this truth table from
the Gray-code polynomial index 103, confirms the full cycle, and recomputes
all three ranks independently.

## 2. Every binary linear gauge of the joined LFSR family

For each primitive six-bit feedback mask in

\[
\{3,27,33,39,45,51\},
\]

let `L(x)=(x div 2)+32 parity(x & taps)`. It has one nonzero orbit of length
63 and fixes zero. Swap the outgoing edges of states zero and one:
`U(0)=32`, `U(1)=0`, and `U(x)=L(x)` otherwise. This gives one 64-cycle.

The search allows every invertible binary linear change of the six label
bits. Changes of basis within either three-dimensional physical label are
local permutations and preserve all relevant ranks. It therefore suffices
to enumerate ordered complementary subspaces `A direct_sum B=F_2^6`, each
of dimension three, using one basis for each.

There are 1,395 possible row subspaces. Once one complement is fixed, all
others are graphs of linear maps from it to the row subspace, giving
`2^9=512` complements. Thus the exact class size checked is

\[
1395\cdot512\cdot6=4,285,440.
\]

The program enumerates each three-dimensional subspace once by its set of
eight vectors, then all 512 complement graphs. It does not restrict the
gauge to subsets of the original bit coordinates.

### Conservative rank prefilter

For the unjoined linear map, let `A_L,D_L` be its same-label three-by-three
diagonal blocks in the chosen decomposition. Its Choi rank is
`2^(rank A_L+rank D_L)`: the corresponding linear relation splits into
disjoint nonzero all-ones rectangles, each of rank one over `F_2`.

The joined-edge correction is the outer product of the vectors supported
on `{0,L(1)}` and `{0,1}`. In a general linear gauge, either physical vector
can have matrix rank two. Therefore its Choi rank can be **four**. The
prefilter uses the exact product of those two physical ranks and retains
every case whose joined Choi rank could be at most five. In particular,
unjoined linear Choi rank eight is included when the correction has rank
four. No coordinate-only rank-two correction bound is assumed.

The conservative prefilter retained 258,552 cases. Exactly 3,360 had actual
leading rank at most five:

| Leading Choi rank | Cases | Minimum exact first-jet floor at that rank |
|---:|---:|---:|
| 3 | 120 | 2,145 |
| 4 | 96 | 2,272 |
| 5 | 3,144 | 824 |

There were no survivors. A best low-leading-rank representative uses taps
`3`, row basis `(1,2,8)`, and column basis `(4,16,32)`, with vectors encoded
as six-bit integers. It has `(r,rank Phi(U^2),k)=(5,9,8)` and floor 824.

## 3. Reproduction

The default independent check is small:

```sh
python tools/check_low_choi_cycle_searches.py --output /tmp/low-choi-cycle-searches.json
```

It verifies the 1,395 subspace count, all six primitive feedbacks, class-size
identities, and both explicit representative permutations and first jets.
The aggregate scan counts in the output record the completed native runs;
this fast check does not independently repeat those exhaustive searches.

The full native reruns are optional:

```sh
c++ -O3 -std=c++20 tools/search_debruijn_q8.cpp -o /tmp/search-debruijn-q8
/tmp/search-debruijn-q8 2
c++ -O3 -std=c++20 tools/search_lfsr_linear_gauges_q8.cpp -o /tmp/search-lfsr-linear-gauges-q8
/tmp/search-lfsr-linear-gauges-q8
```

The second source includes the first source's packed rank helpers, so keep
them in the same directory. Each program prints representative improvements
as JSON lines and exact aggregate counts to standard error. No large scan
log is required by the fast verification target.
