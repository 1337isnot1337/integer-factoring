# A quantitative global transfer from raw rank sums

## Theorem and scope

The [global-block transfer](global-block-field-transfer.md) can use raw
families with exponentially many pieces. Their total rank cost replaces
the separate maximum-width and polynomial-piece-count assumptions.

Fix `1/2<=rho<1` and an integer `h>=1`. At every positive raw binary
length `t`, suppose an effective exact crossing family is supplied over
one of the coefficient fields allowed by the global-block theorem.
The field may depend on `t`. Each retained piece has displayed positive
integer channel dimensions `r_i,c_i<=2^(2t)`. Assume

```math
S_t:=\sum_{i=1}^{k_t}\sqrt{r_i c_i}
 \le 2^{\,2\rho t+C t/\log_2(t+2)}.                    \tag{1}
```

The channel dimensions are actual cardinalities. They need not be powers
of two, and they must not be rounded separately before applying (1).
They may exceed the minimum matrix ranks if explicit factors are easier
to generate. A zero piece may be discarded; no minimality assumption
is needed for the retained list.

Assume complete raw generation, field representation, scale, dimensions,
and scalar factor queries have uniform `2^O(t^h)` bit cost. This includes
the output lengths and the cost of constructing the coefficient field.
Use the explicit number-field presentation in the global-block theorem,
or the effective finite-field presentation in
[characteristic-independent recovery](characteristic-independent-recovery.md).

Then, with `alpha=log_2(1/rho)`, deterministic factoring has simultaneous
bounds

```math
T(n)=\exp\!\left(O\!\left(\frac n{(\log(n+2))^\alpha}\right)\right),
\qquad
S(n)=O\!\left(\frac n{(\log(n+2))^\alpha}\right).       \tag{2}
```

The output variants and their coefficient-field qualifications are the
same as in the two cited transfers. In particular, the number-field
case retains direct exact divisor counts. The conclusion includes the
endpoint `alpha=1` when the hypothesis holds with `rho=1/2`.

No better crossing family is constructed here. The best proved exponent
remains `1/2`.

## 1. Relation to existing rank-sum transfers

The multinomial rank-sum inequality below is the established finite-seed
amplification argument. The
[gap-free score transfer, pinned at `ae5115b3`](https://github.com/1337isnot1337/integer-factoring/blob/ae5115b3d6b4175693eb0899288b71d7b9b4629f/research/gap-free-seed-score-transfer.md)
already attains a computably enumerated seed infimum in the challenge's
limiting-score sense, without a convergence rate or a generation bound.
Its Section 6 explicitly handles the rank-sum infimum.

The present statement supplies the constant-factor bound (2) and its
simultaneous workspace bound under the quantitative assumptions above.
It combines that rank-sum inequality with one global raw block length.
It does not claim a new finite-seed amplification identity or a stronger
score-selection theorem.

## 2. Piece count follows from the rank sum

Every displayed dimension is a positive integer, so

```math
k_t\le S_t\le 2^{O(t)}.                               \tag{3}
```

The second inequality is uniform because `rho` and `C` are fixed.
All four factor tables have at most `2^(4t)` entries per piece.
Thus enumerating the full raw family still costs `2^O(t^h)` bits and
time, using `h>=1`. The normalization, denominator collection, and
finite-ring compilation from the global-block theorem retain that cost.

Choose a fixed integer `C_*` bounding all compiled data, field arithmetic,
raw queries, and `k_t` by

```math
G(t)=2^{C_*(t^h+1)}.
```

For input length `n`, take the same integer-root choice of raw length
`m=Theta((log(n+2))^(1/h))`, with a sufficiently large fixed divisor
in its definition, so that

```math
G(m)=O(n^{1/128}).                                    \tag{4}
```

Compile this one field and this one raw family. The existing early-leaf
cutoff ensures that every executed grid step has
`lambda>=n^(1/4)/8`, while a preliminary envelope gives
`lambda<=2^16(n+2)^3`. Consequently `m<lambda` for all sufficiently
large inputs. A fixed exhaustive procedure covers the finite exceptions.

## 3. Keep exact channel cardinalities through amplification

Write `lambda=bm+q`, with `0<=q<m`. For every type
`nu=(nu_1,...,nu_k)` of length `b`, let

```math
B_\nu=\sum_{\operatorname{type}(w)=\nu}
 A_{w_1}\otimes\cdots\otimes A_{w_b},
\qquad M_\nu=\frac{b!}{\prod_i\nu_i!}.                 \tag{5}
```

This is an ordered sum, with no averaging or division in the coefficient
ring. The explicit channel dimensions in the two cuts are

```math
R_\nu=M_\nu\prod_i r_i^{\nu_i},\qquad
C_\nu=M_\nu\prod_i c_i^{\nu_i}.                       \tag{6}
```

The channel chooses one word of that type and then one raw channel at
every position. Products of the raw factors give exact factorizations.

Put `w_i=sqrt(r_i c_i)`. The multinomial expansion of the nonnegative
real number `S_m^b` contains each individual term, so

```math
\sqrt{R_\nu C_\nu}
 =M_\nu\prod_i w_i^{\nu_i}\le S_m^b.                  \tag{7}
```

Only now pad each complete channel set to a power of two. The combined
width is at most `2b log_2 S_m+2`. Applying (1) gives

```math
\lceil\log_2R_\nu\rceil+\lceil\log_2C_\nu\rceil
 \le4\rho bm+O(bm/\log(m+2))+2.                       \tag{8}
```

No numerical evaluation of the square roots is part of the machine.
The machine computes the integers in (6); equation (7) is used only
to prove their size bound.

For `q>0`, append the single exact crossing tensor on `q` bits. Its
two identity factorizations each use `2^(2q)` channels. The coefficient
field is unchanged. This simpler tail costs `4q` width and introduces
no additional piece index. Summing all types now gives exactly the
crossing on `lambda` bits. The resulting total width is

```math
4\rho\lambda+O\!\left(\frac{\lambda}{\log(m+2)}+m\right).
                                                               \tag{9}
```

If either width exceeds `2lambda`, replace that cut by the identity
matrix times the complete piece. This cap cannot increase (9).

## 4. Enumeration and local workspace

There are at most `(b+1)^k` types. A count vector and a type ordinal
therefore use at most

```math
f_n(\lambda)=O(k_m\log(b+1))
 =O(n^{1/128}\log(n+2))                               \tag{10}
```

bits, by (4) and the preliminary label envelope. This is smaller than
`lambda/(log n)^delta` for every fixed `delta>0` at every active
label length, for sufficiently large `n`.

Compute multinomial counts as products of binomial coefficients.
Every partial product is at most `M_nu<=k^b`. The first-letter
counts `M_nu nu_i/b` permit exact word unranking. Decode the remaining
channel ordinal by mixed-radix division, using the raw dimensions
`r_i` or `c_i`. These divisions concern integer indices, not field
elements. Arbitrary non-power-of-two channel dimensions cause no change
to the factorization identity.

Equation (3) implies `b log_2 k=O(bm)=O(lambda)`. A word counter,
the multinomial arithmetic, and a channel counter therefore use
`O(lambda)` bits. The count vector (10) is smaller. A factor query
evaluates the selected raw factors sequentially, retaining only one
raw workspace and a constant number of coefficient accumulators.

A complete-piece query enumerates all words in `[k]^b`, keeps those
of the requested type, and sums their raw tensor-product entries.
There are

```math
k^b\le 2^{O(bm)}\le2^{O(\lambda)}                     \tag{11}
```

words. A raw piece entry is computed by its supplied factor sum,
which has at most `2^(2m)` terms. All raw work, arithmetic, and tape
indexing from (4) add at most polynomial factors in `n`. Since
`lambda>=n^(1/4)/8`, this still gives `2^O(lambda)` local bit time
and `O(lambda)` workspace. These queries also implement the cut caps.

Enumerating all count vectors takes `2^O(f_n(lambda))` time, within
the same bound. The complete compiled raw data are held once globally;
they are not duplicated at each recursion level.

## 5. Global cost and recovery

Because `log(m+2)=Theta(loglog(n+4))`, the relative error in (9) is
`O(1/loglog n)+O(m/n^(1/4))`. Across `O(loglog n)` levels its
product is bounded by a constant. The geometric volume envelope and
the main parent-query series in the global-block proof are unchanged.

The total piece-index entropy, summed over one worst recursion path,
is at most

```math
O\!\left(n^{3/4+1/128}\log(n+2)\right)
 =o\!\left(\frac n{(\log(n+2))^\alpha}\right).         \tag{12}
```

Indeed, `K_j^2=V_j/lambda_j<=8n^(-1/4)V_j`, and the geometrically
contracting volumes satisfy `sum_j V_j=O(n)`. Combine this with (10).
The same bound pays for stored ordinals and for their enumeration in
logarithmic running time. Width metadata are smaller.

The globally compiled data and `O(loglog n)` active scalar registers
cost only `O(n^(1/128)loglog(n+4))` bits. All remaining grid, top-query,
and recovery costs are those already proved in the global-block and
characteristic-independent transfers. This proves (2), including its
constant-factor endpoint and the stated field-specific output variants.

## Verification

The proof uses the existing field-compilation, tail, and global-budget
checks. The additional checker independently compares mixed-radix
type-channel unranking with explicit word-channel enumeration. It also
checks complete amplified tensors with non-power-of-two raw channels
and verifies the final-padding distinction.

```sh
python tools/check_rank_sum_global_transfer.py \
  --output evidence/rank-sum-global-transfer.json
```

These are finite interface checks. The asymptotic claim follows from
(1)--(12). This work adds no Lean formalization or improved numerical
factoring exponent.
