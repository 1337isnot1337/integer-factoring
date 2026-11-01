# Counting Smith channel words by a finite min-plus automaton

## Result and scope

The degree-sensitive auxiliary-path bound can count each channel word
once. A finite deterministic automaton records its least possible
auxiliary degree cost. This removes multiplicity caused by several paths
realizing the same channel word.

For auxiliary dimension `D`, cap the certified entry-order lower bounds
at a fixed integer `M>=0`, assigning an exact zero entry the finite lower
bound `M` as well. The resulting automaton has at most

```math
(M+1)^D-M^D                                                   \tag{1}
```

states after the first symbol. Its weighted transitions have rational
entries at rational `s>=1`. Two positive rational potential certificates
therefore give the same effective crossing-to-factoring criterion as
the path theorem, with potentially smaller exponential bases.

All parameters in (1) belong to the fixed local certificate. The state
count is not claimed polynomial in a growing auxiliary dimension. No
improving seed or factoring exponent is asserted here.

## 1. The existing ordered factorization

Use the data, regularized boundaries, and exact coefficient family from
[degree-sensitive-smith-paths.md](degree-sensitive-smith-paths.md).
For one cut write the local channels in the uniform form

```math
\sum_{i=1}^{b}\tau^{-e_i}U_i(\tau;\text{left})
                                      v_i(\tau;\text{right}).   \tag{2}
```

At infinity `e_i=d_i`; at zero `e_i=-nu_i`. The outside factors are
regular series. Let `F` and `E` be nonempty supports containing those
of the regular left and right boundary vectors.

Choose certified nonnegative lower bounds for the entry orders of
`U_i`, including a scalar-factor order if desired. Replace every bound
by its minimum with `M`; replace infinity by `M`. Denote the resulting
integers by `kappa_(i,a,b)`. Every entry is now finite and between zero
and `M`, and remains a valid lower bound. Lowering a degree cost can
admit extra channel words, which is safe for an upper bound.

For a word `i=(i_1,...,i_t)`, define its capped minimum path cost

```math
m(i)=\min_{a_0\in F,\ a_t\in E}
       \sum_{j=1}^t\kappa_{i_j,a_{j-1},a_j}.                    \tag{3}
```

The word's exact auxiliary matrix product can have higher order because
of larger true entry orders or cancellation, but never lower order
than (3).

## 2. A bounded normalized cost vector

Start with a cost vector whose coordinates in `F` are zero and all
others are infinity. Reading channel `i` updates a vector `c` by

```math
\widehat c_b=\min_a(c_a+\kappa_{i,a,b}).                         \tag{4}
```

After one symbol every coordinate is finite. If `m=min_a c_a`, each
updated coordinate belongs to `[m,m+M]`: the lower bound follows from
nonnegative edge costs, and an index attaining `m` proves the upper
bound. Thus the normalized vector

```math
\eta=c-(\min_a c_a)\mathbf1
```

has entries in `{0,...,M}` and at least one zero. This proves (1).

For a normalized state `eta` and symbol `i`, compute

```math
g_b=\min_a(\eta_a+\kappa_{i,a,b}),\qquad
\gamma_i(\eta)=\min_b g_b,\qquad
T_i(\eta)=g-\gamma_i(\eta)\mathbf1.                            \tag{5}
```

Each symbol has exactly one next state. The removed minima add along
the word. If `eta_t` is its terminal state, then

```math
m(i)=\sum_{j=1}^t\gamma_j+\min_{b\in E}(\eta_t)_b.              \tag{6}
```

The first update uses the initial zero/infinity vector instead of an
ordinary normalized state. It is finite and contributes only a fixed
initial vector below. States unreachable from this first update may
be omitted.

## 3. Count weighted words exactly

Fix rational `s>=1`. On the reachable normalized states define

```math
\mathcal T(s)_{\eta,\eta'}
 =\sum_{i:T_i(\eta)=\eta'}s^{e_i-\gamma_i(\eta)}.               \tag{7}
```

Let `b(s)` sum the same weights from the first update into each state.
Define the terminal vector

```math
r(s)_\eta=s^{-\min_{a\in E}\eta_a}.                            \tag{8}
```

Equation (6) gives, for every `t>=1`, the exact identity

```math
Q(t,s)=b(s)^T\mathcal T(s)^{t-1}r(s)
 =\sum_{i\in[b]^t}s^{\sum_j e_{i_j}-m(i)}.                     \tag{9}
```

There is one summand per channel word. No auxiliary-path multiplicity
appears. Dropping the nonnegative terminal cost is also safe, but (8)
retains it without enlarging the automaton.

For comparison, the corresponding capped nondeterministic path matrix
counts every path with weight `s^(sum e-sum kappa)`. Its sum for a word
is at least the largest such path weight, which is the word's weight
in (9). Therefore the deterministic weighted count never exceeds that
of the capped path construction. This comparison uses the same finite
costs; capping a previously infinite cost can introduce paths absent
from an uncapped graph.

## 4. Coefficient ranks and finite rational certificates

For the infinity cut, a word can contribute to coefficient `k` only
if `sum d_i-m(i)>=k-H_infinity`. For the zero cut, it can contribute
only if `sum(-nu_i)-m(i)>=-k-H_zero`. The ordinary full-word
coefficient-splitting proof now gives

```math
\begin{aligned}
r(A_{t,k})&\le C_\infty(t+1)s^{H_\infty-k}Q_\infty(t,s),\\
c(A_{t,k})&\le C_0(t+1)s^{H_0+k}Q_0(t,s).
\end{aligned}                                                \tag{10}
```

Unlike the path expansion, the first factor here keeps the entire
ordered auxiliary matrix product for each channel word. One global
coefficient split still supplies one rank-one outer product. Discarding
words that fail the capped minimum-cost test is valid by (3).

For either automaton, supply a strictly positive rational vector `h`
and a positive rational `beta` with

```math
\mathcal T(s)h\le\beta h.                                     \tag{11}
```

Then

```math
Q(t,s)\le
 (b(s)^T h)\max_\eta\frac{r(s)_\eta}{h_\eta}\,
 \beta^{t-1}.                                                 \tag{12}
```

This is an exact rational certificate, with a fixed boundary constant.
Combining the two cuts cancels `k`. If
`B_det=beta_infinity beta_zero`, every coefficient has

```math
r(A_{t,k})c(A_{t,k})\le C(t+1)^2 B_{\rm det}^{t}.              \tag{13}
```

The polynomial number of coefficient pieces gives the same bound for
the rank sum up to a polynomial factor. Thus

```math
\rho=\frac{\log B_{\rm det}}{4\log q}
```

is the resulting leading width parameter. The exact test
`B_det^5<q^14` is sufficient for an exponent above one half. One may
use the ordinary Smith bound or the earlier path bound in either cut
when it is smaller.

The automata have fixed size. Entries and full-word factors still take
`2^O(t)` bit operations to construct over the fixed field. All ordinary
effectiveness and recurrence requirements of the existing transfer
remain in force. Matrix entries, state counts, and coefficient queries
are not treated as free input-dependent oracles.

## 5. Controls and reproducibility

At `M=0` there is one state, every increment is zero, and (7) becomes
the original independent channel sum. Increasing a valid cap can retain
more degree information. A fixed finite cap also makes determinization
safe when the exact zero graph would have unbounded normalized cost
differences.

The checker independently enumerates short channel words and auxiliary
paths, verifies (6)--(9), and compares exact coefficient ranks of actual
polynomial crossing families against (10)--(13). All potential
inequalities use rational arithmetic. Its finite examples support the
implementation; they do not establish an improving infinite family.

```sh
python3 tools/check_tropical_smith_channel_words.py \
  --output evidence/tropical-smith-channel-words.json
```

The new step is bounded min-plus determinization inside the Smith
coefficient budget. It uses the fixed local factors and ordered-product
framework of the linked notes. This work used substantial OpenAI
assistant support.
