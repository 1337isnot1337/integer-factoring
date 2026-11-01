# Finite-state polynomial crossing families

## Status

This note gives two explicit families of polynomial crossing identities that
can depend on the order of the bit positions. Their auxiliary dimensions are
constant, their coefficient counts are polynomial, and their endpoint is the
exact crossing tensor. It also gives an exact rank algorithm exploiting the
characteristic-two parity support of the local tensors.

The accompanying screen covers **88,006 finite models**. Selected coupled
models have been evaluated through twelve bit positions. These computations
do **not** establish an improved factoring exponent. In particular, a small
maximum rank product at a short length is not an asymptotic contraction
certificate. The two metrics below, their Pareto frontiers, and their growth
increments are retained separately.

The source files are `markov_pairing_search.cpp` and
`markov_pairing_profile.cpp`. Independent checks are in
`check_markov_pairing_transfer.py`; the completed evidence is
`markov-pairing-verification.json`. The compact screen record is
`finite-state-pairing-screen.json`, produced by
`audit_finite_state_screen.py` from the complete raw profiles.

## 1. Local tensors and the score interface

Work over a field of characteristic two, initially over `F_2`. For one-bit
labels define

\[
 R(x,y,u,v)=[x=y][u=v],\qquad
 C(x,y,u,v)=[x=u][y=v],\qquad
 E(x,y,u,v)=[x=v][y=u].
\]

Set

\[
 H=R+C+E=[x+y+u+v=0].
\]

Thus `R+H+C=E`. The equality with the parity indicator follows by checking
the eight even-parity inputs; all three pairing tensors vanish on odd parity.
For words of length `t`, the label dimension is `D=2^t`, and the target is
`E_t=E^{\otimes t}`.

For a polynomial tensor

\[
 F_t(z)=\sum_{d=0}^{L_t}z^d A_{t,d},\qquad F_t(1)=E_t,
\]

write `r_{t,d}` and `c_{t,d}` for the ranks in the cuts `(x,y)|(u,v)` and
`(x,u)|(y,v)`. Define

\[
 P_t=\max_d r_{t,d}c_{t,d},\qquad
 S_t=\sum_d\sqrt{r_{t,d}c_{t,d}},
\]

\[
 \rho_{\max}(t)=\frac{\log_2P_t}{4t},\qquad
 \rho_S(t)=\frac{\log_2 S_t}{2t}.
\]

If there are at most `M_t` nonzero coefficients, then

\[
 \sqrt{P_t}\le S_t\le M_t\sqrt{P_t},\qquad
 0\le\rho_S(t)-\rho_{\max}(t)
 \le\frac{\log_2 M_t}{2t}.                 \tag{1}
\]

Consequently these quantities have the same limiting upper exponent whenever
`M_t=t^{O(1)}`. Their finite values can differ substantially. A strict finite
seed certificate for a score above one half would require

\[
 S_t<2^{\sqrt2\,t}.
\]

Alternatively, a uniform family satisfying
`limsup rho_max(t) <= gamma < 1/sqrt(2)` and polynomial coefficient count
would enter the transfer theorem in `vanishing-slack-score-transfer.md` with
width `4 gamma t+o(t)`. Explicit coefficient construction and ordinary exact
rank factorization cost `2^{O(t)}` bit operations for the fixed families here.
Thus the construction overhead meets that theorem's interface. The missing
ingredient is a suitable asymptotic rank bound, not an endpoint identity or
an implicit exponential coefficient list.

## 2. A three-state nearest-neighbor family

Use the symbols `0,1,2` for `R,H,C`. Fix nonnegative integer vertex costs
`a_s` and transition costs `b_{rs}`. For a word `w` put

\[
 e(w)=\sum_{i=1}^t a_{w_i}
      +\sum_{i=1}^{t-1}b_{w_i,w_{i+1}}.
\]

Define

\[
 F_t(z)=\sum_{w\in\{R,H,C\}^t}
       z^{e(w)}\,w_1\otimes\cdots\otimes w_t.              \tag{2}
\]

**Endpoint and complexity.** At `z=1` every word has weight one, so

\[
 F_t(1)=(R+H+C)^{\otimes t}=E_t.
\]

The degree is at most
`t max_s a_s+(t-1) max_{r,s} b_rs`. For fixed costs this gives `O(t)` pieces.
Each entry is obtained by a three-state transfer recurrence, so the auxiliary
dimension is three. This is a direct polynomial-time computation in `t` for
each polynomial entry.

More generally, replace each edge monomial `z^{b_rs}` by a fixed polynomial
`q_rs(z)` with `q_rs(1)=1`. Both the endpoint and the linear degree bound
remain valid. The polynomial-edge screen uses the alphabet

\[
 \{1,z,1+z+z^2\}.
\]

**A genuine order-dependent example.** Give `R,H` the same cost and the same
edge behavior, and merge them into `S=R+H=E+C`. Take vertex costs
`a_S=0,a_C=5`, with transition cost one exactly when the symbol changes
between `S` and `C`. The entry at the local quadruple `(0,0,0,0)` selects `C`,
and the entry at `(0,1,1,0)` selects `S`. The entries selecting `CSS` and
`SCS` are respectively `z^6` and `z^7`. They are related by a permutation of
positions, so at least one coefficient is not position invariant. This
explicitly places the family outside the hypothesis of
`all-field-binary-word-barrier.md`.

Two further selected variants are:

* the run-cost family `a=(0,0,2)`, with cost one only on `C->C`;
* the count/run bigrading, where `a_C=ceil(t/2)+1`, `a_R=a_H=0`, and the
  only positive transition costs are `R->C` and `H->C`.

The latter encodes both the number of `C` symbols and the number of `S->C`
transitions. The transition count is at most `ceil(t/2)`, so the stated
stride makes this encoding injective. Its degree is `O(t^2)`, and (1) still
applies with an `O(log(t)/t)` difference between the two exponents.

## 3. A two-dimensional noncommuting auxiliary family

Fix three binary `2 by 2` matrices `B_R,B_H,B_C`. Put

\[
 T_s(z)=I+(z+1)B_s.
\]

Let `ell` be a linear functional on `2 by 2` matrices with `ell(I)=1`, and
let `a_R,a_H,a_C` be nonnegative integers. Define the tensor entry

\[
 F_t(z;q_1,\ldots,q_t)=
 \ell\!\left(\prod_{i=1}^t
     \left[\sum_{s\in\{R,H,C\}}z^{a_s}S_s(q_i)T_s(z)\right]\right).
                                                               \tag{3}
\]

The product is taken in position order. Since `T_s(1)=I`, each factor at
the endpoint equals `E(q_i)I`; hence (3) is exactly `E_t` at `z=1`.
Its auxiliary dimension is two and its degree is at most
`t(1+max_s a_s)`. Noncommuting `B_s` permit order dependence. Whether an
individual boundary functional actually detects that order dependence must
be checked; noncommutativity alone is not asserted to be sufficient.

The screen exhausts all `16^3=4096` matrix triples for each of the vertex
vectors `(0,0,1)` and `(0,1,2)`, and for both boundary functionals

\[
 \ell_0(M)=M_{00},\qquad
 \ell_1(M)=M_{01}+M_{10}+M_{11}.
\]

These represent all conjugacy classes of boundary functionals satisfying
`ell(I)=1` over `F_2`. Indeed write `ell(M)=Tr(WM)`. A binary matrix of
trace one has characteristic polynomial either `X(X+1)` or `X^2+X+1`,
and is conjugate respectively to `diag(1,0)` or
`[[0,1],[1,1]]`. Simultaneously conjugating the three arbitrary `B_s`
does not change the exhaustiveness of the screen.

The auxiliary family slightly lowers the best finite rank sum in this
screen: at `t=6`, `S` drops from about `1461.92113` to `1364.59717`.
The latter profile has a larger maximum rank product, `397488` instead of
the uncoupled profile's `132496`. Both remain far above the strict finite
seed threshold. This is a reason to retain both Pareto coordinates.

## 4. Exact XOR-circulant rank reduction

Every local tensor in the span of `R,H,C` vanishes on odd parity. For the
row cut, nonzero entries therefore have

\[
 g=x+y=u+v.
\]

Within such a block use `x,u` as free labels and put `h=x+u`. Addition here
is bitwise XOR. At position `i` the values of `R,H,C` are respectively

\[
 (1-g_i,\ 1,\ 1-h_i).
\]

Thus a coefficient of (2) or (3) is described by a binary kernel
`K_d(g,h)`, and its row block at `g` has entries

\[
 K_d(g,x+u).
\]

This is a `2^t by 2^t` XOR-circulant matrix. In the column cut the block
label is `h=x+u=y+v`, and its kernel is the column `K_d(*,h)`. Therefore

\[
 r_d=\sum_g\operatorname{rank}\operatorname{Circ}(K_d(g,*)),\qquad
 c_d=\sum_h\operatorname{rank}\operatorname{Circ}(K_d(*,h)).       \tag{4}
\]

Only `4^t` polynomial entries are generated, compared with `16^t` entries
in the original tensor. Formula (4) is an exact direct-sum decomposition;
it neither approximates ranks nor uses a numerical singular-value test.

### 4.1 Two exact shortcuts

For a binary kernel `(k_h)`, the circulant acts by multiplication by

\[
 f=\sum_h k_h\prod_{i:h_i=1}(1+\epsilon_i)
\]

in the group algebra

\[
 A_t=\mathbb F_2[\epsilon_1,\ldots,\epsilon_t]/(\epsilon_i^2).
\]

The constant coefficient of `f` is `sum_h k_h`. If it is one, `f` is a
unit in this local algebra and its multiplication map has rank `2^t`.
If the constant is zero but some linear coefficient is one, the rank is
exactly `2^{t-1}`. To see this, replace that variable by `f`, leaving the
other variables fixed. Every element of the maximal ideal squares to zero
in characteristic two, so this defines an algebra map. Its linear part is
invertible; the filtration by powers of the maximal ideal then shows the
map is an automorphism. It conjugates multiplication by that variable,
which has rank `2^{t-1}`, to multiplication by `f`.

The linear coefficient at `epsilon_i` is simply
`sum_{h:h_i=1} k_h`. All remaining kernels are ranked by packed exact
Gaussian elimination. The independent checker compares these shortcuts
with unshortened Gaussian elimination for **all 65,536 kernels at `t=4`**.

### 4.2 Independent endpoint and rank checks

`check_markov_pairing_transfer.py` builds the small and wide native engines
from their current source. It independently expands the original tensor
words in Python, without the parity-block reduction, and compares both full
matrix ranks for 24 cost, polynomial-edge, and auxiliary-matrix examples.
It checks the exact endpoint in each case, compares the native engines on
two six-position instances, and tests a degree-100 polynomial to exercise
the wide engine beyond the 64-bit polynomial range. It also verifies the
`CSS` versus `SCS` position-symmetry counterexample above.

The wide engine currently uses 128-bit polynomial masks and packed binary
kernels. Its default compile-time limit is twelve positions. A build with
`MAX_PROFILE_BITS=14` supports the fourteen-position experiment. These
implementation limits are explicit; they are not restrictions on (2),
(3), or the rank identity (4).

## 5. Completed finite screen

For cost screens, the vertex coordinates lie between zero and the stated
maximum, at least one is zero, and `a_R<=a_C`; all 512 binary transition
matrices are tested. Each polynomial-edge screen tests all `3^9` choices.

| Family | Length | Cases | Smallest integer lower bound for `S` |
|---|---:|---:|---:|
| Vertex costs through 2, binary transitions | 5 | 6,144 | 464 |
| Vertex costs through 5, binary transitions | 6 | 26,112 | 1,460 |
| Polynomial edges, vertex `(0,1,2)` | 6 | 19,683 | 2,104 |
| Polynomial edges, vertex `(0,0,1)` | 6 | 19,683 | 1,460 |
| All auxiliary matrices, four grading/boundary pairs | 6 | 16,384 | 1,363 |

The integer lower bound is `sum_d floor(sqrt(r_d c_d))`. Every displayed
minimum exceeds `2^{ceil(3t/2)}`, which exceeds `2^{sqrt(2)t}`. Thus no
candidate in the finite screen satisfies the strict finite seed criterion;
that conclusion uses only exact integer arithmetic. Floating-point square
roots are used to order search candidates, not to certify an improvement
or an exclusion at a close threshold.

The compact JSON retains each screen's two-coordinate Pareto representatives,
all their exact rank profiles, the raw-file hashes, and selected larger
profiles. It also records growth increments rather than extrapolating a
limit from an early low value.

| Coupled model | `t` | `rho_max(t)` | `rho_S(t)` |
|---|---:|---:|---:|
| Domain wall, `a_C=5` | 6 | 0.66945192 | 0.95276895 |
| Domain wall, `a_C=5` | 8 | 0.68120981 | 0.93740424 |
| Domain wall, `a_C=5` | 10 | 0.70969977 | 0.92465339 |
| Domain wall, `a_C=5` | 12 | 0.72739227 | 0.91415891 |
| `C->C` run cost, `a_C=2` | 6 | 0.68411014 | 0.92353504 |
| `C->C` run cost, `a_C=2` | 8 | 0.69866534 | 0.91193607 |
| `C->C` run cost, `a_C=2` | 10 | 0.71835289 | 0.90217749 |
| Count/run bigrading | 10 | 0.71406753 | 0.90174502 |
| Count/run bigrading | 12 | 0.72868877 | 0.89345779 |

For example, the domain-wall maximum product grows from `351266544` at
`t=10` to `32389991200` at `t=12`. Its normalized incremental exponent is

\[
 \frac{\log_2(32389991200/351266544)}{4(12-10)}
 \simeq0.81585476.
\]

This differs substantially from the cumulative value `0.72739227` and
shows why neither number should be promoted to an asymptotic theorem.
The data leave open the behavior of broader finite-state models and more
carefully chosen transition polynomials.

## 6. Reproduction

The independent checks require Python 3 and a C++17 compiler. From the
repository root, run:

```sh
make verify-search-native
```

To regenerate all eight raw screens, all selected larger profiles, and
the compact evidence, run:

```sh
python3 tools/reproduce_finite_state_screen.py --directory /tmp/integer-factoring-state-screen --output /tmp/integer-factoring-state-screen/finite-state-pairing-screen.json
```

The driver lists every exact grading, transition, boundary, and profile
parameter. It builds its executables and writes the large raw screens in
the requested temporary directory. Larger two-state profiles use the
independently checked [run-compressed engine](two-state-run-rank-compression.md).
Their complete integer rank lists are compared with the original engine
on the validation instances.

The compact [screen evidence](../evidence/finite-state-pairing-screen.json)
retains every reported Pareto profile and larger-length comparison.
The [independent verification](../evidence/markov-pairing-verification.json)
records the direct tensor and XOR-circulant comparisons. Raw-file hashes
identify the original screen runs; timings and those raw hashes can
change when a screen is regenerated. Exact rank lists and the stated
finite inequalities are the mathematical assertions to reproduce.
