# Full coefficient-space optimization for auxiliary crossings

The [binary site-symmetry barrier](binary-site-symmetry-barrier.md) does
not apply to a general ordered product of noncommuting auxiliary
matrices. This experiment therefore tests actual recombinations of four
existing constructions, rather than extrapolating their Smith bounds.

The finite optimization finds strict improvements over raw monomial
pieces in several auxiliary examples. None supplies a post-half
rank-sum certificate or a uniform improved factoring construction.
The field here is exactly `F_2`; new combinations available over field
extensions are outside this enumeration.

## The finite optimization problem

For each existing local matrix polynomial `K`, form its ordered product
with the original fixed boundary:

```math
G_m(z)=e_0^T K_1(z)\cdots K_m(z)e_0
      =\sum_{j=0}^{d} z^j A_j,\qquad G_m(1)=E_m.
```

Every tensor in the coefficient span is
`A(w)=sum_j w_j A_j`, for a binary mask `w`. Enumerate **every** mask,
including coefficient dependencies and zero tensors, and compute its
two exact cut ranks. There are two different optimization questions:

1. Minimize the largest product `r(A(w))c(A(w))` in a decomposition
   of `E_m` into coefficient-span tensors.
2. Minimize the total cost `sum_w sqrt(r(A(w))c(A(w)))` in such a
   decomposition.

They need not select the same pieces. A low finite maximum can be
misleading when several pieces must subsequently be combined.

For the first question, sort masks by rank product and add them to a
binary linear span until it contains the all-ones coefficient mask.
The saved positive certificate lists masks whose XOR is that target.
A separating binary covector annihilates **every** mask below the
claimed threshold but evaluates to one on the target, proving the
matching lower bound. Zero-tensor masks are included from threshold
zero, so coefficient dependencies do not invalidate the test.

For the second question, view coefficient masks as vertices of a
Cayley graph. Adding a mask `w` has weight `sqrt(r(w)c(w))`. A path
from zero to the target is exactly a coefficient-span decomposition.
Use the two rational weights

```math
\frac{\lfloor 2^{20}\sqrt{r(w)c(w)}\rfloor}{2^{20}},\qquad
\frac{\lceil 2^{20}\sqrt{r(w)c(w)}\rceil}{2^{20}}.
```

Integer square roots calculate these bounds exactly. Shortest paths
give rational lower and upper bounds on the true optimum. In addition
to the upper path, the artifact saves an integer potential `phi` with
`phi(0)=0` and

```math
\phi(v+w)-\phi(v)
 \le\lfloor 2^{20}\sqrt{r(w)c(w)}\rfloor
```

for **every** pair of masks `v,w`. Telescoping these inequalities
certifies the lower bound without trusting or rerunning Dijkstra's
algorithm. The verifier checks all inequalities and recomputes all
rank pairs using the opposite pivot order.

## Complete finite results

The families are the [elementary pencil](linear-smith-barrier.md),
the [module-cancellation crossing](cancellation-module-crossing-example.md),
the [ordered-pairing crossing](block-smith-crossings.md), and the
mixed-loop construction implemented in
[`check_cancellation_module_automata.py`](../tools/check_cancellation_module_automata.py).
All use their original endpoint and boundary, rather than a free
auxiliary contraction.

| Family | Copies | Raw monomial maximum | Optimal recombined maximum | Optimal rank sum, approximately |
| --- | ---: | ---: | ---: | ---: |
| Elementary | 1 | 4 | 4 | 4 |
| Elementary | 2 | 36 | 36 | 13.797959 |
| Elementary | 3 | 224 | 224 | 44 |
| Elementary | 4 | 2116 | 2116 | 140.893146 |
| Module | 1 | 4 | 4 | 4 |
| Module | 2 | 36 | 36 | 13.797959 |
| Module | 3 | 476 | **224** | 44 |
| Module | 4 | 14336 | **11776** | 159.113723 |
| Ordered pairing | 1 | 4 | 4 | 4 |
| Ordered pairing | 2 | 128 | **96** | 13.797959 |
| Ordered pairing | 3 | 896 | **784** | 50.332021 |
| Ordered pairing | 4 | 8192 | **6720** | 191.271219 |
| Mixed loop | 1 | 6 | **4** | 4 |
| Mixed loop | 2 | 70 | **50** | 13.797959 |
| Mixed loop | 3 | 1176 | **1078** | 50.332021 |

The decimal column only summarizes the saved exact rational intervals.
The complete run covers 15 cases, 2,588 coefficient functionals and
214,080 endpoint entries. The three-copy module construction can match
the elementary optimal maximum by recombination, despite its larger
raw monomial maximum. At four copies its complete coefficient space
has a strictly worse optimum than the elementary space. The result is
finite; no all-length comparison follows from this table.

## Why a favorable finite maximum is insufficient

For four elementary copies, the optimal maximum `2116` even satisfies
`2116^5<2^56`, the finite `rho<7/10` comparison. Yet its optimal rank
sum is greater than `140`, which is too large for the simple repeated
finite-block envelope below. Thus this apparent small-block gain
does not establish a uniform `rho<7/10` family.

Given a fixed finite decomposition `E=sum_(j=1)^s B_j`, take the
multivariate tensor-power coefficients of `sum_j z_j B_j` on `t`
blocks. There are only `binom(t+s-1,s-1)` types, polynomial in `t`.
For a fixed type `(t_1,...,t_s)`, rank subadditivity in both cuts gives

```math
\sqrt{r(C_{\boldsymbol t})c(C_{\boldsymbol t})}
 \le \binom{t}{t_1,\ldots,t_s}
       \prod_j \sqrt{r(B_j)c(B_j)}^{t_j}
 \le \left(\sum_j\sqrt{r(B_j)c(B_j)}\right)^t.
```

For block word length `m`, this envelope would need rank sum below
`2^(sqrt(2)m)` to beat the half-rate contraction. All 15 certified
rank-sum lower bounds exceed even `2^(3m/2)`, a larger threshold.
For example, the four-copy elementary lower bound exceeds `140>64`.
This rules out improvement through **this** finite-block rank-sum
envelope in the enumerated spaces. More refined cancellations after
tensoring, larger blocks, field extensions, new boundaries, and other
auxiliary constructions remain open.

There is also no automatic passage from an `F_2` crossing to an
integer factoring machine. Such a passage needs its own exact-count
or modular reconstruction argument and complete cost analysis.

## Reproduction and evidence

```sh
make verify-auxiliary-coefficient-span
python3 tools/search_auxiliary_coefficient_span.py \
  --output /tmp/auxiliary-coefficient-span-replay.json
cmp evidence/auxiliary-coefficient-span.json /tmp/auxiliary-coefficient-span-replay.json
```

The first command verifies every saved rank pair, both minmax
certificates, and the rank-sum dual potentials and positive paths.
The second reruns the full optimization and writes source hashes.
The [artifact](../evidence/auxiliary-coefficient-span.json) retains every
functional's rank pair and the exact finite scope. The source module's
independent literal coefficient constructor was also compared through
two copies during generation.

Implementation and analysis used substantial AI assistance from Codex.
These are finite exact computations and ordinary arguments, with no
new Lean proof or stronger factoring exponent.

## Extension-field control

A second complete enumeration uses `F_4=F_2[omega]/(omega^2+omega+1)`.
It covers elementary copies one through four, module and ordered-pairing
copies one and two, and the one-copy mixed loop: nine cases and
**3,792** functionals. Every minmax threshold is identical to its
`F_2` counterpart, and every certified rank-sum interval is identical
at the `2^-20` rounding resolution. The latter does not assert that
the two exact irrational optima are equal.

The packed digit basis is `(omega,omega^2)`, so digit three represents
one and XOR still implements vector addition. The binary-span separator
is valid for the full extension-field search because every scalar
multiple of every generator is included. The Cayley-graph certificate
likewise permits all extension coefficients, including combinations
that descend to `F_2` only after summation.

Two-bit-plane elimination computes the full ranks. A separate dense
implementation in basis `(1,omega)` checks 64 random rectangular
matrices and the selected witness ranks; all binary coefficient
embeddings are also compared. Verification of the saved artifact
uses the opposite pivot order for every functional. Reproduce with

```sh
make verify-auxiliary-extension-span
python3 tools/search_auxiliary_extension_span.py \
  --output /tmp/auxiliary-extension-span-replay.json
```

The [extension-field artifact](../evidence/auxiliary-extension-span.json)
contains every rank pair, threshold separator, rank-sum dual potential,
and upper decomposition. This control closes a small field-extension
loophole in those nine finite cases; it makes no claim over the
algebraic closure, other fields, or untested tensor lengths.
