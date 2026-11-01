# Product-subspace automata and joint scalarization

## A finite refinement of the channel-walk construction

Assume the exact matrix-valued crossing identity and scalar cross-half
relations from `twisted-channel-walk-transfer.md`. All data belong to a
fixed finite field `F` and a fixed matrix algebra `Mat_b(F)`. Choose one
field-linear functional `ell` with `ell(I_b)=1`, used for every piece
and every word of a given amplification. A matrix entry is a permissible
choice, but the construction also permits other linear functionals.

For a row channel `s`, let `U_s` be the linear span of all matrices
appearing in its left factor table and let `V_s` be the span of all
matrices in its right factor table. External inputs at different
positions are independent. Thus the spans of the two ordered products
for a word `s_1...s_k` are exactly

```math
U(s_1...s_k)=\operatorname{span}(U_{s_1}\cdots U_{s_k}),\qquad
V(s_1...s_k)=\operatorname{span}(V_{s_1}\cdots V_{s_k}). \tag{1}
```

Its gathered scalar term is identically zero if and only if

```math
\ell(UV)=0\quad\hbox{for all }U\in U(s_1...s_k),\
                              V\in V(s_1...s_k).     \tag{2}
```

The channel-dependent commutation phase is nonzero and independent
of the external inputs, so it does not change this criterion.

## 1. A finite deterministic automaton

Start with the pair `(span I, span I)`. A channel `s` updates a pair
of subspaces by

```math
(U,V)\longmapsto
  (\operatorname{span}(U U_s),\operatorname{span}(V V_s)). \tag{3}
```

There are finitely many subspaces of `Mat_b(F)`, so the reachable pairs
form a finite deterministic automaton. A crude bound is
`|F|^O(b^4)` states. This may be too large for a practical search; it is
a finite certificate bound, not a claim of efficient construction as
`b` grows.

A smaller supplied automaton can use subspace overapproximations.
Each state `q` carries two subspaces `U_q,V_q`. The initial spaces
must contain `I`. For every channel transition `q --s--> q'`, verify

```math
U_q U_s\subseteq U_{q'},\qquad V_q V_s\subseteq V_{q'}. \tag{4}
```

The transition table is deterministic and includes every channel
from every retained state, except transitions sent to an explicit
discard state whose continuations are certified to contribute zero.
Equivalently one may keep a complete transition table first and
remove states from which no nonzero terminal contribution can ever
be reached. All containments and zero claims are checked on finite
bases.

A zero terminal value at a prefix does not by itself permit removing
that state. For example, with `ell(T)=T_11`, one has
`ell(E_12)=0` but `ell(E_12 E_21)=1`. The state after `E_12`
must be retained because a later channel can make it contribute.
The removal rule concerns all future continuations.

The exact reachable-pair construction needs no approximation and
automatically satisfies (4). The column factors have a separate
automaton of the same kind. Their alphabets still record the common
piece index, even though their channel indices differ.

## 2. Pay the actual terminal bilinear rank

For each row state `q`, choose bases of `U_q,V_q` and form the matrix

```math
H_q=(\ell(U_a V_c))_{a,c},\qquad t_q=\operatorname{rank}_F H_q. \tag{5}
```

The scalar function produced by a channel word ending at `q` has a
factorization through at most `t_q` ordinary field channels. Indeed,
express each product in the two state-space bases and apply any
rank factorization of `H_q`. The channel-word phase is absorbed into
the first factor. In particular, `t_q=0` certifies (2).

For the exact product-space automaton, this is the exact row rank
of an individual scalar channel-word function. The evaluations of
the left product span `U_q`, and those of the right product span
`V_q`, independently. Their coordinate matrices therefore have
full column rank; multiplying `H_q` by those matrices preserves
its rank. For an overapproximating state, only the upper bound is
asserted.

For `ell(T)=phi T psi`, the form `ell(UV)` has rank at most `b` on
all of `Mat_b` and hence on these subspaces. Its actual restricted
rank may be smaller. For a general linear functional the bound
`b^2` is sufficient; it is still a fixed certificate constant.

Let `R_nu(k)` now be the sum of terminal weights `t_q` over row
channel words of length `k` and piece type `nu`. Define `C_nu(k)`
with the independently computed column terminal ranks. The actual
scalar type pieces satisfy

```math
S_k\le\sum_{|\nu|=k}\sqrt{R_\nu(k)C_\nu(k)}.       \tag{6}
```

There is no additional global factor `b` in (6); the exact local
bilinear ranks have already been charged.

## 3. The same finite rational spectral certificate

For positive piece weights `z_i`, make a nonnegative transition
matrix by summing `z_i` for each channel of piece `i` that effects a
given state transition. Distinct channels are counted separately.
Let `e_0` be the initial-state vector and let `t` be the vector of
terminal ranks. The weighted count is

```math
P_R(k,z)=e_0^T T_R(z)^k t_R
        =\sum_{|\nu|=k}R_\nu(k)z^\nu.              \tag{7}
```

Use reciprocal weights in the column automaton. Cauchy--Schwarz gives

```math
S_k\le\sqrt{P_R(k,z)P_C(k,z^{-1})}.                \tag{8}
```

If positive rational vectors satisfy
`T_R(z)h_R<=beta_R h_R` and
`T_C(z^-1)h_C<=beta_C h_C`, then

```math
P_R(k,z)\le (h_R)_0
     \max_q\frac{(t_R)_q}{(h_R)_q}\,\beta_R^k,     \tag{9}
```

and similarly for the column count. Thus the asymptotic rank-sum
base is at most `sqrt(beta_R beta_C)`, with all constants and
inequalities finitely checkable. The effective all-width construction
and the score transfer are unchanged.

This provides a precise sense in which scalarization should be
optimized on surviving products. The functional determines the
terminal ranks and may make an entire component unable to reach
any nonzero terminal state. Removing such a component can lower
the certified exponential growth, as well as the fixed prefactor.
Minimizing a global rank bound of `b` alone cannot change the
asymptotic base.

## 4. Scope and unresolved construction

The automata are derived from actual factor tables and exact field
arithmetic. An arbitrary accepted language, or a zero tested at only
one external input, is insufficient. Any supplied transition
overapproximation must satisfy (4), and every excluded terminal
state must have the zero bilinear form (5).

For a fixed finite field and matrix order, all normalized linear
functionals can in principle be enumerated. The product-space
transitions are independent of the functional, so their spaces can
be reused while terminal ranks are recomputed. This is a search
reduction; its crude enumeration cost is large.

The construction strictly extends adjacent zero-product certificates:
it detects arbitrary-length zero products and products annihilated
only after the common scalarization. It also charges exactly the
restricted bilinear ranks instead of an ambient bound. It does not
yet supply a crossing identity whose certified base improves the
known factoring score.

## 5. Exact finite check

`check_scalarized_channel_automata.py` constructs the complete reachable
product-space automaton for the three square-zero matrices in the
channel-walk note. There are eleven states over `Mat_2(F_2)`. It checks
all eight normalized linear functionals on that matrix algebra.
For six functionals, including `ell(T)=T_11`, only seven states can
reach a nonzero terminal form. The remaining two functionals retain
ten such states. Thus scalarization really does remove entire
reachable components in this example.

Every surviving terminal bilinear form has rank one, so no ambient
factor two is needed. For `ell(T)=T_11`, the exact weighted count at
length `k` is

```math
4^k\,\frac{2^{k+1}+(-1)^k}{3}.                    \tag{10}
```

The checker verifies 132 transition-space identities, 88 terminal
rank computations, 312 comparisons with explicit matrix words, and
62 exact rational supervector inequalities. It also checks (10)
through `k=6`. The recorded run in `scalarized-channel-automata.json`
completed in about two hundredths of a second.
The `E_12 E_21` control explicitly checks that a zero terminal prefix
remains in the live automaton when it has a nonzero continuation.

The exponential base remains eight. Direct scalarization of the
single-copy toy already leaves the ordinary crossing with rank sum
four, so this example supplies no competitive construction. It
checks the product-space reduction, the distinction between
normalized functionals, and the terminal-rank accounting.
