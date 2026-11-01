# Reversal-stable binary alphabets and nilpotent braid dressings

The [F4 braid path family](quantum-braid-path-quotient.md) has an explicit
leading coefficient. Its two factors pair a path with its reversal.
Selecting an exact binary alphabet closed under that reversal makes
the normalized leading operator square-zero. This permits polynomial
unimodular conjugations preserving both strict endpoints.

The selection has polynomial-time unranking and uses the same least
even length as the earlier binary alphabet. Exact small conjugations
cancel leading degrees, but no favorable uniform rank-pressure theorem
or higher factoring score is established.

## 1. A word-level reflection identity

Write `P^rev=(-s_n,...,-s_1)` for reversal of a return path. For a
binary word `w`, let `Jw` reverse the bit order and complement every
bit. For each valid q-splitting from the
[paired-path construction](paired-path-crossings.md),

```math
F_P(Jw)=F_{P^{rev}}(w).                          \tag{1}
```

This holds directly on word entries. Define the nonzero potential

```math
\eta(h,r)=q^{(r+1)(h-r)}[r]_q![h-r]_q!.
```

If a valid step takes `(h,r)` to `(h+s,r')` on bit `b`, its reversed
step starts at `(h+s,r')` with direction `-s` and bit `1-b`.
The two local primal coefficients satisfy

```math
c(h+s,-s;r',1-b)
 =c(h,s;r,b)\frac{\eta(h+s,r')}{\eta(h,r)}.       \tag{2}
```

For an up step this follows from the two ratios
`q^(r+1)[h+1-r]_q` and `q^(h-r)[r+1]_q`.
A down step is the reversed identity. In characteristic two all
signs here are positive. The magnetic states are reversed as well:
the new reversed prefix state is exactly the old next state.
Thus invalid words occur in matched pairs, while the potentials in
every valid return path telescope from `eta(0,0)=1` back to one.
This proves (1), including zero entries.

For the order-five specialization, put

```math
w_P=\prod_{\text{down from }h}[h]_q[h+1]_q
 \in\mathbb F_4^*,\qquad \nu_P=q^{n/2}w_P.
```

Every up crossing of a fixed level is matched by a down crossing.
Reversal therefore preserves this product: `w_(P^rev)=w_P`.

## 2. The leading tensor in closed form

Let `A_j` be the unprojected F4 path coefficients and `D=n^2`, as
in the path-quotient note. The top coefficient is

```math
\boxed{\quad
A_D(P,Q,R,S)=\frac{w_P}{w_R}
 [Q=P^{rev}][S=R^{rev}].\quad}                  \tag{3}
```

Here is a physical-matrix proof of the factorization. In the block-braid
word, the degree-D term is the product of all local rank-one `e_i`.
Commuting disjoint gates orders the square grid by anti-diagonals.
Cell `(j,k)` has gate index `n+j-1-k`. Sorting by `j+k` retains
the order of every pair whose gate indices differ by at most one;
only disjoint gates are exchanged.
Its central layer is `e_0 e_2 ... e_(2n-2)`, a tensor product of
`n` cups and caps. Each triangular half turns the neighboring cups
into the nested rainbow. The local cup matrix is
`C=[[0,q],[1,0]]`, and `q^-1 C^2=I`; each reconnection in the
triangle has coefficient one. Induction through its layers gives

```math
\prod e_i=q^{-n}r_n r_n^T,
\qquad r_n(x,y)=q^{n-|x|}[y=Jx].                \tag{4}
```

All encoder words have weight `n/2`, so the factor in the rightmost
expression is `q^(n/2)`. Combining (1), the orthogonal norm duals,
and (4) gives input overlap `q^(n/2)nu_P[Q=P^rev]` and output
overlap `q^(n/2)nu_R^-1[S=R^rev]`. The factors `q` cancel, giving
(3). This is also a direct derivation of both actual cut factors.

The strict F4 coefficient is `B_D=b_D A_D`, where `b_D` is nonzero:
it is one when `5` divides `n`, and `alpha` otherwise.

Let `S` be any selected set of paths, and let
`k=|{P in S:P^rev in S}|`. Formula (3) gives row rank one and
column rank `k^2` when `k>0`. Normalize the row operator by exchanging
its two column labels, so the crossing endpoint becomes the identity.
Its leading operator `P_D` satisfies

```math
P_D^2=b_D k P_D.                               \tag{5}
```

The integer count is evaluated in characteristic two. In particular,
if `S` is closed under reversal and has even size `L`, the two ranks
are exactly `(1,L^2)` and `P_D^2=0`. An even-size set without that
closure need not have a square-zero leading operator.

## 3. Exact binary alphabets closed under reversal

At length `n=2m`, a path is uniquely `A (B^rev)`, where `A,B` are
half-paths of length `m` ending at the same height. Let `M_h` be the
number of those half-paths. Reversal exchanges `A,B`, so there are

```math
\sum_h\binom{M_h}{2}\text{ two-element orbits},
\qquad \sum_h M_h\text{ fixed paths},
\qquad L_n=\sum_h M_h^2.
```

Choose the least even `n>=4` with `L_n>=2^lambda`. Take as many whole
two-element orbits as fit, followed by the required fixed paths.
If there are too few pairs, the inequality `L_n>=2^lambda` guarantees
enough fixed paths to fill the remainder. The target is even, so
the remainder is even. The selected set has exactly `2^lambda`
elements and is closed under reversal.

The [unranking implementation](../tools/quantum_braid_reversal.py)
counts half-path suffixes, unpacks an unordered pair by triangular
counts, and uses lexicographic unranking for the fixed half-paths.
Its integer counts have `O(lambda)` bits and the operations take
polynomial time. Labels in a selected pair reverse by toggling one
bit; fixed labels reverse to themselves. The least length and its
`n=O(lambda)` bound are unchanged.

Applying the existing paired frames on this new alphabet preserves
both strict endpoints. Equation (5) supplies a nonzero rank-one,
square-zero leading row operator at every positive binary width.

## 4. Polynomial conjugation and a pole-cancellation criterion

Let `R(z)` be the normalized row operator of the selected crossing,
with degree `D` and leading coefficient `P=P_D`. For
`f(z)=a z(1+z)` over `F4`, set `U(z)=I+f(z)P`.
Equation (5) gives `U^-1=U`, and `U(0)=U(1)=I`.
Thus

```math
\widetilde R=URU
 =R+f(PR+RP)+f^2PRP                           \tag{6}
```

gives another exact crossing after undoing the column-label exchange.
Its degree is at most `D+4`, and dense coefficient construction still
costs `2^O(lambda)`.

Because `P` has rank one, write `PR(z)P=ell(z)P`. If
`deg(PR+RP)<=D-3`, `deg ell<=D-4`, and `ell_(D-4)!=0`, choose
`a^2=ell_(D-4)^-1`. The degree-D term in (6) cancels and no higher
term is introduced. This is a sufficient algebraic criterion, not
an assertion that every width satisfies it.

The exact four-site, two-label example with `a=1` reduces degree
16 to 14. The six-site, four-label example with `a=alpha+1` reduces
36 to 35. The eight-site example checked does not reduce its degree.
The new leading ranks can increase, and initial unweighted pressure
checks give no stronger factoring bound. Degree cancellation alone
does not settle the simultaneous two-cut cost.

Run `make verify-braid-reversal` for the NumPy replay. The
[receipt](../evidence/quantum-braid-reversal.json) checks the local
reflection potential, word reflection, exact binary orbit unranking,
complete leading tensors through eight sites, the square-zero property,
independent polynomial conjugation products and strict endpoints.
The first eight lexicographic eight-site labels give a negative control:
their leading operator is not square-zero. The complete factoring
score remains `1/2`. These arguments and tools used substantial Codex
AI assistance; no new Lean or independent acceptance is claimed.
