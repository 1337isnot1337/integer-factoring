# Paired path encoders beyond the classical spin formulas

The [coherent spin construction](bounded-spin-crossing.md) uses a
particularly structured primal frame and its norm dual. The crossing
endpoint needs less: any uniformly generated paired frame suffices.
This note gives a constant-state path construction of such frames over
`F16`, with Fibonacci label growth and polynomial-time entries.
It supplies a new family to study, not an improved factoring score.

The starting split-projection identity is already recorded in
[physical split crossings](physical-split-crossings.md). The contribution
here is an all-length construction with small local maps, rather than a
search over unrelated finite-dimensional projections.

## 1. Arbitrary invertible local splittings

Let `H_h` have dimension `h+1`, and let `W` have dimension two.
For each `h=0,1,2,3`, choose an invertible linear map

```math
U_h:H_{h+1}\oplus H_{h-1}\longrightarrow H_h\otimes W,
```

where `H_-1=0`. The dimensions are `2(h+1)` on both sides.
The up summand of `U_3` is defined for invertibility, but is discarded
when forming retained paths. Retain paths starting at height zero,
staying in `0..3`, and ending at zero.

Use the columns of `U_h` for primal extension and the corresponding
rows of `U_h^-1` for dual extension. Start with the one-dimensional
identity pair. For every retained prefix path `P` ending at height `h`,
this gives vectors `F_(P,r)` and covectors `D_(P,r)`, `0<=r<=h`.
They satisfy

```math
D_{(P,r)}F_{(Q,s)}=[P=Q][r=s].                    \tag{1}
```

Proof is induction on the prefix length. Different parent paths are
orthogonal by the preceding instance of (1). For equal parent paths,
the new pairing is the relevant entry of `U_h^-1 U_h=I`.
Discarding a branch selects a subset of the identity rows and columns,
so it preserves the assertion. No invariant form, orthogonal projection,
or classical `SL2` invariance is needed.

At final height zero, write the frames as label-by-word matrices
`F_n,D_n`. Then `F_n D_n^T=I`. There are
`L_(2r)=Fibonacci_(2r-1)` retained labels. The existing
[suffix-count unranking](../tools/bounded_spin_crossing.py) therefore
selects exactly `2^lambda` labels at a least even length `n=O(lambda)`.
Selecting that common subset of frame rows preserves (1).

For general fixed local splittings, a path's word entry is a product
of matrices of dimension at most four. This gives `O(n)` operations
per encoder entry. A four-leg polynomial contraction has at most
`4^4` internal magnetic states and degree at most `n`; dynamic
programming computes its entries in polynomial time.

## 2. An explicit F16 specialization

Use `F16=F2[a]/(a^4+a+1)`. For nonzero `q` with `q!=q^-1`, put

```math
[k]_q=(q^k-q^{-k})/(q-q^{-1}).
```

For `1<=m<=h`, the weight sector has rows
`(parent magnetic m,bit 0)` and `(parent magnetic m-1,bit 1)`.
Its columns are `(up,magnetic m)` and `(down,magnetic m-1)`.
Choose the block

```math
U_{h,m}=
\begin{pmatrix}
1&-[m]_q\\
q^{-h+m-1}&q^m[h-m+1]_q
\end{pmatrix}.                                  \tag{2}
```

Both boundary sectors have only an up column, with scalar one.
The determinant of (2) is

```math
q^m[h-m+1]_q+q^{-h+m-1}[m]_q=[h+1]_q,             \tag{3}
```

by direct substitution and cancellation. Thus all retained splittings
are invertible whenever `[2]_q,[3]_q,[4]_q` are nonzero.
Their inverse blocks are

```math
\frac1{[h+1]_q}
\begin{pmatrix}
q^m[h-m+1]_q&[m]_q\\
-q^{-h+m-1}&1
\end{pmatrix}.
```

There are twelve such parameters in this fixed field: the elements of
orders five and fifteen. The default `q=a^3` has order five. No claim
about a quantum group or a categorical trace is required for this
construction; equations (1)--(3) suffice.

The weight sectors make each word's magnetic state deterministic:
`r_new=r_old+bit-[down]`. A final height-zero word has exactly `n/2`
ones. This property also simplifies the four-leg contraction below.

## 3. Polynomial crossing and transfer interface

Let `K(z)=(1-z)I+zS` on `W tensor W`, where `S` swaps the two factors.
Compress `K(z)^(tensor n)` with primal `F_n` on `X,Y` and dual `D_n`
on `U,V`. Equation (1) gives both endpoints exactly:

```math
G_n(0;P,Q,R,S)=[P=R][Q=S],\qquad
G_n(1;P,Q,R,S)=[P=S][Q=R].                        \tag{4}
```

Only balanced binary words contribute. For fixed input and output
words, let `s` be the number of mixed positions and `b` the number
of swaps among them. Both are even, since each word has weight `n/2`.
Its polynomial is `z^b(1-z)^(s-b)`. In characteristic two this is
a polynomial in `z^2`. Consequently the coefficient span has the
uniform upper bound `n/2+1`.

The implemented entry algorithm uses the finite magnetic states,
without constructing the ambient tensor. The
[independent checker](../tools/check_paired_path_crossing.py) instead
groups complete balanced word pairs by their common and different
positions for its small controls. It compares the resulting polynomials
coefficient by coefficient with the magnetic algorithm.

As in the [spin-family modular interface](bounded-spin-modular-interface.md),
the binary-selected length is `O(lambda)` and the number of coefficients
is `O(lambda)`. Enumerating all fixed-field coefficient combinations
and factoring their two dense cut matrices costs `2^O(lambda)`.
The same permitted selection cost therefore applies here. If a strict
uniform rank gain is proved over `F16`,
[whole-tensor descent](fixed-field-descent.md) to `F2` costs only the
fixed extension degree four per cut after amplification. The existing
[modular factoring transfer](scale-one-modular-transfer.md) then supplies
the factoring interface. Projecting the local tensor at every site
is not the descent operation in that assertion.

## 4. Enlarged span and uniform support factors

The balanced six-vertex construction also applies to these paired
encoders. With local weights `a,f,b,e,c,d` on `0000,1111,0101,1010,0110,1001`,
equal word weights force equal occurrence counts within each pair.
Thus the compressed tensor is

```math
\sum_{A+B+C=n/2}T_{ABC}U^AV^BW^C,\qquad U=af,
\quad V=be,\quad W=cd.
```

The [same disjoint-window exposure](bounded-spin-six-vertex.md),
`a=f=d=1,b=e=1-z,c=z^(n+1)`, gives a univariate polynomial of degree
at most `n(n+1)/2`, with the crossing endpoint at one. Its coefficients
span every `T_ABC`: windows for different `C` are disjoint, and
the remaining polynomials `(1-z)^(2B)` have distinct leading degrees.
This argument remains valid in characteristic two. The magnetic
entry algorithm uses `O(n^2)` monomial counters and polynomial time.

There is a support bound for **every** weight-preserving local splitting
in section 1, on primal and dual words alike. At prefix length `i`,
write `D_i` for the number of down steps and `W_i` for the number of
word ones. Then `h_i=i-2D_i` and `r_i=W_i-D_i`, so

```math
|i-2W_i|=|h_i-2r_i|\le h_i\le3.
```

The [bounded-walk count](bounded-spin-support-ranks.md) therefore bounds
the physical support by `s_n=2a_(n/2)=O(g^n)`, where
`g=sqrt(2+sqrt(2))<phi^(4/3)` and `L_n` grows as `phi^n`.

In fact each `T_ABC` has explicit factors in both cuts with at most

```math
s_n\binom n{2(B+C)}\binom{2(B+C)}{2C}
 =s_n\frac{n!}{(2A)!(2B)!(2C)!}                 \tag{5}
```

channels. A channel chooses a supported word `x`, a set `J` of
`2(B+C)` positions, and a subset `D` of `J` of size `2C`. Its four
words are `x,x^J,x^D,x^(J xor D)`. The channel contributes

```math
F_P(x)F_Q(x^J)D_R(x^D)D_S(x^{J\mathbin\triangle D}).
```

This is a product of a row factor and a column factor in either cut.
All six-vertex configurations occur exactly once. Invalid or unsupported
words have zero encoder entries. The channels are uniformly unranked
by the existing bounded-walk and subset routines, including each word
and its complement as separate channels. No complement covariance of
the frame is assumed.

In particular the three extreme tensors each have both ranks at most
`s_n=O(L_n^(4/3))`. Unlike the classical family, they need not coincide:
an eight-site entry minor for `q=a^3` proves that all three are
independent. These are several uniformly cheap directions. They do
not yet decompose the crossing. By rank subadditivity and (5), even
all three corner bands of width `o(n/log n)` have total row rank at
most `poly(n)g^n exp(o(n))`, smaller than the crossing's `L_n^2`.

The enlarged coefficient dimension is `O(lambda^2)`. The exponential
selection argument of section 3 would now cost `2^O(lambda^2)` and
does not supply the required effectivity. A bound on the raw pieces
or another explicit selection recipe is an additional obligation
when using arbitrary enlarged-space recombinations.

## 5. Evidence and the remaining question

Run `make verify-paired-paths`; it uses only the Python standard library.
The [source-pinned receipt](../evidence/paired-path-crossing.json) records:

- All 72 interior local inverse blocks for the twelve parameters.
- All 2,388 final frame pairings at lengths two, four, six and eight.
- 240 complete polynomial entries compared with the independent word kernel,
  including the two endpoints and vanishing odd coefficients.
- Entry certificates attaining coefficient-span dimensions five, six and
  seven at lengths eight, ten and twelve for `q=a^3`. The parity bound
  above supplies the matching upper bounds.
- Complete support-channel factorizations for the three extremes and two
  first-corner pieces through six sites: 3,208 entry identities.
- The three independent eight-site extreme tensors and full exposure-map
  ranks through twelve sites, plus enlarged crossing endpoint controls.

The open problem is an effective decomposition of the crossing into
coefficient-span pieces with favorable simultaneous cut ranks, uniformly
as the selected label width grows. Arbitrary invertible splittings from
section 1 remain available if the specialization (2) is too rigid.
The [universal finite rank-sum exclusion](geometric-cut-rank-budget.md)
already rules out a strict seed through 95 labels, over every field.
Thus the 13-, 34- and 89-label checks here verify structure and entry
algorithms; they cannot provide a strict finite seed. A binary-label
seed would need at least 128 labels.

There is also a useful necessary condition for the uniform problem.
Suppose `E_L=sum_(j=1)^m H_j`, with `m=poly(n)` and
`r(H_j)c(H_j)<=L^(2 gamma)` for every piece. Rank subadditivity
forces some piece to have `r>=L^2/m`, hence `c<=m L^(2 gamma-2)`;
another piece has the reversed inequalities. At the post-half
threshold `gamma=sqrt(2)`, the small cut must therefore grow at
most as `L^(2 sqrt(2)-2+o(1))`, whose exponent is about 0.828.
Several pieces with both ranks subquadratic do not replace these
asymmetric directions. Classifying low-opposite-cut directions in
the enlarged span, together with their compatible combinations,
is a more specific next question than finding further cheap corners.

The construction and finite span certificates establish neither such
a rank theorem nor a higher score. The strongest complete score remains
`1/2`. These arguments and tools used substantial Codex AI assistance;
no new Lean work or independent acceptance is claimed.
