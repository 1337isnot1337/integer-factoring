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

## 4. Evidence and the remaining question

Run `make verify-paired-paths`; it uses only the Python standard library.
The [source-pinned receipt](../evidence/paired-path-crossing.json) records:

- All 72 interior local inverse blocks for the twelve parameters.
- All 2,388 final frame pairings at lengths two, four, six and eight.
- 240 complete polynomial entries compared with the independent word kernel,
  including the two endpoints and vanishing odd coefficients.
- Entry certificates attaining coefficient-span dimensions five, six and
  seven at lengths eight, ten and twelve for `q=a^3`. The parity bound
  above supplies the matching upper bounds.

The open problem is an effective decomposition of the crossing into
coefficient-span pieces with favorable simultaneous cut ranks, uniformly
as the selected label width grows. Arbitrary invertible splittings from
section 1 remain available if the specialization (2) is too rigid.
The construction and finite span certificates establish neither such
a rank theorem nor a higher score. The strongest complete score remains
`1/2`. These arguments and tools used substantial Codex AI assistance;
no new Lean work or independent acceptance is claimed.
