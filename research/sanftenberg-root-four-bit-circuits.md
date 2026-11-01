# Root-four bit circuits for the semiprime transform

## Scope and source

This ordinary arithmetic analysis uses the pinned
[`RiemannGaussian` source](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635).
It studies the exact word operations before claiming any global circuit-width
or factoring result. The source defines the transform ring as
`Z/(2^T+1)`, uses root four, and stores canonical residues in `T+1` bits:
[`SemiprimeMacroFermatNetwork.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroFermatNetwork.lean#L21-L49),
[`SemiprimeMacroBitNetwork.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroBitNetwork.lean#L21-L45).
The source's explicit twiddle stage multiplies by a power of four
([`MacroBitNetwork`, lines 52–78](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroBitNetwork.lean#L52-L78)).

## Twiddles are signed cyclic rotations

Let `M=2^T+1`. Then `2^T=-1` in `Z/MZ`, so `4^T=1`. Reduce a twiddle
exponent to `0<=a<T` and write `2a=epsilon*T+s`, where
`epsilon` is 0 or 1 and `0<=s<T`. A canonical residue has the form

```math
x=h2^T+\sum_{i=0}^{T-1}b_i2^i,
\qquad h\in\{0,1\},\quad h=1\Longrightarrow b_i=0\text{ for all }i.
```

For `4^a x=2^{2a}x`, the low bits move by the cyclic shift `i -> i+s mod T`.
Each moved bit receives sign
`(-1)^(epsilon+floor((i+s)/T))`; the top bit contributes
`(-1)^(1+epsilon) h 2^s`. At destination `s`, the two contributions combine
to `(-1)^epsilon (b_0-h)2^s`. Canonicality makes `b_0-h` one of `-1,0,1`.
Thus the result before canonical reduction is a signed binary sum

```math
S=\sum_{j=0}^{T-1}c_j2^j,\qquad c_j\in\{-1,0,1\},\qquad
|S|\le 2^T-1,
```

congruent to `4^a x` modulo `M`. Its canonical representative is `S` when
`S>=0`, and `S+M` when `S<0`.

This gives a fixed `O(T)`-gate circuit. Read the signed digits from low to high
with carry state `k_j`:

```math
c_j+k_j=y_j+2k_{j+1},\qquad y_j\in\{0,1\},\quad k_j\in\{-1,0,1\},\quad k_0=0.
```

The terminal carry is `-1` exactly when `S<0`. If it is `-1`, increment the
low `T` output bits once; the increment carry is the canonical top bit. This
computes `S+2^T+1` and handles `S=-1` by producing the valid top-only word
`2^T`. Every position uses a constant-size signed-carry and increment cell.
The rotation wires are a permutation, not a generic constant multiplier.

## Modular additions and pointwise products

For canonical `x,y` in `[0,M-1]`, addition computes `x+y<2M` and conditionally
subtracts `M` once. Subtraction computes the signed difference in
`[-(M-1),M-1]` and conditionally adds `M` when negative. Ripple carry, one
sign/comparison chain, and per-bit selection use `O(T)` Boolean gates for
either operation. Therefore one radix-2 butterfly—two modular add/subtract
outputs and its root-four twiddle—has an `O(T)`-gate implementation.

For a pointwise product, schoolbook multiplication computes an integer
`P=x*y` with `O(T^2)` gates. Write `P=L+H2^T` with `0<=L<2^T` and
`0<=H<=2^T`. Since `2^T=-1`, the residue is `L-H mod M`; one signed
correction by `M` returns its canonical `T+1`-bit word. The usual row-wise
schoolbook circuit keeps one `O(T)`-bit accumulator and the multiplier word,
so an individual product gadget has a path decomposition of width `O(T)`.

With `T` transform positions and `log_2 T` butterfly stages, the transform
circuits have `O(T^2 log T)` gates. The `T` pointwise product gadgets have
`O(T^3)` gates in total. These are polynomial in the original input length
because `T^2=Theta(ell)`. This gate-count result does not bound the treewidth
of their composition.

## The global separator is still the question

The previous
[dense-table analysis](sanftenberg-transform-treewidth-floor.md) gives a
`Q_d` clique-blow-up minor and a width floor
`Omega(ell/sqrt(log ell))` for extensional four-word butterfly tables.
The circuits above change that local representation: add/subtract use carry
chains and the root-four twiddle is a signed rotation. But the complete
arithmetic network also connects `T` variable pointwise products, both forward
transforms, the inverse transform, and the coefficient/carry recovery. A
local `O(T)`-width product gadget does not by itself give a sub-`ell/sqrt(log
ell)` separator for this connected network.

The next proof obligation is an explicit tree decomposition or elimination
order for the *entire* bit circuit, including canonical-word checks, all
pointwise multipliers, the inverse transform, carry recovery, and every
interface wire. A bound below `ell/sqrt(log ell)` could change the evaluator
cost; gate count alone cannot. The source's semiprime promise and complete
primitive-machine cost also remain separate obligations. No stronger
factoring exponent or Lean result is claimed here.

No upstream code was copied. Substantial OpenAI assistant support was used in
the source and circuit analysis.
