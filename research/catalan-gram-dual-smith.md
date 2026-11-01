# Catalan cup Gram poles near characteristic two

## Result and scope

The full cup frame has a natural Laurent-field dual

\[
D^{\rm G}=(F^TF)^{-1}F^T
\]

when its Gram matrix is invertible over the rational function field. Near
`q=1` in characteristic two, this dual has poles. An unrestricted
Clebsch--Gordan path basis gives a uniform `O(n^2)` pole bound for this
decoder, so any polynomial cutoff uses only polynomially many Laurent
coefficients. The exact finite Smith profiles below show arithmetic jumps
through sixteen physical strands. This does **not** improve the crossing
rank sum or establish a stronger factoring exponent. The regular
unitriangular coordinate dual remains the current constructive default.

## 1. Loop Gram matrix

Let `n=2r`, and let `P,Q` range over the `C_r` noncrossing perfect matchings
of `n` ordered positions. Their cup vectors have coefficients `1` on `(1,0)`
and `q` on `(0,1)`. Overlaying `P` and `Q` gives `ell(P,Q)` alternating
loops. Summing the two alternating bit orientations on each loop gives

\[
 (F^TF)_{P,Q}=q^r\delta^{\ell(P,Q)},\qquad
 \delta=q+q^{-1}.                                      \tag{1}
\]

The scalar `q^r` is a unit, so it does not affect the Smith valuations. Every
overlay has at least one loop. Write

\[
F^TF=q^r\delta M_r(\delta),\qquad
[M_r(\delta)]_{P,Q}=\delta^{\ell(P,Q)-1}.             \tag{2}
\]

For `q=1+u` over `F_2`,

\[
\delta=q+q^{-1}=\frac{u^2}{1+u}.                       \tag{3}
\]

Thus a `delta`-adic Smith exponent `s` for `F^TF` corresponds to a pole of
order `2s` in `u` for its inverse.

## 2. Exact truncated-series certificate

For each `k`, multiplication by `M_r(delta)` on
`F_2[delta]/(delta^k)` is represented by a `k C_r` square lower-triangular
block Toeplitz matrix. If the Smith exponents of `M_r` are `s_i`, the kernel
dimension is

\[
\dim_{\mathbb F_2}\ker(M_r\bmod\delta^k)
    =\sum_i\min(k,s_i).                               \tag{4}
\]

First differences recover the number of exponents at least `k`; the factor
`delta` in (2) shifts every exponent up by one. The checker performs exact
packed-bit elimination on these Toeplitz matrices. A zero tail at consecutive
precisions certifies that all exponents have been recovered for each listed
finite size.

| Physical strands `n` | Catalan dimension `C_(n/2)` | Gram Smith exponents over `F_2[[delta]]` |
|---:|---:|---|
| 2 | 1 | `1^1` |
| 4 | 2 | `1^2` |
| 6 | 5 | `1^4, 3^1` |
| 8 | 14 | `1^8, 3^6` |
| 10 | 42 | `1^16, 3^26` |
| 12 | 132 | `1^32, 3^100` |
| 14 | 429 | `1^64, 3^364, 7^1` |
| 16 | 1430 | `1^128, 3^1288, 7^14` |

Here `s^m` means exponent `s` with multiplicity `m`. The one-loop residue
has rank `2^(r-1)` in all eight tested sizes. The next exponent jumps from
three to seven at `r=7`; this is why a uniform order-six pole claim based on
the smaller cases would be false. The finite largest exponent follows the
pattern `2^floor(log2(r+1))-1` through `r=8`, but that pattern is only a
conjecture from this screen.

The checker also independently verifies (1) against physical word
contractions over `F_16` for all Gram entries through eight strands.

## 3. A uniform pole budget from spin paths

Work over `F_2(q)` and let `V` be the two-dimensional fundamental module.
The unrestricted spin-path splitting `H_h tensor V = H_(h+1) plus H_(h-1)`
has, in each interior magnetic sector `1<=m<=h`, the block

\[
U_{h,m}=\begin{pmatrix}
1 &[m]_q\\
q^{-h+m-1}&q^m[h-m+1]_q
\end{pmatrix},\qquad \det U_{h,m}=[h+1]_q.              \tag{5}
\]

Boundary sectors are one-dimensional with scalar one. Each determinant is
nonzero in `F_2(q)`, so iterating these explicit splittings gives a direct
sum over all paths. The number of paths ending at height zero after `2r`
steps is `C_r`.

For the ordinary word pairing, the magnetic components have norms

\[
g_{h,m}=q^{-m(h-m)}{[h]_q!\over[m]_q![h-m]_q!}.
\]

Substitution of (5), followed by the quantum-binomial recurrence, gives
the local diagonal identity

\[
U_h^T(\operatorname{diag}(g_{h,m})\otimes I)U_h
 =\operatorname{diag}(g_{h+1,m},q[h]_q[h+1]_qg_{h-1,m}).\tag{6}
\]

Thus the return-path vectors `V_P` are orthogonal and have norms

\[
\nu_P=\prod_{\text{down from }h}q[h]_q[h+1]_q.       \tag{7}
\]

They are regular Laurent-polynomial physical vectors. The cup tensors are
invariant under the same splitting and form an independent frame by the
unitriangular minor in
[`catalan-cup-crossings.md`](catalan-cup-crossings.md). Since both frames
have `C_r` vectors, they span the same invariant subspace. Write `V=FC`;
the coordinate decoder of that note computes `C` without poles. Then

\[
D^{\rm G}=C\,\operatorname{diag}(\nu_P^{-1})V^T,
\qquad D^{\rm G}F=I.                                  \tag{8}
\]

Indeed, `F=VC^{-1}` and `V^TV=diag(nu)`. The matrices `C` and `V` are
regular at `q=1`, so all poles in (8) come from the path norms.

In characteristic two the quantum-integer recurrence gives

\[
v_\delta([k]_q)=2^{v_2(k)}-1.                          \tag{9}
\]

The identity follows from `[2m]_q=delta[m]_q^2` and the fact that odd
quantum integers are units at `delta=0`. A Dyck path of semilength `r` has
`r` down-steps. At each one, exactly one of `h,h+1` is even; it is at most
`r+1`, so its valuation is at most `r`. Therefore

\[
v_\delta(\nu_P)\le r^2,
\qquad \operatorname{pole}_u(D^{\rm G})\le 2r^2,
\qquad \operatorname{pole}_u(D^{\rm G}\otimes D^{\rm G})\le4r^2.
                                                               \tag{10}
\]

Here `delta=u^2/(1+u)`. The braid word and crossing normalization are
regular at `q=1`; hence the full strict crossing has no Laurent powers below
`u^(-4r^2)`. Every nonzero Laurent grade other than zero vanishes at both
strict endpoints. This gives an effective polynomial-size jet window for
rank recombination. It does not bound the two-cut ranks of the pieces.

The checker verifies the unrestricted path vectors, their nonzero
orthogonal norms, and their coordinates in the full cup frame through eight
physical strands. The uniform pole estimate (10) follows from the displayed
all-height splitting and valuation argument; no factoring improvement is
claimed.

## 4. Reproduction

Run the finite exact screen with:

```sh
make verify-catalan-gram-smith
```

The recorded profiles are in
[`catalan-gram-dual-smith.json`](../evidence/catalan-gram-dual-smith.json).
