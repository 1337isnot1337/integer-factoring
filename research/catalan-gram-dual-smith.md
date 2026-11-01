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
quantum integers are units at `delta=0`. At a down-step from height `h`,
exactly one of `h,h+1` is even. Its contribution is
`w(h)=2^(v_2(k))-1`, where `k` is that even index; equivalently,
`w(2m-1)=w(2m)=2^(1+v_2(m))-1`.

A Dyck path of semilength `r` is the contour walk of a plane rooted tree with
`r` edges. The number `c_h` of down-steps from height `h` is the number of
tree vertices at depth `h`. If the maximum height is `H`, then
`c_1,...,c_H >= 1` and `sum c_h=r`. Conversely, realize any such profile
by attaching `c_1` children to the root, then attaching all `c_(h+1)` next-
level children to one vertex at depth `h`. Thus for fixed `H`, the maximum
path-norm valuation is

\[
S_H+(r-H)W_H,\qquad S_H=\sum_{h=1}^H w(h),\quad
W_H=\max_{1\le h\le H}w(h).
\]

The prefix records of `w(h)` occur at `H_k=2^k-1`, with `W_{H_k}=H_k`.
Between consecutive records `W_H` is constant and each added `w(h)` is at
most that record, so the objective cannot increase by extending `H` past a
record. Pairing heights `2m-1,2m` gives

\[
S_{2^k-1}=2\sum_{m=1}^{2^{k-1}-1}(2^{1+v_2(m)}-1)+(2^k-1).
\]

For `k=1`, this is `S_1=1`. For `k>=2`, set `d=k-1`. For each
`j=0,...,d-1`, there are `2^(d-j-1)` integers `m` with `v_2(m)=j` in
`1,...,2^d-1`, so

\[
\sum_{m=1}^{2^d-1}2^{v_2(m)}=d2^{d-1}.
\]

Substitution gives `S_{H_k}=(k-1)2^k+1`. Therefore the exact maximum
path-norm valuation is

\[
M_r=\max_{k\ge1,\ 2^k-1\le r}
\left((k-1)2^k+1+(r-2^k+1)(2^k-1)\right).
\]

Since `C` and `V` in (8) are regular at `q=1`, this yields

\[
v_\delta(\nu_P)\le M_r,\quad
\operatorname{pole}_u(D^{\rm G})\le2M_r,\quad
\operatorname{pole}_u(D^{\rm G}\otimes D^{\rm G})\le4M_r.
\tag{10}
\]

A simpler envelope is `M_r <= r^2/4+(r+1)log_2(r+1)`, since
`S_{H_k}<=k2^k<=(r+1)log_2(r+1)` and `H_k(r-H_k)<=r^2/4`. In particular
the paired pole is at most `r^2+4(r+1)log_2(r+1)`, improving its leading
quadratic coefficient by half over the height-only estimate.

The checker computes the exact path-norm valuation distribution by a
height/step dynamic program. Through semilength 32, it confirms the Catalan
path counts and that the observed maximum equals the closed form; an
independent max-plus recurrence matches the closed form through semilength
512. At semilength 32, `M_32=304` versus the coarse height bound `528`. The braid
word and crossing normalization are regular at `q=1`, so the full strict
crossing has no Laurent powers below `u^(-4M_r)`. Every nonzero Laurent grade
other than zero vanishes at both strict endpoints. The window is polynomial
in the strand count, but it does not bound the two-cut ranks of the pieces.

### Mean pole order and path-stratified support

The worst path is much more singular than a typical path. Let `T_(r,h)` be
the total number of down-steps from height `h`, summed over all Dyck paths of
semilength `r`. Marking a depth-`h` vertex in a plane rooted tree gives the
generating function `z^h C(z)^(2h+1)`, where `C=1+zC^2` is the Catalan
generating function. Lagrange inversion therefore gives the exact count

```math
T_{r,h} = [z^r]z^h C(z)^{2h+1}
        = \frac{2h+1}{2r+1}\binom{2r+1}{r-h}.
```

Consequently the exact mean path-norm valuation is

```math
\mathbb E W_r
 = \frac1{C_r}\sum_{h=1}^r
    w(h)\frac{2h+1}{2r+1}\binom{2r+1}{r-h}.             \tag{11}
```

It has order `Theta(r log r)`. To see the upper bound, choose a uniformly
random down-step from a uniformly random Dyck path. Its height distribution
is `p_(r,h)=T_(r,h)/(r C_r)`. The binomial ratio formula gives

```math
p_{r,h}\le \frac{2h+1}{r}\exp\!\left(-\frac{h^2}{2r}\right).
```

Write `w(2m-1)=w(2m)=sum_(a>=0)2^a [2^a divides m]`. At dyadic level `a`,
the contributing heights lie in two residue classes of spacing
`d=2^(a+1)`. Summing the displayed Gaussian envelope on those classes gives
probability `O(1/d)` for `d<=sqrt(r)` and
`O((d/r) exp(-c d^2/r))` for `d>sqrt(r)`. Multiplying by the level weight
`2^a=d/2` costs `O(1)` per level, and there are `O(log r)` levels. For the
lower bound, throughout `sqrt(r)<=h<=2sqrt(r)` the same binomial ratio is
bounded below by a positive constant, so `p_(r,h)=Theta(1/sqrt(r))`. Each
dyadic spacing `d<=sqrt(r)/2` then contributes `Theta(1)`; there are
`Theta(log r)` such levels. Thus the mean weight per down-step is
`Theta(log r)`, proving the claim.

This yields a path-dependent precision invoice. The diagonal factors
`nu_P^(-1)` have `2W(P)` possible negative Laurent grades in `u`; across all
paths their total principal-part slots are
`sum_P 2W(P)=2 C_r E[W_r]=Theta(C_r r log r)`, instead of the uniform
`O(C_r r^2)` padding. For a pair of paths, the corresponding sum is
`Theta(C_r^2 r log r)`. These are scalar support counts, not cut ranks: the
missing step is a two-cut factorization that preserves the path strata and
turns their lower support into a favorable rank sum. No factoring score
improvement follows from (11) alone.

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
The same checker cross-checks (11) against the independently enumerated path
valuation distributions at semilengths `1,2,4,8,16,24,32` and records exact
means at powers of two through 512.
