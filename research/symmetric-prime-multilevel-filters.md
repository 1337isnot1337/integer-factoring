# Symmetric multilevel filters in every prime characteristic

## Result and scope

The occupation-weight reduction extends to every prime characteristic.
An arbitrary symmetric polynomial filter on `r` selected `p`-state
identity modes has only `r(p-1)+1` scalar polynomial coordinates. The
compression retains every binomial coefficient, every endpoint sign,
every exact beta order, and every full infinity degree.

This is a constructive search interface. It permits nonhomogeneous and
multilevel filters and combines exactly across disjoint mode blocks.
The statement does not supply a new factoring exponent or an absolute
pressure floor for the whole family.

## 1. Vandermonde grouping retains the complete digit weights

Let `p` be prime, `n=p^a`, `P` a positive odd integer, and `m=nP`.
Select any `r<=n` physical identity indices. For occupation vectors
`k in {0,...,p-1}^r`, write

\[
f_k(x)=\alpha_{\mathbf{p-1}-k}(x)=f_{|k|}(x),\qquad C=1-x^m.
\]

Reindexing the exact [prime-state joint transform](joint-binary-filter-synthesis.md#5-the-same-construction-in-every-prime-characteristic)
by `ell=k-h` gives

\[
\beta_k(x)=\sum_{\ell\le k}
 (-1)^{|\ell|}\left(\prod_{i=1}^r\binom{k_i}{\ell_i}\right)
 C^{|k|-|\ell|}f_{|\ell|}(x).
\tag{1}
\]

For every integer `h`, ordinary Vandermonde convolution gives

\[
\sum_{\substack{0\le\ell_i\le k_i\\\sum_i\ell_i=h}}
 \prod_i\binom{k_i}{\ell_i}
=\binom{|k|}{h}.
\tag{2}
\]

This identity holds over the integers before reduction modulo `p`.
It accounts for the full product of digit binomial coefficients.
For `w=|k|`, equations (1)--(2) yield

\[
\boxed{\beta_w(x)=\sum_{h=0}^w
 (-1)^h\binom wh(1-x^m)^{w-h}f_h(x),
 \qquad 0\le w\le r(p-1).}
\tag{3}
\]

No additional factorial or multiplicity is discarded. The `f_h` are
the coefficient polynomials of the same diagonal occupation filter as
in (1); the grouped binomial coefficient is precisely their full
state multiplicity in that transform.

The matrix in (3) is an involution. Indeed, its `(w,j)` entry after
squaring is

\[
\binom wj C^{w-j}(-1)^j
 \sum_{h=j}^w(-1)^h\binom{w-j}{h-j},
\]

which is one when `w=j` and zero otherwise. In particular, the natural
filter is recovered exactly:

\[
f_h=1\text{ for all }h
\quad\Longrightarrow\quad
\beta_w=(C-1)^w=(-x^m)^w.
\tag{4}
\]

Thus its order is `mw`, with endpoint sign `(-1)^w`. For every filter,

\[
\beta_w(1)=(-1)^w f_w(1).
\tag{5}
\]

These signs matter in odd characteristic. All admissible endpoint
values remain nonzero, and every beta polynomial has finite order.

## 2. Exact feasibility in the reduced coordinates

Put `N=r(p-1)` and prescribe degree/order bounds `(D_w,B_w)` for
`0<=w<=N`. The coefficients in (3) lie in `F_p[x^m]`, so the
residue decomposition from the binary symmetric construction applies
unchanged. There are at most `2(N+1)+1` residue patterns. Within each
pattern, the coefficient system has `N+1` polynomial coordinates,
with degree and order bounds measured in units of `m`.

Let `E subset F_p^(N+1)` be the sum of the residue endpoint images.
Choose a fixed extension field of order `Q>N+1`. A vector in the
extended image has every coordinate nonzero exactly when no coordinate
projection of `E` is zero: the coordinate-zero conditions are at most
`N+1` proper hyperplanes, whose union cannot cover the space over
`F_Q`. The smallest power of `p` above `N+1` suffices.

The [synthesis routine](../tools/symmetric_prime_filter_synthesis.py)
computes these images over the prime field and retains complete
coefficient lifts. Returned polynomials are replayed in (3) before
their actual orders and degrees are reported. The constant physical
endpoint repair uses the inverse of `f_w(1)` in the chosen fixed field
and contributes no polynomial degree or Smith order.

## 3. The physical measure and block cost

For the actual selected energies `lambda_i=2n-1-2j_i`, define

\[
\prod_{i=1}^r\left(\sum_{d=0}^{p-1}e^{-t\lambda_i d}y^d\right)
=\sum_{w=0}^{r(p-1)}c_w(t)y^w,\qquad t=P\log s.
\tag{6}
\]

The weight distributions are

\[
p_w=\frac{c_w(t)}{\sum_hc_h(t)},\qquad
q_w=\frac{c_w(t)e^{-tJw}}{\sum_hc_h(t)e^{-tJh}}.
\tag{7}
\]

Using the exact `B_w=ord_x beta_w` and full `D_w=deg f_w`, the
complete relative two-cut cost remains

\[
\boxed{\mathcal R=
\left(\sum_wq_w e^{2t(nw-B_w/P)}\right)
\left(\sum_wp_w e^{2tD_w/P}\right).}
\tag{8}
\]

Equation (6) uses each actual digit state once and keeps its individual
physical energy. Its coefficients are generally not the binomial or
uniform digit-sum multiplicities, because the energies need not agree.

For disjoint blocks, the original state transform (1) factors exactly.
Their beta polynomials multiply, their zero orders and infinity degrees
add, and both occupation measures factor. Consequently
`mathcal R_total=product_blocks mathcal R_block`, including for different
filters and different physical energy lists. Untouched modes remain
in the full Gaussian and physical partitions.

A source filter at any prime-power interval `n_0>=r` and odd
ramification `P_0` transports to an ambient `n,P` by substituting
`x -> x^(nP/(n_0P_0))`, whenever the exponent is an integer. One may
take a common odd multiple of the source ramifications. This allows
two- and four-mode qutrit blocks even though neither block size is
itself a power of three: their source intervals are `n_0=3` and
`n_0=9`, respectively. In an odd ambient interval, unfilled modes are
retained as explicit Gaussian spectators.

## 4. Original-state controls

The [independent checker](../tools/check_symmetric_prime_filter_synthesis.py)
compares (3) with the original multivariate transform (1), using the
existing digit-state coefficient engine. It includes characteristic
two, qutrit two- and four-mode blocks, and characteristic-five controls.
Its [record](../evidence/symmetric-prime-filter-synthesis.json) contains:

- 18 exact natural-filter controls, including every endpoint sign;
- 655 nonconstant original digit-state polynomial comparisons;
- 48 feasibility comparisons with full, uncompressed digit-state
  coefficient equations;
- 26 synthesized original-state polynomial controls;
- a nine-state characteristic-three product of two different blocks,
  including every product order, full infinity degree, and endpoint;
- six exact rational cost factorizations at distinct physical energies
  5 and 1 and at positive, negative, and zero charge.

In the different-block control, the two univariate digit filters are

\[
(1+x^3+2x^6,\ 1,\ 1+x+2x^4)
\]

and

\[
(1,\ 1+2x^3+x^7,\ 1).
\]

Every endpoint is one in `F_3`, and the beta endpoints are the required
signs `(-1)^w`. The full-state replay retains the original digit
binomial coefficients throughout.

```sh
python3 tools/check_symmetric_prime_filter_synthesis.py --output /tmp/symmetric-prime-filter-replay.json
```

## 5. Search interface and limits

The [qutrit tiled driver](../tools/search_symmetric_prime_tiled_filters.py)
builds finite exact degree/order catalogues, records their source
bounds and literal polynomial hashes, and chooses a profile separately
for every block's actual energies. It searches contiguous, interleaved,
reflected, and deterministically shuffled mode partitions. Its manifests
record every untouched spectator and every completed optimization.

The numerical objective is the positive root of the full paired-type
pressure, with the complete `p`-state physical partition. It does not
optimize a nearly zero absolute pressure as a substitute for the root.
Continuous searches remain exploratory, and a strict gain would require
an independent literal and outward-pressure replay. A finite catalogue
or unsuccessful optimizer cannot exclude the complete symmetric family.
