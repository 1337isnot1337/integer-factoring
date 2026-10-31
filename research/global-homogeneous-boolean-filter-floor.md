# Every homogeneous one-level Boolean filter has positive critical pressure

## Result and scope

Every [homogeneous one-level Boolean occupation filter](homogeneous-boolean-filter-bounds.md)
has strictly positive complete critical pressure, at every real aggregate
charge and every positive scale. The result permits arbitrary overlapping
monomials, arbitrary degree up to the number of selected modes, and every
admissible power-of-two interval length. No read-once restriction remains.

More precisely, let `n` be a power of two, let `P` be positive and odd, put
`m=nP`, and select any `r<=n` distinct physical identity indices. Let `g`
be a nonzero homogeneous multilinear polynomial over `F_2`, of degree
`1<=d<=r`. For any integer `L>md`, the occupation filter is

\[
f_k(x)=1+(x^{md}+x^L)g(k).
\tag{1}
\]

At the critical parameter `rho=1/sqrt(2)`, its full pressure is strictly
positive for every `J in R` and `t=P log s>0`. The endpoint remains one,
and the full degree at infinity is included. This is an ordinary theorem
with exact compact certificates, without Lean formalization.

This closes the homogeneous one-level family. Arbitrary multilevel
occupation filters, nonhomogeneous cancellation patterns, and unrelated
constructions remain outside the theorem. There is no improved factoring
exponent or strict crossing witness in this result.

## 1. Reduce the filter family to a Gaussian support floor

Write

\[
R=\sqrt2,\qquad f(u)=\log(1+e^{-u}).
\]

Let `G_n(J,t)` denote the complete binary Gaussian critical pressure from
the [global Gaussian theorem](global-charged-gaussian-floor.md). The
filtered pressure is `G_n(J,t)+log mathcal R`, where the exact cost ratio is

\[
\mathcal R=(1-q+qa)(1-p+pH/a),\qquad
H=e^{2ntd},\quad 0<a<1.
\tag{2}
\]

Here `p` and `q` are the support probabilities of `g` under the physical
and charge-shifted product occupation measures. The previous
[charge-band theorem](homogeneous-boolean-charge-bands.md) proves

\[
\log\mathcal R>-n f(2nt)
\tag{3}
\]

uniformly over every choice in (1) and every real charge. The
[large-scale comparison](homogeneous-boolean-filter-bounds.md#4-uniform-comparison-over-every-aggregate-charge)
also proves `mathcal R>1` whenever `t>=t_n^*`, where

\[
e^{2nt_n^*}(1-e^{-t_n^*})=1.
\tag{4}
\]

The left side of (4) is strictly increasing. Since `G_n(J,t)>0` at all
charges and scales, it suffices to prove

\[
\boxed{G_n(J,t)>n f(2nt)\quad(0<t\le t_n^*).}
\tag{5}
\]

We prove (5) at every admissible power-of-two length. We do not need or
assert this stronger Gaussian bound above `t_n^*`.

The charge-band theorem has already proved the bound in (5), at every
scale, for `|J|>=2n`, and for `|J|<=2n-2-2R` when `n>=17`. By evenness,
the new work concerns the positive interior edge and finitely many small
lengths.

## 2. Two uniform estimates near the exterior boundary

The full pressure is

\[
\begin{aligned}
G_n(J,t)={}&\sum_{i<n}[f((2i+1-J)t)+f((2i+1+J)t)]\\
&+\sum_{i<2n}f((2i+1-2n)t)\\
&-R\sum_{i<n}[f((4i-2n+2-J)t/R)+f((4i-2n+2+J)t/R)].
\end{aligned}
\tag{6}
\]

For `J>=0`, each charge derivative in the first line lies in `[0,t]`.
Each paired derivative in the last line lies in `[-t,0]`. This follows
directly from `-1<f'<0` and the monotonicity of `f'`. Therefore

\[
|\partial_JG_n(J,t)|\le nt.
\tag{7}
\]

The existing exterior comparison gives

\[
G_n(2n,t)\ge\sum_{i<n}k((2i+1)t),\qquad
k(u)=2[f(u)+f(2u)-Rf(Ru)].
\tag{8}
\]

Its previously proved scalar inequality says `k(u)>f(2u)` for all
`u>0`. We need one additional compact strengthening:

\[
\boxed{k(u)-f(2u)>1/10\qquad(0\le u\le2).}
\tag{9}
\]

The outward-integer checker proves (9) by a complete subdivision of
`[0,2]`. Every leaf has a strictly positive lower bound for
`2f(u)+f(2u)-2R f(Ru)-1/10`; the proof record stores the full cover.
This is the only new scalar inequality used below.

Put `z=nt` and let `M` count the indices with `(2i+1)t<=2`. Combining
(8)--(9) with Jensen's inequality gives

\[
G_n(2n,t)>n f(2nt)+M/10.
\tag{10}
\]

If `z<=B`, where `B>=1`, the elementary count gives
`M>=n/B-1/2`. For `z<=1`, all `n` indices are counted. For `z>=1`,
use `M=min(n,floor(1/t+1/2))>=n/z-1/2`.

## 3. The edge at every power-of-two length at least 512

Write `J=2n+c` in the remaining edge band. Then

\[
-(2+2R)\le c\le0,\qquad |c|<5.
\tag{11}
\]

We need only `t<=t_n^*`, by Section 1.

### 3.1 Bounded values of `z=nt`

For `n>=1024` and `z<=4`, equations (7) and (10) give

\[
\begin{aligned}
G_n(2n+c,t)-n f(2nt)
&>n/40-1/20-5z\\
&\ge n/40-1/20-20>0.
\end{aligned}
\tag{12}
\]

The final inequality already holds at `n=1024` and increases with `n`.

At `n=512`, the entire required scale range satisfies `z<3`. Indeed,
`e^3>20` implies `e^6>400`, while
`1-e^{-x}>=x/(1+x)` for `x>=0`. Thus

\[
e^6(1-e^{-3/512})>400\cdot\frac3{515}>1,
\]

so (4) gives `512 t_512^*<3`. Apply (10) with `B=3`:

\[
G_{512}(1024+c,t)-512f(1024t)
>512/30-1/20-15>0.
\tag{13}
\]

### 3.2 Values `z>=4`

The exact identity from the global Gaussian theorem holds for every
real `c`, including the negative range in (11):

\[
G_n(2n+c,t)=\sum_{i<n}g_i(c,t),
\tag{14}
\]

where, with `b=2i+1`,

\[
g_i(c,t)=2f(bt)+f((c+2b-1)t)+f((c+2b+1)t)
          -2R f((c+2b)t/R).
\tag{15}
\]

The first five summands equal `G_5(10+c,t)>0`. For `i>=5`, set
`u=bt` and `alpha=2+c/b`. Since `b>=11` and `c>-5`,

\[
\alpha>2-5/11>3/2.
\]

Convexity and the increasing dependence on `alpha` give

\[
g_i(c,t)\ge2[f(u)+f(\alpha u)-Rf(\alpha u/R)]
\ge2h(u)>0,
\tag{16}
\]

where the previously proved scalar kernel is

\[
h(u)=f(u)+f(3u/2)-R f(3u/(2R))>0.
\]

Its existing small-argument proof gives, for `0<u<=1/4`,

\[
h(u)>8/21-5/16=23/336>1/16.
\tag{17}
\]

Thus every index `i>=5` with `(2i+1)t<=1/4` contributes more than
`1/8`. When `z>=4`, the upper index cutoff lies below `n`; the number
of these indices is at least `1/(8t)-6`. Discard all other positive
terms in (14) to obtain

\[
G_n(2n+c,t)>1/(64t)-3/4=n/(64z)-3/4.
\tag{18}
\]

For `z>=4`, the decreasing function `z e^{-2z}` and `e^8>2500>512`
give `e^{-2z}<1/(128z)`. Therefore

\[
n f(2z)<n e^{-2z}<n/(128z).
\tag{19}
\]

Finally write `n=2^k`, `k>=10`. The established threshold bound is

\[
t_n^*<\frac{\log(2n)}{2n-1}.
\]

Since `n/(2n-1)<=2/3` and `log 2<1`, our scale range gives
`z<(2/3)(k+1)`. Hence

\[
\frac n{128z}>\frac{3\cdot2^k}{256(k+1)}
\ge\frac{3\cdot1024}{256\cdot11}=\frac{12}{11}>\frac34.
\tag{20}
\]

The sequence `2^k/(k+1)` increases for `k>=1`. Equations (18)--(20)
prove (5) on this branch. The branch is empty at `n=512` by (13).
This completes every admissible length `n>=512`.

## 4. The nine finite lengths

It remains to treat `n in {1,2,4,8,16,32,64,128,256}`. Define

\[
F_n(J,t)=G_n(J,t)-n f(2nt).
\]

### 4.1 An analytic neighborhood of zero scale

At `t=0`,

\[
F_n(J,0)=(3-2R)n\log2>2n/21.
\tag{21}
\]

For `0<=J<=2n`, the positive zero-cut energy sum is at most `4n^2`,
and the positive infinity-cut energy sum is `n^2`. The total absolute
negative physical energy is

\[
\sum_{i<n}|4i-2n+2+J|\le n^2+nJ\le3n^2.
\]

Using `-1<f'<0` in (6), and noting that the derivative of the subtracted
floor is nonnegative, gives `partial_t F_n>=-8n^2`. Consequently

\[
F_n(J,t)>n(2/21-8/128)>0
\quad(0<t\le1/(128n),\ 0\le J\le2n).
\tag{22}
\]

### 4.2 Compact domains and the large-scale cutoff

The checker gives complete positive covers for these nine rectangles:

| Lengths | Scale interval | Charge interval |
|---|---|---|
| `1,2,4,8,16` | `[1/(128n),3/(2n)]` | `[0,2n]` |
| `32,64,128,256` | `[1/(128n),3/n]` | `[2n-5,2n]` |

For the second row, the existing interior charge band covers every
smaller nonnegative charge: the two regions overlap because
`2+2R<5`. The exterior theorem covers `J>=2n` in both rows.

These compact scale intervals reach past `t_n^*`. For `n<=16`,

\[
e^3(1-e^{-3/(2n)})>20\cdot\frac3{2n+3}
\ge\frac{60}{35}>1.
\tag{23}
\]

For `n<=256`,

\[
e^6(1-e^{-3/n})>400\cdot\frac3{n+3}
\ge\frac{1200}{259}>1.
\tag{24}
\]

Equations (22)--(24), the compact certificates, and evenness prove
(5) at all nine lengths. Together with Section 3 and the established
charge bands, this proves (5) for every power of two. Equations
(2)--(4) then prove the theorem stated at the start.

## 5. What the exact certificate establishes

The [checker](../tools/check_global_homogeneous_boolean_filter_floor.py)
uses the existing 80-bit outward-integer exponential and logarithm
bounds. It evaluates the complete pressure after the reciprocal
identity has removed its exact linear part. For `J>=0`, that part is

\[
t\left[n^2+\sum_{i<n}(J-(2i+1))_+
             -\sum_{i<n}|4i-2n+2+J|\right].
\tag{25}
\]

All knots in the bracket are integers. Evaluating its interval endpoints
and every interior integer therefore gives its exact range. All remaining
thermal arguments are nonnegative. Every pressure enclosure includes the
complete paired physical term and the subtracted floor `n f(2nt)`.

The [evidence](../evidence/global-homogeneous-boolean-filter-floor.json)
stores every subdivision path and positive lower numerator. An independent
tree traversal checks coverage and recomputes every leaf. Negative
controls require rejection of a deleted leaf and an inflated bound.
The completed record contains 13,323 positive leaves and 26,636 independently
audited tree nodes. The scalar gap uses 100 leaves; the nine length bases
use 526, 844, 1,585, 2,917, 5,556, 475, 452, 435, and 433 leaves,
respectively. It also checks 535 exact signed profiles and 84 independent
complete-expression interval overlaps.
Exact signed-energy multiset controls and interval comparisons against
the independent complete-pressure expression check the cancellation
formula. The compact certificates, together with the analytic arguments
above and the explicitly cited earlier theorems, cover all real parameters
in the claimed family.

The stronger interior assertion

\[
G_n(J,t)>\sum_{i<n}f((4i+2)t)
\]

is false. At `n=4,J=0,t=3/5`, its difference has a strictly negative
outward interval. The checker records that counterexample exactly.
The weaker support floor (5) avoids this failure.

```sh
python3 tools/check_global_homogeneous_boolean_filter_floor.py --output /tmp/global-homogeneous-filter-replay.json
```
