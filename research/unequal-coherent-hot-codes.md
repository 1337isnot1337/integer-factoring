# Unequal paired alphabets at a coherent energy boundary

This finite screen uses different physical alphabet sizes on the two
paired directions. It includes complete tensors, both natural Smith
profiles, separate optimized physical types, and the asymmetric cost.
None of its 18 retained rational tests passes the required support
threshold 1. This is a finite diagnostic, not an exclusion of unequal
alphabets or other gauges.

## 1. Exact tensors

Use the interval-16 stabilized reflection tensor and the retractions in
[the coherent-boundary construction](coherent-hot-boundary-crossings.md).
At energy 32 the lower code has seven letters and the boundary has two.
At energy 48 the lower code has fifteen letters and the boundary has three.
The X direction retains one coherent boundary vector, or its coordinate
control. The Y direction retains no boundary vectors, the full boundary,
or one of five specified two-dimensional boundary spaces.

The dimension pairs are (8,7), (8,9), (16,15), (16,17), and (16,18).
The report specifies every left/right matrix. Each pair satisfies KJ=I,
so the four maps preserve the strict crossing endpoint. The constructor
retains the complete tensor and removes its exact monomial content; the
vacuum entry proves there is no further scalar factor. The two complete
native profiles have opposite scalar shifts. Removing them gives
nonnegative orders sigma_infinity and sigma_zero, each of length q_X q_Y.

## 2. Full asymmetric cost

Let a,b>0 satisfy a^2+b^2=1, with shared intrinsic scale u>0. Set

\[
L=a\log\sum_i e^{-u\sigma_{\infty,i}/a}
 +b\log\sum_i e^{-u\sigma_{0,i}/b}.
\]

The opposite scalar shifts cancel. For paired types p_X,p_Y the full cost is

\[
\frac{L+u(\mathbb E_{p_X}E_X+\mathbb E_{p_Y}E_Y)}
 {H(p_X)+H(p_Y)}.
\tag{1}
\]

The physical means include both orientations. The Gibbs variational identity
on the two alphabets minimizes (1) at the positive root R of

\[
L=R\left(\log\sum_x e^{-uE_X(x)/R}
             +\log\sum_y e^{-uE_Y(y)/R}\right).
\tag{2}
\]

The [asymmetric transfer](asymmetric-modular-transfer.md) requires R<1.
For an exact test, choose a=p/h,b=q/h with p^2+q^2=h^2, a rational
0<z<1, and u=pq log(1/z). Put

\[
P_\infty=\sum_i z^{qh\sigma_{\infty,i}},\quad
P_0=\sum_i z^{ph\sigma_{0,i}},\quad
Q_X=\sum_x z^{pqE_X(x)},\quad Q_Y=\sum_y z^{pqE_Y(y)}.
\]

At exact rational types proportional to z^(pq E_X) and z^(pq E_Y),
the test for (1)<1 becomes

\[
P_\infty^pP_0^q<(Q_XQ_Y)^h.
\tag{3}
\]

The checker clears all denominators and compares positive integers.
Failure of (3) concerns this retained rational scale, not the complete
one-dimensional scale search.

## 3. Completed finite screen

The six branch triples are (3,4,5), (4,3,5), (5,12,13), (12,5,13),
(20,21,29), and (21,20,29). All 108 numerical scale searches completed
across 18 tensors in approximately 24.393 seconds. Each tensor retains a
rational approximation to its best sampled scale and repeats (3) exactly.
All 18 exact tests fail the strict threshold.

The best numerical support is approximately 1.009628782 for dimensions
16 by 18, with coherent X and full Y boundary. This is above 1; the
numerical optimizer does not certify a global scale minimum.

The report is
[`unequal-coherent-hot-screen.json`](../evidence/unequal-coherent-hot-screen.json).
The [screen source](../tools/screen_unequal_coherent_hot_codes.py) retains
the maps, complete profiles, costs, rational comparisons, and source hashes.
The [replay](../tools/check_unequal_coherent_hot_codes.py) reconstructs every
tensor, endpoint and complete profile, repeats saved costs, and repeats
the exact comparisons. Its result is
[`unequal-coherent-hot-screen-replay.json`](../evidence/unequal-coherent-hot-screen-replay.json).

```sh
OPENBLAS_NUM_THREADS=1 python tools/check_unequal_coherent_hot_codes.py
```
