# Complete binary Gaussian floors at the three small residual bases

## Result

For every scale \(t>0\), the complete binary Gaussian critical
pressure is strictly positive in each of the following ranges:

| Length \(n\) | Charge interval |
|---|---|
| 1 | \(0\le J\le2\) |
| 2 | \(0\le J\le4\) |
| 3 | \(0\le J\le2\) |

Evenness supplies the corresponding negative charges. The pressure
has the full two-cut and paired-physical normalization in the
[charged Gaussian theorem](charged-gaussian-pressure-floor.md).
These bases are analytic pressure statements; no prime-power
restriction on an algorithmic realization is needed to evaluate them.

The proof uses elementary bounds outside a compact scale interval and
a complete outward-integer certificate inside it. There are no
floating-point sign decisions.

## 1. Exact cancelled profiles

Write \(R=\sqrt2\) and \(f(x)=\log(1+e^{-x})\). The positive
thermal contribution at both cuts is always

\[
Z_n(J,t)=\sum_{i=0}^{n-1}\left[
 f(|2i+1-J|t)+f((2i+1+J)t)+2f((2i+1)t)
\right].
\tag{1}
\]

The full pressure is

\[
G_n(J,t)=g(J)t+Z_n(J,t)-2R\sum_{e\in\mathcal E}f(et/R),
\tag{2}
\]

with the following complete data:

| Length and charge range | Ground coefficient \(g(J)\) | Physical positive energies \(\mathcal E\) |
|---|---|---|
| \(n=1,\ 0\le J\le2\) | \((1-J)_+\) | \(J\) |
| \(n=2,\ 0\le J\le2\) | \((J-1)_+\) | \(2-J,\ 2+J\) |
| \(n=2,\ 2\le J\le4\) | \((3-J)_+\) | \(J-2,\ J+2\) |
| \(n=3,\ 0\le J\le2\) | \((1-J)_+\) | \(J,\ 4-J,\ 4+J\) |

Here \(x_+=\max(x,0)\). The table follows by applying
\(f(-x)=f(x)+x\) to every original energy. Zero physical energies
are retained, with their full \(f(0)=\log2\) contribution.

## 2. Small scales

At scale zero the continuous limiting value is

\[
G_n(J,0)=2n(2-R)\log2>16n/21,
\tag{3}
\]

using \(R<10/7\) and \(\log2>2/3\).
For all real \(x\), \(-1<f'(x)<0\). Differentiate the original
complete energy expression, before cancellation. A lower derivative
bound is the negative of the sum of positive cut energies and the
absolute negative physical energies. Throughout the stated charge
ranges that sum is at most \(4n^2+8\).

One can check this directly from the short energy lists: its actual
maxima for \(n=1,2,3\) are 6, 24, and 38, respectively. Each sum
is convex and piecewise linear in \(J\), so only the charge endpoints
are needed for its maximum. Thus for \(0<t\le1/32\),

\[
G_n(J,t)>\frac{16n}{21}-\frac{4n^2+8}{32}>0
\qquad(n=1,2,3).
\tag{4}
\]

All final comparisons in (4) are rational and are replayed exactly.

## 3. Large scales

Use \(y=J\) except in the third row of the table, where
\(y=J-2\). Then \(0\le y\le2\), a positive cut term is
\(f(|1-y|t)\), and the negative physical contribution has at most
three terms, each with coefficient \(2R\).

There are two cases. In the half of each row where \(g(J)=0\),
put \(a=|1-y|\) and let \(b\) be the smallest physical energy
divided by \(R\). The table gives

\[
b-a\ge R-1.
\tag{5}
\]

For the even row this is
\((2-y)/R-(1-y)\ge R-1\) on \(0\le y\le1\).
For each other row it is
\(y/R-(y-1)\ge R-1\) on \(1\le y\le2\).
Concavity of \(\log(1+z)\) on \([0,1]\) and
\(f(x)\le e^{-x}\) show that the ratio of the entire negative
part to the retained positive term is at most

\[
\frac{6R}{\log2}e^{-(R-1)t}
<\frac{27}{2}e^{-2t/5}<1\qquad(t\ge9).
\tag{6}
\]

Here \(R>7/5\), \(R<3/2\), \(\log2>2/3\), and
\(e^3>20\) suffice.

In the other half of each row, put \(v=|1-y|\in[0,1]\).
The ground and the retained positive term give

\[
vt+f(vt)\ge\max(vt,\log2),
\]

while the smallest physical argument is \((1-v)t/R\).
For \(t\ge6\), the elementary bound from the
[periodic-continuation proof](charged-gaussian-pressure-floor.md) is

\[
\max(vt,\log2)e^{(1-v)t/R}\ge t.
\tag{7}
\]

For clarity, on \(v\le(\log2)/t\) the left side decreases; on
the other branch its logarithm is concave. Its minimum is therefore
\(\min\{t,(\log2)e^{(t-\log2)/R}\}=t\) for \(t\ge6\).
The complete negative contribution is at most
\(6R e^{-(1-v)t/R}\), with \(6R<9\le t\).
This proves strict positivity for every row at \(t\ge9\).

## 4. Compact certificate and replay

The remaining domains have scale interval \([1/32,9]\), split
into the four table rows. The
[checker](../tools/check_small_charged_gaussian_bases.py) evaluates (2)
with the integer-outward arithmetic from
[`charged_gaussian_floor.py`](../tools/charged_gaussian_floor.py).
It keeps the cancelled ground term and every complete thermal factor;
no energy is omitted. Absolute values and positive parts use their
exact interval images.

The full subdivision trees contain 1,627 positive rectangles:

| Length and charge interval | Positive rectangles | Smallest lower numerator |
|---|---:|---:|
| \(n=1,\ [0,2]\) | 189 | 129261008993479454648 |
| \(n=2,\ [0,2]\) | 411 | 95082161405874356377 |
| \(n=2,\ [2,4]\) | 228 | 129334181005874972918 |
| \(n=3,\ [0,2]\) | 799 | 48264974374160096274 |

Every lower numerator has denominator \(2^{80}\). The record stores
the complete partition paths and all lower bounds. Thirty-five exact
profile controls independently compare the ground coefficients and
physical thermal multisets against the original energies. The same
controls check the derivative estimate used in (4).

```sh
python tools/check_small_charged_gaussian_bases.py --output /tmp/small-bases-replay.json
```

The command regenerates all rectangle trees with integer arithmetic.
The portable record is
[`small-charged-gaussian-bases.json`](../evidence/small-charged-gaussian-bases.json).
Together, (4), (6)--(7), and the complete compact certificate prove the
three base ranges at every positive scale.
