# Centered divisor windows combine low-bit classes

The [quarter-bit endpoint](circuit-box-ratio-small-roots.md) can use fewer
low bits when it retains the ordinary interval width. Centering its lattice
variable also removes both extra guard bits from the full balanced band.
The exact implementation uses the existing integer LLL certificates and
Sturm root isolation; no new lattice-reduction assumption is introduced.

For an odd \(n\)-bit integer \(N\), let

\[
I=[L,U)\subseteq[2^k,2^{k+1}),\qquad
1\le L<U\le N+1,\qquad U-L\le2^\ell.
\]

Take the least such nonnegative integer \(\ell\), and put

\[
\boxed{h=\max\!\left(0,\ell-
             \left\lfloor\frac{k^2}{n}\right\rfloor\right).}
\tag{1}
\]

**Theorem.** For any fixed residue modulo \(2^h\), every divisor of \(N\)
in that residue class and in \(I\) can be listed in polynomial bit time and
space. The polynomial is uniform in \(N,I\), and the residue. A larger
precision may also be used. The local output has at most \(2n-1\) entries
unless the routine takes its constant-size direct branch.

This is a centered, explicitly charged application of the classical linear
small-root method. It does not improve the repository's factoring score.
Streaming all required classes still gives
\(2^{n/4}\mathrm{poly}(n)\) full factoring time, which is weaker than the
repository's existing logarithmic-saving bound. Its constructive benefit
is that a single query now combines many of the preceding endpoint's
low-bit classes while preserving the ordinary interval selector.

## 1. Center the actual candidate progression

Write \(R=2^h\). Let \(a_*\) be the first integer in \(I\) with the chosen
residue, and let \(C\) be the number of such integers. Empty classes are
discarded. When \(C\le3\), direct division takes at most three tests.
Even residue classes are also discarded when \(h>0\), because \(N\) is odd.

For \(C\ge4\), set

\[
j_0=\lfloor C/2\rfloor,\qquad a_0=a_*+Rj_0,
\qquad a=a_0+Rx.
\]

The integer roots to inspect are exactly in the half-open interval

\[
-j_0\le x<C-j_0.
\tag{2}
\]

Thus \(|x|\le j_0\). Since \(h\le\ell\) and the original interval has
width at most \(2^\ell\),

\[
C\le2^{\ell-h},\qquad
j_0\le2^{\ell-h-1}.
\tag{3}
\]

Let \(X=2^s\) be the least power of two at least \(j_0\). Equations
(1)--(3) give

\[
1\le s\le\ell-h-1
\le\left\lfloor\frac{k^2}{n}\right\rfloor-1.
\tag{4}
\]

The lower bound follows from \(C\ge4\). It also ensures that the shared
lattice engine's nontrivial parameter range applies.

The exact monic congruence is

\[
f(x)=x+c,\qquad c\equiv a_0R^{-1}\pmod N.
\tag{5}
\]

For any divisor \(a=a_0+Rx\) of \(N\), equation (5) implies
\(a\mid f(x)\), because \(a\) is odd. No primality assumption on \(a\)
is used.

## 2. One uniform lattice bound

For integers \(1\le m<d\), take the standard monic-linear lattice generated
by

\[
N^{m-i}f(x)^i\quad(0\le i<m),\qquad
x^j f(x)^m\quad(0\le j<d-m).
\tag{6}
\]

Scale the coefficient of \(x^j\) by \(X^j\). Its determinant is

\[
\mathcal D=N^{m(m+1)/2}X^{d(d-1)/2}.
\tag{7}
\]

Every integer combination of (6), evaluated at a valid root, is divisible
by \(a^m\). The usual \(\delta=3/4\) LLL guarantee and Cauchy's inequality
give a nonzero combination \(g\) with

\[
\sum_j |g_j|X^j
\le\sqrt d\,2^{(d-1)/4}\mathcal D^{1/d}.
\tag{8}
\]

It is enough to make the right side strictly less than \(2^{km}\).
Then \(|g(x)|<2^{km}\le a^m\), so divisibility forces \(g(x)=0\).
This also applies to negative centered roots and roots with \(|x|=X\).

Choose \(d=2n\) and \(m=2k\). Using \(N<2^n\) and (4), the binary
logarithmic margin in (8) is at least

\[
\begin{aligned}
&km-\frac12\log_2 d-\frac{d-1}4
       -\frac{nm(m+1)}{2d}-\frac{s(d-1)}2\\
&\qquad\ge
\frac n2-\frac k2+\frac{k^2}{2n}
       -\frac14-\frac12\log_2(2n)\\
&\qquad=
\frac{(k-n/2)^2}{2n}+\frac{3n}{8}
       -\frac14-\frac12\log_2(2n)>0
       \qquad(n\ge8).
\end{aligned}
\tag{9}
\]

For \(n<8\), a nontrivial radius in (4) is possible only for

\[
(n,k)=(4,3),(5,4),(6,4),(6,5),(7,4),(7,5),(7,6).
\]

All seven pass the exact integer inequality below. Every remaining small
case has at most three candidates and uses direct division. This closes
the theorem for every input size.

The implementation chooses the first smaller lattice dimension that meets
the same sufficient test. For its chosen \(d,m,s\), define

\[
E=4kmd-d(d-1)-2nm(m+1)-2sd(d-1).
\]

It verifies \(d^{2d}<2^E\), using integer bit lengths. The proved choice
\(d=2n,m=2k\) guarantees termination. This test is independent of numerical
logarithms and floating-point approximations.

The lattice has \(O(n)\) rows and columns with \(O(n^2)\)-bit integer
entries. LLL and exact root isolation therefore have polynomial bit cost.
The routine retains an integer row-combination certificate and verifies
the strict weighted coefficient bound. It isolates every integer root of
the resulting nonzero polynomial in (2), then checks divisibility by \(N\).
The polynomial has degree less than \(2n\), which bounds the local output.

The underlying small-root construction and LLL norm estimate are classical;
see Cohn and Heninger,
[Ideal forms of Coppersmith's theorem and Guruswami--Sudan list decoding](https://arxiv.org/abs/1008.1284),
Theorems 1.1 and 2.1 and Section 3. The preceding calculation supplies the
exact finite constants and centered interval convention used here.

## 3. How the existing engine is reused

The shared engine normally uses a radius \(2^{k+1-H}\) at a formal
precision \(H\). Set \(H=k+1-s\), and supply the formal residue
\(a_0 2^{H-h}\). Its linear constant is still

\[
a_0 2^{H-h}(2^H)^{-1}\equiv a_0(2^h)^{-1}\pmod N,
\]

and its radius is exactly \(2^s\). Here \(H-h\ge0\), since the window lies
within one dyadic band. This is a change of lattice coordinates, not a
claim that additional low bits of the unknown divisor are known.

The root interval remains the precise half-open interval (2). The existing
Sturm routine includes negative roots, excludes its upper endpoint, and
tests integer candidates exactly. Neither ordinary interval endpoints nor
individual residues are silently expanded.

## 4. Consequences for complete queries

For a full dyadic band, \(\ell=k\), and

\[
h=k-\left\lfloor\frac{k^2}{n}\right\rfloor
 =\left\lceil k-\frac{k^2}{n}\right\rceil
 \le\left\lceil\frac n4\right\rceil.
\tag{10}
\]

The previous uniform endpoint used \(\lceil n/4\rceil+2\) bits. Thus the
new full-band schedule uses at most one quarter as many residue classes;
unbalanced bands may use much fewer. A shorter ordinary interval reduces
the required low precision by the same number of bits until it reaches
zero. In particular, if \(\ell\le\lfloor k^2/n\rfloor\), the whole window
is one polynomial-cost query with no low-residue condition.

More explicitly, an unbalanced full band with \(k\) near \(\beta n\)
requires \(2^{\beta(1-\beta)n+O(1)}\) streamed classes. This is the familiar
classical small-root tradeoff, with the width and constants charged here.
It does not give a stronger all-input logarithmic-saving exponent.

For the smallest prime factor, scan dyadic bands in increasing order up to
\(\lfloor\sqrt N\rfloor\). Inspect every required odd residue of a band,
retain its least divisor, and stop at the first nonempty band. That divisor
is prime because any composite divisor has a smaller proper divisor. If
all bands are empty, return \(N\); even inputs return two. Equation (10)
gives \(2^{n/4}\mathrm{poly}(n)\) time and polynomial bit space.

The rectangle wrapper first caps the positive row interval at \(N\), splits
it into dyadic windows, and checks the partner \(N/a\) against the column
interval. It can return the count and smallest row using polynomial space.
An explicit `collect_pairs=True` option retains every output pair and
therefore additionally uses space proportional to that output. There are
no modular wraparound candidates: each reported pair satisfies \(ab=N\)
as integers. Arbitrary positive rectangle endpoints are allowed, with their
written bit lengths charged in the interval comparisons.

## 5. Exact controls

Run:

```bash
python tools/check_circuit_box_centered_windows.py \
  --output evidence/circuit-box-centered-windows.json
```

The completed checker verifies 120,169 universal integer budgets through
512-bit parameters and selected parameters through 16,384 bits, including
all seven exceptional small cases. It compares 160 local windows, 80
arbitrary positive rectangles, and 999 smallest-factor outputs with literal
integer arithmetic. Its 37 independent certificate checks rebuild the
lattice rows and retained short polynomial and verify divisibility at every
known factor. Five invalid inputs are rejected.

The large local controls use independently supplied factorizations
\(N=3^u5^v\). They include a 64-bit full dyadic band, a 255-bit input with a
window of width \(2^{32}\), and a 1,019-bit input with a window of width
\(2^{150}\). The two short-window queries require no low bits and include
negative centered roots. Their expected complete factor lists are computed
from all \((u+1)(v+1)\) supplied prime-power products, not from scanning
the queried windows. These are local checks, not complete factorizations
of unknown large inputs.

| Input bits | Low bits used | Old nonempty odd classes combined | Lattice dimension | Sturm queries |
| --- | ---: | ---: | ---: | ---: |
| 64 | 15 | 8 | 19 | 15 |
| 255 | 0 | \(2^{31}\) | 3 | 34 |
| 1,019 | 0 | \(2^{149}\) | 3 | 152 |

The 17-bit complete example \(67591=257\cdot263\) takes 48 residue queries
and 25 lattice queries. The 19-bit example \(271421=269\cdot1009\) takes
71 residue queries and 32 lattice queries. The preceding fixed-precision
wrapper used 512 residue queries on each example. These are exact query
counts; no wall-clock speedup is inferred from them.

The [receipt](../evidence/circuit-box-centered-windows.json) records every
large window, its exact candidate count, the number of old nonempty odd
residue classes combined, the complete lattice certificate, and the
charged root-isolation work. It separately records the two completed
17-bit and 19-bit smallest-factor examples. Source hashes are checked
before and after the suite. No numerical frequency sum or uncharged
Fourier enumeration is used by this endpoint.
