# Shared-minor budgets for graded cube alphabets

The [homogeneous binary floor](cube-homogeneous-profile-budget.md)
extends to a computable necessary profile budget for every finite
graded alphabet. The budget uses actual minors shared between the
three cuts. It does not assume a Gaussian tensor, independent modes,
or a particular physical basis inside an energy level.

Two exact consequences give a useful boundary for further construction:

- Every homogeneous cube tensor with physical heights \(0,1,3,4\)
  has common-scale, unweighted profile pressure strictly above
  \(7/12>1/\sqrt3\), for every common physical type.
- The same necessary budgets allow an explicit integer abstract profile
  on \(256\) letters with pressure below
  \(72/125<1/\sqrt3\). This is a target for construction, not an actual
  tensor or a new factoring exponent.

The [exact envelope module](../tools/cube_graded_minor_budget.py)
stores energy multiplicities and convex-hull segments. It can therefore
handle the \(256^4=4294967296\) profile dimension in the second statement
without enumerating that many channels.

## 1. A general shared-prefix theorem

Let \(A\) be a physical alphabet with \(q\ge2\) letters. Give letter
\(a\) a nonnegative integer energy \(\varepsilon(a)\), and choose a
zero-energy letter \(0\). Suppose a common polynomial cube tensor

\[
 T_{a_{000},\ldots,a_{111}}(x,y,z)
\]

is homogeneous of degree \(\sum_v\varepsilon(a_v)\), and at
\((1,1,1)\) equals the product of the four antipodal equality tensors.
All three axis flattenings have dimension \(q^4\) and identity
endpoint, after a column permutation. Their ascending Smith orders
are \(s_{a,1},\ldots,s_{a,q^4}\), for axis \(a\). Define

\[
 \gamma_{a,k}=\sum_{i\le k}s_{a,i},\qquad
 \bar\gamma_k=\frac13\sum_a\gamma_{a,k}.
\]

For \(r\in\{1,2,4\}\), let \(S_r(k)\) be the sum of the \(k\)
smallest energies among words of length \(r\), with multiplicity.
Then

\[
 \boxed{
 \begin{aligned}
 \bar\gamma_k&\le\frac23S_1(k) &&(k\le q),\\
 \bar\gamma_k&\le S_2(k)      &&(k\le q^2),\\
 \bar\gamma_k&\le2S_4(k)     &&(k\le q^4).
 \end{aligned}} \tag{1}
\]

**Proof.** Choose a collection of antipodal pairs, allow their
representative letters to range over a set \(X\), and use the paired
letters as the column configurations. Fix every remaining letter to
\(0\). The selected minor has identity endpoint, so its determinant
is nonzero. Every determinant term has degree

\[
 2\sum_{w\in X}\varepsilon(w).
\]

One antipodal pair is shared by all three cuts. The sum of the three
coordinate valuations of its determinant is at most its homogeneous
degree. This gives \(\sum_a\gamma_{a,k}\le2S_1(k)\).

For each pair of axes, two antipodal pairs have representatives on
the zero side of both cuts. The same selected minor occurs in both
cuts. Thus \(\gamma_{a,k}+\gamma_{b,k}\le2S_2(k)\). Add the three
axis-pair inequalities and divide by six.

On a single cut all four antipodal pairs are available, giving
\(\gamma_{a,k}\le2S_4(k)\) directly. This proves (1).

The choice of \(X\) is arbitrary. Restricting it to full Cartesian
products loses useful bounds at intermediate dimensions.

## 2. The exact convex envelope

An ascending order vector has a convex sequence of cumulative sums.
Take the lower convex hull of all the points on the right of (1),
including \((0,0)\). Let its successive slopes be \(\sigma_i\),
repeated for the integer length of each hull segment.

Every feasible \(\bar\gamma_k\) lies below this hull. Indeed, convexity
bounds it by the chord between any two applicable endpoint bounds.
The slopes are therefore optimistic orders: for every \(0<t<1\),
weak majorization and convexity of the logarithmic exponential sum
give

\[
 \frac13\sum_a\log\left(\sum_i t^{s_{a,i}}\right)
 \ge\log\left(\sum_i t^{\bar s_i}\right)
 \ge\log Z_*(t),
 \qquad Z_*(t)=\sum_i t^{\sigma_i}. \tag{2}
\]

Here \(\bar s_i=\frac13\sum_a s_{a,i}\). The hull need not be attained
by three actual order vectors, much less by three cuts of one tensor.

For efficient computation, let \(P(z)=\sum_a z^{\varepsilon(a)}\).
The coefficient multiplicities of \(P,P^2,P^4\) give all needed
energy counts. Within a run of equal energy, \(S_r(k)\) is linear.
It is therefore enough to retain the two ends of each run: convexity
supplies every interior bound. The work depends polynomially on the
number of distinct input energies and their binary encodings, with
the tensor power fixed at four. Large alphabet multiplicities enter
only as integers.

## 3. A four-letter exclusion beyond Gaussian kernels

Take physical heights \(h\in\{0,1,3,4\}\), and tensor degree
\(3\sum_v h_v\). This includes the full occupation alphabet of two
physical modes with heights one and three, as well as arbitrary
non-Gaussian tensors with the same grading.

The exact envelope has the following partition polynomial:

\[
\begin{aligned}
 Z_*(t)={}&1+t^2+t^4+t^6+2t^9+4t^{12}+2t^{15}+t^{18}
           +6t^{19}\\
 &+17t^{24}+24t^{30}+22t^{36}+28t^{42}+36t^{48}\\
 &+28t^{54}+22t^{60}+24t^{66}+17t^{72}
           +8t^{78}+6t^{84}+4t^{90}+t^{96}.
\end{aligned} \tag{3}
\]

Its coefficients sum to \(4^4=256\). Let a common physical type have
probabilities \(p_h\), mean height \(E=\sum_h hp_h\), and entropy
\(H(p)=-\sum_h p_h\log_2p_h>0\). Homogeneity gives total common
coefficient degree \(24mE\) in a word tensor power. The profile
pressure consequently obeys

\[
 \rho_{\mathrm{prof}}\ge
 \frac{\log_2 Z_*(t)-8E\log_2t}{4H(p)}. \tag{4}
\]

Set \(W(v)=1+v+v^3+v^4\). The checker forms the integer polynomial

\[
 P(s)=Z_*(s^7)^3-W(s^{24})^7.
\]

Its first nonzero degree is \(14\). After dividing by \(s^{14}\),
the quotient has degree \(2002\). Every coefficient of

\[
 (1+u)^{2002}
 \frac{P\bigl(u/(1+u)\bigr)}{\bigl(u/(1+u)\bigr)^{14}}
\]

is a positive integer. There are \(2003\) coefficients and their
minimum is \(3\). The receipt records the coefficient digest; the
checker reconstructs and verifies every coefficient exactly.

Therefore

\[
 Z_*(t)>W(t^{24/7})^{7/3}\qquad(0<t<1).
\]

The entropy identity

\[
 H(p)+E\log_2v\le\log_2W(v)
\]

now gives \(\log_2 Z_*(t)>\frac73H(p)+8E\log_2t\).
Substitution into (4) proves \(\rho_{\mathrm{prof}}>7/12\).
The threshold comparison is exact: \(3\cdot7^2>12^2\).

This rules out the half-exponent target for the stated four-letter
grading and profile method, even after adding homogeneous interactions
between the two modes. It does not rule out other gradings or larger
alphabets.

## 4. An exact larger-alphabet target remains feasible

Now take one vacuum letter and \(255\) distinct nonzero letters, all
of height one. As before, entries have degree three times total
physical height. The following ascending integer abstract profile
obeys every inequality in (1):

| Order | Multiplicity |
|---:|---:|
| 0 | 1 |
| 2 | 255 |
| 4 | 255 |
| 6 | 65025 |
| 13 | 262650 |
| 14 | 62985 |
| 18 | 66325500 |
| 24 | 4228250625 |

Its total multiplicity is \(256^4\). To obtain it, round each constant
slope segment of the convex envelope to its two adjacent integers,
preserving the segment's total order, then sort the complete vector.
The checker verifies every compressed raw prefix constraint after
rounding. Convexity supplies all intervening prefixes. This is an
integer feasible profile for the relaxation, not a claim that integer
rounding gives a universal profile lower bound.

Write \(Z(t)\) for the polynomial specified by the table. Choose

\[
 t=\frac6{13},\qquad p=\frac{19}{20},\qquad R=\frac{72}{125},
\]

where \(p\) is the total nonzero probability and each nonzero letter
has probability \(p/255\). Its entropy is
\(H_2(p)+p\log_2 255\). The checker verifies exactly that

\[
 Z(t)^{2500}t^{-19000}
 <\left(\frac{20^{20}255^{19}}{19^{19}}\right)^{288}. \tag{5}
\]

Taking logarithms and dividing by the positive entropy gives

\[
 \frac{\log_2Z(t)-8p\log_2t}
      {4\bigl(H_2(p)+p\log_2 255\bigr)}
 <\frac{72}{125}<\frac1{\sqrt3}.
\]

The second strict inequality is the integer comparison
\(3\cdot72^2<125^2\). Equation (5) is checked by clearing rational
denominators and comparing integers; no numerical optimization is
part of the certificate.

**No tensor realizing this profile is known here.** The point of the
witness is that these shared-minor and degree constraints alone do
not close the large, equal-height alphabet route. A successful
construction must still produce one common multivariate tensor,
certify all three generic cut profiles, extract common coefficient
pieces, and pay the complete evaluation and factorization costs.

## 5. Exact replay

```sh
python3 tools/check_cube_graded_minor_budget.py \
  --output /tmp/cube-graded-minor-budget.json
```

The checker compares compressed envelopes with every literal prefix
for forty small alphabets, verifies all \(2003\) positive coefficients
in the four-letter exclusion, and replays the complete integer
large-alphabet witness. The
[receipt](../evidence/cube-graded-minor-budget.json) records the exact
profiles, parameters, comparison hashes, and construction status.
The established factoring exponent remains \(1/2\).
