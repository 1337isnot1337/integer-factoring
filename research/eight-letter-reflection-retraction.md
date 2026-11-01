# An eight-letter nonorthogonal reflection crossing

## Result and scope

An explicit homogeneous linear retraction of the sixteen-letter reflection
crossing gives a polynomial eight-letter crossing with a strict common-type
certificate

\[
\rho<\frac{18}{25}=0.72.
\tag{1}
\]

The certificate uses a specified rational type and the exact scale
`s=(73/71)^18`. At that same scale, the coordinate restriction with the same
physical heights fails the trial value `rho=18/25`. More generally, the new
crossing has strictly smaller two-cut pressure than that coordinate
restriction for every finite `s>=3/2`.

This is a constructive improvement within an eight-letter family. It does
not improve the factoring exponent: the certified trial rate in (1) remains
above `1/sqrt(2)`. No global optimality is asserted for either restriction,
and the unrestricted sixteen-letter reflection seed remains numerically
stronger. The construction and exact profiles work over every field of
characteristic two.

## 1. The ambient reflection crossing

Use the `M=8` crossing from
[reflection-pair compression](reflection-dual-rail-compression.md). All its
polynomial exponents are even. Remove that common ramification, so the raw
polynomial tensor `T(z)` has degree 48 and physical dimension sixteen.
Write a physical word as four bits and set

\[
H(a)=7a_0+5a_1+3a_2+a_3.
\tag{2}
\]

Its two regular endpoint tensors are

\[
Z_{xyuv}(z)=z^{H(x)+H(y)-16}T_{xyuv}(z),
\qquad
W_{xyuv}(w)=w^{H(u)+H(v)+16}T_{xyuv}(1/w).
\tag{3}
\]

The normalized kernels satisfy the exact physical identity

\[
W_{xyuv}=Z_{xuyv}.
\tag{4}
\]

The constructor forms every Gaussian coefficient as the Pfaffian of the
corresponding principal kernel submatrix. It verifies the kernel version of
(4), all 65,536 literal ambient swap entries at `z=1`, and every polynomial
restoration used below. Thus the common scalar shifts in (3) are retained
through the computation.

## 2. A nonorthogonal retraction on one equal-height pair

Retain the ordered coordinate labels

\[
C=(0,1,2,4,9,8,10,12),
\qquad h=(0,7,5,3,8,1,6,4).
\tag{5}
\]

Let `P` have columns `e_c` in this order. The coordinate retraction is
`L_0=P^T`. Change only its row corresponding to the retained label 9:

\[
L_4=e_9^{\mathsf T}+e_6^{\mathsf T},
\qquad L_j=e_{C_j}^{\mathsf T}\quad(j\ne4).
\tag{6}
\]

The omitted label 6 has height eight, as does label 9. Hence both maps
preserve the physical height. Since 6 is absent from `C`, they satisfy
`LP=I_8`. Define

\[
S=(L\otimes L)T(P\otimes P).
\tag{7}
\]

The same retraction on both physical wires preserves the literal swap:

\[
S(1)=(L\otimes L)\mathrm{swap}(P\otimes P)
     =\mathrm{swap}.
\tag{8}
\]

The complete compressed tensor is polynomial of actual degree 32. Its
regular tensors are obtained by the same homogeneous contraction of `Z`
and `W` in (3). The checker also reconstructs the raw tensor from its zero
normalization, reverses it at its actual degree, and checks all 4,096
infinity entries independently.

The two normalized Smith profiles need not agree: realignment changes how
`L` and `P` act on the cut. Compute both complete cuts. If their orders are
`a_i` and `b_i`, the cancelling scalar shifts in (3) give

\[
B(s)=\left(\sum_i s^{-a_i}\right)
     \left(\sum_i s^{-b_i}\right),
\qquad G_j=4h_j.
\tag{9}
\]

The coordinate control uses `L_0`, the same embedding, and the same height
multiset. Its two normalized profiles coincide.

## 3. Exact profiles and a uniform pressure gain

Let `r=1/s` and let `c(r)` be the coordinate counting polynomial. Its
coefficients, in ascending degree from 0 through 32, are

```text
[1,1,0,1,1,1,1,1,2,2,2,2,3,3,3,4,5,4,4,5,3,3,3,2,1,1,1,0,1,1,0,1,1].
```

For the nonorthogonal restriction, the zero and infinity counting
polynomials are exactly

\[
\begin{aligned}
z(r)={}&c(r)-r^{16}-r^{17}-r^{19}-r^{20}
              +r^{24}+r^{25}+r^{27}+r^{28},\\
i(r)={}&c(r)+r^{20}+2r^{21}+r^{22}+r^{23}+r^{24}\\
       &\hspace{12mm}-r^{25}-r^{26}-r^{28}-r^{29}-r^{31}-r^{32}.
\end{aligned}
\tag{10}
\]

Each polynomial has 64 terms counted with multiplicity. The zero cut
improves some orders, while the infinity cut loses some; the product is
the relevant quantity.

Direct integer multiplication gives

\[
c(r)^2-z(r)i(r)=r^{16}(1-r)Q(r).
\tag{11}
\]

Here `Q(0)=1`, its first negative coefficient has degree 23, and the sum
of the absolute values of all its negative coefficients is 749. Therefore

\[
Q(r)\ge1-749(2/3)^{23}>0
\qquad(0\le r\le2/3).
\tag{12}
\]

The final strict inequality is checked as `749*2^23<3^23`. Equations
(11)--(12) prove the stated strict pressure gain for every finite
`s>=3/2`. This interval contains the exact certificate scale below.
The full coefficient list of `Q` is included in the evidence.

## 4. An explicit rational common type

Put

\[
s=(73/71)^{18},\qquad \rho_0=18/25,
\qquad
A=\sum_{j=0}^{7}(71/73)^{25h_j}.
\tag{13}
\]

The common type is the rational probability vector

\[
p_j=\frac{(71/73)^{25h_j}}{A}.
\tag{14}
\]

Its entries are positive, so its entropy is nonzero. Because all entries
are rational constants, a common type block can realize them exactly by
taking a multiple of their common denominator. The asymptotic type
overhead is unchanged.

For the common-type pressure criterion, (14) is the Gibbs type at this
scale and trial rate. Indeed, writing `H(p)` for entropy in natural units,

\[
4\rho_0H(p)=4\rho_0\log A+4\log(s)\sum_jp_jh_j.
\tag{15}
\]

Thus it suffices to establish `B(s)<A^(72/25)`. No numerical logarithms
are needed: the checker verifies

\[
\bigl[z((71/73)^{18})i((71/73)^{18})\bigr]^{25}<A^{72}
\tag{16}
\]

after clearing all denominators. Every operand is an explicitly
constructed positive integer. The same comparison with `c(r)^2` has the
opposite strict sign. Consequently this exact type and scale certify
(1) for the mixed crossing, while the coordinate control fails that
same trial rate at that scale. This is not an assertion that the
coordinate control fails the trial rate at all possible scales or
under different gauges.

The exact rational vector, integer-comparison recipe, bit lengths, and
hashes of the compared integers are stored with the profiles. Floating
optimization guided the choice of scale during discovery; it plays no
part in the certificate.

## 5. Replay and limits

The [constructor](../tools/eight_letter_reflection_retraction.py) accepts
general homogeneous binary embeddings and retractions of this ambient
seed, checks `LP=I`, and returns both complete profiles. The
[checker](../tools/check_eight_letter_reflection_retraction.py) evaluates
the coordinate control and (6), reconstructs all endpoint entries, and
checks 28 independent truncated Toeplitz ranks, including precision 33
beyond every displayed Smith order.

```sh
python3 tools/check_eight_letter_reflection_retraction.py \
  --output evidence/eight-letter-reflection-retraction.json
```

The [evidence](../evidence/eight-letter-reflection-retraction.json)
contains all 256 orders, all maps, and every exact comparison. The
coefficients and minors lie in the prime field, so extending the
coefficient field within characteristic two changes none of these
profiles or identities.

The useful construction change is the nonorthogonal retraction in a
degenerate physical-height space. Its benefit is established by the
full product of two differently changing cuts. No conclusion is drawn
about all linear restrictions, other heights, or the global optimum of
the sixteen-letter seed, and the proved factoring exponent stays at
`alpha=1/2`.
