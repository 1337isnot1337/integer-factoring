# A strictly better 96-letter coherent crossing

## Certified result

An explicit 96-letter physical retraction has the exact strict pressure
bound

\[
\boxed{\rho<\frac{7100058652501}{10^{13}}
             =0.7100058652501.}
\tag{1}
\]

This improves the earlier coordinate 96-letter certificate
`rho < 0.710005874881191`. Both complete cuts have 9,216 Smith orders.
The two coherent profiles differ, and both are used in the certificate.

The new map improves the coordinate product pressure for every
`0 < r <= 4/5`. It does not dominate at every scale: the coordinate map
is better sufficiently near `r=1`. This is a small constructive gain for
the stated alphabet, not a factoring exponent above one half. The bound
in (1) remains above `1/sqrt(2)`.

## 1. The physical map

Use twelve physical modes with heights `1,3,...,23`. A physical occupation
word `A` has height `h(A)=sum_(a in A)(2a+1)`. Retain all 92 coordinate
states of height less than 25, together with the five-particle word `31`,
which has height 25.

Order the ten three-particle words of height 25 by their integer masks.
On the first four of them, write the standard coordinate vectors as
`e_0,e_1,e_2,e_3`. Append three embedding columns

\[
P_1=e_0+e_1,\qquad P_2=e_0+e_2,\qquad P_3=e_0+e_3,
\]

and the paired retraction rows

\[
L_1=e_0+e_2+e_3,\qquad
L_2=e_0+e_1+e_3,\qquad L_3=e_3.
\tag{2}
\]

Their dot products give `LP=I`. Every derivative of a height-25 boundary
state lies in the fully retained lower-height part. Hence both map spaces
are closed under every bit deletion. All paired labels preserve height
and occupation. The complete compound reduction therefore applies to
both cuts, with their distinct compressions.

The coordinate comparison retains the first three of those boundary
triples instead. Its physical height multiset is identical. The
[constructor](../tools/compute_coherent_q96_candidate.py) records the exact
boundary order, embedding, retraction, source hashes, and all complete
block profiles. The [coherent record](../evidence/coherent-q96-four-boundary-candidate.json)
contains the finished calculation. The independent generic constructor
produced the [coordinate reference](../evidence/coherent-q96-coordinate-reference.json).

## 2. Exact rational pressure witness

Put

\[
p=7100058652501,\qquad q=10^{13},\qquad
a=27855864849641,\qquad x=\frac{a-1}{a},\qquad r=x^p.
\]

For the complete coherent profiles `nu_i` and `mu_i`, define

\[
Z_0(r)=\sum_i r^{\nu_i},\qquad
Z_\infty(r)=\sum_i r^{\mu_i},\qquad
A(y)=\sum_j y^{h_j}.
\]

The exact certificate proves

\[
\frac{(Z_0(r)Z_\infty(r))^q}{A(x^q)^{4p}}
 <\frac7{20000}<1.
\tag{3}
\]

This is precisely the strict balanced pressure inequality at `rho=p/q`.
Every physical probability is rational:

\[
\pi_j=\frac{x^{q h_j}}{\sum_k x^{q h_k}}.
\tag{4}
\]

Finite word types can realize (4) after a common denominator is cleared.
These large fixed denominators affect constants, not the asymptotic
exponent of the existing finite-alphabet transfer. No irrational type
or unspecified limiting choice is used.

The coordinate profile fails the same rational rate and scale: its
powered ratio is strictly larger than one. This is a direct comparison
of actual complete profiles, in addition to the comparison with the
previously displayed coordinate bound.

The large exponents in (3) do not require materializing large integers.
The existing outward interval engine stores a 256-bit mantissa and a
separate binary exponent. Repeated squaring needs a logarithmic number
of interval multiplications in each displayed exponent. The final
coherent ratio lies below `7/20000` with a large strict margin. All
decisions use integer arithmetic; numerical scale optimization only
suggested the explicit integers above.

## 3. A whole interval of strict improvement

Let `C_0,C_infinity` denote the coordinate counting polynomials. Exact
integer multiplication gives

\[
C_0(r)C_\infty(r)-Z_0(r)Z_\infty(r)
 =r^{63}(1-r)Q(r),\qquad \deg Q=176.
\tag{5}
\]

The scaled polynomial `5^176 Q(4x/5)` has a nonnegative expansion

\[
5^{176}Q(4x/5)
=\sum_{i=0}^{176}b_i x^i(1-x)^{176-i},
\qquad b_i\ge0,
\]

with both endpoint coefficients positive. The checker independently
expands this integer identity. Thus (5) is strictly positive throughout
`0 < r <= 4/5`. The rational witness in section 2 lies inside that
interval.

At the other endpoint,

\[
Q(1)=-8441856<0.
\]

Continuity therefore makes (5) negative sufficiently near one. The
construction supplies no all-scale domination claim. This distinction
matters if one tries to extend the coherent map to larger alphabets or
to a regime where the favorable counting variable `r` tends to one.

## 4. Completeness and replay

The coherent calculation resolves all eleven occupation blocks at each
cut, with series precision 128. The largest resulting zero and infinity
orders are 116 and 109. The exact polynomial endpoint proves the generic
rank of each block independently of its precision window. Missing pivots
would cause a retry or failure in the constructor.

The [pressure checker](../tools/check_coherent_q96_pressure.py) verifies
the exact maps, 4,608 deletion conditions, all block labels, and the
coverage of all 36,864 orders across the two maps and two cuts. It checks
the source hashes of the completed computations, fourteen independent
small rational interval controls, (3), and the Bernstein identity in
section 3. Its [receipt](../evidence/coherent-q96-pressure.json) retains
the full integer coefficients and interval mantissas.

The default checker reads the already completed native profiles. It does
not claim to recompute them. A separate explicit option repeats both
expensive native calculations:

```sh
python3 tools/check_coherent_q96_pressure.py \
  --output /tmp/coherent-q96-pressure.json
python3 tools/check_coherent_q96_pressure.py --recompute-profiles \
  --output /tmp/coherent-q96-pressure-full-replay.json
```

No global optimization assertion or new complete factoring score follows
from the numerical search values. The exact gain is the finite physical
crossing certificate (1), with the precise scale interval (5).
