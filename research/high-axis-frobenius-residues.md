# Frobenius residue coefficients change the higher axis jets

## Registered structural construction

Base: `5b20436`. The earlier degree-64
[Plucker kernel](high-axis-plucker-kernel.md) changed the quadratic rank,
but a higher-jet minor still excluded its rate. Keep the same nine axes,
degree 64, antipodal endpoint and quadratic coefficient. Change the
five- and six-letter coefficients using Frobenius identities, so that a
specific higher exterior flag has smaller generic rank. This is a new
coefficient relation rather than a degree, field or sampling extension
of the closed source.

For each nonempty subset S of size k, write `64=k*e+r`. Choose the largest
power of two f dividing r such that `binom(k,r/f)` is odd, and put

```
c_S=(product x_a)^e * elementary_(r/f)((x_a^f)_(a in S)).
```

When r=0 use the monomial `(product x_a)^e`. The choice exists for every
k=1,...,9. Only the following coefficients change from the preceding kernel:

```
k=5: c_S=(product x_a)^12 * sum x_a^4,
k=6: c_S=(product x_a)^10 * elementary_2((x_a^2)).
```

Both are homogeneous of degree 64 and equal one at the all-one point,
since `binom(5,1)` and `binom(6,2)` are odd. The same Cayley construction,
with `c_empty=sum_(S nonempty)c_S`, therefore has its exact endpoint and
zero diagonal. Every axis matrix squares to `x_a^128 I`.

## The changed generic quartic flag

At the a-axis work in the exterior algebra on the other eight symbols
epsilon_b over the rational function field. Write

```
u_j=sum_b x_b^j epsilon_b,
E_l(j)=elementary_l((x_b^j epsilon_b)).
```

The five-letter leading quartic coefficient becomes

```
sum_(|B|=4) (product_(b in B)x_b^12)*(sum_(b in B)x_b^4)*epsilon_B
    = u_16 E_3(12) = u_16 u_12 E_2(12).
```

The six-letter leading quintic coefficient becomes

```
E_2(12) E_3(10) = E_2(12) u_10 E_2(10).
```

These are coefficient identities in characteristic two, with no numerical
specialization. In particular the quartic now contains two independent
linear factors. After scaling epsilon_b by x_b^12, the quadratic E_2(12)
is the nondegenerate eight-variable alternating form with all off-diagonal
entries one. The symplectic pairing of u_12 and u_16 is
`sum_b x_b^4`, nonzero generically. A symplectic change of basis therefore
reduces the quartic multiplication to two fixed exterior factors times a
nondegenerate quadratic in six remaining variables. Its rank is
`(2^6-2^3)/2=28`: the three two-variable complexes have two-dimensional
homology each, and the tensor differential has homology dimension eight.
The old quartic E_4(13) has rank 46, as an exact binary constant matrix
after invertible variable scaling. The quadratic jet of the preceding
construction remains unchanged. Smaller rank in this higher flag is a
mechanism, not yet a complete pressure improvement.

## Diagnostic and stopping rule

Construct complete specialized profiles for at most the three previously
used nonzero F4 patterns, starting with all ones. A favorable specialization
still cannot certify the generic multivariate profile. After each query,
test whether some prefix contains `m>256/3` orders with average at most
`64/3`. If so, independently compute its literal binary Toeplitz ranks
(using scalar restriction for a nonbinary specialization), verify every
rank against the proposed complete profile, and stop: the same retained
channel/Jensen argument excludes this new kernel at every scale and common
occupied frequency. A single nonzero specialized minor suffices.

Only if that obstruction is absent in all three patterns, and a specialized
fit is below `1/3`, proceed to symbolic generic divisibility or low-rank
factorization and rational outward pressure. Keep failed/unresolved queries.
No neighboring coefficient, field, degree or copy-count sweep follows.
The complete factoring score remains `1/2`; no stronger bound, Lean work
or independent mathematical review is claimed. Substantial OpenAI assistant
support, retaining the preceding Cayley, transfer and Smith-tool credits.

## Exact outcome and all-scale stop

The first all-one specialization already supplies the registered stopping
certificate, so the two later patterns are not run. Its complete native
profile has an exploratory minimum near `0.4016937058`, still above `1/3`.
That numerical fit is not used in the exclusion. The
[diagnostic archive](../evidence/high-axis-frobenius-residues.json) retains
all 256 specialized orders and the explicit early-stop status.

The [independent checker](../tools/check_high_axis_frobenius_residues.py)
computes the literal binary Toeplitz ranks

```
R(28)=840,       R(29)=960.
```

Thus exactly `m=120` orders are at most 28, with total
`29*120-960=2520`. Specialization gives the generic bound
`Delta_120<=2520`, simultaneously at every axis by symmetry. Since
`120>256/3` and `2520/120=21<64/3`, every positive scale satisfies

```
log Z_a(t) >= 120 log(1+t^21)
           >= (256/3) log(1+t^(64/3)),        0<t<1.
```

The same binary entropy inequality as in the preceding kernel gives
`R_9(p,t)>=1/3` for every scale and common occupied frequency. The binary
pivots persist over every characteristic-two extension field. The other
136 channels are retained as nonnegative extra terms; no generic full
profile or favorable tail is asserted.

The [receipt](../evidence/high-axis-frobenius-residue-controls.json) also
checks 280 symbolic quartic and 560 symbolic quintic monomials against the
two factor identities, old quartic rank 46 and new rank 28, all 13,122
nonzero-pattern entries at two field fixtures, and 96 complete untruncated
small profiles. The rank-28 lower fixture has a nonzero symplectic pairing;
the preceding algebra gives the matching generic upper bound. All 29
Toeplitz ranks are independent of native Smith arithmetic.

The higher flag change is real but insufficient. This particular kernel
stops; no coefficient or degree sweep follows. Run
`make verify-high-axis-frobenius-residues`. The complete score remains
`1/2`. The broader question is which cross-axis compatibility condition
forces the low-prefix obstruction, rather than another adjustment to a
specialized profile.
