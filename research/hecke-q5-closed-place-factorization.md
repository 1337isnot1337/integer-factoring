# Exact Hecke determinants also exclude two closed places

## Result

The five-letter three-row Hecke tensor has matrix order 25 and primitive
degree 36. Exact factorization of both F8 determinant polynomials gives:

```text
cut    determinant degree    irreducible factors    largest degree*multiplicity
 0             814                  17                          200
 1             900                  16                          332
```

Infinity is included when taking each maximum, with degree one.
The resulting degree-weighted budget has positive margin
`25*36-200-332=368`.

Consequently, the unweighted uniform-type supporting criterion cannot
improve the bound using **one irreducible closed place per cut**, at any
positive supporting angle. This includes the full conjugate orbit and
its reconstruction-degree weight. It strengthens the
[two-geometric-point result](hecke-all-algebraic-place-barrier.md) for this
five-letter tensor.

Multiple closed places allocated within one cut, physical-height gauges,
and changes to the tensor remain outside this result. The complete
factoring exponent is unchanged.

## 1. Exact polynomial computation

The checker reconstructs the literal tensor and computes the determinants
of its two 25-by-25 polynomial matrices over
`F8=F2[x]/(x^3+x+1)`. It does not compute over integers modulo eight.

Fraction-free Bareiss elimination uses ordinary polynomial arithmetic.
Every division by the previous pivot is checked exactly. The native
representation has a precision strictly larger than `2*n*D`, which bounds
every intermediate numerator degree. Thus no series truncation is used
in this determinant calculation. The implementation also checks the
minor-degree bounds during elimination.

The finite determinant degrees are 814 and 900. Their homogeneous
infinity multiplicities are therefore 86 and zero. Every F8-rational
root multiplicity agrees with the independently checked local Smith data.

The factorization separates repeated factors in characteristic two,
then uses distinct-degree and equal-degree splitting. Every returned
factor passes the exact finite-field Rabin irreducibility test.
Their powered product, including the leading unit, reconstructs the
entire original determinant coefficient by coefficient. Random choices
are used only to find splits; the resulting identities and irreducibility
conditions are checked exactly.

The complete polynomials and factors are in the
[receipt](../evidence/hecke-q5-polynomial-determinants.json).
The largest nonrational factor degree is 92 in the first cut and 65 in
the second. Every nonrational factor has multiplicity two or four.

## 2. Why the degree-weighted budget is decisive

At an irreducible polynomial `f` of degree `e`, let the local Smith orders
be `nu_j`. Their sum is the multiplicity of `f` in the determinant.
A closed-place CRT condition of order `k` contributes degree `e*k`.
Thus the corresponding unweighted supporting expression, for one closed
place in each cut, is

```math
\Psi(s)=\frac{Ds+a\log\sum_j e^{-s e_0\nu_{0,j}/a}
                    +b\log\sum_j e^{-s e_1\nu_{1,j}/b}}{2\log5},
\qquad a,b>0,\quad a^2+b^2=1.
```

The factorization proves `e_0*sum(nu_0)<=200` and
`e_1*sum(nu_1)<=332`. Each exponentially tilted mean is at most its
uniform mean. Therefore the derivative numerator is at least

```math
36-\frac{200+332}{25}=\frac{368}{25}>0.
```

The minimum occurs as `s` tends to zero, where `Psi(0)=a+b>=1`.
No individual local profile at the large-degree factors is needed for
this conclusion. Their determinant multiplicities suffice.

This uses the degree-weighted single-group convention of the repository's
[mixed-place supporting pressure](../tools/asymmetric_mixed_place_pressure.py).
It does not treat a higher-degree closed place as a free linear place,
and it does not extend this derivative argument to arbitrary allocations
among several closed places in the same cut.

## 3. Independent controls and reproduction

Run:

```sh
python3 tools/check_hecke_q5_determinants.py
```

The checker includes:

- All 4,096 constant two-by-two F8 determinants.
- 48 polynomial determinants of orders one through four, compared with
  full permutation expansions.
- 36 determinant evaluations at independently chosen F64 points outside
  the embedded F8 subfield, using direct field Gaussian elimination.
- Every rational-place multiplicity and both infinity orders, compared
  with the prior Smith/Toeplitz records.
- Exact reconstruction and Rabin tests for all 33 irreducible factors.

The native determinant engine, Python finite-field factorizer, complete
factor lists and source fingerprints are retained. These are ordinary
exact computational certificates, not Lean verification or independent
mathematical acceptance.
