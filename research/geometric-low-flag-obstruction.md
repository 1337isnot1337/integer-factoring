# The low image flag is rank-one-free in every characteristic-two field

The eleven-dimensional low image flag of the retained five-letter
crossing contains no nonzero rank-one matrix over **any field of
characteristic two**, including an algebraic closure. This strengthens
the earlier exhaustive `F_4` check. It follows from 25 explicit
polynomial identities over `F_2`, each of total degree at most three.

Consequently the [single-pure-term two-copy floor](pure-two-copy-flag-floor.md)
extends to all characteristic-two coefficient fields:

\[
 \rho\ge\frac57>\frac1{\sqrt2}
\]

for one physically pure rational correction to the retained tensor
square at its fixed natural physical heights. The correction may
have arbitrary Laurent leading orders and may correlate the two
copies within each physical vector.

This extension concerns the rank-one flag obstruction. It does not
extend the finite `F_4` rank-twelve computation used for the
[eight-term floor](low-rank-two-copy-flag-floors.md) to other fields.
Different physical heights and general higher-rank corrections
remain separate questions.

## 1. The exact flag and its annihilator

Let `U` be the span of the eleven constant physical image vectors
with relative Smith order at most nine, reshaped as five-by-five
matrices. Its basis and reference digest are retained in the
[certificate](../evidence/geometric-low-flag-obstruction.json).
The two image flags and two coimage flags have this same span.

Fourteen independent binary linear forms annihilate `U`. Encode
a five-by-five form by assigning bit `5i+j` to matrix entry `(i,j)`.
In the order used below, the forms are

```
16916481, 4210722, 16388, 13156360, 1187856, 64, 1056896,
1450240, 524800, 4195328, 1058816, 4096, 1376256, 10485760
```

The checker verifies their rank is fourteen, the image basis has
rank eleven, and all 154 pairings are zero. Their common kernel
is therefore exactly `U`. These rank statements remain true over
every characteristic-two field: the same nonzero binary minors
remain nonzero after extending scalars.

A rank-one matrix is `a b^T`. For each annihilating form `H_i`,
write its corresponding bilinear polynomial as

\[
 f_i(a,b)=\sum_{r,c=0}^4 (H_i)_{rc}a_r b_c.
\]

Thus a nonzero rank-one matrix belongs to `U` exactly when two
nonzero vectors `a,b` satisfy all fourteen equations `f_i=0`.

## 2. Twenty-five ordinary polynomial identities

Every nonzero vector over a field has a first nonzero coordinate.
Let those indices be `p` for `a` and `q` for `b`. Scale each vector
independently so that

\[
 a_0=\cdots=a_{p-1}=0,\quad a_p=1,
 \qquad
 b_0=\cdots=b_{q-1}=0,\quad b_q=1.                       \tag{1}
\]

The 25 choices of `(p,q)` cover every possible nonzero pair.
On every such chart, the certificate gives explicit polynomials
`g_0,...,g_13` such that

\[
 \boxed{\sum_{i=0}^{13}g_i f_i=1}                       \tag{2}
\]

after the substitutions (1). The total degree of every displayed
identity is at most three. Across all charts, the multipliers
contain only sixty monomials.

For example, on the chart `a_0=b_0=1`, the identity is

\[
 \begin{aligned}
 1={}&f_0+b_3f_2+b_2f_3+b_4f_6+b_2f_8\\
     &+(b_2+b_3)f_9+(1+b_4)f_{11}.
 \end{aligned}                                        \tag{3}
\]

Here each `f_i` has already undergone the chart substitutions.
Several other charts have a constant equation `f_i=1`, so their
certificate consists of a single term.

These are identities in an **ordinary polynomial ring over `F_2`**.
Variables retain their full exponents. No equations such as
`x^2=x`, `x^4=x`, or `x^16=x` are imposed. If all `f_i` vanished
at a point over any characteristic-two field, (2) would evaluate
to `0=1`. Hence no such point exists. This proves the geometric
rank-one exclusion, without inferring it from a finite-field scan.

The constructor finds (2) by exact binary elimination among
monomial multiples of the chart equations. The independent checker
does not rerun that search: it expands the retained multipliers
and equations literally in ten named variables and verifies that
their sum is the constant polynomial one. It also rejects a
deliberately damaged identity.

## 3. Finite-field and boundary controls

As a separate control, the checker also fixes each projective
left vector `a` and solves the fourteen homogeneous equations
for the five coordinates of `b`.

| Coefficient field | Projective left factors checked | Rank of every right-factor system |
|---|---:|---:|
| `F_2` | 31 | 5 |
| `F_4` | 341 | 5 |
| `F_16` | 69,905 | 5 |

The `F_16` arithmetic uses the quotient by `x^4+x+1`; the checker
constructs multiplication and verifies that every nonzero element
has an inverse. These finite checks corroborate the polynomial
identities but are not needed to infer the result over an
algebraic closure.

The order cutoff is sharp for rank one. At relative order ten,
the twelve-dimensional flag already contains the binary outer product

\[
 (0,1,1,1,0)^T(0,0,1,0,0).
\]

The checker verifies that this matrix lies in the larger flag and
does not lie in `U`. Thus the preceding argument has not accidentally
excluded all rank-one matrices from every flag.

## 4. The pressure consequence

The single-pure-term theorem uses the rank-one exclusion in `U`,
its preservation under the grouped tensor square, a regular-row
Toeplitz argument, and rank-one valuation interlacing. Every one
of those arguments is valid over an arbitrary field. The only
previous restriction to `F_4` came from the finite rank-one test.

The grouped tensor square `U tensor U` still contains no nonzero
rank-one matrix: a nonzero rank-one element would have a nonzero
rank-one coefficient contraction belonging to `U`. The higher-jet
argument then gives, at each two-copy place,

\[
 Z_{\rm new}(s)\ge P(s)^2-s^{10}.
\]

The reference profiles also remain unchanged under extension of
the binary coefficient field. The exact integer pressure
certificate from the single-term theorem therefore applies
without alteration. Its sixty dyadic intervals and analytic small
interval give the same `5/7` floor over every characteristic-two
coefficient field.

## 5. Reproduction

```sh
python3 tools/geometric_low_flag_obstruction.py
python3 tools/check_geometric_low_flag_obstruction.py
```

The [constructor](../tools/geometric_low_flag_obstruction.py),
[literal checker](../tools/check_geometric_low_flag_obstruction.py),
[polynomial identities](../evidence/geometric-low-flag-obstruction.json),
and [independent controls](../evidence/geometric-low-flag-checks.json)
are all retained. The factoring exponent claim is unchanged.

This research and its implementation used substantial OpenAI assistant support.
