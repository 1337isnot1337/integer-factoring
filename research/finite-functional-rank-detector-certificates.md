# Finite-field functional certificates for observable transfers

Before enumerating auxiliary matrices, one can test whether the
scalar entry space already has a strong linear rank detector.
The detector may be a diagonal matrix whose entries are repeated
linear coefficient functionals. Over a specified finite field,
the required inequalities have a finite exact certificate.

This is a certificate method for the observable word-rank spectral
bound. It does not lower-bound arbitrary grouped-piece true ranks,
and a certificate checked over one finite field does not
automatically apply over its extensions.

## 1. The finite certificate

Fix a field `F_q`, a linearly independent tensor basis
`B_0,...,B_(s-1)`, and the crossing's coefficient vector `e` in
that basis. For a coefficient vector `z`, write `t(z)=sum z_i B_i`
and `sigma(z)=sqrt(r(t(z))c(t(z)))`.

Supply linear functionals `ell_j`, nonnegative integer weights
`w_j`, and a positive integer denominator `D`. Check, for every
nonzero coefficient vector up to scalar multiplication,

```math
\boxed{\left(\sum_{j:\ell_j(z)\ne0}w_j\right)^2
 \le D^2r(t(z))c(t(z)).}                         \tag{1}
```

The check uses only integer ranks and arithmetic. Scalar
multiplication changes neither the two ranks nor which functional
values are nonzero, so one representative per projective point
is sufficient.

Form the diagonal linear map having `w_j` copies of `ell_j` on
its diagonal. Its rank at `z` is the sum in (1). The
[linear rank-detector lemma](linear-rank-detectors-for-observable-tables.md)
therefore gives the exact lower bound

```math
\boxed{\operatorname{spr}(G_{\rm live})
 \ge\frac1D\sum_{j:\ell_j(e)\ne0}w_j.}          \tag{2}
```

This holds for every observable matrix table over `F_q` whose
entries lie in the specified scalar space, with any finite
auxiliary dimension. The size of the diagonal detector need not
be constructed in memory to verify the certificate.

### Finding candidates

One may choose candidate functionals, impose
`sum_j w_j [ell_j(z)!=0] <= sigma(z)` as linear inequalities,
and maximize the objective in (2). A numerical linear program
only proposes weights. Scaling, rounding weights downward, and
checking (1) yields the actual certificate.

Functionals vanishing on the crossing have zero objective and
can be omitted. Every other functional can be normalized so
`ell_j(e)=1`. This reduces the candidate list to `q^(s-1)`
functionals when the crossing is a basis vector. There are
`(q^s-1)/(q-1)` projective tensor cases.

## 2. A subalphabet entry space

Use five external labels. Let `E,R,C` be the three standard
pairings, and let

```math
H=(E-R-C)\,[x,y,u,v\in\{0,1,2\}].               \tag{3}
```

In characteristic two the minus signs are plus signs. The four
tensors are independent: a pair of distinct labels, one outside
the subalphabet, isolates each of `E,R,C`, while `H` vanishes
there and is nonzero inside the subalphabet.

Write an entry as `aE+bR+cC+hH`. This space passes the simplest
coefficient-of-`E` test: in characteristic two, `H` has cost four
and `E+R+C+H` has cost nine. The latter is below `5^(sqrt(2))`.
The full functional certificate is substantially stronger.

### Certificates over F2 and F4

Over each of these two fields, the following integer weights
satisfy (1) with `D=1`:

| Functional | Weight |
| --- | ---: |
| `a` | 7 |
| `a+b` | 5 |
| `a+c` | 3 |
| `a+h` | 2 |
| `a+c+h` | 1 |

Every functional is one on the crossing. Hence every observable
word-rank spectral transfer in this space has base at least
`18`, separately over F2 and F4. The checker verifies all fifteen
F2 projective entries and all eighty-five F4 projective entries.
The F4 conclusion comes from the additional check, not from the
F2 check alone.

### A certificate over F3

For the signed tensor (3), the following weights satisfy (1)
over F3, again with `D=1`:

| Functional | Weight |
| --- | ---: |
| `a` | 6 |
| `a+h` | 1 |
| `a+c` | 3 |
| `a+2c+2h` | 1 |
| `a+b` | 3 |
| `a+2b+h` | 1 |

They give base at least `15`. The checker covers all forty
projective entries over F3. These bounds close the entry spaces
for this spectral method before any larger auxiliary-state
search is needed.

## 3. Why field scope matters

Let

```math
M(a,b)=a\begin{pmatrix}1&0\\0&1\end{pmatrix}
       +b\begin{pmatrix}0&1\\1&1\end{pmatrix},
\qquad t(a,b)=M(a,b)_{x,u}[y=0][v=0].            \tag{4}
```

Its cut ranks are `(rank M,1)`. Over F2, all three nonzero
coefficient vectors give invertible matrices, so the functional
`ell(a,b)=a` obeys the rational detector inequality with weight
seven and denominator five: `7^2<=5^2*2`.

Over F4, choose `a=1,b=omega`, where
`omega^2+omega+1=0`. The matrix in (4) has rank one. The same
claimed inequality would require `49<=25`, which is false.
This exact control demonstrates why a finite-field certificate
must state its coefficient field.

## 4. Reproducible verification

`check_finite_functional_detectors.py` constructs all three fields'
projective entry tables independently, verifies basis independence,
checks the stated integer certificates, and checks the field-scope
counterexample. It uses exact field elimination and no numerical
optimization.
