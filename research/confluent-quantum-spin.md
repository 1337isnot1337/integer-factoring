# A confluent seven-letter quantum-spin crossing

## Question, result and stopping condition

Can collision of the quantum-spin Casimir eigenvalues at a primitive seventh
root produce the short flag needed by the surviving
[asymmetric seven-letter profile](seven-letter-asymmetric-profile.md)?
Unlike the previous rational-parameter screen, this tests a nonsemisimple
root-of-unity specialization. It changes the representation, rather than
increasing a field or parameter cutoff. The stopping test is the actual
polynomial tensor, its common endpoint and both local Smith profiles.

The confluent tensor exists over `F_8`, has degree six and has the flip
endpoint. One specialization has first-cut rank one. Nevertheless, its
**bare two-place pressure is at least `7^4`, over every finite coefficient
extension and every two places**. Its best bare coefficient is therefore
`rho=1`. This fixed degeneration cannot cross `1/sqrt(2)` and is stopped.
No stronger complete factoring bound follows. Physical gauges, collective
constructions, multipoint criteria and other root-of-unity representations
are outside this statement.

The construction regularizes the projector formulas in
[quantum-crossing-initial-forms.md](quantum-crossing-initial-forms.md),
which credits Bytsko's [higher-spin R-matrix paper](https://arxiv.org/abs/math/0412482).
The definitions below and the exact replay specify the tensor; a
Yang--Baxter identity is not needed for its endpoint or for this exclusion.

## 1. Paired Hermite regularization

Use `F_8=F_2[a]/(a^3+a+1)`, write `zeta=a`, and put `Q=zeta+epsilon`.
For `k=6`, take the Laurent-polynomial Casimir `C(Q)`, eigenvalues
`lambda_j=[j]_Q[j+1]_Q`, scalar numerators `f_j(u,Q)` and normalization
`D(Q)` from equations (1) and (4) of the quantum-spin note. The integer
weights `H=3-i` require no square-root choice. Every entry of `C(Q)` is
regular at `Q=zeta`, and `D(zeta)=1`.

The eigenvalue triples through `epsilon^2`, with field elements encoded
in the polynomial basis `1,a,a^2`, are

```text
j       0       1       2       3       4       5       6
lambda (0,0,0) (7,6,6) (5,7,2) (4,2,7) (5,7,5) (7,6,4) (0,0,1).
```

Thus the pairs `(0,6),(1,5),(2,4)` collide to order exactly two, with
leading differences `1,2,7`. Their distinct limiting eigenvalues are
`0,7,5`; the singleton `j=3` limits to `4`. For each pair `(j,l)`, the
coefficients of `f_j-f_l` vanish through order one. Define `r_u(X)` of
degree at most six by these seven conditions:

```math
r_u(lambda_j(zeta))=f_j(u,zeta)/D(zeta),
```

```math
r'_u(lambda_j(zeta))=
 \frac{[epsilon^2](f_j(u,Q)-f_l(u,Q))}
 {D(zeta)[epsilon^2](lambda_j(Q)-lambda_l(Q))}
 \quad ((j,l)=(0,6),(1,5),(2,4)),
```

and `r_u(lambda_3(zeta))=f_3(u,zeta)/D(zeta)`. The prime denotes the
ordinary derivative in `X`, including in characteristic two. Hermite
interpolation at these distinct nodes, of multiplicities `2,2,2,1`,
is invertible over `F_8`.

This is a genuine regularization of the scalar spectral interpolant.
Replace the second evaluation row for each colliding pair by its divided
difference. Its coefficient row tends to the derivative row above; its
right-hand side is regular because the numerator difference has order
at least two and the eigenvalue difference has order exactly two. The
resulting interpolation matrix has an invertible limiting Hermite matrix,
so its inverse and the interpolating coefficients are regular. Since
`C(Q)` is regular, the matrix interpolant has the indicated specialization.

Define the tensor as `R(u)=flip*r_u(C(zeta))`. Its entries have degree
at most six in `u`. At `u=0` all value conditions are one and derivative
conditions zero, so `r_0(X)=1` and `R(0)=flip`. The replay checks all 2401
entries of this endpoint. Integer total-weight conservation gives both
cuts exact square blocks of order at most seven.

## 2. Local profiles and determinant budgets

The replay computes every profile at all nine points of `P^1(F_8)`.
Infinity uses `t^6 R(t^-1)`, with the global degree six. Both determinants
are nonzero homogeneous forms of degree `49*6=294`, since the endpoint
is invertible. Let `m_c(P)` be the sum of the Smith orders in cut `c`.
The exact totals are:

| Place | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | Infinity | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| First cut | 0 | 0 | 28 | 42 | 56 | 84 | 14 | 70 | 0 | 294 |
| Second cut | 0 | 56 | 0 | 66 | 0 | 0 | 6 | 6 | 56 | 190 |

At `u=5`, the first profile is `0^1,1^12,2^36`: rank one does not
supply a deep flag. At `u=3`, the second is `0^28,2^12,4^6,6^3`.
All profiles are in the [exact receipt](../evidence/confluent-quantum-spin.json).

The first determinant has spent its entire degree at the displayed
rational points. It has no other zeros, so `m_0(P)<=84` at every point
of the algebraic closure. For the second determinant, 104 degrees remain.
Any point outside `P^1(F_8)` has at least two distinct Frobenius conjugates
with equal multiplicity. Consequently its multiplicity is at most 52.
Together with the displayed points this gives `m_1(P)<=66` everywhere.
In particular,

```math
m_0(P)+m_1(Q)<=150<294
```

for any two places, in either cut orientation, after any finite field
extension. This uses the same determinant section throughout; the degree
budget is not restarted at another place.

## 3. Exact pressure exclusion

For the bare two-place criterion, write the first local orders as
`mu_1,...,mu_49`, the other cut's orders as `nu_1,...,nu_49`, and put

```math
B(s)=\left(\sum_i s^{6-mu_i}\right)
     \left(\sum_i s^{-nu_i}\right),\qquad s>=1.
```

Jensen's inequality gives

```math
B(s)>=49^2 s^{6-(sum_i mu_i+sum_i nu_i)/49}
     >=2401 s^{144/49}>=2401.
```

Equality holds at `s=1`. Thus `log B/(4 log 7)>=1`, excluding the target
`rho<1/sqrt(2)` without numerical optimization. The root-of-unity
collision is legitimate, but this tensor distributes its determinant
zeros too broadly to help this criterion. No larger-field replay or
finer optimization of this fixed tensor is warranted.

## Reproduction and review scope

```sh
make verify-confluent-quantum-spin
```

The [checker](../tools/check_confluent_quantum_spin.py) reconstructs the
jet conditions, Hermite interpolant, tensor and profiles. It compares the
scalar-restriction native Smith engine with 440 independent block-Toeplitz
nullity calculations over `F_8`, and verifies the endpoint and determinant
ledger. The receipt contains exact field elements and integer orders;
no floating-point estimate is used for the exclusion. These are finite
exact checks supporting the displayed ordinary proof, not independent
mathematical review or a new Lean formalization.

The [clean full default replay](../evidence/confluent-quantum-spin-full-verification.json)
passes at exact source `570f6166360eeca50243e822465358fe1ffd8ce0`
in 731.40 seconds. The subsequent archive commit changes only documentation
and verification records.
