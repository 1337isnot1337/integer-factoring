# A thirteen-letter global negation sector

## Registered realization diagnostic

After the [combined thirteen-letter profile target](joint-bare-profile-feasibility.md),
test an actual crossing rather than further optimizing necessary profiles.
Base: `7f39198`. The complete factoring score is `1/2`.

Use three copies of the [three-gate qutrit cell](cyclic-pauli-crossing-diagnostic.md),
lifted to `F7` with cubic root `omega=2`. Normalize each Weyl projector by
`1/3`, and write the gate as

```
P_U = (I+U+U^2)/3,
F_U,lambda(z) = z I + (1-lambda)(1-z) P_U.
```

Use the same three Pauli vectors and top eigenphases
`(lambda1,lambda2,lambda3)=(omega,omega^2,omega)` as the earlier cell.
Verify the full polynomial endpoint and both leading rank-one assertions;
the field lift is not presumed to preserve them without checking.

The new physical image is the global negation-odd subspace of functions
on `F3^3`, of dimension `(27-1)/2=13`. For each pair `{a,-a}` choose
the representative with first nonzero coordinate one. The lift column
is `L_a=e_a-e_(-a)` and the projection is `R=L^T/2`, so `RL=I13`.
Apply `R` on the two output legs and `L` on the two input legs. This
preserves the literal crossing endpoint, with the scalar field cost and
all physical labels explicit. `F7` supplies both `1/2` and `1/3` and a
cubic root; this is a changed physical representation, not a field sweep.

The uncoupled bare cell floor applies before this noncoordinate projection.
The exact Frobenius identity for the earlier four-to-two `F4` CSS code
concerns a different code and field. No extension of either obstruction
to this global involution sector is assumed. The full degree-nine tensor
must be constructed, including every projected entry.

For each cut compute exact block-Toeplitz ranks `R4,R5`. Their difference
is the number of Smith orders at most four. If both differences are at
least 38, the retained channels imply, at every `t>0`,

```
Z_zero(t) >= 38 exp(-4t),
Z_infinity_normalized(t) >= 38 exp(-4t),
B(t)=exp(9t) Z_zero(t) Z_infinity_normalized(t) >= 38^2 exp(t).
```

Rational logarithm bounds can certify
`log38/(2 log13)>1/sqrt(2)`, ending the bare diagnostic without computing
unneeded higher orders. Unknown channels are retained as additional
nonnegative contributions, not assigned guessed orders. If this test does
not exclude the target, compute complete untruncated orders before any
pressure claim. No neighboring-code, field, copy-count or parameter sweep
is registered. A failed lift or a proved bare floor stops this proposal.

No actual gain or independent mathematical review is claimed. No Lean.
Substantial OpenAI assistant support.

## Exact outcome and stopping certificate

The lift succeeds: the complete `F7` cell has the stated endpoint and
both leading ranks one. The global odd split has `RL=I13`, and all
28,561 projected endpoint entries agree with the literal crossing.
The degree-nine polynomial is constructed in full. Independent literal
summation agrees on 128 selected complete polynomials; the four-letter
two-copy control agrees at all 256 entries and has both complete order
profiles independently reconstructed.

For the thirteen-letter tensor the exact Toeplitz ranks are:

| Cut | `R4` | `R5` | Orders at most four |
| --- | ---: | ---: | ---: |
| Zero, `(x,y)` versus `(v,u)` | 52 | 122 | 70 |
| Normalized infinity, `(x,u)` versus `(v,y)` | 235 | 404 | 169 |

The [checker](../tools/check_negation_odd_qutrit.py) validates every rank
with a nonsingular selected minor and a full exact factorization of the
Toeplitz matrix. It explicitly inverts that minor, checks the inverse,
and verifies all entries of `A=A[:,J] G^(-1) A[I,:]`. Thus the minor
gives the lower rank bound and the factorization the upper rank bound.
The [receipt](../evidence/negation-odd-qutrit.json) retains the pivot sets,
matrix/factor hashes, polynomial hashes, controls and rational logarithm
bounds. The smaller control also compares with the independent literal
Toeplitz and polynomial-series algorithms.

Both counts exceed the registered threshold 38. The preceding retained-channel
argument therefore excludes this explicit tensor's **bare balanced** rate
at every scale. The rational logarithm bounds verify
`log(38)^2>2 log(13)^2`. Higher orders are unnecessary for that conclusion;
no complete thirteen-letter profile is claimed. Free physical gauges and
other constructions are not covered by this stopping certificate.

Run `make verify-negation-odd-qutrit`. The registered bare realization
diagnostic stops here. The combined thirteen-letter necessary-profile
target remains unrealized, and the complete score remains `1/2`.
