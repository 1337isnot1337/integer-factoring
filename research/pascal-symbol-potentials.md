# Correlated symbol potentials on restricted Pascal crossings

## Constructive extension and boundary

The [free mode-gauge search](permuted-pascal-free-gauges.md) extends to
arbitrary potentials on complete occupation words. These potentials can
correlate different occupied modes. An arbitrary coordinate subalphabet
can also be retained. The evaluator computes both complete profiles from
the original polynomial tensor, so it requires neither a Gaussian
factorization after restriction nor the earlier sufficient gauge cone.

The [implementation](../tools/pascal_letter_gauges.py) and
[search driver](../tools/search_pascal_letter_gauges.py) retain both
independent pair costs. This is a further search interface, with exact
finite validation. There is no strict improving pressure certificate;
the established factoring score remains **1/2**.

## 1. Every coordinate subalphabet retains the typed endpoint

Let `P[i,j]=[j subset i]` on the `L` Pascal modes, and let `S` be
any set of distinct occupation words. Take the principal physical
restriction of the existing polynomial Gaussian tensor to `S` on
each of its four legs. At the endpoint its entries are

```math
 T_1(b,x,v,y)=[b=y]\det P[v,x],
```

where unequal-size minors are zero. For an occupation word `s`, put
`W(s)=sum_(i in s) i`. A nonzero term of `det P[v,x]` matches every
input mode to a containing output mode. Consequently it either fixes
every mode, giving `v=x`, or strictly increases `W`. The diagonal
minors equal one.

Prescribing the same word type on `X,V` fixes their total `W`. In a
tensor power no strictly increasing endpoint entry can then occur.
The endpoint is exactly `[B=Y][X=V]` on these type classes. A separate
common type may be prescribed on `B,Y`. The argument permits every
coordinate subset `S`, without a deletion-closure assumption.

Rational symbol heights at either place add only a linear-in-word-length
Laurent degree range after fixed ramification, and their value at the
endpoint is one. The [paired-type transfer](paired-type-crossing-pressure.md)
therefore applies. It supplies a uniform crossing family whenever its
strict pressure condition is actually certified.

## 2. Complete original-polynomial pressure

Write the eight symbol potentials as `h_(c,leg)(s)`, where `c` is
zero or infinity. For a common positive scale `t`, compute the complete
orders of the actual two cut matrices with those row and column shifts:
`nu_(c,i)`. Infinity reverses the polynomial variable and charges its
raw negative degrees. Define

```math
 B=\left(\sum_i e^{-\nu_{0,i}}\right)
   \left(\sum_i e^{-\nu_{\infty,i}}\right),
 \qquad
 G_{BY}(s)=\sum_c(h_{c,B}(s)+h_{c,Y}(s)),
 \qquad
 G_{XV}(s)=\sum_c(h_{c,X}(s)+h_{c,V}(s)).
```

The paired Gibbs root is

```math
 \log B=2\rho\left(
   \log\sum_{s\in S}e^{-G_{BY}(s)/(2\rho)}+
   \log\sum_{s\in S}e^{-G_{XV}(s)/(2\rho)}\right).
```

Both physical pair energies are recomputed at every objective query.
For additive symbol potentials, the original-Fock calculation recovers
the earlier Clifford pressure, including its scalar normalization.
The search permits nonadditive symbol potentials and can lift a saved
mode-gauge start, preserving that start's exponent permutation and flag.

The finite implementation currently constructs the literal Fock tensor
on at most four modes. Every generic cut is full rank: at the endpoint
it is a permutation of copies of a principal unipotent triangular
compound of `P`. A missing promised native pivot is an error, not an
assumed zero. Floating searches do not replace an exact rational
pressure certificate.

## 3. Independent controls

The [checker](../tools/check_pascal_letter_gauges.py) verifies all endpoint
entries for four specified kernels. Five complete integer-profile cases
include arbitrary symbol heights, negative heights, a full 16-letter
alphabet, its eight-letter even-occupation subalphabet, and a nonmonotone
exponent permutation. Pure Python polynomial elimination agrees with
**740 original-polynomial Smith orders** from the native engine.

Eight separate numerical controls recover the earlier additive-mode
pressure. Forty-eight central differences check the nonadditive symbol
gradients. These derivative comparisons are numerical and are recorded
separately from the exact endpoint and integer-profile checks.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_pascal_letter_gauges.py \
  --output /tmp/pascal-letter-gauge-controls.json
```

The [receipt](../evidence/pascal-letter-gauge-controls.json) retains the
actual heights and integer profiles, numerical errors, and source hashes.
The search and numerical controls require NumPy and SciPy in addition
to the native engine.

## 4. Registered coordinate-alphabet screen

The [compressed packet](../evidence/pascal-coordinate-alphabet-screen.json.xz)
retains 339 registered cases: the full alphabet, every single omission,
both occupation parities, and forty seeded samples at each size seven
through fourteen. The OS-random seed and the complete candidate list
were saved before running the profiles. There are 329 distinct alphabets;
the sampled list intentionally retains any repetitions.

Every case has both complete integer profiles at the inherited additive
heights. Subsequent scale minimization is numerical. The best retained
fit is the full 16-letter reference, approximately `0.731080935974`.
This screen does not optimize every symbol potential, prove a global
scale minimum, or exclude other alphabets or kernels. The archive
contains the full driver and all actual profiles, not only fitted values.

The constructive next step is to optimize the larger symbol-potential
space and use occupation sectors to compute larger physical restrictions
without expanding their entire ambient Fock tensor.

This work reuses and credits the earlier weighted Pascal, arbitrary
permutation, paired-type, native Smith, and original-Fock constructors.
Substantial OpenAI assistant support was used for the implementation,
proof, computation, and text. No new Lean formalization is claimed.
