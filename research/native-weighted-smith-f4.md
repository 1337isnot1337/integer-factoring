# Optional native weighted Smith search over F4

## Scope

The [C++ engine](../tools/native_weighted_smith_f4.cpp) and its
[Python wrapper](../tools/native_weighted_smith_f4.py) preserve the exact
finite-field arithmetic used in the residual-code searches. They support
both binary tensors and tensors over
\(\mathbf F_4=\mathbf F_2[\omega]/(\omega^2+\omega+1)\).
The physical alphabet can have size 2 through 32.

**This is an optional search tool, not a numerical proof certificate.**
Polynomial arithmetic is exact, but physical heights, pivot comparisons,
pressures, gradients, and type optimization use floating point numbers.
Any claimed improvement requires a separate exact replay, such as the
rational weighted-minor or exact pressure certificates elsewhere in this
repository. No new pressure bound is asserted by this package.

The engine uses a 512-coefficient series window. It returns an error if
this window does not reveal every pivot. It does not report a partial
profile as a full profile. This protection addresses truncation; it does
not remove the separate need to certify floating pivot choices.

## Arithmetic and matrix decomposition

A polynomial is represented by two binary bit masks, its coefficients of
\(1\) and \(\omega\). Multiplication uses three binary carryless products.
If these products are
\(p=a_0b_0\), \(q=a_1b_1\), and
\(r=(a_0+a_1)(b_0+b_1)\), the resulting planes are
\((p+q,r+p)\). A binary fast path uses just one product. The compiler uses
the processor's carryless multiplication instruction when available; a
portable bit loop supplies the same operation otherwise.

Unit-series inversion starts with the inverse of the actual leading
field coefficient. In particular, \(\omega^{-1}=\omega^2\), so an
initial inverse of one would be incorrect for a general F4 pivot.

For each flattening, the constructor computes the connected components
of its bipartite nonzero support. Balanced components are profiled
separately, and their Smith orders are merged. This decomposition is
exact and applies to arbitrary support patterns. It can substantially
reduce the work for parity- or charge-conserving tensors.

The optional `PairedOracle` constructor forms two conjugate copies with
a common logical core bit and independent bulk bits. It supports one
or two bulk bits per copy. Its four parity components avoid constructing
a dense 1024 by 1024 matrix for the larger case. Coefficients violating
these parity components are rejected. The constructor does not itself
prove a crossing endpoint or a common-type polynomiality statement.

## Use

The wrapper requires Python's standard library and an available C++17
compiler. Set `CXX` if the compiler is not named `c++`. NumPy and SciPy
are optional dependencies of optimization programs, not of this wrapper
or its finite control checker. The first call compiles the adjacent
source into a temporary cache indexed by its contents and build flags.
Compilation finishes before the completed library is installed atomically.
No generated binary is stored in the repository.

```python
from native_weighted_smith_f4 import Oracle

# values[x + q*y + q*q*u + q**3*v] = [low_plane, high_plane]
oracle = Oracle(values, q, endpoint_check=True)
params = [0.0] * oracle.parameters
params[-1] = 1.0  # log(s)
cost, gradient, spectra = oracle.evaluate(
    params, [1 / q] * q, cut=1, want_spectra=True
)
rho, probabilities, spectra, weights = oracle.types(params, cut=1)
```

The parameter vector contains infinity heights, zero heights, and then
\(\log s>0\). Each set of heights lists the four legs in order
\((x,y,u,v)\), with physical labels \(1,\ldots,q-1\); label zero has
height zero. The supplied heights already include their factor of
\(\log s\). Cut 1 places the column flattening at infinity and the row
flattening at zero. Cut 0 exchanges them. Returned spectra list infinity
first and zero second.

`evaluate` returns the logarithmic pressure after restoration for the
given common type. `raw_evaluate` returns the logarithm of the product
of the two partition sums before physical-type restoration; it supports
search programs that impose pair types or other mean-height constraints.
`prefixes` returns selected attaining prefix-minor masks for \(q\le8\).

For a common rational denominator that is a unit at zero, the optional
`denominator_degree=d` adds \(d\) to every infinity Smith order. The
wrapper includes the corresponding pressure and gradient correction.
This option is arithmetic bookkeeping. The caller must separately prove
that the rational tensor is admissible after the intended restriction.

The native library stores one active seed per process. Constructing a
second oracle invalidates the first instance's active status. Reusing it
raises an error until `first_oracle.activate()` reloads its seed. A lock
prevents simultaneous native calls from observing an inconsistent seed.

## Independent finite controls

Run this optional check from the repository root:

```sh
python tools/check_native_weighted_smith_f4.py
```

The [recorded evidence](../evidence/native-weighted-smith-f4.json) compares
192 weighted orders with the independent exact F4 Bareiss oracle using
rational weights. It also checks every parameter derivative of a small
example, selected prefix masks, denominator restoration, rejection of an
inactive seed, and rejection of a profile hidden beyond the series window.
Finally, all 4096 entries of a conjugate paired seed are checked directly,
and the paired and generic constructors give identical profiles.

These controls are reproducible without NumPy or SciPy. They validate
finite implementation cases and the intended failure behavior; they do
not establish a universal floating point error bound.
