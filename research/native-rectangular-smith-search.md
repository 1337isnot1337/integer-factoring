# Optional native search for rectangular Smith pressure

The optional native engine evaluates both weighted Smith profiles of an
\(A\times B\times B\times A\) polynomial tensor over \(\mathbf F_2\).
It also optimizes the two independent physical distributions through the
[paired-type Gibbs identity](paired-type-crossing-pressure.md).
The physical tensor may have unequal crossing-pair alphabets; both cost
flattenings are square matrices of size \(AB\).

This is an exploratory search tool. Its polynomial arithmetic is exact,
but its pivot weights, pressure values, and type optimization use floating
point. The mathematical transfer theorem and the exact endpoint proofs
are separate from this implementation. A numerical value below a target
must be replayed with exact Smith or rank bounds before it is a certificate.

## Build and interface

Python uses only the standard library. A C++17 compiler is required for
the optional native engine. Set `CXX` to select a compiler. The wrapper
compiles the adjacent source into a temporary cache keyed by its contents,
compiler command, and platform. On supported x86 processors it uses
carryless multiplication instructions; a portable integer implementation
is available otherwise. NumPy and SciPy are not required for the wrapper,
checker, or command-line interface. Separate optimization scripts may use
them.

The binary polynomial mask \(p\) represents
\(\sum_{j:\,p_j=1}w^j\). The dense tensor index is

\[
 x+A\bigl(y+B(u+Bv)\bigr).
\]

The parameter vector lists infinity heights, then zero heights, then
\(\log s>0\). At each endpoint, list the four legs in order
\(x,y,u,v\), with each leg's labels \(1,\ldots,q_j-1\). Label zero
has height zero. Each height is already multiplied by \(\log s\).
There are \(4(A+B-2)+1\) parameters.

```python
from native_rectangular_smith import RectangularOracle

oracle = RectangularOracle(values, A, B,
                           denominator_degree=0,
                           endpoint_check=True)
cost, gradient, spectra = oracle.evaluate(
    params, pi_x_v, pi_y_u, cut=1, want_spectra=True)
rho, pi_x_v, pi_y_u, spectra, G_X, G_Y = oracle.types(params, cut=1)
```

Cut 1 uses \((x,u)\mid(y,v)\) at infinity and
\((x,y)\mid(u,v)\) at zero. Cut 0 reverses these choices. Spectra list
infinity first and zero second, with orders multiplied by \(\log s\).
The returned `cost` is the logarithm of the product profile bound after
both physical type means have been restored.

The optional `denominator_degree=d` restores division by a scalar
polynomial of degree \(d\) that is a unit at zero. Every infinity order
increases by \(d\), and the logarithmic pressure decreases by
\(d\log s\). This restoration is applied exactly once. Supplying this
number does not prove that the rational tensor is an admissible crossing.
For the rational rectangular constructions, use `endpoint_check=False`
and retain their separate polynomiality and endpoint proof.

Only one seed is active inside a loaded native library. A generation guard
rejects queries on an older instance after another seed has been loaded.
Call `oracle.activate()` to restore its retained seed. Failed seed loads
also invalidate the preceding generation. A lock serializes access to
the native seed.

The internal power-series window contains 512 coefficients. If it cannot
resolve a complete profile, the query raises an error; it does not return
a shorter profile. Input masks whose degree exceeds 511 are rejected.
The current engine accepts \(2\le A,B\le64\); large dense matrices can
still be expensive in time and memory.

## Reproduction

Run the independent exact controls:

```sh
python3 tools/check_native_rectangular_smith.py \
  --output evidence/native-rectangular-smith.json
```

The checker compares eight weighted profile pairs to the independent
fraction-free Python oracle at exact rational weights. It covers both
orientations, square and rectangular alphabets, and two scalar denominator
restorations. It checks the pressure reconstructed from those exact orders,
the paired Gibbs entropy identity, 24 finite-difference gradient coordinates,
active-seed isolation, failed loads,
oversized polynomial rejection, and unresolved singular profiles.

To generate the residual-code seeds and inspect one pressure query:

```sh
python3 tools/check_rectangular_residual_repair.py \
  --seed-directory /tmp/rectangular-seeds
python3 tools/native_rectangular_smith.py \
  /tmp/rectangular-seeds/rectangular-partial-012.json \
  --log-s 0.1 --cut 1
```

An optional `--params PATH` accepts a JSON list containing the full weighted
parameter vector. An optional `--output PATH` saves the numerical result.
These commands do not claim a new factoring exponent.

## Four-element-field extension

The companion `native_rectangular_smith_f4.py` has the same rectangular
interface and paired-type calculation. Import its `RectangularOracle` and
supply each polynomial as a pair `(low, high)` of binary coefficient masks,
representing `low + omega*high` with `omega^2=omega+1`. Its adjacent C++ source
uses exact field operations and separates balanced support components before
elimination. Each leg retains its own alphabet and height offsets.

The [mixed-field repair](rectangular-f4-repair.md) supplies admissible rational
seeds. Pass the recorded `denominator_degree` and `endpoint_check=False`;
the paired-type endpoint proof supplies admissibility. The engine restores
this scalar denominator once. A finite profile alone does not prove the
endpoint or an improved asymptotic factoring rate.

```sh
python3 tools/check_native_rectangular_smith_f4.py \
  --output evidence/native-rectangular-smith-f4.json
python3 tools/check_rectangular_f4_repair.py --seed-directory /tmp/f4-seeds
python3 tools/native_rectangular_smith_f4.py \
  /tmp/f4-seeds/rectangular-f4-12-0-1-2-M4.json --log-s 0.1
```

The F4 checker compares both cuts at exact rational weights for four shapes,
including unequal alphabets and nonzero denominator corrections. It checks
24 gradient coordinates, the paired Gibbs formula, and failed-load,
oversized-input, endpoint, and singular-profile controls. An additional
independent fraction-free replay on the actual mixed `(4,3)` seed agreed
in both cuts to less than `3e-16` in the native floating outputs.
