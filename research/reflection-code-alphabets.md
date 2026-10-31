# Restricting the physical alphabet of a reflection crossing

This checkpoint tests a direct way to leave the Gaussian class while retaining
an exact crossing endpoint. It does **not** improve the factoring bound.
The best value found in the recorded 236-alphabet screen is the unrestricted
reflection baseline, approximately `rho = 0.7166536597281878`.

## Construction

Let `K_M(z)` be the strict reflection crossing from
[`reflection-dual-rail-compression.md`](reflection-dual-rail-compression.md), with
physical alphabet `{0,...,2^(M/2)-1}`. Choose a subset `C` of this alphabet and
restrict all four physical legs to `C`. If `q = |C|`, then

\[
 K_{M,C}(1;x,y,u,v)=[x=v][y=u]\qquad(x,y,u,v\in C).
\]

Thus every restriction remains an exact `q`-letter crossing. A non-Cartesian
subset, such as an even-parity code, generally gives a non-Gaussian tensor.
This restriction acts on the **physical alphabet**, rather than on a conserved
virtual charge or a subset of one-particle Smith channels.

The local Smith profiles must therefore be recomputed from the complete
restricted flattenings. Selecting entries from the original one-particle
profile would give the wrong object.

The screen uses the original natural diagonal gauges. For `M=8`, write a
physical state as four bits and assign bit energies

\[
 (28,20,12,4),\qquad
 E(x)=28x_0+20x_1+12x_2+4x_3.
\]

At zero, the `(x,y)|(u,v)` flattening has row height
`(E(x)+E(y))/2` and column height zero. At infinity, the
`(x,u)|(y,v)` flattening has row height `E(u)/2` and column height `E(v)/2`.
The native exact polynomial elimination returns the restored profiles
`nu_C` and `mu_C`. These profiles can contain negative orders at infinity;
the scalar shifts must be retained consistently.

## Pressure calculation

For `t = log(s) > 0`, put

\[
 L_C(t)=\log\sum_i e^{-t\nu_{C,i}}
       +\log\sum_i e^{-t\mu_{C,i}}.
\]

The paired-type variational formula from
[`paired-type-crossing-pressure.md`](paired-type-crossing-pressure.md)
reduces the physical type optimization to the scalar equation

\[
 L_C(t)=2R\log\sum_{x\in C}e^{-tE(x)/R},
 \qquad \rho=R/2.
\]

The search solves this equation numerically, then minimizes over
`10^-5 <= t <= 2`. It optimizes the actual quotient defining `rho`, not a
pressure with a frozen denominator. The latter can misleadingly approach zero
when the physical type collapses to one symbol.

## Recorded screen

The deterministic screen contains 236 distinct subsets for `M=8`: the full
alphabet, both parity classes, three Hamming balls, twelve additional energy
balls, fourteen additional binary hyperplanes, and 204 fixed-seed random
subsets. The categories are reported after deduplication. The random seed is
`20261031`.

| Alphabet | Numerical `rho` |
| --- | ---: |
| All sixteen symbols | 0.7166536597281878 |
| All symbols except `15` | 0.7166714030963784 |
| All symbols except `7` | 0.7166860282214897 |
| All symbols except `7,15` | 0.7167032790564378 |
| `{0,4,8,12}` | 0.7310809359736002 |
| Even-parity symbols | 0.8157520646957259 |
| Odd-parity symbols | 0.8227089855423472 |

These are exploratory numerical minima for the listed restrictions and fixed
gauge directions. They do not exclude other subsets, nonuniform gauges, other
reflection intervals, or more general linear combinations of physical states.
No strict pressure certificate below `1/sqrt(2)` was found.

## Exact controls and reproducibility

[`check_reflection_code_alphabets.py`](../tools/check_reflection_code_alphabets.py)
checks all eleven subsets of size at least two for `M=4`, and three selected
`M=8` subsets: `{0,4,8}`, `{0,4,8,12}`, and `{0,3,12,15}`. It checks 1,269
complete endpoint entries, recomputes 28 local Smith profiles with an
independent Python polynomial implementation, and performs 252 block-Toeplitz
rank checks. All checks pass. This validates the restriction and profile
pipeline on the stated cases; it does not turn the numerical pressure search
into an exhaustive optimization certificate.

Run from the repository root:

```sh
python3 tools/check_reflection_code_alphabets.py \
  --output evidence/reflection-code-alphabets.json
python3 tools/search_reflection_code_alphabets.py --random 20 \
  --out evidence/reflection-code-alphabet-search.json
```

The search requires NumPy, SciPy, and the existing native F4 polynomial engine.
The evidence files contain the exact profile multiplicities and source hashes.

The useful outcome is a checked method for testing correlated physical
alphabets. These small-code and subset controls show that removing physical
states can worsen the entropy-normalized contraction despite decreasing some
individual ranks.
