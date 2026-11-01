# Complete affine equations in F4 reference frames

The finite-window Smith backend now accepts genuinely F4-valued reference
tensors and frames. It preserves the complete degree-bounded strict-endpoint
coefficient universe. The historical binary-frame solver already allowed F4
coefficient choices; this extension allows its resulting nonbinary tensors
to serve as new reference frames. No improving contraction is established.

## Exact construction

The [frame implementation](../tools/windowed_f4_smith.py) uses the
[proved coefficient guard](windowed-smith-frames.md): an internal precision
`P > 2 M + b`, with maximum resolved order `M` and permitted row increase
`b`. The determinant degree bound supplies a finite retry ceiling. Field
polynomials use two packed coefficient bits and exact arithmetic in
`F4=F2[omega]/(omega^2+omega+1)`.

The [affine constructor](../tools/windowed_f4_affine_lattices.py) imposes
every coefficient of each row condition on `T=T0+(1+t)Q`, with
`deg Q < D`. A regular right frame can be omitted: multiplication by it
and its regular inverse preserves each row's divisibility ideal. Left
frames still retain every physical column equation.

For an F4 equation with coefficient planes `a,b` and unknown planes `x,y`,
the two binary equations have coefficients `(a,b)` and `(b,a+b)`. These
couple the planes when the reference frame is nonbinary. The source has
`D*q^4` F4 variables and twice as many binary variables. The compressed
fiber retains this coupled binary system. Witness checks reconstruct the
literal endpoint, full left-frame products, and complete Smith prefixes.

## Independent controls and actual five-letter pilot

Run `make verify-windowed-f4-affine`. The
[checker](../tools/check_windowed_f4_affine.py) verifies 48 augmented binary
row-space equivalences against the original solver, 561 independently
computed complete fraction-free orders, and 948 literal endpoint entries.
Controls include nonbinary references, compressed fibers, field-valued
frame units, known nonzero Smith orders, and an explicit precision retry.
The [receipt](../evidence/windowed-f4-affine.json) retains the inputs and
source hashes. The original source and historical receipts are preserved.

An [actual degree-42 five-letter pilot](../evidence/windowed-f4-q5-pilot.json.xz)
built 112,100 binary equations in about 2.65 seconds, representing 26,250
F4 coefficient variables. Its row-floor fiber has 3,700 binary free
coordinates. A sampled witness has 620 nonbinary entries. All fifty
orders of its two complete profiles were independently reconstructed by
the fraction-free oracle; the binary-frame backend also checked its
literal row conditions. The complete pilot took about 83.5 seconds.
These are measured host timings, not a runtime theorem. This pilot used a
binary reference and a nonbinary witness; small controls also rebuild
frames from genuinely nonbinary references.

The predecessor binary adaptive run has now completed its four registered
frames and all 800 walk attempts. Its
[complete packet](../evidence/adaptive-q5-affine-v3-complete.json.xz) retains
the manifest, driver, witnesses and frame-level counters. Direct target
tests failed in all four frames; the best measured ratio stayed about
`0.89459228`. The `completed_walks` diagnostic was only updated on selected
events and is not the coverage counter: `walks_attempted` is 200 per frame.
This is finite frame coverage. The subsequent
[subalphabet cofactor argument](subalphabet-adjugate-caps.md) excludes the
target independently of these experiments.

The backend is ready for targets surviving stronger necessary constraints.
The target exclusions and numerical pilot do not improve the complete
factoring exponent **1/2**. Substantial OpenAI assistant support was used;
no new Lean formalization was performed.
