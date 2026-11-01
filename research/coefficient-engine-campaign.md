# Coefficient-engine campaign

Owner instruction: work indefinitely toward a complete stronger deterministic
factoring bound. The strongest complete exponent is 1/2. Work directly on
main, coordinate in Discussion #51, preserve current evidence and attribution,
and do not spawn subagents without an explicit user request.

## Completed checkpoints

- `e51a9c1`, verification `5135578`: explicit inverse for supplied one-bit
  character coordinates. The full default suite passed at e51a9c1. This is a
  decoder, not a coefficient provider. See
  [the proof](cyclotomic-character-decoder.md) and
  [verification receipt](../evidence/character-decoder-verification.json).
- [Theta sections](lambert-theta-sections.md): exact coefficient provider via
  DTheta/Theta and characteristic-two sections, with full polynomial state
  cost charged. Literal divisor controls pass; cost remains exponential.
  The common linear section space already needs dimension at least
  floor(2^k/36) at phase order three. The cheap square test there demonstrates
  that this is a representation restriction.
  The full default suite passed at `4d73ddd`: 337 Python command blocks,
  63 unit tests, five crossing tests and 919 Markdown files. See the
  [frozen receipt](../evidence/theta-section-verification.json).
- [Tate coordinates](tate-coordinate-coefficients.md): construct the finite
  etale algebra A_q of degree (q-1)/2 directly from the nodal division
  polynomial, lift its branch X_q(v), and compile internal coefficients in
  poly(q,log(n+1)) cost. At odd original indices the target is
  `[t^N] X_q(Delta(t))`. This composition is not solved. Its scalar algebra
  trace has the cheap identity `T^2+T=Delta(t)+Delta(t^q)`; the remaining
  coefficient space has trace zero and dimension (q-3)/2.

## Active next construction

Seek a nonlinear representation of the trace-zero Delta composition. Do not
repeat the supplied-coordinate decoder, the failed common linear numerator
basis, or merely compute more internal coefficients in v.

The current positive model is the uniform scalar trace, with an explicit
trace-one eta and decomposition `X=eta*Tr(X)+(X-eta*Tr(X))`. At q=3 this
solves everything; at q=5 one binary coordinate remains, and at q=7 two.
First concrete target: derive the q=5 trace-zero scalar's exact equations
and test whether it has a cheap arithmetic predicate, a useful promise
case, or an explicit remaining convolution cost. Keep the unrestricted
engine goal separate from squarefree or semiprime promises.

Other useful representation to examine: weighted algebraic coefficient
automata composed with odd-square theta series, using parameterized
quadratic predicates rather than an exponentially large common linear
basis. The full cost of any lattice/coset/state family must be charged.

The goal remains a uniform original-index evaluator strong enough to improve
the complete exponent beyond 1/2. Internal algebraic coefficient queries and
one cheap trace do not establish that result.

Primary leads: Meier--Ozornova, *Rings of modular forms and a splitting of
TMF_0(7)*, Appendix A.3 (arXiv:1812.04425); Bostan--Christol--Dumas,
*Fast computation of the Nth term of an algebraic series over a finite prime
field* (arXiv:1602.00545). Check the exact statements before using them.
