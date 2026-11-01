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
- [Phase five](q5-quadratic-coefficients.md): the group law gives
  `X^2+B X+B^3=0` and `X=B(c+V)`, where
  `V=sum_(n>=1)(t^(n^2)+t^(5n^2))` and `B=1+V+V^2`.
  Four positive binary quadratic-form parities give the full odd-index
  value in sqrt(N) poly(log N) time. The alternative tau-mod-four query is
  explicitly charged; squarefree/semiprime conclusions are promises only.

## Active next construction

Seek a nonlinear representation of the trace-zero Delta composition. Do not
repeat the supplied-coordinate decoder, the failed common linear numerator
basis, or merely compute more internal coefficients in v.

The current positive model is the uniform scalar trace, with an explicit
trace-one eta and decomposition `X=eta*Tr(X)+(X-eta*Tr(X))`. At q=3 this
solves everything; at q=5 one binary coordinate remains, and at q=7 two.
The q=5 scalar is now reduced to four binary quadratic predicates with its
cost charged. Do not repeat that derivation or replace its divisor-count
query with an unproved squarefree-part oracle.

Next concrete target: phase seven. Check a cubic parametrization of its
canonical kernel and test finite polynomial expressions using several
unary-square residue series. Use the known q3/q5 identities as controls,
hold out larger prefixes, and require an ordinary proof before asserting
any identity at all indices. A finite prefix fit is only a candidate.

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
