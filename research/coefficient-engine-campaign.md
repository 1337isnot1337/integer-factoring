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

## Active next construction

Explore a change of arithmetic variable instead of a common linear section
basis. For odd coefficient indices, the signed additive divisor coefficient
is the nonconstant coefficient of the x-coordinate of a q-torsion point on
the Tate curve in characteristic two. Translating y by a4 gives
`y^2+xy=x^3+Delta(t)` with x unchanged. Here
`Delta(t)=sum_(r>=0) t^((2r+1)^2)` modulo two.

Proposed route, not yet a repository theorem: construct the q-division
polynomial on `y^2+xy=x^3+v`, lift the branch with constant
`x0=Z/(1+Z)^2`, and compile its coefficients in v to a finite algebraic
section machine. The original index problem then remains coefficient
extraction after v=Delta(t). A fast coefficient of the algebraic series in
v must not be confused with a fast coefficient of its composition in t.

First checks: division-polynomial recurrences, simple branch derivative at
v=0, exact Tate/Delta coefficient comparisons at odd indices, then a
Furstenberg rational-diagonal compiler with its complete construction cost.
The potentially expensive quadratic composition must remain explicit.

Primary leads: Meier--Ozornova, *Rings of modular forms and a splitting of
TMF_0(7)*, Appendix A.3 (arXiv:1812.04425); Bostan--Christol--Dumas,
*Fast computation of the Nth term of an algebraic series over a finite prime
field* (arXiv:1602.00545). Check the exact statements before using them.
