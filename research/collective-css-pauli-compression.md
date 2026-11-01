# Collective CSS compression of a polynomial crossing

## Research question and stopping test

Can an error-detecting physical code remove low-order Pauli channels from
both cuts of a tensor-power crossing, with enough gain to compensate for
its smaller physical alphabet? The concrete first test uses four copies
of the three-letter cyclic cell and a nine-dimensional ternary CSS code.
This changes the physical image across copies. It is not a rank-one
correction to the retained five-letter seed, a coordinate downset, or a
boundary-only modification preserving a protected natural core.

The all-copy pure-correction floor assumes one pure correction and the
fixed five-letter natural gauge. The homogeneous middle-fiber floor
covers its saved 120 restrictions. Neither proves a floor for this
collective code. Coverage by a broader code obstruction is not assumed.
The code's error detection may improve jets, but the rate loss and
complete high-order tails must be charged. No favorable rate is presumed.

First establish exact split maps and the one-site annihilation identity,
construct the entire compressed polynomial crossing without truncation,
and compute complete orders at both places with independent controls.
Compare normalized rank-pressure against the uncoded tensor power.
Continue only if the code yields a useful structural profile improvement
or a threshold candidate. If the loss of alphabet outweighs the detected
errors, do not enumerate nearby checks or larger code lengths absent a
new rate-versus-distance mechanism. The complete factoring target is
`rho<1/sqrt(2)` with effective transfer and every cost retained; the current
score is still `1/2`.

## Explicit algebraic code

Work over `F4`, with `omega^2+omega+1=0`, and use ternary coordinate
labels `s=(s0,s1,s2,s3)`. Put `a=(1,1,1,1)` and `b=(1,1,2,2)` in `F3^4`.
They satisfy `a.b=0`. Partition the hyperplane `b.s=0` into the nine
cosets of the line `<a>`. For each coset C, define `P[s,C]=[s in C]`.
The supports are disjoint and each has three entries, so in characteristic
two **`P^T P=I9`**. The left inverse is `L=P^T`; no factorial division or
complex inner product is used.

For every nonidentity one-site ternary Weyl operator W, `P^T W P=0`.
A nonzero shift changes the b-check because every b-coordinate is nonzero.
A pure nonzero phase sums to `1+omega+omega^2=0` along each a-orbit.
These two cases also cover a simultaneous phase and shift. This is an
ordinary identity over F4, not an inference from complex quantum states.

The stabilizer/CSS viewpoint is prior coding theory; see Ketkar,
Klappenecker, Kumar and Sarvepalli, [Nonbinary stabilizer codes over finite
fields](https://arxiv.org/abs/quant-ph/0508070). The finite-field split,
polynomial contraction and pressure question here are given explicitly.
No claim of a new error-detecting code or external novelty is made.

Apply P to all four external legs of four independent cyclic cells.
Because L=P^T and LP=I, the endpoint remains the exact nine-letter
crossing. All 6,561 polynomial entries and the full degree twelve remain.
Neither single-site error detection nor a leading-rank calculation alone
constitutes a rank or pressure bound for the compressed tensor.

## Exact outcome: the code reproduces a Frobenius-rescaled cell

Let `T(z;x,y,v,u)` denote the cell in matrix order `(x,y)|(v,u)`.
Identify each nine-letter label with `(i,j)` in `F3^2`, and put
`M(i,j)=(i,i+j)`. Let `sigma` square the coefficients in F4. The full
contracted polynomial satisfies the exact identity

```
T_code(z;x,y,v,u)
  = sigma((T tensor T)(z^2; Mx,My,Mv,Mu)).
```

This is an equality of all **6,561 complete polynomials**, including their
zero entries. It is not only a coincidence of orders. Each side has
degree twelve. In particular all odd powers in the coded tensor vanish.
The checker expands all 531,441 literal coset summands, with packed
polynomial arithmetic separate from the vectorized constructor, and
compares the result with two independently multiplied cell entries at
every external quadruple. The source cell is separately checked by
independent Weyl actions on all 81 entries. The exact finite coefficient
identity proves the displayed statement over F4 and all of its extensions;
it does not depend on numerical tolerances.

The zero orders are
`(0^1,2^4,4^8,6^16,8^20,10^16,12^16)`.
At infinity the raw orders are
`(-12^1,-10^12,-8^40,-6^24,-4^4)`.
They are twice the two-copy cell orders. The first nonzero relative
order doubles from one to two at both places, while the four physical
cells encode only two logical cells. Thus this gain does not improve the
normalized bare pressure. More strongly, coefficient automorphism and
matched symbol permutation preserve every weighted cut order, and
`z -> z^2` just doubles the positive valuation scale. Arbitrary physical
heights and scales on this coded source therefore give exactly the same
support family as unrestricted heights on the uncoded two-cell source.
There is no new free-height optimization problem to run here.

The [receipt](../evidence/collective-css-pauli-controls.json) also retains
all 32 one-site nonidentity Weyl annihilations and **324 independent
untruncated orders**, covering both places at zero heights and unequal
integer row/column weights. It is an exact control of this specified
code and identity, not an obstruction to every collective code. Repeating
this four-to-two block code, or composing it with itself, cannot change
the source family: tensor products, coefficient Frobenius, relabeling and
monomial parameter rescaling preserve the identity.

The stopping condition is met. Do not sweep neighboring checks or merely
increase code length on this evidence. A new coding proposal must explain
how it avoids the exact rescaled-cell identity and improves the charged
rate-versus-profile tradeoff. The complete factoring score remains `1/2`.

```
make verify-collective-css-pauli
```

Directed by Cinna Davis / 1337isnot1337 with substantial OpenAI assistant
support. Coding theory, the credited cyclic source, and existing exact
Smith engines are reused; no novelty claim is made for those ingredients.
