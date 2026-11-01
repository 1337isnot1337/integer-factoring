# Arbitrary binary codes and monomial gauges retain a three-quarter floor

## Result and scope

Let `C` be any finite set of binary words of a common length, with
`q=|C|>1`. Restrict all four physical legs of a tensor power of the
elementary crossing to `C`. The set may be an arbitrary code inside
a Hamming layer, a union of layers, or a code with no weight symmetry.

For every diagonal physical monomial gauge at both endpoints, the
common-type Smith-pressure parameter satisfies

```math
\boxed{\rho\ge\frac34>\frac1{\sqrt2}.}             \tag{1}
```

The opposite assignment of the two physical cuts has the stronger
floor `rho>=1`. These bounds hold over any coefficient field. The
exact finite replay uses characteristic two, as the factoring search
does.

This extends [the full-Hamming-type theorem](direct-type-elementary-pressure-floor.md).
It also covers another common-type selection on the codeword alphabet
and arbitrary codeword-dependent monomial potentials. Section 7 also
covers different codes on the two paired sets of legs. It does not
cover non-diagonal physical frames, independently chosen unpaired codes,
new scalar seeds, auxiliary extensions, or a direct coefficient-rank
method that bypasses these Smith pressures. It does not lower-bound
integer factoring.

## 1. The restricted crossing and its gauges

For binary physical labels define

```math
C_{xyuv}=[x=u][y=v],\qquad E_{xyuv}=[x=v][y=u],
\qquad T(z)=(1-z)C+zE.                            \tag{2}
```

Let `T_C` be `T^tensor(t)` restricted to the same code `C` on all
four legs. Its endpoint is the exact swap, so both relevant cut
matrices are nonsingular. The row cut is `(x,y)|(u,v)`; the column
cut is `(x,u)|(y,v)`.

Choose arbitrary integer functions `h_l^0,h_l^infinity:C->Z`, one
for each physical leg `l`. At zero multiply an entry by
`z^(sum_l h_l^0)`. At infinity use the Laurent matrix
`tau^(sum_l h_l^infinity) T_C(tau^-1)`.
Any uniform shifts used to make the matrices regular are restored
when their weighted Smith sums are computed.

Write `Phi_s(M)=sum_i s^(-nu_i(M))` for the finite Laurent Smith
orders at scale `s>1`, and put

```math
B(s)=\Phi_s(M_{\rm row}^{\infty})
      \Phi_s(M_{\rm column}^{0}),\qquad
G(c)=\sum_{l=1}^4(h_l^0(c)+h_l^\infty(c)).       \tag{3}
```

The standard common-type criterion at parameter `rho` would require

```math
B(s)<\left(\sum_{c\in C}s^{-G(c)/(4\rho)}\right)^{4\rho}.
                                                        \tag{4}
```

We prove the reverse inequality at `rho=3/4`. The function
`R log sum_c exp(-G(c) log(s)/R)` increases with `R>0`, since its
derivative is the entropy of a positive Gibbs distribution. Thus
the same certificate cannot attain any smaller `rho`.

## 2. Sum and difference blocks retain one cheap Smith factor

Treat words as vectors in `{0,1}^t` with integer addition. The row
matrix splits into blocks indexed by the integer sum

```math
S=x+y=u+v.
```

In a block, all pairs agree outside the positions where `S_i=1`.
If their number is `d`, the entry from the pair `(x,y)` to `(u,v)`
is `z^(d_H(x,u))(1-z)^(d-d_H(x,u))`. In particular the crossing
entry `(u,v)=(y,x)` is exactly `z^d`.

The column matrix splits into blocks indexed by the integer
difference

```math
D=x-u=v-y.
```

At its `d` nonzero coordinates the crossing is forced and contributes
`z^d`. Outside those coordinates, the entry is a power of `1-z`.
The same crossing entry `(x,y,u,v)=(x,y,y,x)` is again exactly
`z^(d_H(x,y))`.

Physical monomial gauges preserve both block decompositions.
The smallest Laurent Smith order of a nonzero block is its minimum
entry valuation. It is therefore at most the valuation of each
of the crossing entries just described. Consequently each sum block
contributes a term to the infinity pressure at least

```math
s^{d_H(x,y)-h_x^\infty(x)-h_y^\infty(y)
                 -h_u^\infty(y)-h_v^\infty(x)},   \tag{5}
```

for every pair in that block. Each difference block contributes a
term to the zero pressure at least

```math
s^{-d_H(x,y)-h_x^0(x)-h_y^0(y)-h_u^0(y)-h_v^0(x)}.
                                                        \tag{6}
```

We use only one Smith factor per nonempty block. Other factors and
cancellations inside a block can only increase the resulting bound.

## 3. An entropy inequality on the Boolean cube

For independent Boolean random vectors `X,Y`, with no assumption
of independent coordinates within either vector,

```math
H(X+Y)\ge\tfrac34(H(X)+H(Y)),\qquad
H(X-Y)\ge\tfrac34(H(X)+H(Y)).                     \tag{7}
```

Here addition and subtraction are coordinatewise integer operations.
All entropy logarithms may use any one fixed base.

First consider independent Bernoulli variables of parameters `p,q`.
Let `h` be binary entropy in bits and set
`a=p(1-q)`, `b=(1-p)q`. The only information lost by forming their
integer sum occurs when the sum is one. Hence

```math
h(p)+h(q)-H(B_p+B_q)=(a+b)h\left(\frac{a}{a+b}\right).
```

The expression is zero when `a+b=0`. Otherwise the elementary bounds

```math
4u(1-u)\le h(u)\le2\sqrt{u(1-u)}                 \tag{8}
```

give

```math
\begin{aligned}
(a+b)h\left(\frac{a}{a+b}\right)
&\le2\sqrt{ab}\\
&\le p(1-p)+q(1-q)\\
&\le\tfrac14(h(p)+h(q)).
\end{aligned}                                    \tag{9}
```

This proves the one-coordinate version of (7). For vectors, use
the entropy chain rule and condition each output coordinate on
both complete input prefixes. This is stronger conditioning than
knowing just their sum prefix, and the two next input bits remain
independent under it. Summing the one-coordinate inequality gives
the first claim in (7). Replace `Y` by `1-Y` for the second claim;
this changes neither its entropy nor independence.

For completeness, both bounds in (8) have short analytic proofs.
Write `x=|1-2u|`. The entropy expansion is

```math
h((1-x)/2)=1-\sum_{m\ge1}
 \frac{x^{2m}}{\log(2)\,2m(2m-1)}.
```

The positive coefficients sum to one. Since `x^(2m)<=x^2`, this
gives `h>=1-x^2=4u(1-u)`. For the upper bound put
`f(x)=sqrt(1-x^2)-h((1-x)/2)`. Both endpoint values are zero, and

```math
f'(x)=\frac{x}{\sqrt{1-x^2}}
\left(\frac{g(x)}{\log 2}-1\right),\qquad
g(x)=\frac{\mathrm{arctanh}(x)\sqrt{1-x^2}}x.
```

The function `g` decreases strictly from one to zero because

```math
g'(x)=\frac{x-\mathrm{arctanh}(x)}
             {x^2\sqrt{1-x^2}}<0.
```

Thus `f` first increases and then decreases to zero, so `f>=0`.

## 4. The gauge charges cancel in the entropy bound

Choose any probability distribution `p` on the code, and take
independent `X,Y` with that law. Let `S=X+Y`, `D=X-Y`.
The variational formula for a finite weighted sum, applied to
the block contributions (5) and (6), gives

```math
\log B(s)\ge H(S)+H(D)-\mathbb E[G(X)]\log s.
                                                        \tag{10}
```

The two expected Hamming-distance terms cancel. Each leg charge
has expectation under the same law `p`, which gives exactly the
symbol energy in (3). Inequality (7) now yields

```math
\log B(s)\ge3H(p)-\sum_{c\in C}p(c)G(c)\log s.
```

Maximizing over `p` is another finite Gibbs variational identity:

```math
\boxed{B(s)\ge
 \left(\sum_{c\in C}s^{-G(c)/3}\right)^3.}        \tag{11}
```

This proves (1) with all the stated physical gauges and codeword
potentials. Without gauges, it says simply `B(s)>=q^3` at every
scale, independently of the code's length or geometry.

## 5. The opposite orientation has a floor of one

We use a small Smith-majorization fact. Suppose a square Laurent
matrix has, for every index subset `I`, a principal minor of
valuation at most `sum_(i in I) a_i`. Then the sum of its `k`
smallest Smith orders is at most the sum of the `k` smallest `a_i`.
The determinantal characterization of Smith orders proves this
immediately by choosing that index subset. It follows that

```math
\Phi_s(M)\ge\sum_i s^{-a_i}.                     \tag{12}
```

One direct justification for the last step is that the partial-sum
inequalities imply
`sum_i (r-nu_i)_+ >= sum_i (r-a_i)_+` at every integer `r`.
The geometric summation of these inequalities gives (12), after
one common shift if any orders are negative.

At zero the ungauged row matrix is the identity. Every principal
minor therefore has valuation zero before applying its row and
column gauges. Apply (12) with the straight pairing `(u,v)=(x,y)`.
At infinity, permute the columns of the column-cut matrix to make
its endpoint at `z=1` the identity. Every principal minor is a
nonzero polynomial, so its degree is nonnegative and its infinity
valuation is nonpositive before the gauges. Apply (12) with the
crossing pairing `(y,v)=(u,x)`.

The product of the two resulting bounds is

```math
\begin{aligned}
B_{\rm opposite}(s)\ge{}&
 \left(\sum_c s^{-h_x^0(c)-h_u^0(c)}\right)
 \left(\sum_c s^{-h_y^0(c)-h_v^0(c)}\right)\\
&\cdot
 \left(\sum_c s^{-h_x^\infty(c)-h_v^\infty(c)}\right)
 \left(\sum_c s^{-h_y^\infty(c)-h_u^\infty(c)}\right)\\
\ge{}&\left(\sum_c s^{-G(c)/4}\right)^4,
\end{aligned}                                    \tag{13}
```

where the last inequality is Holder's inequality. This is the
common-type threshold at `rho=1`, proving the stronger claim.

## 6. Exact replay

The [checker](../tools/check_arbitrary_binary_code_pressure.py) verifies
the sum and difference block formulas directly against local tensor
products, including every nontrivial binary code of length at most
three. It also computes full Laurent Smith profiles for several
irregular codes and nonuniform physical gauges. Evaluation scales
that are twelfth powers of rational numbers make both Gibbs
thresholds rational, so the comparisons require no floating point.

The analytic entropy proof above supplies the all-code and all-scale
theorem. The finite checks are independent algebra and indexing
controls, rather than a numerical substitute for it. Evidence is
saved in [arbitrary-binary-code-pressure.json](../evidence/arbitrary-binary-code-pressure.json).

Run:

```sh
python3 tools/check_arbitrary_binary_code_pressure.py
```

## 7. Different paired code embeddings have the same entropy floor

Allow two distinct ordered codes
`A=(a_1,...,a_q)` and `B=(b_1,...,b_q)` in the same Boolean cube.
Embed the symbols on legs `x,v` through `A`, and on legs `y,u`
through `B`. The endpoint is still the exact swap of the abstract
`q` symbols. The two code sets may overlap, be disjoint, or have
unrelated weight distributions.

The sum and difference block decompositions still hold. Choose
independent abstract indices with any common law `p` and put
`X=a_i`, `Y=b_j`. Each random word has entropy `H(p)`. The crossing
entry is still `z^(d_H(X,Y))` in both blocks. Thus (10)--(11) give
the same `rho>=3/4` floor, with each leg's height evaluated on its
own embedding.

For the opposite orientation, use the same crossing entries again.
Their distance scores are now `-d_H(X,Y)` in the row-zero block
and `+d_H(X,Y)` in the column-infinity block. These cancel just
as before. Consequently **both cut assignments** satisfy

```math
B(s)\ge\left(\sum_{i=1}^q s^{-G(i)/3}\right)^3    \tag{14}
```

for every paired embedding and every stated monomial gauge.
The second code need not use the same ordering or the same words
as the first one.

There is a sharper opposite-orientation result for complement-aligned
codes, where `b_i=1-a_i` for every abstract symbol `i`. The leading
coefficient `E-C` of the binary column-cut matrix is a signed
permutation: its nonzero entry in row `(x,u)` has
`(y,v)=(1-x,1-u)`. Tensoring and restricting to the paired codes
preserve this signed permutation. Thus the degree-`t` reversal has
unit principal minors after that column permutation.

For the row-zero matrix, use the swap column permutation at `z=1`.
Every principal minor is a nonzero polynomial of degree at most
`t|I|`, so its zero valuation is at most `t|I|`. The minor lemma
from Section 5 gives two partition products with opposite scalar
factors `s^-t` and `s^t`. They cancel. Holder then gives

```math
B_{\rm opposite}(s)\ge
\left(\sum_i s^{-G(i)/4}\right)^4.                \tag{15}
```

The complement alignment matters for the gauged formula: it makes
the same abstract index collect the four marginal charges. An
arbitrary reordering of the complement code still has (14), but
(15) would use a correspondingly permuted type constraint.

The replay includes four nonuniform paired-code controls, including
disjoint codes, overlapping codes, and one complement-aligned pair.
It checks all forty full Laurent profiles and thirty exact scale
comparisons across the common and paired controls.
