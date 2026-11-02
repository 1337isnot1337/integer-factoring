# A counterexample to the trace-character selection-entropy bound

## Status and result

Base: `dd84ac3`. This note records a counterexample to the attempted
selection-entropy bound in an uncommitted draft, and withdraws the
corresponding claim. The draft asserted that an exact trace-character
solution set `S` with `r` predicates of maximum total degree `D` satisfies

~~~
v*log2(D+1) >= v*w - log2|S|,                              (1)
~~~

and concluded that exact pointwise selection of the factor or prefix graph
costs `2^Omega(b)` independently of `r`. That conclusion is not established.
The inequality and its degree lemma are false. A later note,
[Weil bounds and the native-word pointwise-selection budget](trace-selection-weil-budget.md),
proves a correct rank/degree tradeoff for thin pointwise trace sets and
retains this counterexample in its receipt. The present note preserves the
explicit counterexample and the valid `F_2`-ANF part of the failed draft.
Contraction/runtime and weighted cancellation remain open as stated there;
the complete score remains **1/2**. The correction was posted at the
[coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338124).

The problem is a field mismatch. Over `F_2`, the trace character has a
small algebraic-normal-form (ANF) degree, but the Fourier inversion divides
by `2^r` and loses the indicator. Over an odd field, division is valid, but
the odd-field multilinear degree of a trace bit is not its `F_2`-ANF
degree. A trace bit can be an XOR of up to `v*w` coordinate bits (and
already up to `w` when `v=1`). Its odd-field degree can therefore be as
large as `v*w`, even when the written predicate has degree one.

## 1. The posted counterexample

Let `q=8`, represented as `F_2[U]/(U^3+U+1)`, so `w=3` and `v=1`. Use the
two linear predicates

~~~
f_1(U) = 4U,      f_2(U) = 7U,      b=(1,0).               (2)
~~~

Both predicates have total degree `D=1`. The trace bits and the selected
set are

~~~
U          0  1  2  3  4  5  6  7
Tr(4U)     0  0  1  1  0  0  1  1
Tr(7U)     0  1  1  0  1  0  0  1
1_S(U)     0  0  0  1  0  0  1  0
~~~

Thus `S={3,6}` and `|S|=2`. The left side of (1) is

~~~
v*log2(D+1) = 1,
~~~

while the right side is

~~~
v*w - log2|S| = 3 - 1 = 2.
~~~

So (1) predicts `1>=2`, which is false. Direct Mobius inversion over `F_3`
gives multilinear degree `3` for `1_S`, above both the claimed degree bound
`1` and the Reed--Muller lower bound `2`. Over `F_2` the same indicator has
degree `2`.

The failure of the degree lemma is visible with one predicate. Take
`F_2[U]/(U^3+U^2+1)` and `f(U)=U`, with target `1`. Then
`Tr(U)` is the XOR of all three coordinate bits, so its `F_3`-multilinear
degree is `3`, while `v*log2(D+1)=1`. Its solution set has size `4`, so the
one-predicate instance of (1) is `1>=1` and does not by itself refute the
inequality; it refutes the claimed degree upper bound. The two-predicate
case above is what refutes (1) itself.

## 2. Where the proof fails

The draft's Lemma bounded the ANF degree of `Tr(f)` by a Hamming-weight
count of the exponents in `f`. That count is valid for the `F_2`-ANF, not
for the `F_p`-multilinear ANF with `p` odd.

For Boolean variables `X_1,...,X_s`, the XOR is

~~~
X_1 xor ... xor X_s
  = (1 - product_(j=1)^s (1-2X_j)) / 2.                    (3)
~~~

Over an odd field the coefficient of `X_1...X_s` is nonzero, so the
odd-field degree is exactly `s`. Over `F_2` the same function has degree
one. A trace bit is an `F_2`-linear functional of the coordinate bits;
after choosing a basis it is an XOR of some subset of those bits. Its
odd-field degree can therefore be as large as `v*w`. In the one-predicate
example above, `Tr(U)` is the XOR of all three bits.

The Fourier identity used in the draft is correct over an odd field:

~~~
sum_(lambda in {0,1}^r) (-1)^(lambda.b)
  product_i (-1)^(lambda_i Tr(f_i(U)))
  = 2^r * 1_S(U).                                         (4)
~~~

Since `p` is odd, division by `2^r` recovers `1_S`. But the sign
`(-1)^(Tr(f_i(U))) = 1-2*Tr(f_i(U))` has the odd-field degree of the trace
bit, not its `F_2`-ANF degree. The Hamming-weight estimate does not bound
that degree. Over `F_2`, (4) reduces to `0=0` and the division is
impossible. Neither field supplies the claimed pointwise bound (1).

## 3. What survives: an `F_2` necessary condition

The `F_2`-ANF part of the draft is valid. Let `f` have total degree at most
`D`. In the `F_2`-ANF,

~~~
deg_2(Tr(f)) <= v*log2(D+1).                              (5)
~~~

Indeed, write a monomial as `c*product_j U_j^(e_j)`. Each Frobenius power
`U_j^(2^t)` is `F_2`-linear in the `w` coordinate bits of `U_j`. Hence
`U_j^(e_j)` has `F_2`-ANF degree at most the number of one-bits of `e_j`;
the trace adds Frobenius conjugates without increasing degree. Summing
over exponents gives

~~~
sum_j wt(e_j) <= sum_j log2(e_j+1)
              <= v*log2(1+D/v)
              <= v*log2(D+1).
~~~

Consequently, for an exact solution set cut out by `r` trace bits,

~~~
deg_2(1_S) <= r*v*log2(D+1).                              (6)
~~~

The Reed--Muller minimum-distance bound over `F_2` says that a nonzero
multilinear polynomial of degree at most `d` is nonzero on at least
`2^(v*w-d)` points. For `0<|S|<q^v`,

~~~
r*v*log2(D+1) >= v*w - log2|S|.                           (7)
~~~

This is only a necessary condition on the predicate count and degree. It
does not imply the `r`-independent bound (1). The Artin--Hasse invoice is
`2^r (K+2)^O(v)` with `K` at least the written degree parameter `D`.
With `E=v*w-log2|S|`, choosing `r` and `log2(D+1)` at the crossover
`r*v*log2(D+1) ~ E` gives a formal cost balance `2^O(sqrt(E))`. That is
compatible with the necessary condition and is far below `2^b` for the
factor graph. It is not an achieved representation: no evaluation circuit,
coefficient construction, or uniform selection algorithm is supplied.
The later Weil-budget note proves a stronger thin-set tradeoff,
`r+log2(D(D+1))>=h` for density at most `2^-h`; that argument, not (7),
supplies the valid pointwise selection budget and rules out the sublinear
parameter for the native thin sets.

## 4. Consequence for the factor and prefix graphs

The draft's conclusion that exact pointwise trace-character selection of the
factor graph, or of a dyadic prefix graph, costs `2^Omega(b)` independently
of `r` is withdrawn. The Reed--Muller lower bound on the intrinsic `F_2`
degree of the indicator still applies, but the trace-character interface
does not control that degree over the odd field used for Fourier inversion.
The later [Weil-budget note](trace-selection-weil-budget.md) supplies the
valid pointwise selection bound. What remains open there is contraction,
runtime, weighted cancellation, and representations outside that interface.
No complete factoring representation and no stronger complete score are
claimed; the score remains **1/2**.

## 5. Finite controls

The [checker](../tools/check_trace_anf_counterexample.py) verifies the two
explicit counterexamples, the full family of one- and two-predicate linear
trace systems over the two cubic fields, the XOR degree gap
`deg_2=1` versus `deg_3=s`, the `F_2` loss and `F_3` recovery of the
Fourier identity, and representative cases of the corrected `F_2` bound (5).
Run `make verify-trace-anf-counterexample`. These are finite controls for
the stated identities and arithmetic. They are not an independent proof of
the ordinary lemmas, an Artin--Hasse representation, or a factoring
improvement.

## 6. Credit and scope

The Artin--Hasse interface, its Fourier projection, and its cost invoice are
from [the Artin--Hasse note](artin-hasse-trace-projection.md) and the owner's
`dd84ac3` checkpoint. The correct thin-set tradeoff and its native-word
application are from the [Weil-budget note](trace-selection-weil-budget.md)
and its contributors. The failed selection argument and this correction were
developed by worker-1 with assistant support. The Reed--Muller and
Hamming-weight facts are standard. This note has no independent mathematical
review and no Lean formalization. It is a retraction and a negative result,
not a lower bound for factoring and not a score improvement.
