# Enlarging the coherent spin crossing without losing uniformity

## Construction and open target

The original [bounded-spin crossing](bounded-spin-crossing.md) compresses
the pencil `(1-z)I+zS`. Its full eight-site coefficient span is
[closed by an exact certificate](bounded-spin-rank-algebra.md).
The following six-vertex construction has a larger span: seven rather
than six dimensions at eight sites, with genuinely noncommuting cut
operators. The earlier certificate therefore does not apply to it.

This is an all-length strict crossing with polynomial entry generation,
degree `O(n^2)`, the same exact binary-width encoders, and an integer
companion of scale `12^n`. A favorable uniform two-cut rank bound remains
open. The construction does not establish an improved factoring score.

## 1. Six local weights and three balanced products

Allow exactly the following nonzero entries of a binary four-leg tensor:

| Bits `x,y,u,v` | Weight |
| --- | --- |
| `0000` | `a` |
| `1111` | `f` |
| `0101` | `b` |
| `1010` | `e` |
| `0110` | `c` |
| `1001` | `d` |

Every path encoder ends at spin zero and is supported on words of weight
`n/2`. In any surviving four-word term, the four weight equations force
equal numbers of the `a/f`, `b/e`, and `c/d` vertices, respectively.
For example, subtracting the `u` and `v` weights from the `x` and `y`
weights first gives equality of the two crossing orientations; subtracting
individual weights then gives equality of the identity orientations.
The total number of zeros and ones gives equality of `a/f`.

Consequently the compressed tensor depends only on

```math
U=af,\qquad V=be,\qquad W=cd,
```

and is homogeneous of degree `r=n/2` in those products:

```math
\mathcal G_n(U,V,W)=\sum_{A+B+C=r}T_{A,B,C}U^A V^B W^C.       \tag{1}
```

The same primal and paired dual encoders as before are used. No square
root of any product is required to define the coefficient tensors.
At `(U,V,W)=(1,0,1)` the underlying local tensor is the swap, so

```math
\sum_{A+C=r}T_{A,0,C}=E.                              \tag{2}
```

## 2. A univariate crossing exposing the entire span

Set `M=n+1` and choose actual local polynomial weights

```math
a=f=d=1,\qquad b=e=1-z,\qquad c=z^M.
```

At `z=1` these are exactly the swap entries. The resulting complete
crossing is

```math
\boxed{\widetilde G_n(z)=
 \sum_{A+B+C=r} T_{A,B,C}(1-z)^{2B}z^{MC},\qquad
 \widetilde G_n(1)=E.}                              \tag{3}
```

Its degree is at most `r(n+1)`. More importantly, its coefficient span
is **exactly** the full monomial tensor span in (1), over every allowed
coefficient field.

To prove this, fix `C`. Its terms have exponents in
`[MC,MC+2(r-C)]`. Windows for different `C` are disjoint because
`M>2r`. Inside one window the polynomials `(1-z)^(2B)` have distinct
degrees and leading coefficient one. Descending elimination recovers
every `T_(A,B,C)` from the univariate coefficients. Conversely every
univariate coefficient is a linear combination of those tensors.
This triangular argument works in characteristic five even when some
intermediate binomial coefficients vanish.

The old pencil is the more restrictive substitution
`U=1,V=(1-z)^2,W=z^2`. Its overlapping windows can collapse tensor
directions; (3) removes that collapse with only a polynomial degree cost.
This is a new coefficient space, not another basis choice within the old
six-dimensional eight-site space.

## 3. Entry algorithm and complete construction costs

The entry constructor tracks the four magnetic states from the original
encoder. At each site it visits only the six vertices in the table.
In addition it tracks two nonnegative counters: the number `ell` of
`b/e` vertices, and the number `C` of `c` vertices. Give both `b/e`
vertices a temporary formal weight `t`, and the `c` vertex weight `W`.
The other three weights are one.

At the terminal state, weight balance makes `ell` even and the result
`t^ell W^C` becomes `T_(r-ell/2-C,ell/2,C)`. Thus at most `4^4`
magnetic states and `O(n^2)` counter pairs suffice. There are `n`
steps and six local vertices, giving `O(n^3)` field operations per full
entry. Expanding (3) also takes `O(n^3)` elementary operations.
The implementation is
[`bounded_spin_six_vertex.py`](../tools/bounded_spin_six_vertex.py).

For exactly `2^lambda` labels, use the original path unranking and
selection with `n=O(lambda)`. Degree, coefficient count and entry
construction cost are polynomial in `lambda`.

The same norm clearing gives the integer identity

```math
\widetilde H_n=12^n\widetilde G_n,\qquad
\widetilde H_n(1)=12^n E.                            \tag{4}
```

Although the exponents now reach `O(n^2)`, each word contribution has
at most `n` factors `1-z`; its binomial coefficient magnitudes are
bounded by `2^n`. The former bound `31104^n` on every cleared coefficient
therefore still applies. Coefficient bit lengths remain `O(n)` and the
integer entry construction has polynomial bit cost.

These costs establish an effective crossing construction. Effective
low-rank factors, a favorable uniform rank sum, and any required
fixed-characteristic factoring transfer remain separate obligations.

## 4. Exact finite structure and a smaller rank problem

The optional NumPy checker independently compresses coordinate signature
blocks. Fix `(x&y,x xor y)`, let `s=|x xor y|`, and for a crossed column
`(v,u)` put `h=|x xor u|`. Both are even. The entry of that block is
the single monomial

```math
U^{r-s/2}V^{(s-h)/2}W^{h/2}.                        \tag{5}
```

This differs from the magnetic-state algorithm, providing an independent
entry check. Complete tensors at lengths four, six and eight have:

| Length | Path labels | Formal monomials | Actual tensor span | Degree bound |
| --- | ---: | ---: | ---: | ---: |
| 4 | 2 | 6 | 3 | 10 |
| 6 | 5 | 10 | 4 | 21 |
| 8 | 13 | 15 | 7 | 36 |

At length eight take the six old basis tensors
`(C0,C2,C3,C5,C4,C6)` and the additional tensor `T_(0,0,4)`.
They are independent. The latter has rank pair `(34,34)`.
Both normalized cut spaces contain noncommuting operators.

There is nevertheless an exact smaller module: in each 169-dimensional
cut space, the images of all commutators span a 28-dimensional subspace,
and that subspace is invariant under every coefficient operator. Hence
the quotient has dimension 141 and the induced coefficient operators
commute. This is a finite module identity, not an asserted all-length
decomposition or a favorable asymptotic rank bound.

For a basis making this invariant subspace first, each coefficient
combination has block form `[[A,B],[0,C]]`, with `A` only 28-by-28.
Elementary elimination gives

```math
\operatorname{rank}\begin{pmatrix}A&B\\0&C\end{pmatrix}
 =\operatorname{rank}C+\operatorname{rank}[A\;B K],   \tag{6}
```

where the columns of `K` form any basis for `ker C`. Thus a primary
decomposition of the commuting quotient can reduce exact rank queries
to small matrices. Neither commutativity nor categorical dimension alone
would justify replacing those ordinary rank calculations.

The next structural question is whether the invariant commutator modules,
their quotients, and the two-cut ranks in (6) admit useful recurrences as
`n` grows. The eighth-length extra direction is substantial evidence
that the original pencil's primary-algebra certificate cannot simply be
reused; it is not evidence of a new factoring exponent.

## Replay

Run `make verify-bounded-spin-six-vertex` (NumPy required only for this
checker). It compares forty complete monomial entries with the independent
magnetic constructor, checks sixteen integer companion entries, verifies
the complete crossing endpoints, and checks full-rank univariate
coefficient recovery. It recovers the old eight-site coefficient tensor
entry by entry and verifies the seven-dimensional span and invariant
commutator modules. The
[receipt](../evidence/bounded-spin-six-vertex.json) pins the sources.
The construction and checks used substantial Codex AI assistance.
The complete factoring score remains `1/2`.
