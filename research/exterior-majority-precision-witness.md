# A majority-state precision witness for the exterior Cauchy lift

## Result and scope

Replacing the Boolean multilinear Cauchy deformation by its exterior-power
lift does not, by itself, remove the large precision requirement of local
jet recovery. There are explicit row and column tables with entries only
zero and one, and with coefficient array exactly of the form `A(Y)B(X)`,
for which the exterior deformation has value one at zero but vanishes to
the required high two-adic precision at every active interpolation center.

The construction uses a full packet's existing ports. It requires no
auxiliary selector ports, signed table entries, or arbitrary linear
combination of minors.

**Exterior witness theorem.** Let `L=2^a`, let `u>=1`, and let `v`
be the least odd integer at least `u`. Assume `v<=a`. For any positive
integer `g`, put `M=gL` and use `g` full packets. There are `0/1`
tables `A,B` for the exterior lift such that

```math
F(0)=1,\qquad \deg F=\frac34ML,
```

and, for every nonzero `t` in `{-L,...,L-1}` with `i=v_2(t)<u`,

```math
v_2\bigl(F(t+2Lz)\bigr)
\ge \frac{ML}{8v\,2^i}
\ge \frac{M D_t}{16u},
\qquad D_t=L/2^i.
\tag{1}
```

For a polynomial, `v_2` in (1) means the minimum valuation of all
its coefficients. Thus the statement concerns its entire local
polynomial, not only its constant value or a short initial jet.
If `u` is odd, the first bound has denominator `8u`.

The theorem concerns the **proposed exterior lift** defined below.
It does not establish this witness for the original Boolean lift,
for the special multiplication/divisor tensors, or for all factoring
algorithms. It is a limitation on recovering `F(0)` from only the
specified local jets. The pinned input states themselves can be
recognized and evaluated directly.

The Cauchy matrix and active-center geometry come from
[OpenAI's quantitative contraction source](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex)
at revision `4eaee87548250b2682affda7d05d947c05554bae`.
The exterior lift is an alternative investigated here; it is not
the deformation asserted by that source. Its actual deformation
uses products of the Boolean linear functions `1-z` and `z`.

## 1. The proposed exterior lift and its actual input states

For one packet, use slots `0,...,L-1` and set

```math
S=\{-L,-L+1,\ldots,L-1\},\quad
C(T)=\prod_{t\in S}(T-t),\quad
W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)}.
\tag{2}
```

The denominator's linear factor cancels a numerator factor.
The source proves `W(0)=I` and coefficientwise integrality of
`W(t+2Lz)` at every `t in S`.

Identify a Boolean port assignment with its occupied subset.
Define

```math
\Gamma(W)_{Y,X}=
\begin{cases}
 \det W_{Y,X},&|Y|=|X|,\\
 0,&|Y|\ne|X|.
\end{cases}
```

The empty determinant is one. Since `Gamma(I)=I`, the expression

```math
F_{\rm ext}(T)=\sum_{X,Y}A(Y)B(X)\,\Gamma(W(T))_{Y,X}
\tag{3}
```

has the original contraction as its value at zero. Equation (3)
is the proposed alternative deformation studied by this note.

Let `I` be the slots whose lowest `v` binary bits contain strictly
more than `v/2` ones. Because `v` is odd, complementing those bits
interchanges the two outcomes, so `|I|=L/2`. Choose

```math
A(Y)=[Y=I],\qquad B(X)=[X=I].
\tag{4}
```

Then (3) is the single principal minor

```math
F_I(T)=\det W_{I,I}(T),\qquad F_I(0)=1.
\tag{5}
```

Both tables are literal `0/1` delta tables. Their joint coefficient
array is rank one. They can also be written as products of local
bit-pinning predicates at the row and column vertices. All `L`
packet slots remain real ports; the bound counts occupied and
unoccupied ports alike.

For several full packets, use the same state in each packet and
the product delta tables. The matrix is block diagonal, so the
deformation is `F_I(T)^g`. Its local contents multiply by `g`.
It remains to prove the one-packet bound.

## 2. Exact local content of an arbitrary principal minor

This section applies to any subset `I`, of size `h`. Put

```math
V(I)=\prod_{\substack{x,y\in I\\x<y}}(y-x).
```

The Cauchy determinant formula gives

```math
F_I(T)=
\frac{(-1)^{h(h-1)/2}V(I)^2}{C'(0)^h}
\frac{C(T)^h}{\prod_{x,y\in I}(T+x-y)}.
\tag{6}
```

Every difference `y-x` lies in `S` and occurs at most `h` times
in the denominator. Thus (6) is a polynomial. Its degree is

```math
h(2L-h).
\tag{7}
```

The factorization can be verified by ordinary Cauchy elimination;
the same determinant identity underlies
[the exact Cauchy channel note](cauchy-two-adic-smith.md).

Fix `t in S` and write `alpha=2L=2^(a+1)`. The grid is a complete
residue system modulo `alpha`, so

```math
v_2(C'(t))=v_2(C'(0)),\qquad
v_2\bigl(C(t+\alpha z)\bigr)=a+1+v_2(C'(0)).
```

For the second equality, the factor indexed by `t` contributes
`alpha z`. Every other factor has content `v_2(t-q)<=a`.
The content of a product equals the sum of the contents: after
removing powers of two, the product of two nonzero polynomials
over `F_2` remains nonzero.

Applying that fact to the canceled factors in (6) gives the exact
integer formula

```math
\beta_t:=v_2\bigl(F_I(t+\alpha z)\bigr)
=h(a+1)+2v_2(V(I))
 -\sum_{x,y\in I}\min\bigl(v_2(t+x-y),a+1\bigr),
\tag{8}
```

where `v_2(0)=infinity`. This equality accounts for the normalized
Vandermonde numerator as well as every denominator. Omitting that
normalization would give an incorrect content estimate.

## 3. Rewrite the content as residue discrepancies

Let

```math
c_{k,r}=\#\{x\in I:x\equiv r\pmod{2^k}\},\qquad
C_k(t)=\sum_{r\bmod2^k}c_{k,r}c_{k,r+t}.
```

Counting each power of two dividing a difference gives

```math
2v_2(V(I))
=\sum_{k=1}^{a-1}\left(\sum_r c_{k,r}^2-h\right),
\qquad
\sum_{x,y\in I}\min(v_2(t+x-y),a+1)
=\sum_{k=1}^{a+1}C_k(t).
```

At level `a`, all residue counts are zero or one, so
`sum_r c_(a,r)^2=h`. Substitution in (8) yields

```math
\boxed{
\beta_t=
h-C_{a+1}(t)
 +\sum_{k=1}^{a}
  \left(\sum_r c_{k,r}^2-C_k(t)\right).
}
\tag{9}
```

Every term is nonnegative. The boundary term is nonnegative
because `C_(a+1)(t)<=h`. The other terms satisfy

```math
\sum_r c_{k,r}^2-C_k(t)
=\frac12\sum_r(c_{k,r}-c_{k,r+t})^2.
\tag{10}
```

If `v_2(t)=i`, its residue modulo `2^(i+1)` is `2^i`.
The term at level `k=i+1` is therefore

```math
S_i=\sum_{r=0}^{2^i-1}
 (c_{i+1,r}-c_{i+1,r+2^i})^2,
\qquad \beta_t\ge S_i.
\tag{11}
```

Only this one term is needed for the lower bound. Equation (9)
also explains why a multiscale occupancy pattern is useful:
its imbalances at different bit levels contribute to different
active centers.

## 4. Majority gives a large discrepancy at every active level

Return to the majority subset in Section 1. For `i<v`, define

```math
\Delta_r=c_{i+1,r}-c_{i+1,r+2^i}.
```

Let `b(x)=[x in I]`. Summing these differences gives

```math
\sum_{r=0}^{2^i-1}\Delta_r
=\sum_{x=0}^{L-1}b(x)(-1)^{\operatorname{bit}_i(x)}.
\tag{12}
```

Write `v=2m+1`. Changing bit `i` from zero to one changes the
majority outcome exactly when the other `2m` bits contain `m`
ones. Higher slot bits do not affect the predicate. Thus

```math
\frac1L\left|\sum_r\Delta_r\right|
=\frac12\frac{\binom{2m}{m}}{4^m}.
\tag{13}
```

An elementary bound suffices. Put `a_m=binom(2m,m)/4^m`.
Then `a_0=1`, `a_(m+1)=a_m(2m+1)/(2m+2)`, and induction gives

```math
a_m^2\ge\frac1{4m+1}.
```

The induction step follows from the exact identity

```math
(2m+1)^2(4m+5)-4(m+1)^2(4m+1)=1.
```

Consequently the square of (13) is at least
`1/(16m+4)`, which is larger than `1/(8v)`. Cauchy--Schwarz on
the `2^i` differences gives

```math
S_i=\sum_r\Delta_r^2
\ge\frac{(\sum_r\Delta_r)^2}{2^i}
\ge\frac{L^2}{8v\,2^i}.
\tag{14}
```

For every active center, `i<u<=v`, so (11)--(14) prove the
first inequality in (1) for one packet. Since `v<=u+1<=2u`,
the second follows. The degree from (7), with `h=L/2`, is
`3L^2/4`. Multiplying `g` packet factors proves the theorem.

The source schedule has `a>=2u` at every active call, so it
allows the required odd `v` without changing the packet size.
The witness therefore fits the same relation between the
packet parameter and the active-center cutoff.

## 5. What this rules out, and what remains open

Suppose a recovery rule receives only the local polynomials

```math
F(t+2Lz)\pmod{2^{P_t}},\qquad v_2(t)<u.
```

If all budgets satisfy

```math
P_t\le\left\lfloor\frac{M D_t}{16u}\right\rfloor,
```

the witness above supplies zero at every one of those inputs.
The zero polynomial supplies exactly the same data, but its
value at zero has different parity. Thus such local data cannot
determine `F(0) mod 2` uniformly over even the rank-one `0/1`
state subclass of the exterior lift. At some center, a uniform
recovery method must exceed a budget of this order, or use
additional information.

This conclusion does not rely on a small bound for the output:
the witness already has `F(0)=1`. It also avoids the realizability
gap of an arbitrary linear combination of Cauchy minors.
The single minor in (5) has exactly the input-state form required
by (3), and uses the existing full packet ports.

The original Boolean deformation is a different polynomial.
For the same delta tables, it is a product of row sums and their
complements, rather than a determinant. Equations (6)--(14) do
not apply to that polynomial. Nor does this construction prove
that the source's divisor-prefix tensors generate such a minor.
Those restrictions remain separate questions.

This is also not a runtime lower bound. An evaluator can inspect
the pinned states and return one directly, and another recovery
method might use information beyond these local jets. The result
closes the specific suggestion that exterior lifting alone makes
the prescribed uniform local recovery lose only a much smaller
precision, while preserving the rest of that approach.

## 6. Exact finite verification

Run

```sh
python3 tools/check_exterior_precision_witness.py \
  --output evidence/exterior-precision-witness.json
```

The checker verifies majority counts, exact residue discrepancies,
the content formula, and the integer form of (1) for packet
lengths through 256. For smaller packets it independently expands
the principal determinant by permutations, compares it with the
Cauchy product formula, expands every active local polynomial,
and measures its coefficient valuations directly. It also checks
that repeated packets multiply the local content exactly.

The completed run checks 5,496 active local contents, 100 majority
discrepancy bounds, 366 direct Cauchy-content calculations, five
expanded permutation determinants, 40 expanded local polynomials,
and 12 repeated-packet cases.

The evidence records the minimum normalized content
`beta_t*u*2^i/L^2` in each case and the exact checker hash.
It contains no floating-point arithmetic. These finite controls
verify signs, normalization, and parameter endpoints; the
all-size statement follows from the proof above.
