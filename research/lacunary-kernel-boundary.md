# What a lacunary characteristic-two kernel would still need

## Scope

Power-of-two mode weights produce a real graded state space with much
slower growth than ordinary partitions. This observation alone does not
produce a crossing identity. This note separates that construction from
two exact limitations of the most direct scalar-kernel replacement.

The weighted reciprocal pencil is investigated separately in
[weighted-energy-reciprocal-pencil.md](weighted-energy-reciprocal-pencil.md).
Its exact projected matrices have an asymptotic `rho=1` obstruction,
despite the small graded state count. None of these statements rules out
a different lacunary exchange relation or a multicomponent construction.

## 1. Divided powers preserve the binary-partition state count

Give a variable `x_j` weight `2^j`, for each `j>=0`. Over
characteristic two, its divided-power algebra has basis `x_j^[n]` and
multiplication

```math
x_j^{[a]}x_j^{[b]}={a+b\choose a}x_j^{[a+b]}.
```

The binomial coefficient is odd exactly when the binary addition of
`a,b` has no carry. Consequently the algebra is isomorphic, as a graded
algebra, to the tensor product of square-zero variables indexed by
`(j,k)`, where the variable for `x_j^[2^k]` has weight `2^(j+k)`.
Every nonzero binary-digit product has coefficient one.

There are `s+1` bit variables of weight `2^s`, and the degree-`E`
dimension is

```math
b(E)=[q^E]\prod_{s\ge0}(1+q^{2^s})^{s+1}
    =[q^E]\prod_{j\ge0}(1-q^{2^j})^{-1}.              (1)
```

The second equality is an identity of ordinary counting series, using
the unique binary expansion of each occupation number. At `E=2^R`,
the elementary estimates in the weighted-energy note give
`log b(E)=Theta((log E)^2)`.

This construction avoids two incorrect shortcuts. Using ordinary
factorial vacuum pairings in characteristic two makes occupations of
size at least two degenerate. Keeping only one square-zero variable at
each mode leaves only one basis state at each integer energy, by unique
binary expansion. The full divided-power digit presentation has neither
defect. It still needs creation and annihilation operators with a proved
exchange identity before its dimensions can serve as contraction widths.

## 2. The Artin--Schreier kernel and its exact scalar fixed field

In any characteristic-two constant field `F`, define the formal series

```math
K(z)=\sum_{j\ge0}z^{2^j},\qquad K(z)^2+K(z)=z.         (2)
```

It is not rational over `F(z)`. If a rational function `f` satisfied
`f^2+f=z`, any finite pole of `f` would remain a pole of twice the
order, which the right side lacks. Thus `f` would be a polynomial.
A nonconstant polynomial of degree `d` makes `f^2+f` have even degree
`2d`, whereas `z` has degree one. A constant also cannot work.

The nontrivial automorphism of `F(z,K)` sends `K` to `K+1`.
For a rational function `R(T)` with coefficients in `F`,

```math
R(K(z))\in F(z)
\quad\Longleftrightarrow\quad
R(T)=R(T+1)
\quad\Longleftrightarrow\quad
R(T)\in F(T^2+T).                                    (3)
```

For polynomials, write uniquely
`R(T)=A(T^2+T)+T B(T^2+T)`. Translation by one changes the
second term by `B`, so invariance is equivalent to `B=0`.
The rational-function statement follows from the same quadratic
extension and its involution.

Therefore, if `R` is defined at both zero and one, scalar
rationalization forces

```math
R(0)=R(1).                                           (4)
```

In particular, an additive polynomial
`L(T)=sum_i c_i T^(2^i)` makes `L(K)` rational exactly when
`sum_i c_i=0`, which also makes `L(1)=0`. The scaled version uses
the involution `T -> T+a` and fixed field `F(T^2+aT)`.

Equation (4) is a precise obstruction only when a proposed identity
needs to distinguish these two conjugate kernel values. It is not a
general obstruction to interpolation in the base variable `z`:
`K^2+K=z` itself is a rationalization that retains distinct base
values at `z=0` and `z=1`. Explicit base-variable dependence, extra
components, or different evaluation points require their own analysis.

## 3. There is no finite scalar locality multiplier

Consider the two-sided formal distribution

```math
\mathcal C(z)=\sum_{j\ge0}(z^{2^j}+z^{-2^j}).          (5)
```

No nonzero Laurent polynomial annihilates this distribution. Let
`L(z)=sum_(k=a)^b l_k z^k`, with `l_b != 0`. Choose `j` with
`2^(j-1)>b-a`. The coefficient of `z^(2^j+b)` in `L C`
is exactly `l_b`: the next smaller positive power is too far away,
larger powers require a shift outside `[a,b]`, and negative powers
cannot reach this exponent. Hence `L C` is nonzero.

Multiplication by a Laurent polynomial is legitimate coefficient by
coefficient because it uses only finitely many summands. No product of
two infinite two-sided distributions is used in this argument.

This closes the direct replacement in which the lacunary distribution
is the commutator kernel and a finite scalar polynomial is supposed to
restore ordinary locality. It does not exclude a finite-window identity,
a different operator algebra, or a matrix-valued locality relation.

For a truncated series `K_J=sum_(j=0)^J z^(2^j)`, the exact identity is

```math
K_J^2+K_J=z+z^{2^{J+1}}.                              (6)
```

The boundary term must be accounted for. A truncation argument may be
valid if its complete coefficient window proves the term irrelevant;
dropping it from a global exchange identity is not justified.

## 4. Remaining constructive target

A useful positive result would supply explicit operators on the
divided-power state space from Section 1, prove their crossing/exchange
identity, and bound the actual two-cut factors of its restricted
coefficients. One would then count the number of pieces, generation
cost, field cost, and repeated contraction loss. The small value of
`b(E)` supports investigating this target; it does not establish any
of those missing properties.
