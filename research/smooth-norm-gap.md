# Polynomial smooth-norm conversion with a radius gap

The [bounded smooth-norm interface](bounded-smooth-norm.md) converts a
supplied primitive smooth represented norm into a quotient of two small
smooth forms. Its exact balancing step costs `2^O(b/log b)` in the worst
case. A positive gap between an inner and an outer radius makes this step
polynomial, without the `sqrt(Y)` radius loss of greedy balancing.

This is a supplied-witness conversion theorem. It neither locates a smooth
represented norm nor proves a smoothness-density estimate at a smaller
radius. The strongest complete factoring score remains `1/2`.

For relation rows alone, [smooth-norm word extraction](smooth-norm-relations.md)
now removes the balancing requirement altogether, at the original radius.
The present theorem addresses applications that require two explicit small
forms as output.

## 1. An integer-only approximate divisor algorithm

Suppose a positive integer `n` is supplied as a list of prime factors
`p_1,...,p_k`, with repetitions. For every integer `z >= 1`, a deterministic
algorithm finds a divisor `d` such that

```math
d\le d_*\le(1+1/z)d,
\qquad d_*=\max\{e:e\mid n,\ e^2\le n\}.
```

Its bit cost is polynomial in `bits(n)` and `z`. It uses no floating-point
logarithms, real-number comparison, or enumeration of the complete divisor
set. This is the standard list-trimming approximation idea for subset
selection, applied directly to integer products; no new general
approximation-method claim is made. See the classical
[Ibarra--Kim paper](https://doi.org/10.1145/321906.321909) for the subset-sum
approximation background. The proof here is self-contained.

For `k=0`, return `1`. Otherwise put

```math
D=2kz,\qquad q=(D+1)/D,\qquad T=\lfloor\sqrt n\rfloor.
```

Start with the increasing list `L_0=[1]`. At step `j`, form and sort

```math
C_j=L_{j-1}\ \cup\ \{p_j a:a\in L_{j-1},\ p_j a\le T\}.
```

Remove duplicates. Scan the candidates in increasing order. Keep the first
candidate and then keep `v` exactly when `D*v > (D+1)*u`, where `u` is the
most recently kept candidate. The result is `L_j`. Return its largest
element after the last step. Every retained integer is a divisor of `n`.

### Approximation proof

After step `j`, every subset product `a <= T` of the first `j` factors has
a retained representative `a'` satisfying

```math
a'\le a\le q^j a'.
```

This holds at `j=0`. If a product at step `j` uses `p_j`, multiply the
representative of its preceding product by `p_j`; otherwise keep the
preceding representative. The resulting candidate is no greater than `a`,
so it is not removed by the cutoff. Trimming replaces it by a no-larger
retained value within a factor `q`. This proves the induction. Truncating
at `T` loses no feasible product because every item is at least two.

Apply the invariant to `d_*`. The returned maximum `d` is no smaller than
its representative, and it is itself feasible. Finally, entirely rational
inequalities give

```math
q^k
=\sum_{j=0}^k {k\choose j}(2kz)^{-j}
\le\sum_{j=0}^{\infty}(2z)^{-j}
=\frac{2z}{2z-1}
\le1+\frac1z.
```

The middle comparison uses `{k choose j} <= k^j`; the last uses `z >= 1`.
Thus `d <= d_* <= (1+1/z)d`.

### Bit cost

Successive retained values differ by a factor greater than `q`. Bernoulli's
inequality gives `q^D >= 2`. Therefore every `D` positions multiply the
retained value by more than two. If `b=bits(n)`, every retained list has at
most `D*b+1` entries. There are at most `b` prime occurrences, all candidate
values have at most `b` bits, and `D=2kz`. Sorting, multiplication and the
integer trimming comparisons therefore have polynomial bit cost in `b,z`.
The provided implementation retains the actual divisors, so no separate
solution-reconstruction oracle is needed.

The same argument works for arbitrary integer items at least two, but then
the optimum is over their subset products. Primality is what makes every
divisor available in the factoring application.

## 2. Inner/outer radius guarantee

Let `R_in,R_out,z` be positive integers satisfying

```math
(z+1)R_{\rm in}\le zR_{\rm out}.
```

If `n` has a factorization `n=a*b` with `a,b <= R_in`, the preceding
algorithm returns a split with both parts at most `R_out`.

To see this, interchange the two inner factors if necessary so that
`a <= sqrt(n)`. Then `d_* >= a` and `n/d_* <= R_in`. Hence

```math
d\le\sqrt n\le R_{\rm in}\le R_{\rm out},
\qquad
\frac nd\le(1+1/z)\frac n{d_*}
\le(1+1/z)R_{\rm in}\le R_{\rm out}.
```

The algorithm need not be told `R_in`. It can generate its approximate
divisor and check exactly whether the two factors fit `R_out`. Failure of
that check is not an exact negative certificate at the outer radius.

For example, `220=11*20` fits exact radius `20`. With the sorted prime list
`[2,2,5,11]` and `z=1`, trimming returns `10`, giving `10*22`. A rejection at
radius `20` is therefore a rejection of this approximate output only. The
same example occurs for the genuine primitive form `(11,1,220)`, of
discriminant `-9679`, at the vector `(0,1)`.

## 3. Constructing the two forms directly

Let `Q=(A,B,C)` be primitive and positive definite, with discriminant
`Delta < 0`. Suppose `(x,y)` is primitive,

```math
n=Q(x,y)>0,\qquad \gcd(n,\Delta)=1,
\qquad n=de,\quad d,e\le R_{\rm out},
\qquad 4R_{\rm out}^2<|\Delta|.
```

Choose integers `u,v` with `xv-yu=1`. In the properly equivalent form whose
first column is `(x,y)`, the middle coefficient is

```math
B'=2Axu+B(xv+yu)+2Cyv.
```

Thus `B'^2 == Delta mod 4n`. Center `B'` modulo `2d` to a residue
`b_d` in `(-d,d]`, and center `-B'` modulo `2e` to `b_e` in `(-e,e]`.
Define

```math
Q_d=\left(d,b_d,\frac{b_d^2-\Delta}{4d}\right),
\qquad
Q_e=\left(e,b_e,\frac{b_e^2-\Delta}{4e}\right).
```

Both forms are integral and primitive. Their trailing coefficients exceed
their leading coefficients because `4R_out^2 < |Delta|`, so the displayed
forms are already reduced. Their classes satisfy

```math
[Q_d][Q_e]^{-1}=[Q].
```

One way to verify the class identity is to use the proper ideal represented
by the equivalent form `(n,B',C')`. At each prime dividing `n`, the common
root `B'` fixes one unramified prime-ideal orientation, since `n` is
coprime to `Delta`. Dividing its norm into `d` and `e` splits the exponents
in that same orientation. The two resulting ideals multiply to the norm
`n` ideal, even when `gcd(d,e)>1`. The second displayed form uses the
conjugate orientation, so their quotient gives the original class. This
is the prime-ideal factorization used in
[the existing conversion proof](bounded-smooth-norm.md), expressed directly
through the two known divisors.

No prime-form enumeration, modular square-root search or class-group
composition is required by this constructor. Extended gcd and the displayed
integer formulas have polynomial bit cost. Independent finite checks below
use class-group composition to test the identity; that composition is not
part of the conversion algorithm.

## 4. Consequence for smooth-pair witnesses

Let `S_R` denote the reduced forms with leading coefficient at most `R`,
coprime to `Delta`, and `Y`-smooth. Suppose

```math
(z+1)R_{\rm in}\le zR_{\rm out},\qquad
4R_{\rm out}^2<|\Delta|.
```

Every supplied primitive norm witness for membership in
`S_(R_in)/S_(R_in)` can now be converted to an explicit pair in
`S_(R_out)/S_(R_out)` in

```math
Y\,\mathrm{poly}(b)+\mathrm{poly}(b,z)
```

bit operations, where the arithmetic inputs have `O(b)` bits. Trial division
through `Y` supplies the smooth factorization if it is not already given.
The inner pair supplies a balanced divisor by the existing bounded-norm
equivalence. Section 2 yields an outer split, and Section 3 constructs the
pair. For bounded primitive representations of a reduced form with
`n <= R_in^2 < |Delta|/4`, the coordinates and all formulas have
`O(log|Delta|)` bits.

The factor-base words of both output forms still have length at most
`log_2(R_out)`: their prime factors come directly from the represented
norm. Thus the [positive relation kernel](positive-relation-kernel.md)
can use an inner smooth subset to prove the existence of a pair, while an
extractor returns an outer pair and keeps the same logarithmic word bound.
The index-transfer argument also only needs both pairs to lie in the
factor-base subgroup.

This removes the converter's `2^O(b/log b)` term when the inverse radius gap
`z` is polynomial in `b`. In particular it would preserve faster prospective
arithmetic search bounds, rather than imposing a score-1 bottleneck itself.
It does not make any such search bound available.

The density premise must concern `S_(R_in)`. A density estimate for
`S_(R_out)` cannot simply be reused after shrinking the radius. This note
does not assert that the exact Lenstra--Pomerance density statement quoted
elsewhere in the repository already supplies every requested inner radius.
Fast norm search, the appropriate inner-density input, and deterministic
coverage of the missing ambiguous directions remain separate obligations.

## 5. Exact evidence

Run:

```sh
python3 tools/check_smooth_norm_gap.py
```

The [retained record](../evidence/smooth-norm-gap.json) contains:

- `16,384` exact approximation controls covering all integers through
  `4,096` and four values of `z`, with complete divisor enumeration used
  only as the independent finite oracle;
- `1,024` controls with reversed item order;
- `5,092` promised smooth represented norms over ten discriminants and
  three radius gaps, all converted to the correct class quotient;
- `5,092` nonprimitive scaled controls after coordinate normalization;
- `9,148` direct-form identity checks where the two norm factors share a
  prime;
- a `337`-bit square with `3^32 = 1,853,020,188,851,841` divisors, whose
  exact optimum is its known square root. The algorithm processes
  `6,578,976` candidates and retains at most `104,258` values; the full
  divisor set is not enumerated.

The last control is an exact algorithmic example, not a factoring benchmark
or an asymptotic proof by experiment. The all-input guarantee is Sections
1--4. The implementation is
[smooth_norm_gap.py](../tools/smooth_norm_gap.py); the independent checker is
[check_smooth_norm_gap.py](../tools/check_smooth_norm_gap.py). The work starts
from repository commit `a5e05feb756b4ad35ae1877d54a5f2e4eb387c19`, uses substantial
OpenAI assistant support, and adds no Lean formalization.
