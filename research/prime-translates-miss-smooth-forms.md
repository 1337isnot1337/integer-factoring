# Small prime translates can miss every recognizable smooth form

The [prime-form sample theorem](polylog-prime-form-generators.md) supplies
small character bias. It does not permit the smooth-pair extractor to be
replaced by testing whether a sampled translate is itself a recognizable
smooth reduced form. The following family disproves that replacement for
the actual arithmetic set, rather than for an arbitrary black-box subset.

For every odd prime `A`, let

$$\Delta=1-4A^4,\qquad Q=(A,1,A^3),\qquad R=A^2-1.$$

Let `S_R` be the set of reduced primitive forms of discriminant `Delta`
whose leading coefficient is at most `R`, coprime to `Delta`, and
`Y`-smooth. If `X<A` and `Y<A`, then

$$[Q]u\notin S_R$$

for every prime form `u` of norm at most `X` coprime to `Delta`, in
either orientation. In fact, for a prime form of norm `p<A`, the
translated reduced form has leading coefficient exactly `Ap`. Thus no
choice of probability weights on these
prime forms repairs the single-set hitting claim.

Since `bits(abs(Delta)) = 4 log_2 A + O(1)`, the conclusion applies to
every prime cutoff and smoothness bound of size `2^o(b)` along this
family, once `A` is large enough. It does not contradict the positive
Fourier guarantee for `S_R/S_R`.

This is an obstruction to that specific extraction shortcut. It is not
a factoring lower bound. In particular, the displayed discriminants
were selected for a transparent form calculation, not as hard factoring
instances.

## 1. The first two primitive norm scales are separated

More generally, suppose `C >= A >= 2` and

$$Q_C(x,y)=Ax^2+xy+Cy^2.$$

If `y=1` or `y=-1`, the integer term `Ax^2+xy` is nonnegative, so
`Q_C(x,y) >= C`. If `abs(y)>=2`, completing the square gives

$$Q_C(x,y)=A\left(x+\frac y{2A}\right)^2
          +\left(C-\frac1{4A}\right)y^2
 \ge4C-1/A>C.$$

Therefore every primitive representation of a positive integer below
`C` has `y=0` and `x=1` or `x=-1`, and its value is exactly `A`.
For the displayed family, `C=A^3`, so

$$0<Q(x,y)<A^3,\quad\gcd(x,y)=1
 \quad\Longrightarrow\quad Q(x,y)=A.\tag{1}$$

The target form is reduced and primitive. Also

$$4R^2<4A^4-1=|\Delta|,$$

so `S_R` lies in the same small-leading-coefficient range as the earlier
[smooth-pair interface](bounded-smooth-norm.md).

## 2. A prime translate would force an impossible smooth norm

There is a direct composition proof. Since `A` and `p` are distinct,
the product ideal has primitive norm `Ap`. Choose its centered middle
coefficient `B`, with `-Ap < B <= Ap`. The associated form is

$$\left(Ap,B,\frac{B^2+4A^4-1}{4Ap}\right).$$

For `p<A`, its trailing coefficient is strictly larger than `Ap`.
It is already reduced, so its leading coefficient is `Ap`, which is
coprime to `Delta`, at most `R`, and not `Y`-smooth. Every translate
fails exactly the smoothness requirement.

The represented-norm argument also explains the word extension below.

Suppose `[Q]u` lies in `S_R`, with `u` a prime form of norm `p<=X`.
Choose its small primitive ideal representative of norm `a<=R`.
Composing with the inverse prime ideal for `u` gives the class `[Q]`
and an integral ideal of norm `ap`. Removing its integer content gives
a primitive ideal, still in `[Q]`, of norm

$$n=ap/c^2\le ap\le XR<A^3.$$

Its prime factors are all at most `max(X,Y)<A`. Thus the reduced form
`Q` primitively represents a positive `max(X,Y)`-smooth integer below
`A^3`. Equation (1) says that this integer must be `A`, which is prime
and larger than that smoothness bound. This is a contradiction.

The argument covers both prime-ideal orientations, possible common
content, and nonmaximal quadratic orders. All ideals used have norms
coprime to the discriminant. It uses the same primitive-ideal-to-form
correspondence as the supplied-witness converter.

## 3. A necessary word length for the same shortcut

The same proof applies if `u` is supplied as a word of at most `ell`
oriented prime forms of norms at most `X<A`. Multiplying the
corresponding integral ideals and removing content gives a primitive
representative for `u` of norm at most `X^ell`, supported on primes
below `A`.

Consequently

$$R X^\ell<A^3
 \quad\Longrightarrow\quad [Q]u\notin S_R.\tag{2}$$

A word that does hit `S_R` must therefore satisfy

$$\ell\log X\ge\log(A^3/R)>\log A.$$

For `X=2^O(sqrt(b log b))`, this forces
`ell=Omega(sqrt(b/log b))`. The statement is about a norm-bounded
prime-word sampler and the recognizable set `S_R`; it does not exclude
other representations, pair extraction, larger norms, or arithmetic
methods that find a word without enumerating all words.

In particular, fixed-length convolution of subexponential prime samples
does not establish a universal single-set hitter in this family.

## 4. What remains valid

The small-bias theorem and the
[positive correlation theorem](positive-relation-kernel.md) continue to
give a translate in `S_R/S_R` when their numerical hypotheses hold.
That quotient allows a represented norm as large as `R^2`, near
`abs(Delta)/4`, which crosses the second norm scale `A^3` in this
example. There is no contradiction.

The [direct short-word interface](smooth-norm-short-words.md) also
remains valid: it accepts represented smooth norms beyond the small
radius and does not claim that the original reduced representative
must itself be smooth.

The missing arithmetic search cannot be filled simply by enumerating
the new short prime sample and smoothness-testing its reduced forms.
That precise search rule fails on every sufficiently large member of
the displayed family.

## 5. Exact finite controls

Run:

```sh
python3 tools/check_prime_translate_obstruction.py
```

The [checker](../tools/check_prime_translate_obstruction.py) constructs
the actual prime forms, composes them with the target, and tests the
resulting reduced form against the full definition of `S_R`. It also
checks all two-prime words satisfying the displayed product-norm
condition in its finite range, and primitive-coordinate controls for
(1). It does not infer a uniform theorem from those tests; Sections
1--3 provide the proof.

The [evidence record](../evidence/prime-translate-obstruction.json)
checks ten prime parameters, 242 prime translates, 248 two-prime words,
and 32,956 primitive coordinate controls. The largest parameter is
`A=1009`, with discriminant `-4145955690243`; all 148 sampled prime
translates in that case fail smoothness. The checker also verifies the
exact leading coefficient `Ap`, and its corresponding product formula
for the bounded two-prime words. No class-group enumeration or unknown factorization
is used by the counterexample checker. The existing small-form
composition routine is used only for finite controls, with no claimed
polynomial-bit implementation bound.

Substantial OpenAI assistant support was used. No Lean formalization
or stronger factoring score is claimed.
