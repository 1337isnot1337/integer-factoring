# The missing quantitative short-box theorem

The arithmetic framework now constructs and replays finite factoring
certificates, but it has no improved uniform runtime. This note specifies
the exact target family and a sufficient localized smooth-norm statement.
It proves the conditional transfer and keeps the two unproved inputs
separate: localization of witnesses, and coverage of factoring directions.

The [supplied-norm word theorem](smooth-norm-relations.md) is used throughout.
No balanced-divisor hypothesis is required by the relation construction.
Balanced pairs enter only the already-proved existence argument.

## 1. Effective parameters and the target family

Let N be an odd composite non-prime-power integer, and put b=bits(N).
Handle primality, perfect powers and the finite small-input range first.
Trial division through 2^ceil(sqrt(b)) either finds a factor or leaves at
most sqrt(b)+O(1) distinct prime factors. These steps have at most the
smaller arithmetic cost used below.

Consider the constant list of d<=12 with gcd(d,N)=1 and Delta=-dN a valid
negative discriminant. A nontrivial gcd with a candidate d is already a
factor. This list contains the multipliers in Lenstra--Pomerance Algorithm
10.1. Call a multiplier **good** when the source's prime-density hypotheses
and smooth-form density conclusion apply. Its multiplier theorem supplies
at least one good choice with the prescribed source parameters; the
algorithm need not identify which choice is good in advance.

Write ell_b=ceil(log2(b+2)), s_b=ceil(sqrt(b*ell_b)), and fix

```math
Y_b=2^{\lceil s_b/2\rceil},\qquad \delta_b=2^{-2s_b}.
```

These are explicit dyadic parameters, fixed before witness search. They
fall in the range of Lenstra--Pomerance Theorem 8.1 for all sufficiently
large b. Its x is (1/2)sqrt(N), independently of the multiplier. If
u=ln(x)/ln(Y_b), then u<=b/s_b and ln(b)<=ell_b*ln(2). Consequently

```math
u(\ln u+13\sqrt{\ln u})
\le (\ln2)s_b\left(1+\frac{13}{\sqrt{\ln b}}\right)
\le2(\ln2)s_b
```

once ln(b)>=169. Thus the source's density bound is at least delta_b.
The upper range condition follows from

```math
\frac{\ln Y_b}{\sqrt{\ln x\,\ln\ln x}}\longrightarrow1/\sqrt2<1.
```

Its lower range condition also holds. In the source's notation,
ln(v) is of order sqrt(b*ln b) and ln(w) is of order sqrt(b), so
v,w eventually dominate every fixed multiple of ln(N), while ln(v)
dominates every fixed multiple of ln(u). The remaining logarithmic ratio
in Theorem 9.1 grows as sqrt(ln u), and an auxiliary z=ceil(exp(b^2))
satisfies its lower bound. Thus that theorem's inequalities hold eventually.
Source
constants and all thresholds are effective; enlarge one fixed b_0 to cover
them. This is an asymptotic parameter choice, not a practical claim about
the size of that threshold. All finite smaller inputs can be handled
directly. No smoothness heuristic is introduced.

For each candidate discriminant:

1. Construct all prime forms at unramified splitting primes <=Y_b. Denote
   them g_1,...,g_m and let H be the subgroup they actually generate. Thus
   m<=Y_b. Full generation of the class group G is not assumed.
2. Let R=floor(sqrt(|Delta|-1)/2), X=R^2 and L=floor(log2 X). Let S_R be
   the reduced forms of smooth leading coefficient <=R, coprime to Delta.
   For good multipliers the source's counted subset gives
   |S_R|>=delta_b*|G|, and hence density at least delta_b within H.
3. Fix one uniform [biased-word construction](arithmetic-biased-words.md)
   with epsilon=delta_b/2, applied to these generators and the class-size
   bound |G|<|Delta|<2^(b+4). Let its indexed multiset have size J and words
   w_u. Put W=max ||w_u||_1. The existing construction gives
   J=2^O(s_b) and log2(1+W)=O(b+s_b). Retain multiplicities.
4. Put K=max(W+L+1,2^(2b)). Search the fixed list
   T=supp(mu) union {identity}, with the known words; the extra identity
   is only an additional search option. For each generator and sample form

```math
Q_{i,u}=\operatorname{red}(g_i^K u),\qquad u\in T.
```

This is the exact family in the conditional statement below. The sampler
and ordering are fixed before norm search. There is no input-dependent
advice, unknown group decomposition or choice of a favorable target after
search. Binary powering costs depend on bits(K), not its numerical value.

## 2. What is already proved about these targets

The [positive kernel](positive-relation-kernel.md) supplies, for every i
and good multiplier, some sample u in the original support with
g_i^K u=a b^-1 for a,b in S_R. The
[bounded-norm equivalence](bounded-smooth-norm.md) therefore supplies a
primitive vector for Q_{i,u} with

```math
0<n=Q_{i,u}(x,y)\le X,\qquad
\gcd(n,\Delta)=1,\qquad P^+(n)\le Y_b.
```

That witness happens to have a balanced norm split. The word theorem
accepts every smooth represented norm satisfying the displayed conditions,
including those without a balanced split.

None of these results gives a short coordinate bound. For a reduced form
Q=(A,B,C), completing the square gives

```math
(2Ax+By)^2+|\Delta|y^2=4An.
```

The existing norm bound implies only

```math
|y|\le\sqrt A,\qquad
|x|\le R/\sqrt A+\sqrt A/2.
```

These can be exponential in b. The associated full ellipse can contain
exponentially many lattice points. Global density of smooth classes does
not make those points local in the chosen reduced-form coordinates.

## 3. The precise unproved localization input

Fix a rational beta>1/2. For an effective constant C define a dyadic
height bound of order

```math
h_b=2^{\lceil Cb/\ell_b^\beta\rceil}.
```

Fixed rational powers can be rounded with integer root computations; an
equivalent bound within a fixed factor in the exponent suffices.

**Uniform short-box hypothesis U_beta (unproved).** There are effective
constants C,b_0 such that for every b>=b_0, every good discriminant above
and every generator index i, at least one target Q_{i,u}, u in T, has a
primitive integer vector satisfying

```math
|x|,|y|\le h_b,\quad
0<Q_{i,u}(x,y)\le X,\quad
\gcd(Q_{i,u}(x,y),\Delta)=1,\quad
P^+(Q_{i,u}(x,y))\le Y_b.
```

This is an analytic/arithmetic hypothesis, not a consequence claimed from
the class-group density theorem. No balanced-divisor condition is hidden
in U_beta, and no promise is made for a particular individual sample.

A stronger quantitative route can be stated as a localized counting
estimate. For each target let M_{i,u}(h) count primitive vectors in the
displayed box with the norm and coprimality conditions, and let F_{i,u}(h)
count those additionally Y_b-smooth. A sufficient estimate is

```math
\sum_u\mu(u)F_{i,u}(h_b)
\ge\eta_b\sum_u\mu(u)M_{i,u}(h_b)>0,
\qquad \eta_b\ge2^{-C_1s_b}.
```

This weighted local smoothness estimate is also unproved. M counts the
actual clipped ellipse, so this formulation does not incorrectly assume
every target has h_b^2 admissible points. The global smooth-class density
delta_b cannot simply be substituted for this localized eta_b.

## 4. Conditional algorithm and complete cost

Under U_beta, enumerate the fixed target list and the signed coordinate
box for each i until finding a valid norm. Trial division through Y_b
recognizes smoothness and supplies its factorization. The word theorem
gives v_n with ||v_n||_1<=L. Hence

```math
r_i=K e_i+w_u-v_n
```

is a relation with correction norm <=W+L<K. All rows are strictly
diagonally dominant, so their determinant E_H is nonzero and is a multiple
of |H|. This uses no relation-span growth or parity-distribution premise.

There are m(J+1)=2^O(s_b) targets. Norm evaluation and coprimality tests
have polynomial bit cost because both the coefficients and coordinate
bit lengths are O(b). The complete search costs at most

```math
m(J+1)(2h_b+1)^2\,Y_b\,\operatorname{poly}(b)
=2^{O(b/\ell_b^\beta)}.
```

Factor-base construction and sampler generation are at 2^O(s_b) cost.
Every relation entry has O(b+s_b) bits, and bits(E_H) is
O(m(b+s_b))=2^O(s_b). Exact determinant computation has polynomial cost
in this encoded input, still 2^O(s_b). The index correction
lcm(1,...,ceil(1/delta_b)) has that same scale and yields a global
class-number multiple for each good discriminant.

If good multipliers are not individually certified, running all candidates
produces a list containing a valid corrected multiple; it does not label
every candidate's correction as valid. That distinction does not obstruct
checking an actual recovered integer divisor.

## 5. The separate splitting hypothesis

U_beta supplies relation matrices, not missing ambiguous directions.
The additional sufficient condition is:

**Coverage hypothesis C (unproved).** For every preprocessed input, at
least one good multiplier has a nontrivial image of H[2] in the coprime
factorization group of N.

There is a proved, computable sufficient test once the subgroup basis has
been extracted. Write d_odd for the odd part of d and define

```math
\kappa(\Delta)=
\begin{cases}
1,&\Delta\text{ odd or }\Delta\equiv4\pmod{16},\\
4,&\Delta\equiv0\pmod{32},\\
2,&\text{otherwise}.
\end{cases}
```

Lenstra--Pomerance Theorem 2.5 gives kernel size kappa for the homomorphism
from G[2] to coprime factorizations of the odd part of Delta. Forgetting
the known prime-power blocks of d_odd contributes kernel size
2^omega(d_odd), because gcd(d,N)=1. Hence the kernel dimension of the
projection to factorizations of N is exactly

```math
k_d=\log_2\kappa(\Delta)+\omega(d_{\rm odd}).
```

If dim_F2 H[2]>k_d, the restriction cannot have trivial image. This
criterion requires only d, the discriminant residue and the extracted
subgroup rank; it does not require factoring N. It is sufficient, not
necessary. For d=4 it says that one nontrivial torsion direction suffices
when N==3 mod4, while rank at least two suffices when N==1 mod4. The
stored 60-bit and 80-bit certificates respectively have ranks one and two
in exactly these two cases. They therefore also satisfy this coverage
test independently of a claim that their factor bases generate G.

Containing all of G[2] is sufficient but stronger than C. For a multiplier
satisfying C, the deterministic two-group algorithm recovers a basis of
H[2]; some basis element then yields a proper divisor. Trial division above
bounds its rank by O(sqrt b), so extraction costs 2^O(s_b) with the stated
relation-data sizes. The Delta=-348 example in the
[density/index note](arithmetic-density-lift.md) shows why C cannot be
silently omitted.

Consequently U_beta together with C gives a uniform deterministic factoring
algorithm at 2^O(b/(log b)^beta) bit cost. Recursing on proper divisors adds
only a polynomial factor: the recursion has O(b) nodes, and the displayed
cost is increasing for all sufficiently large input lengths. All finite
exceptions can be handled directly. This is a **conditional transfer**;
neither U_beta nor C has been proved here. Removing balancing also removes
the exact divisor-enumeration cost from this relation route, so the
transfer itself imposes no beta<=1 ceiling.

## 6. Exact localization-gap controls

The [fixed-box records](../evidence/arithmetic-norm-search-exhaustion.json)
contain two concrete targets where a small searched box is exhausted but
a supplied smooth witness exists outside it. For the target at generator
prime 197 of N=1208925818180866012192151, no accepted witness is found in
0<=x<=4096, 1<=y<=64 (nor by the leading-coefficient test), while
(x,y)=(225,89) gives smooth norm 37449200274818585. The separated-factor
input's generator-prime-127 target similarly exhausts height 256 and has
the outside witness (204,359), of norm 139743385530235837.

Both norms have prime-form word length nine. The replay checks the failed
boxes and the outside norm words exactly. It also verifies Y_b*n<=R^2
for an upper bound on every norm in each failed box, so any smooth value
there would already pass greedy balancing. These are failures of localized
norm search, rather than rejections caused by balanced-divisor output.
These finite examples do not
disprove U_beta: they do not instantiate its asymptotic constants, sampler
or signed box. They do demonstrate the distinction between having a
bounded smooth represented norm and locating one in a specified box.

Substantial OpenAI Codex assistance was used. This note gives an ordinary
conditional algorithm proof and an explicit open analytic target, with no
new Lean formalization or improved complete factoring score.
