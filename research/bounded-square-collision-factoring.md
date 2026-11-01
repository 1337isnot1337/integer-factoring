# Factoring from bounded-frontier square collisions

## Status

This deterministic finite search factors the recorded 51-, 60-, and
80-bit inputs directly from two class-group roots with the same square.
It needs no relation determinant or supplied annihilator. The arithmetic
certificates and all frontier costs are explicit.

The search uses a bounded, pruned frontier. No uniform success bound at
a useful width is proved, and the complete factoring score remains
unchanged. The square-collision principle uses classical ambiguous-form
factoring; the contribution is the executable frontier search, its
careful handling of inverse-key ties, and the finite certificates.

## 1. Square keys and helpful ambiguity

Let \(N>1\) be odd and let \(G\) be the proper class group of
discriminant \(-4N\). For a root class \(h\), use the key

\[
\kappa(h)=\min\{\operatorname{red}(h^2),
                    \operatorname{red}(h^{-2})\},
\]

where the minimum uses canonical form order only for exact lookup.
Store a sign \(s\in\{1,-1\}\) with \(h^{2s}=\kappa(h)\). If two
roots have the same key, their oriented quotient

\[
t=h_1^{s_1}h_2^{-s_2}
\]

satisfies \(t^2=1\). The classification in
[early-stop extraction](early-stop-class-factoring.md) supplies coprime
factors of \(N\) unless \(t\) is the identity or the exceptional class

\[
U=(2,2,(N+1)/2),\qquad N\equiv1\pmod4.
\]

These unhelpful classes form a group of size at most two. A collision
with one of them can be deduplicated. The root words and their orientations
are retained for every positive factoring certificate.

### The inverse-key tie must be checked

When \(h^2=h^{-2}\), the square key is itself ambiguous. It can already
be helpful even though identifying \(h\) with \(h^{-1}\) would merge
their inverse orbits. The implementation therefore tests every such key
for a factor **before** deduplication. If it is helpful, it supplies the
two explicit roots \(h,h^{-1}\), with words \(e,-e\); their quotient is
the helpful class \(h^2\).

If a self-inverse key is instead the identity or \(U\), the two possible
quotient orientations differ only by an unhelpful class. Choosing one
orientation is then safe. This check is needed for the complete-radius
guarantee below, not just for an optimization.

## 2. A deterministic bounded frontier

The inputs are \(N\), a factor-base size \(m\), a width \(W\), and a
depth bound \(D\). The driver constructs the first \(m\) odd split
unramified prime forms. It receives no input factors, class number,
annihilator, group decomposition, or random sample.

Starting from the identity, generate the \(2m\) signed-prime neighbors
of each retained root. Every new root has an explicit signed word and
a path depth. Compute its square key, check a helpful self-inverse key,
then check collisions against the retained index. Also check equal keys
among the newly generated candidates before pruning that layer.

Unhelpful duplicates are removed. Among the remaining keys, retain at
most \(W\) root representatives using the stratified leading-coefficient
policy from [bounded word navigation](bounded-class-word-search.md).
Keep the retained keys from earlier layers in the collision index.
Sorted arrays, sorting, and binary search make the lookup deterministic.
No group-compatible ordering assumption is used.

A helpful collision returns a proper divisor and both root words. If
the budget expires, the result is inconclusive. It does not certify
primality, full generation, or absence of helpful two-torsion.

## 3. Fully charged costs

There are at most \(2mWD\) generated neighbors. Each uses one class
composition for the move and one for squaring. Each generated neighbor
can cause at most one retained-index or same-layer collision comparison,
which uses one more composition. Thus the three recorded loop counts
satisfy

\[
C_{\rm move}=C_{\rm square}\le2mWD,\qquad
C_{\rm compare}\le C_{\rm move},\qquad
C_{\rm total}\le6mWD.
\tag{1}
\]

Testing a square key for ambiguity uses its existing coefficients and
gcd arithmetic; it does not require another class composition. Setup,
inversion, sorting, integer operations, and independent certificate
replay have their own polynomial costs. The counts in (1) are exact
frontier-loop counts, not complete wall-clock benchmarks.

There are at most \(1+WD\) retained keys and \(2mW\) temporary
candidates in one layer. A word has \(m\) coordinates of magnitude at
most \(D\). Given the factor-base forms, total bit time and memory are
polynomial in \(b,m,W,D\), with \(b=\lceil\log_2(4N)\rceil\), using
the exact polynomial-bit class operations. The depth and width are
charged by value. Constructing the first split primes is additionally
charged in the largest enumerated prime; the finite records use small
primes. No unproved polynomial bound for that prime enumeration is
silently inserted.

## 4. Completed finite constructions

| Quantity | 51-bit input | 60-bit input | 80-bit input |
| --- | ---: | ---: | ---: |
| \(N\) | 1621982935400141 | 1000000016000000063 | 906694351050364540708897 |
| \(m,W,D\) | \(9,1024,32\) | \(12,2048,40\) | \(12,4096,128\) |
| Stopping depth | 9 | 13 | 123 |
| Retained keys | 5700 | 19745 | 488737 |
| Move compositions | 102600 | 436308 | 11677262 |
| Squaring compositions | 102600 | 436308 | 11677262 |
| Comparison compositions | 39807 | 126455 | 2855737 |
| Total recorded compositions | 245007 | 999071 | 26210261 |
| Root-word lengths | 9, 9 | 5, 11 | 25, 31 |
| Returned factors | 32452843, 49979687 | 1000000007, 1000000009 | 824633708473, 1099511627689 |

Each successful quotient is the reduced form \((p,0,q)\) for the two
returned factors. The 80-bit factors differ by 274877919216. The factor
base in that run is

\[
7,11,19,23,29,31,37,53,61,79,83,89.
\]

The 80-bit roots have path depths 75 and 123 but shorter collected words
of lengths 25 and 31 because opposite generator moves can cancel. Their
canonical forms are

\[
(550575495617,204103641814,1665727829338)
\]

and

\[
(43695262369,-25027861526,20753987969114).
\]

Their common inverse-square key is
\((146381215778,-140749370886,6227895721157)\). The recorded orientations
give quotient \((824633708473,0,1099511627689)\).

These are finite certificates for this particular deterministic search.
They do not certify full generation by the small bases or establish a
new practical factoring benchmark. In particular, the 80-bit run charges
all 26210261 recorded compositions; success is not described merely by
the two short final words.

## 5. A complete-radius guarantee, with its width charge

Let \(H\) be the subgroup generated by the supplied prime forms. Suppose
\(H\) contains a helpful ambiguous class with a signed word of length at
most \(2r\). Split that word into two root words of length at most
\(r\), so that their quotient is helpful and their squares agree.
Therefore a complete radius-\(r\) search finds a factor.

The deduplication does not invalidate this statement. Before a factor is
found, an unhelpful collision identifies roots only under inversion and
multiplication by the unhelpful group \(H\cap\{1,U\}\). Both signed
generator moves are available, so this equivalence respects the neighbor
graph and the path-depth bound. If inversion would identify the two
helpful roots, their shared square is itself helpful; the explicit
self-inverse-key check handles exactly that case.

For a simple sufficient width, let

\[
J=\#\{\{z,z^{-1}\}:z\in H^2\}
  =\frac{|H^2|+|H^2[2]|}{2}.
\]

Width \(W\ge J\) prevents pruning before the search finds a factor.
If a diameter bound \(T\) for \(H\) is known and \(H\) contains a
helpful class, depth \(\lceil T/2\rceil\) consequently suffices at this
charged width. Under the explicitly conditional generator and diameter
theorem in [the fixed-strip note](fixed-strip-smooth-square.md), one may
take \(H=G\) and polynomial \(T\), after the usual prime-power
preprocessing. The quantity \(J\) can still be exponential in the input
length. The conditional short diameter is not a small-width guarantee.

### Why complete small-prime frontiers can already be exponential

There is a direct arithmetic injectivity window. For a nonzero signed
prime word \(v\), put \(n=\prod p^{|v_p|}\). If
\(\prod g_p^{2v_p}=1\), the oriented ideal of norm \(n^2\) is primitive
and principal. Thus

\[
n^2=x^2+Ny^2,\qquad \gcd(x,y)=1.
\]

Here \(y\ne0\), since otherwise primitivity forces \(n=1\). Also
\(x\ne0\), since \(n\) is coprime to \(N\). The integers \(n,N\) are
odd, so \(x,y\) have opposite parity. Hence

\[
n^2\ge\min(N+4,4N+1)=N+4.
\tag{2}
\]

For primes at most \(Y\), two radius-\(r\) words with the same square
key therefore differ only by the literal sign change of their exponent
vectors whenever

\[
Y^{4r}<N+4.
\tag{3}
\]

Indeed, their difference or sum has \(\ell_1\) norm at most \(2r\),
and a nonzero such vector would contradict (2). In particular, the
complete radius-\(r\) index contains at least

\[
2^{r-1}\binom mr
\tag{4}
\]

distinct square keys, counting just words with \(r\) distinct nonzero
coordinates of magnitude one, modulo simultaneous sign change.

For \(b=\log_2N\), \(Y=b^B\), \(m\asymp Y/\log Y\), and fixed
\(B>1\), choose \(r\sim b/(4B\log_2b)\) satisfying (3). The lower
bound (4) is

\[
2^{((B-1)/(4B)-o(1))b}.
\]

The maximum frontier width in a complete search through that radius
must be at least the total index size divided by \(r+1\), so it too is
exponential. This is a bound on complete word-ball enumeration, not a
lower bound for arithmetic algorithms in general. A useful pruned
frontier theorem would need to prove that its selection rule preserves
a helpful relation; expansion and short diameter alone do not do that.

The scale of (2) is attained by
\(N=3^{4r}-4\). The generator at 3 has a helpful square key at radius
\(r\), namely \((3^{2r},4,3^{2r})\), yielding the factors
\(3^{2r}-2\) and \(3^{2r}+2\). The displayed form has discriminant
\(-4N\), norm \(3^{2r}\), and middle coefficient \(-2\) modulo 6;
its oriented prime-power ideal therefore gives the class
\(g_3^{-2r}\). It is ambiguous. These inputs are elementary
difference-of-squares examples, used only to check the exact boundary
and the inverse-key tie. The width-one controls for \(r=1,\ldots,8\)
find that key exactly at depth \(r\), with no pruning.

## 6. Independent replay and reproduction

The checker independently reconstructs the root words, their signed
squares, the common key, and the ambiguous quotient. It verifies the
coprime-factor formulas directly from the quotient coefficients. The
three returned prime pairs are trial-prime certified only after
recovery. The checker supplies no factors or annihilator to the search.

It also converts the helpful signed word \(v\) to an ordinary integer
certificate. The supplied-word norm constructor receives only the word
and prime forms and returns primitive integers satisfying

\[
x^2+Ny^2=n^2,\qquad n=\prod p^{|v_p|},\qquad \gcd(n,N)=1.
\]

The two gcds \(\gcd(n-x,N)\) and \(\gcd(n+x,N)\) recover the same
coprime factors. This conversion is polynomial in the expanded word-norm
length, which is at most \(2D\log_2Y\), apart from the polynomial
uniform-cap constructor overhead. It is separate from the counted
frontier loops. For the 51-bit case one such certificate is

\[
n=1275148838046875,\quad x=-14573386450661,\quad y=31659888.
\]

Once these three integers are supplied, the norm identity and proper
factor gcds can be verified without class-group operations. The replay
record includes the corresponding integers for all three inputs.

Finite controls cover 128 odd inputs from 3 through 257. They recover
coprime splits for all 65 inputs with more than one distinct prime
factor, and report no helpful collision on the 63 prime powers. These
small outcomes are not promoted to a uniform theorem. The proper
subgroup generated by the prime form at 7 for \(N=87\) supplies an
explicit unsuccessful composite control. Width one also misses the
51-bit input within depth 32 and correctly reports an inconclusive
result.

The inverse-key regression uses \(N=6557\), width one, and depth two.
It must detect the helpful key \((81,4,81)\) directly. A separate
injectivity control checks all 560 sign-normalized four-prime words
from eight split primes at a discriminant just beyond the bound (3).

The records are
[51 bits](../evidence/bounded-square-collision-51bit.json),
[60 bits](../evidence/bounded-square-collision-60bit.json),
[80 bits](../evidence/bounded-square-collision-80bit.json), and
[the independent replay and controls](../evidence/bounded-square-collision-replay.json).

```sh
python3 tools/bounded_square_collision_factoring.py --n 1621982935400141 --primes 9 --width 1024 --depth 32
python3 tools/bounded_square_collision_factoring.py --n 1000000016000000063 --primes 12 --width 2048 --depth 40
python3 tools/bounded_square_collision_factoring.py --n 906694351050364540708897 --primes 12 --width 4096 --depth 128
python3 tools/check_bounded_square_collision.py
```

The code is
[`bounded_square_collision_factoring.py`](../tools/bounded_square_collision_factoring.py).
The 80-bit search is substantially larger than certificate replay; its
millions of generated candidates are fully included in the cost record.
