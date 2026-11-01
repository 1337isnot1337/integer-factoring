# Linear supplied-data charges for oriented ideal refinement

## Results and scope

Two elementary charging arguments sharpen the earlier
[oriented basis construction](oriented-coprime-basis.md).

First, a supplied pairwise-coprime norm factorization can acquire all
required quadratic orientations using at most
`sum_i floor(log_2 R_i)` orientation gcds. Combining this step with
Bernstein's classical coprime-basis and factorization algorithms gives
a deterministic construction in

$$S\,(n+\log(S+2))^{O(1)}\quad\text{bit operations},\tag{1}$$

where `S=sum_i bit_length(R_i)` counts every supplied key, including
repeated norms, and `n` bounds the discriminant and individual norm bit
lengths. The output consists of sparse exact oriented ideal words.

Second, the unchanged incremental implementation makes at most

$$2\sum_i\lfloor\log_2 R_i\rfloor
  -\#\{i:R_i>1\}\tag{2}$$

sparse usage writes over its entire insertion history. Thus repeated
rewriting of old usages is not a hidden quadratic cost. Large product
arithmetic and downstream presentation work remain separately charged.

Both results process supplied keys. Neither proves that a deterministic
factoring search obtains enough useful relations within a new bound.
The complete score exponent remains `1/2`.

## 1. Input and classical norm factorization

Fix a negative quadratic-order discriminant `Delta`. Each input key
`(R_i,B_i)` has positive odd norm, is coprime to `Delta`, and satisfies

$$B_i^2\equiv\Delta\pmod{4R_i}.$$

Normalize each root modulo `2R_i`. The corresponding invertible ideal
is `I(R_i,B_i)=(R_i,(-B_i+sqrt(Delta))/2)`.

A norm factorization consists of integers `A_j>1`, pairwise coprime,
and sparse positive exponent entries `e_ij` satisfying

$$R_i=\prod_j A_j^{e_{ij}}.\tag{3}$$

Unused basis entries are discarded. No `A_j` is assumed prime or free
of powers. Bernstein's [2005 paper](https://cr.yp.to/lineartime/dcba-20040404.pdf),
Theorems 18.2 and 21.3, constructs a coprime basis and the factorizations
over it in essentially linear integer input length. The construction
uses multiplication, exact division, gcd, and equality tests. It is an
algorithmic ingredient, not a supplied prime-factorization oracle.

Norm values may be deduplicated for that classical stage. Each occurrence
of a key must still be counted and restored before orientations are
processed. In particular, many roots sharing one norm do not become
free input.

Equation (3) gives the weighted incidence inequality

$$\sum_{i,j:e_{ij}>0}\log_2 A_j
 \le\sum_i\log_2 R_i.\tag{4}$$

It also bounds the number of exponent entries and their encoding size.

## 2. One orientation pass per norm component

For a fixed `A=A_j`, choose the root `B_0` of any incident key and reduce
it modulo `2A`. Start with the one-piece list `[A]`. Process each other
incident root `B_i` once. For every current piece `d`, split it into

$$d_+=\gcd\left(d,\frac{B_i-B_0}{2}\right),\qquad
d_-=d/d_+,\tag{5}$$

discarding unit pieces. All root differences are even because they have
the parity prescribed by the discriminant.

To justify (5), let `p^v` be any prime-power component of `d`. The two
quantities `(B_i-B_0)/2` and `(B_i+B_0)/2` have product divisible by
`p^v`. They cannot both be divisible by `p`: otherwise `p` would divide
`B_0` and hence `Delta`. Therefore exactly one quantity contains the
whole `p^v` factor. It follows that

$$d_+d_-=d,\qquad\gcd(d_+,d_-)=1,$$

and the root is uniformly plus or minus on each resulting piece.
Splitting later into divisors preserves every previously established
orientation agreement. A restart over earlier roots is unnecessary.

The piece list is always a pairwise-coprime factorization of `A` into
integers greater than one. It has at most `floor(log_2 A)` entries.
Consequently the total number of gcds in all passes is at most

$$\sum_j\#\{i:e_{ij}>0\}\lfloor\log_2 A_j\rfloor
\le\sum_i\lfloor\log_2 R_i\rfloor.\tag{6}$$

The last inequality applies (3) to each input row and uses positive
integer exponents. Every gcd operand has at most `n+O(1)` bits after
root normalization.

For each final piece `d`, set `C` to the smaller of the two residues
`B_0` and `-B_0` modulo `2d`. Every incident key has root `C` or `-C`
there, so assign its old exponent the corresponding sign. Distinct
original components remain coprime. The resulting sparse word has the
literal identity

$$I(R_i,B_i)=\prod_{w_{ij}>0}I(d_j,C_j)^{w_{ij}}
 \prod_{w_{ij}<0}\overline{I(d_j,C_j)}^{-w_{ij}}.\tag{7}$$

For prime-power exponents in (3), this follows from uniqueness of the
unramified lift of a chosen quadratic root. Agreement at the component
selects the same prime ideal at every larger incident exponent. This is
the same literal-ideal argument used by the original oriented basis
theorem; negative entries denote conjugate integral factors in (7).

Each final output incidence consumes at least one bit of its input
norm. There are at most `S` nonzero signed entries, and each exponent
has magnitude at most `n`.

## 3. Charging validation and output

Pairwise coprimality need not be validated by repeatedly taking a gcd
against an ever larger running product. Let `P=prod_j A_j`. Compute
`P mod A_j^2` for every `j` using product and remainder trees. Since
`A_j` divides `P`,

$$\gcd\left(A_j,\frac{P\bmod A_j^2}{A_j}\right)
 =\gcd(A_j,P/A_j).\tag{8}$$

All norms are pairwise coprime exactly when every value in (8) is one.
The total tree input length is `O(S)`. Standard fast integer arithmetic
therefore performs this check within the bound (1).

Validate each sparse norm word by exponentiation and multiplication,
stopping if a partial product exceeds its claimed norm. All valid
intermediate products have at most `n` bits. Key validation, orientation
gcds, sorting, signed-word output, and reduction of the quadratic forms
each take `S` times a polynomial in `n+log(S+2)`. Literal ideal HNF
verification is also polynomial per sparse entry: its ideals have
bounded dimension and norm at most that of the input being checked.

Together with the cited constructive norm stage, these facts prove (1).
The output and temporary oriented lists use
`O(S(n+log(S+2)))` bits, apart from the classical norm algorithm's own
charged workspace. No dense key-by-atom matrix is required.

The [executable wrapper](../tools/oriented_norm_basis_lift.py) implements
(5), (8), sparse factorization validation, and independent literal HNF
verification. It accepts an explicit norm basis and norm words. It does
not implement or time Bernstein's norm-basis engine. The complete
complexity theorem uses that published engine; the executable receipts
test the orientation and validation steps precisely as described.

## 4. Amortized usage writes in the streaming implementation

The [incremental implementation](../tools/incremental_oriented_coprime_basis.py)
already reuses an old usage dictionary when its atom norm and canonical
root stay unchanged. Consider the history of one supplied key. Its
initial word has `r` entries, with atom norms `a_1,...,a_r`.

Give an occurrence of an atom of norm `a>1` the integer potential

$$\Phi(a)=\lfloor\log_2 a\rfloor-1\ge0.$$

When this occurrence is rewritten to `k` child occurrences, the norm
identity has the form `a=prod_j d_j^{f_j}`, with every `f_j>=1`.
If `k>=2`, then

$$\Phi(a)-\sum_j\Phi(d_j)\ge k-1.$$

The `k` new usage writes are at most twice this potential decrease.
If `k=1`, an actual replacement has `d<a`: a same-norm atom with its
canonical root would have been reused. Thus `a=d^f` for `f>=2`, and
the potential decreases by at least one, which pays for its one write.

Each occurrence has a unique old parent because distinct old atoms
have coprime norms. There is no merger of two old usages into the same
child for a fixed key. All successive replacement trees can therefore
be charged independently.

Their initial total potential is at most
`floor(log_2 R)-r`. The `r` initial writes plus every later copied write
are bounded by

$$r+2(\lfloor\log_2 R\rfloor-r)
 =2\lfloor\log_2 R\rfloor-r.$$

A nonunit key has `r>=1`, and a unit key has no usages. Summing proves
(2). The bound concerns actual dictionary-entry writes. It does not
count exports of all sparse words, quadratic-form reduction, product
tree arithmetic, or downstream relation presentation operations as free.

## 5. Reproducible checks

The [receipt](../evidence/oriented-norm-basis-lift.json) records 800
independent coprimality comparisons and 15 ordered small ideal batches.
A mixed-root control starts from one composite norm supported on eight
primes. Its 256 keys use all root-sign patterns and exponents through
five, reaching 186-bit norms. The orientation pass makes seven splits
and returns the eight separate prime supports; all literal ideals pass.

The generated controls verify 1,193 literal ideals and check (2) after
all 1,193 corresponding insertions. Seven malformed factorizations or
corrupted ideal certificates are rejected.

The complete 80-bit and 96-bit fixtures add 58,962 literal ideal checks.
Supplying their retained norm bases reconstructs every original atom and
signed word exactly. Coarsening equal unsigned incidence columns is also
checked; these particular fixtures have no mergeable columns. Their
retained 15,799 and 28,993 usage writes obey (2).

```sh
python3 tools/check_oriented_norm_basis_lift.py \
  --output /tmp/oriented-norm-basis-lift-replay.json
```

This implementation and analysis used substantial OpenAI assistant support.
