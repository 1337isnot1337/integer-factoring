# What one list of candidate norms can cover

The [fixed-norm orientation decoder](factored-norm-orientation-search.md)
removes coordinates from the search output. This note records a
counting limitation on a particular way of using that decoder.
It is not a lower bound for target-dependent norm selection or
integer factoring.

Fix a negative quadratic-order discriminant `Delta`, and consider
positive integers `n` coprime to `Delta`, with all prime factors
split. A single integer

$$n=\prod_{i=1}^{r}p_i^{e_i}$$

is primitively represented by at most `2^r` proper ideal classes.
Indeed, a primitive integral ideal of this norm must choose one
whole orientation at each rational prime, giving exactly the
`2^r` candidate ideals

$$\prod_{i=1}^r P_{i,\pm}^{e_i}.$$

Different ideals may have the same class, so this is an upper
bound on represented classes. Inert prime factors only decrease
the bound, since they preclude a primitive representation under
the coprimality hypothesis.

Consequently, if a list `A` of integers, each with `Omega(n)<=L`,
contains a primitively represented norm for **every** class in
`C_Delta`, then

$$|A|\ge\frac{h(\Delta)}{2^L}.\tag{1}$$

The list may depend arbitrarily on `Delta`. The restriction is that
the same list must work for all target classes at that discriminant.

## Consequence at the new short-norm parameters

Under the pinned analytic input, the squarefree-product injection
in [the powered-prime note](powered-prime-smooth-norms.md),
equations (5)--(7), gives

$$h(\Delta)\ge
 \left(\frac{2^v}{32vk}\right)^k,$$

where `v=(1/sqrt(2)+o(1))sqrt(b)`,
`k=(1/sqrt(2)+o(1))sqrt(b)`, and `b=bits(abs(Delta))`.
Since `kv=b/2+O(sqrt(b))`, this implies

$$h(\Delta)\ge2^{b/2-O(\sqrt b\log b)}.\tag{2}$$

For the multiplicity `L=(sqrt(2)+o(1))sqrt(b)`, equations (1)--(2)
therefore require

$$|A|\ge2^{b/2-O(\sqrt b\log b)}.\tag{3}$$

Thus the cheap orientation test does not turn one short universal
list of integers into an extractor for every class. Enumerating
a list satisfying this stronger universal specification is
exponential in the discriminant bit length.

There are two important limits to this statement. A norm selector
may depend on the target form, and the relation program only
needs the particular generator-power targets, not all classes.
Those targets number at most the size of the factor base. The
counting argument does not rule out a short list tailored to
that smaller family, or a different fast search for each target.
It also places no restriction on algorithms that choose later
candidates using results of earlier representation tests.

## The height-only and imprimitive interfaces

The stronger
[arbitrary-norm decoder](arbitrary-factored-norm-representation.md)
allows all representations, including imprimitive ones. For a fixed
coprime norm `n`, each split prime power `p^e` now gives `e+1` signed
class contributions. An odd inert prime of even exponent gives only
forced integer content; an odd exponent prevents representation.
Consequently a fixed odd norm is represented by at most

$$\prod_{p\mid n,\ p\ \mathrm{split}}(e_p+1)\le\tau(n)$$

proper classes. Class collisions can only reduce this count.

Suppose every candidate has `log_2 n<=(C+o(1))b`, with fixed `C>0`.
The elementary divisor bound proved in the decoder note gives

$$\tau(n)\le2^{(C+o(1))b/\log_2b}.$$

Combining this with (2), a single list that covers **every** class must
therefore satisfy

$$|A|\ge2^{b/2-(C+o(1))b/\log_2b}. \tag{4}$$

Thus removing both the multiplicity and primitivity promises does not
make a subexponential universal list sufficient. The target-dependent
and restricted-target limitations stated above still apply in full.
In particular, (4) does not rule out the desired selector for the
polynomially many weighted generator-power targets.

The orientation counting is already implemented and checked in
the fixed-norm decoder. Equation (2) uses the proved finite
squarefree injection and its existing parameter checker. No
new computational assumption or finite experiment is needed
for this deduction.

This note used substantial OpenAI assistant support. Its analytic
dependency is the same pinned theorem used by the powered-prime
note.
