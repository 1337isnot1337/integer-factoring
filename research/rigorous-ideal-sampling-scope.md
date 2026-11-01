# Scope of the recent rigorous ideal sampler

The deterministic selector in
[the height-only relation interface](height-only-smooth-norm-relations.md)
remains open. A recent primary source supplies a strong randomized
sampler, but its stated theorem does not supply that selector.

## Primary-source statement

The source is Koen de Boer, Alice Pellet-Mary, and Benjamin Wesolowski,
*Rigorous methods for computational number theory*,
[arXiv:2512.01588v2](https://arxiv.org/html/2512.01588v2).
The checked sections are Theorem 6.4, Corollary 6.5, Theorem 7.1,
Proposition 8.12, and Theorem 9.5.

Under ERH, the sampler multiplies its input ideal by independent random
small prime ideals, applies an archimedean distortion, and samples a
lattice point uniformly from a large region. The resulting probability
of a smooth output is bounded below using the density of smooth ideals.
The random walk is close in total variation to uniform on an Arakelov
class-group coset. Theorems 1.1 and 9.5 are randomized statements.

For an imaginary quadratic field, the continuous distortion space has
dimension zero. The discrete random walk remains. With the full class
group and trivial modulus, Corollary 6.5 chooses

$$N=\left\lceil16+2\log(1/\varepsilon)+\log h_K\right\rceil$$

independent prime ideals. Proposition 8.12 still uses a randomized
uniform lattice-point sampler. Expected polynomial sampling cost is
not a bound on the number of possible samples.

## A support-size calculation

Let `G` be a finite set of size `h`, let `U` be its uniform probability
distribution, and let `mu` be supported on `S`. With the convention

$$\|\mu-U\|_1=\sum_{g\in G}|\mu(g)-1/h|,$$

one has

$$\|\mu-U\|_1\ge2\left(1-\frac{|S|}{h}\right). \tag{1}$$

Indeed, outside `S` the sum is `(h-|S|)/h`. On `S`, the sum of the
signed differences is the same number; the sum of their absolute values
is at least their signed sum. Adding proves (1).

Consequently, an `epsilon`-close sampler must satisfy

$$|S|\ge(1-\varepsilon/2)h. \tag{2}$$

If a sampler is implemented by a deterministic map of at most `r`
random bits, then `|S|<=2^r`, so

$$r\ge\log_2h+\log_2(1-\varepsilon/2). \tag{3}$$

This conclusion does not depend on how the samples are encoded, or on
whether different seeds produce the same class. It also applies to a
nonuniform distribution over a prescribed finite seed set: that seed
set still needs at least the number of elements in (2).

For any family of quadratic discriminants with

$$h\ge2^{b/2-o(b)},\qquad b=\mathrm{bits}(|\Delta|),$$

enumerating the full support of a constant-error near-uniform class
sampler costs `2^{b/2-o(b)}` visits. The class-number lower bound used in
[the norm-coverage note](norm-selector-coverage.md) gives such a bound
under this repository's pinned all-number-field analytic input.
Thus exhaustive support enumeration does not turn this particular
near-uniform sampling theorem into the required deterministic running
time.

The observation also explains why increasing the prime cutoff does not
resolve that conversion. A stronger spectral estimate can reduce the
number of walk steps, but any distribution that is close to uniform
still has the support-size requirement (2).

## What the calculation does and does not leave open

The argument concerns reproducing a near-uniform distribution, or
enumerating all of its possible classes. It is not a lower bound for
finding a successful smooth norm for one target.

A deterministic algorithm could choose a much smaller, target-dependent
set of samples that meets the success set. The density lower bound alone
does not identify such a set. A conditional-expectation construction
would need an efficiently computable estimate of the success probability
after each partial choice. Neither the support calculation nor the
existing orientation decoder supplies that estimate.

Likewise, the class-group relation problem only needs the designated
generator-power targets, not every ideal class. A short list tailored
to those targets is not excluded. The existing
[arbitrary-norm decoder](arbitrary-factored-norm-representation.md)
would process a supplied successful norm within the proved
subexponential overhead.

The recent paper therefore provides useful rigorous probability and
sampling results. An additional deterministic selection theorem is
still needed to use those results for a larger score here. The paper's
ERH hypothesis also cannot simply be replaced by the repository's
fixed zero-free strip without redoing its analytic estimates.

This note used substantial OpenAI assistant support. It contains an
elementary support calculation and a scope check, not a new factoring
algorithm or a lower bound for integer factoring.
