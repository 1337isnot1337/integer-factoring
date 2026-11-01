# A fixed smooth norm needs only a short orientation search

The remaining norm finder does not have to return representation
coordinates. It may first choose an integer norm. When the norm has
few prime factors, the missing orientations and primitive coordinates
can be recovered within the arithmetic time budget.

Let `Q` be a primitive reduced positive definite form of discriminant `Delta<0`.
Suppose

$$n=\prod_{i=1}^r p_i^{e_i},\qquad e_i\ge1,\qquad
 \gcd(n,\Delta)=1,$$

with distinct primes `p_i`, is supplied in factored form. Also supply
one prime-form orientation at each split prime, using a root
representative in `[0,2p_i)`. There is a
deterministic algorithm that either finds a primitive representation
`Q(x,y)=n` or proves that none exists, in

$$2^{\lceil r/2\rceil}
  \mathrm{poly}(\mathrm{bits}(|\Delta|),\mathrm{bits}(n))\tag{1}$$

bit operations and space. It uses canonical reduced forms and does
not compute the class-group order or its invariant factors.

This is a meet-in-the-middle specialization of the classical ideal
orientation criterion. A recent primary account is
[Mamah, arXiv:2502.11402v2, Theorem 5.4 and Algorithm 2](https://arxiv.org/abs/2502.11402v2),
stated there for a squarefree input norm. The proof below includes
arbitrary prime powers under the primitive and coprime hypotheses,
and explicitly charges the prime-form input. No new generic
subset-product algorithm is claimed.

## 1. Primitive ideals have one orientation per rational prime

If an inert prime divides `n`, no primitive representation is
possible. The prime ideal above an inert rational prime `p` is the
rational ideal `(p)`, so any integral ideal containing that factor
has nontrivial rational-integer content.

For a split prime, write `(p_i)=P_i conjugate(P_i)`. Unique ideal
factorization away from the conductor applies because
`gcd(n,Delta)=1`. A primitive integral ideal of norm `n` must choose
exactly one orientation at every prime:

$$I=\prod_{i=1}^r P_{i,s_i}^{e_i},\qquad
 P_{i,+1}=P_i,\quad P_{i,-1}=\overline{P_i},\quad
 s_i\in\{+1,-1\}.\tag{2}$$

where negative orientation means the conjugate ideal, not a
fractional ideal. It cannot contain both orientations above the
same prime, because their product is the rational ideal `(p_i)`.
Conversely, every product in (2) is primitive and has norm exactly
`n`.

Under the proper ideal/form correspondence, `Q` primitively
represents `n` exactly when one of these products has class `[Q]`.
Therefore, writing `g_i=[P_i]^(e_i)`, the decision is precisely

$$\prod_{i=1}^r g_i^{s_i}=[Q],\qquad s_i\in\{+1,-1\}.\tag{3}$$

Prime powers remain essential. For example, `(2,-1,3)` of
discriminant `-23` primitively represents `4` at `(1,1)`, but does
not represent `1`. Replacing a norm by its squarefree part would
lose valid inputs in this application.

## 2. Match the two halves without a group decomposition

Split the primes into two sets of sizes `floor(r/2)` and
`ceil(r/2)`. Enumerate every signed product in the first half,
reducing after each multiplication and recording its sign vector.
Sort the resulting canonical forms lexicographically.

For each signed product `v` in the second half, compute the
canonical form of `[Q]v^-1` and look for it in the first table by
binary search. A match gives all signs in (3). If every second-half
product fails, there is no primitive representation of this norm.

Building each product tree requires a number of group operations
linear in its number of leaves. Sorting and lookup add only a
polynomial factor in `r` and the coefficient bit lengths. Prime
powers are computed by binary powering. The
[two-dimensional HNF operations](direct-discriminant-short-norms.md)
give polynomial-bit class-group composition and reduction, proving
the cost in (1). Deterministic sorting and binary search avoid an
unproved worst-case bound for hash-table lookup.

The method never expands the full class group. Its tables contain
only the possible orientations of the supplied norm. Distinct
orientation choices may have the same class; retaining either
matching sign vector is sufficient.

## 3. Recover the actual coordinates

The successful signs give an explicit integral ideal product in
(2). Use the existing HNF word-to-norm constructor to multiply its
prime-power ideals and retain the proper changes of variables during
Gauss reduction. The resulting primitive norm is exactly `n`:
there are no opposite orientations at a common rational prime, so
no rational-integer content is removed.

The retained determinant-one matrix gives primitive coordinates for
the reduced form of `Q`. If the input form was not reduced, apply
its retained reduction matrix to these coordinates. All integers
have polynomial bit length in the written form and `log n`.
For this optional unreduced-input version, the complexity bound also
includes the bit lengths of the supplied coefficients.

Thus a successful group match produces a verifiable arithmetic
witness, rather than just a class-membership answer. The output can
also be supplied directly as a signed prime word.

## 4. Cost of supplying the norm and prime forms

Equation (1) assumes the prime factors and one root orientation at
each split prime are supplied. It does not give a polynomial-time
factorization algorithm for an arbitrary integer `n`.

In the arithmetic relation program, the relevant primes all lie
below the already charged cutoff `Z`. Trial division through that
cutoff can factor a candidate `n`, reject a rough candidate, and
check its total multiplicity. The factor-base construction supplies
the prime forms. Both tasks fit `Z^O(1) poly(b,log n)` deterministic
time, including elementary prime enumeration and root search.

An absent factor-base prime either fails the split or coprimality
condition and can be rejected before orientation matching. The
implementation accepts the certified split-prime root records as
input; it does not silently assume a fast deterministic square-root
routine for arbitrary large primes.

For the refined multiplicity budget

$$\Omega(n)\le L=O\left(\sqrt{b/\log b}\right),\qquad
 Z=2^{O(\sqrt{b\log b})},$$

we have `r=omega(n)<=L`, so orientation matching costs only

$$2^{O(\sqrt{b/\log b})}\mathrm{poly}(b).$$

Moreover `log n<=L log Z=O(b)`. Factoring over the bounded base and
the orientation search therefore remain within
`2^O(sqrt(b log b))` time. A future deterministic norm-selection
procedure may output just such an integer `n`; it need not also
find `x,y`, a prime orientation vector, or a balanced smooth pair.

What remains open is selecting a norm that passes this complete
test. Enumerating every smooth integer of the allowed size is not
charged as a subexponential search by this result.

## 5. Exact implementation and controls

The implementation is
[factored_norm_representation.py](../tools/factored_norm_representation.py).
Run:

```sh
python3 tools/check_factored_norm_representation.py
```

The [checker](../tools/check_factored_norm_representation.py) compares
the meet-in-the-middle answer with a separate complete integer
search. For every permitted integer `y`, the identity

$$4an=(2ax+by)^2+|\Delta|y^2$$

reduces the independent test to an integer square test and a
divisibility check. This includes both positive and negative
representation decisions, repeated prime powers, split primes above
`2`, nonfundamental discriminants, and transformed unreduced forms.

A larger known-word case checks many distinct primes without
enumerating the large representation ellipse. Its explicit
coordinates are still verified directly. The
[evidence record](../evidence/factored-norm-representation.json)
separates this witness check from the exhaustive small controls.
It records 6,092 complete comparisons over 21 discriminants: 1,478
primitive representations and 4,614 negative decisions. Of these,
3,786 cases retain a prime power with exponent greater than one.
Another 350 controls use unreduced equivalent forms. The larger
case has an 18-prime, 214-bit norm and uses 512 left-table entries
and 171 second-half probes.

Substantial OpenAI assistant support was used. No Lean formalization
or stronger complete factoring score is claimed.
