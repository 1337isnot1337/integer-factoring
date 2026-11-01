# An explicit integer cutoff for the prime-form construction

This note makes the effective constants in
[the prime-form generator and biased-sample theorem](polylog-prime-form-generators.md)
concrete. It uses the same pinned Cureton--Thorner zero-free theorem and
the classical degree-two Hecke functional equation. It supplies no new
factoring exponent and no independent verification of that external theorem.

Let `Delta < 0` be a quadratic-order discriminant, put
`D = abs(Delta)`, and let `b = bits(D)`. The following integer cutoff is
sufficient for full generation:

$$X=\bigl(10000(b+1)\bigr)^{32}.\tag{1}$$

Both orientations of every split odd prime below `X` and coprime to
`Delta` generate the complete proper ideal class group. This is a
polynomial asymptotic bound with an enormous constant and degree; it is
not a practical prime enumeration prescription.

For a rational `0 < epsilon <= 1/2`, put

$$k=\lceil\log_2(1/\epsilon)\rceil,\quad B=b+k+2,\quad
 y=\left\lceil\frac{10000B}{\epsilon}\right\rceil,\quad X=y^{32}.
 \tag{2}$$

The prime-form distribution constructed in the preceding note, using
dyadic logarithm precision `j` with `2^j >= 48/epsilon`, has all
nonprincipal character biases at most `epsilon`. All cutoffs and
precision parameters in (1)--(2) are computed by integer arithmetic.

## 1. A numerical form of the smoothed explicit formula

Use `theta=31/32`, and let `chi` be a primitive finite-order Hecke
character of an imaginary quadratic field, with analytic conductor
`Q >= 3`. For

$$W_\chi(x)=\sum_{\mathrm N\mathfrak a\le x}
 \Lambda_K(\mathfrak a)\chi(\mathfrak a)
 \frac{\mathrm N\mathfrak a}{x}
 \left(1-\frac{\mathrm N\mathfrak a}{x}\right),$$

we will prove the deliberately coarse uniform bound

$$\left|W_\chi(x)-\delta_\chi x/6\right|
 \le 100x^{31/32}\log(Q+3),\qquad x\ge2.\tag{3}$$

The shifted contour identity in the preceding note is

$$W_\chi(x)=\delta_\chi x/6
 -\sum_\rho\frac{x^\rho}{(\rho+1)(\rho+2)}
 -r_\chi/2+I_\chi(x).\tag{4}$$

Here `r_chi <= 2`. We bound the zero sum and final integral separately.

### The zero sum

For `sigma >= 3/2`, the degree-two Euler product gives

$$|L'/L(\sigma+it,\chi)|
 \le 2\sum_{n\ge2}\frac{\log n}{n^\sigma}.$$

The function `log(t)/t^sigma` is decreasing for `t >= 2`. Use

$$\sum_{n\ge2}\frac{\log n}{n^\sigma}
 \le\frac{\log2}{2^\sigma}
       +\int_2^\infty\frac{\log t}{t^\sigma}\,dt
 \le\frac{\log2}{2^\sigma}+\frac1{(\sigma-1)^2}.$$

The last quantity is less than two at `sigma=2` and less than five at
`sigma=3/2`. This gives

$$|L'/L(2,\chi)|<4,\qquad
 |L'/L(3/2+it,\chi)|<10.\tag{5}$$

The Hadamard logarithmic-derivative identity at `s=2` now gives

$$\sum_\rho\mathrm{Re}\frac1{2-\rho}
 \le\tfrac12\log Q+7.$$

Indeed, the two possible pole terms sum to at most `3/2`, and the
gamma logarithmic derivative is `psi(2)-log(2pi) < 1`.
For `0 < Re(rho) < 1`,

$$\frac1{1+(\mathrm{Im}\rho)^2}
 \le4\,\mathrm{Re}\frac1{2-\rho}.$$

Consequently the zero-free strip and
`abs((rho+1)(rho+2)) >= 1+(Im(rho))^2` imply

$$\left|\sum_\rho\frac{x^\rho}{(\rho+1)(\rho+2)}\right|
 \le x^{31/32}(2\log Q+28).\tag{6}$$

### The left contour

For `Re(z)>0`, the digamma integral identity
([DLMF 5.9.13](https://dlmf.nist.gov/5.9.E13)) is

$$\psi(z)=\log z-\int_0^\infty e^{-zt}
 \left(\frac1{1-e^{-t}}-\frac1t\right)dt.$$

The parenthesized function lies between zero and one. Therefore
`abs(psi(z)-log(z)) <= 1/Re(z)`. Applying this at `3/2+it`, and using
`psi(z+2)=psi(z)+1/z+1/(z+1)` at `z=-1/2+it`, gives

$$|\psi(3/2+it)|\le\log(|t|+2)+3,$$

$$|\psi(-1/2+it)|\le\log(|t|+2)+7.$$

The Hecke functional equation and (5) thus give, on `s=-1/2+it`,

$$|L'/L(s,\chi)|\le\log Q+24+2\log(|t|+2).$$

Also

$$\frac1{|(s+1)(s+2)|}\le\frac2{1+t^2}.$$

Finally,

$$\int_{-\infty}^{\infty}\frac{\log(|t|+2)}{1+t^2}\,dt
 \le\tfrac52\pi\log2,$$

because `abs(t)+2 <= 2(1+abs(t))` and
`1+abs(t) <= sqrt(2(1+t^2))`. The elementary integral
of `log(1+t^2)/(1+t^2)` is `2pi log2`. Inserting these estimates in
the vertical integral yields

$$|I_\chi(x)|\le x^{-1/2}(\log Q+30).\tag{7}$$

Combining (4), (6), and (7) gives the bound
`x^(31/32)(3 log Q+59)`. Since `Q >= 3`, this is at most the
right-hand side of (3). No constant in this derivation depends on the
quadratic field, class number, conductor, or character order.

## 2. Full generation at the explicit cutoff

The conductor inequality is `Q_chi <= D`, and
`log(D+3) <= b+1`. If the prime-form subgroup at cutoff `X` is proper,
the principal and an appropriate nonprincipal character therefore satisfy

$$|W_1(X)-W_\chi(X)|\ge X/6-200X^{31/32}(b+1).$$

For (1) this is `11X/75`, greater than `X/8`.
On the other hand, all terms away from primes dividing `2 Delta` cancel,
so the missing-prime bound is

$$|W_1(X)-W_\chi(X)|\le(b+1)\log X.$$

Writing `z=10000(b+1)`, we have `log X=32 log z <= 32z`, while
`X/8=z^32/8 > 32(b+1)z`. This is a contradiction for every `b >= 2`.
There is no unspecified small-discriminant branch.

## 3. Exact bias budget at the explicit cutoff

Let `A_chi(X)` be the sum over the two orientations of split odd primes
below `X`, excluding primes dividing `Delta`, with weight
`log(p) (p/X)(1-p/X)`.

For a rational prime at most `sqrt(X)`, the contribution of all higher
prime-ideal powers and inert prime ideals is at most `(log X)/2`.
The excluded prime ideals contribute at most `(b+1)(log X)/2`.
Thus (3) gives the uniform error

$$\left|A_\chi(X)-\delta_\chi X/6\right|
 \le 100X^{31/32}(b+1)
 +\tfrac12\sqrt X\log X+\tfrac12(b+1)\log X.$$

For `X=y^32`, use `log y <= y` to bound this error divided by `X` by
the exact rational expression

$$e=\frac{100(b+1)}y+\frac{16}{y^{15}}
       +\frac{16(b+1)}{y^{31}}.\tag{8}$$

At (2), each summand in (8) is at most `epsilon/96`. Hence
`e <= epsilon/32`, and

$$A_1(X)/X\ge1/6-e\ge29/192>1/12.$$

The unrounded prime distribution has nonprincipal bias at most

$$\frac e{1/6-e}\le\frac{6\epsilon}{29}<\frac\epsilon4.$$

The exact dyadic logarithm rounding changes the distribution in `l1`
norm by at most `24/2^j <= epsilon/2`. The final bias is therefore
strictly below `epsilon`. The constants leave unused slack deliberately.

## 4. Reproducible arithmetic checks

The companion [integer-cutoff module](../tools/prime_form_cutoffs.py)
computes (1)--(2) and the rational envelope (8), without constructing
the enormous prime list. Run:

```sh
python3 tools/check_prime_form_cutoffs.py
```

The [record](../evidence/prime-form-cutoffs.json) contains 4,096 exact
ceiling-logarithm checks, seven generator-cutoff checks, and 42 sample
budget checks, including `b=2^20` and `epsilon=2^-1024`. The checker
validates the exact budget across this parameter grid and the large-bit
controls. It checks the integer arithmetic implementing
the proof, not the analytic theorem or infinitely many L-functions.
The [prime-form sample checker](../tools/check_prime_form_bias.py)
continues to test the constructor itself at much smaller finite cutoffs.

This note used substantial OpenAI assistant support. Its analytic
dependency and source credit are the same as in the preceding note.
