# A finite modular higher-spin crossing screen

## Result and scope

There are twenty-two explicit polynomial swap crossings in this
checkpoint: characteristic `p=3,5,7,11`, alphabet size `p`, and every
nonzero spectral scale in the corresponding prime field. They come
from a rational Casimir construction, followed by exact coefficient
reduction modulo `p`.

For each crossing, use the two natural gauge parameters
`a,b` in `{-2,-1,0,1,2}` defined below. An exact Bernstein certificate
proves that the resulting common-type pressure parameter is at least
`3/4`, at every real evaluation scale `s>1`. This gives **550 all-scale
exclusions**. It rules out an improvement past `1/sqrt(2)` from this
finite family with those gauges.

The result concerns these pressure certificates. It is not a lower
bound for factoring, for actual coefficient ranks, for other spins or
characteristics, or for arbitrary physical gauges and perturbations.

## 1. Exact polynomial construction

The family was motivated by rational higher-spin `sl_2` R-matrices.
A primary reference is Bykov and Zinn-Justin,
[*Higher spin sl2 R-matrix from equivariant (co)homology*](https://arxiv.org/abs/1904.11107).
The finite construction below is specified completely by rational
matrix arithmetic. Its crossing endpoint is checked directly; no
Yang--Baxter identity is needed for the pressure calculation.

Let `V_n` have basis `e_0,...,e_n`, with

$$E e_i=i e_{i-1},\qquad F e_i=(n-i)e_{i+1},\qquad
  H e_i=(n-2i)e_i.$$

On `V_n` tensor `V_n`, use the total operators and define

$$C=H^2/4+H/2+FE.$$

The subspace with `i+j=k` is invariant. In the basis `(i,j)` of that
subspace, the diagonal entry of `C` is

$$ (n-k)(n-k+1)+i(n-i+1)+j(n-j+1), $$

and its two possible off-diagonal entries send `(i,j)` to
`(i-1,j+1)` with coefficient `i(n-j)`, and to `(i+1,j-1)` with
coefficient `j(n-i)`. Thus every block is given explicitly by integers.

Over the rationals, this block has the distinct Casimir values
`J(J+1)` for `J=|n-k|,...,n`. Set

$$P_J=\prod_{H\ne J}
 \frac{C-H(H+1)I}{J(J+1)-H(H+1)},$$

where the product ranges over the spins present in the block. Define

$$\lambda_J(u)=\frac1{n!}
 \prod_{h=1}^{J}(u+h)\prod_{h=J+1}^{n}(u-h),\qquad
 R_n(u)=\sum_J\lambda_J(u)P_J. \tag{1}$$

The spin-`J` summand has swap eigenvalue `(-1)^(n-J)`. Since
`lambda_J(0)=(-1)^(n-J)`, equation (1) gives

$$R_n(0)=\mathrm{swap},\qquad
 [u^n]R_n(u)=I/n!. \tag{2}$$

The checker also verifies both equalities coefficient by coefficient
over the rationals.

Now take `n=p-1`, for one of the four listed primes. In these finite
cases, the least common multiple of all coefficient denominators is
exactly `n!`, which is coprime to `p`. Every combined coefficient of
`R_n` can therefore be reduced modulo `p`. For every `c` in the
nonzero elements of `F_p`, put

$$T_c(z)=R_n(c(z-1))\pmod p. \tag{3}$$

This is an explicitly listed polynomial tensor of degree at most
`p-1`, and

$$T_c(1;x,y,u,v)=[x=v][y=u]. \tag{4}$$

The modular reduction must be performed **after** the rational
projectors in (1) have been summed. The separate projector
denominators need not be invertible modulo `p`: the Casimir values
for `J` and `p-1-J` coincide in that field. Reducing those projectors
separately would not define the construction used here. This checkpoint
does not assert a denominator theorem for untested spins.

## 2. The two actual cuts

The zero cut is `(x,y)|(u,v)`. At infinity, the cut is the
realignment `(x,u)|(y,v)`. The latter is not a partial transpose.

The tensor preserves total occupation `x+y=u+v`. Consequently, the
zero matrix splits into blocks of fixed `x+y`, while the realigned
matrix splits into blocks of fixed `x-u=v-y`. The checker computes
the complete Smith orders of these small blocks by exact truncated
power-series elimination over `F_p`.

As an independent indexing control, it also constructs both full
literal cut matrices for every scale `c` in characteristics three
and five. All sixty gauged full profiles agree with the block
calculations. Equation (4) is separately checked on all 163,478
physical entries across the twenty-two seeds.

The natural zero gauges are

$$ (h_x^0,h_y^0,h_u^0,h_v^0)=(ax,ay,0,0), $$

and the natural infinity gauges are

$$ (h_x^\infty,h_y^\infty,h_u^\infty,h_v^\infty)
   =(0,0,bu,bv). \tag{5}$$

Let their complete Laurent Smith order lists be `nu^0` and
`nu^infinity`. The implementation restores every uniform monomial
shift introduced while computing these lists. In particular, the
infinity profile already includes the negative degree charge.

The infinity orders are independent of the nonzero spectral scale
`c`. Indeed,

$$T_c(z)=T_1(cz+1-c),$$

and this is an invertible change of local parameter at infinity.
Under that change, each monomial gauge acquires only a regular unit
on its row or column. The checker additionally compares the ungauged
infinity profiles for every `c`.

For example, `p=3`, `c=1`, and `a=b=0` give

$$\nu^0=(0,1,1,1,1,2,2,2,2),$$

$$\nu^\infty=(-2,-1,-1,-1,0,0,0,0,2).$$

The resulting rank-cost product is

$$B(s)=(1+4s^{-1}+4s^{-2})(s^2+3s+4+s^{-2}).$$

## 3. Exact all-scale exclusion

For any of the listed gauge pairs, write

$$B_{a,b,c}(s)=
 \left(\sum_i s^{-\nu_i^0}\right)
 \left(\sum_i s^{-\nu_i^\infty}\right).$$

For a common physical type, the total energy of symbol `j` is

$$G_j=2(a+b)j,\qquad 0\le j<p.$$

The common-type criterion at `rho=3/4` would require

$$B_{a,b,c}(s)<
 \left(\sum_{j=0}^{p-1}s^{-G_j/3}\right)^3. \tag{6}$$

For all 550 cases, the checker proves the strict reverse inequality
for every `s>1`. Put `x=s^(-1/3)`. The difference between the left
and right sides of (6) is an integer Laurent polynomial in `x`.
After removing its smallest monomial exponent, it becomes an ordinary
polynomial `P(x)` with a nonzero constant term.

Write `P(x)=sum_j c_j x^j`. Its degree-`N` Bernstein coefficients
have the same signs as the exact integers

$$d_k=\sum_{j\le k}c_j\binom{N-j}{k-j},
  \qquad 0\le k\le N,$$

because `d_k=binom(N,k)b_k`. Every coefficient `d_k` is strictly
positive in each recorded certificate. Thus `P(x)>0` on `[0,1]`,
which proves the required inequality on `s>1`. No sampled real scales
or floating optimization enter this conclusion. The largest polynomial
and Bernstein degree in the replay is 300.

The function

$$R\log\sum_j \exp(-G_j\log(s)/R)$$

strictly increases with `R>0`: its derivative is the entropy of the
corresponding positive Gibbs distribution. Therefore excluding
`R=3`, or `rho=3/4`, also excludes every smaller `rho` under the same
criterion. In particular, this family cannot cross `1/sqrt(2)` with
the tested gauges.

## 4. Replay and remaining scope

Run:

```sh
python3 tools/check_modular_spin_crossings.py \
  --output evidence/modular-spin-crossings.json
```

The [checker](../tools/check_modular_spin_crossings.py) uses exact
rational arithmetic before reduction and exact finite-field arithmetic
afterward. The [evidence](../evidence/modular-spin-crossings.json)
records all block profiles and all 550 positivity certificates.

Other physical gauges, larger spins, different modular constructions,
auxiliary-state extensions, and direct type-restricted rank calculations
remain outside this finite exclusion. The useful conclusion for the
current search is that the tested modular degenerations and their
natural gauge box do not provide a cheaper pressure witness.
