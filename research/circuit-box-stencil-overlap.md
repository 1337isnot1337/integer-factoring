# Proper stencil roots force conjugate overlaps

The [equivalent-norm stencil](equivalent-norm-stencils.md) has an exact
correlation that holds through whole prime powers. Under a simple determinant
condition, a common residual gcd always has opposite proper root signs.
For the standard four-vector stencil, that condition is automatic for odd
residual norms. Its residual overlaps therefore form a graph, with no prime
occurring at three vertices.

This gives a gcd-only formula for the full rational integer content of the
product of the residual ideals. It does not prove that enough overlaps,
smooth norms, or independent relations occur. It is not a general
performance barrier or an improvement to the factoring exponent.

## 1. Proper roots and the exact polarization identity

Let

\[
Q(x,y)=Ax^2+Bxy+Cy^2,\qquad D=B^2-4AC<0,
\]

be a primitive positive integral form. Write

\[
\mathcal B(v,w)=Q(v+w)-Q(v)-Q(w)
\]

for its integral polarization. For a primitive integer vector \(v\), choose
an integer complement \(u\) with \(\det(v,u)=1\). The properly transformed
form is

\[
\bigl(Q(v),b_v,Q(u)\bigr),\qquad b_v=\mathcal B(v,u).
\tag{1}
\]

Changing \(u\) to \(u+jv\) changes \(b_v\) by \(2jQ(v)\). Thus the
proper root \(b_v\bmod2Q(v)\) is independent of the chosen complement.
It retains a sign; it is not automatically the smaller representative
of the two roots \(\pm b_v\).

Take two primitive vectors \(v,w\), with proper complements \(u,t\), and
put \(\delta=\det(v,w)\ne0\). There are integers \(\alpha,\beta\) with

\[
w=\alpha v+\delta u,\qquad v=\beta w-\delta t.
\]

Polarizing these two equalities gives the exact integer identity

\[
\boxed{\delta(b_v+b_w)=2\beta Q(w)-2\alpha Q(v).}
\tag{2}
\]

Let \(g\) be an odd common divisor of \(Q(v),Q(w)\), with
\(\gcd(g,\delta)=1\). The two middle coefficients have the parity of
\(D\), so \((b_v+b_w)/2\) is integral. Divide (2) by two and invert
\(\delta\) modulo \(g\) to obtain

\[
\boxed{b_v\equiv-b_w\pmod{2g}.}
\tag{3}
\]

This proves the full gcd statement, including unequal prime-power
valuations in the two norms. If also \(\gcd(g,D)=1\), the corresponding
odd unramified ideals satisfy

\[
I(g,b_v)I(g,b_w)=(g),\qquad
I(r,b)=\left(r,\frac{-b+\sqrt D}{2}\right).
\tag{4}
\]

Equation (4) is a literal integral-ideal identity. Its right side has
norm \(g^2\). No prime factorization of \(g\) is needed to compute or
check (3).

If roots are stored canonically up to sign, record
\(b_v\equiv\epsilon_v c_v\pmod{2Q(v)}\), with
\(\epsilon_v\in\{1,-1\}\). Then the correct stored-root test is

\[
\epsilon_v c_v\equiv-\epsilon_w c_w\pmod{2g}.
\tag{5}
\]

The canonical roots themselves can agree or be opposite on the common
component. Forgetting the proper signs would give an incorrect rule.

## 2. The four-vector stencil needs no odd determinant exclusions

For the four vectors used in the preceding experiment, proper complements
give the following table:

| Vector | Norm | Proper middle coefficient |
| --- | --- | --- |
| \((1,0)\) | \(A\) | \(B\) |
| \((0,1)\) | \(C\) | \(-B\) |
| \((1,1)\) | \(A+B+C\) | \(B+2C\) |
| \((1,-1)\) | \(A-B+C\) | \(B-2C\) |

Every pair determinant is \(\pm1\) or \(\pm2\). Consequently every odd
common residual gcd satisfies the determinant hypothesis in (3), without
removing any additional odd primes.

For a larger stencil, a sufficient condition is that every residual norm
be coprime to every pair determinant. In particular, if all residual
primes exceed a bound \(Y\) and all nonzero pair determinants have absolute
value at most \(Y\), the condition holds. The weaker pairwise condition
in the next section is sufficient as well. These hypotheses can be checked
with ordinary gcds; primality of a residual is not required.

The determinant condition cannot be omitted. Let \(p\equiv1\pmod4\) be
prime and \(a=p^e\). Choose an even root \(b\) of \(-4\) modulo \(4a\), and set
\(Q=(a,b,(b^2+4)/(4a))\). The vectors \(v=(1,0)\) and \(w=(1,a)\)
have determinant \(a\). Their norm gcd is \(a\), and their proper roots
agree modulo \(2a\). They are not opposite, since \(p\nmid b\).
The checker retains this counterfamily through whole prime powers.

## 3. Exact graph cancellation for residual ideals

Take pairwise nonparallel primitive vectors \(v_1,\ldots,v_t\). Let
\(R_i\mid Q(v_i)\) be odd unramified residual norms, carrying their proper
roots \(b_i\), and suppose

\[
\gcd\!\left(\gcd(R_i,R_j),\det(v_i,v_j)\right)=1
\qquad(i\ne j).
\tag{6}
\]

**No prime divides three residual norms.** If an odd prime \(p\) divided
three, (3) would make each pair of their roots opposite modulo \(p\).
The first two such equalities make the first and third roots agree; the
third equality then forces that root to be zero. This contradicts
\(b_i^2\equiv D\not\equiv0\pmod p\).

Put

\[
g_{ij}=\gcd(R_i,R_j),\qquad
G=\prod_{i<j}g_{ij},\qquad
T_i=R_i\bigg/\prod_{j\ne i}g_{ij}.
\tag{7}
\]

The edge gcds greater than one are mutually coprime: a prime in two
different edges would occur at three or more vertices. Each quotient in
(7) is therefore integral. The \(T_i\) are pairwise coprime, because a
prime occurring at two vertices is removed to its minimum exponent from
both. An unequal exponent leaves a power at only one vertex.

The resulting exact ideal identity is

\[
\boxed{\prod_i I(R_i,b_i)=(G)\prod_i I(T_i,b_i).}
\tag{8}
\]

For a prime shared by two vertices, (3) assigns opposite branches. Their
minimum common power contributes a rational principal ideal, and their
remaining power stays on one branch. A prime occurring once is unchanged.
This proves (8) prime by prime. The residual product on the right is
primitive, so \(G\) is the **largest positive rational integer content**
of the ideal product, not merely a divisor of that content.

For example, \(Q=(65,3,7)\) has discriminant \(-1811\). Its four norms
are \(65,7,75,69\). The two nontrivial edge gcds are 5 and 3, and (7)
gives

\[
G=15,\qquad(T_1,T_2,T_3,T_4)=(13,7,5,23).
\]

The surviving factor 5 records the unequal valuations in 65 and 75.
The checker verifies (8) as an equality of integer HNF lattices.

The prescribed reduction uses \(\binom t2\) pair gcds, integer products,
and exact divisions. All these costs remain charged. For a fixed small
stencil they are polynomial in the written norm and root lengths. If
known prime-base parts were removed first, (8) can be combined with those
known words. It supplies a relation only when the remaining \(T_i\) permit
one; it gives no bound on how often that happens. Every stencil form is
in the same class, so a combined relation also retains its actual
multiplier \(t\) on that class.

## 4. Independent controls and scope

Run:

```bash
python tools/check_circuit_box_stencil_overlap.py \
  --output evidence/circuit-box-stencil-overlap.json
```

The [receipt](../evidence/circuit-box-stencil-overlap.json) records 78
primitive reduced forms at 13 discriminants. Its catalogue contains 1,248
proper keys, 4,992 changes of proper complement, 9,360 exact polarization
identities, and 1,442 nontrivial unramified pair controls. All four possible
pairs of canonical signs occur. Each accepted common component is also
checked by independent literal ideal multiplication.

The four-vector stencil and a 16-vector stencil each have 78 complete
graph-identity controls. The larger stencil removes determinant-supported
primes before applying (6). Independent trial factorizations of these
small residuals check 598 prime-incidence profiles. Trial factorization is
used only for this independent finite control, not in equations (7)--(8).

There are 185 separately constructed prime-power pairs with exact prescribed
valuations, each paired with a counterexample when determinant coprimality
is removed. The largest shared component is \(5^{128}\), of 298 bits; its
second norm has exactly 257 factors of 5 and has 1,493 bits in total. The
remaining part of that norm is never factored. Source hashes bind the
checker and the independent quadratic-form and HNF implementations.

These results explain a definite source of correlation between stencil
norms. They do not turn the preceding finite negative comparison into a
general lower bound on larger stencils or other relation searches.
