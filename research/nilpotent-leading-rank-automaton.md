# Exact leading auxiliary ranks from a finite fusion automaton

The [complete occupancy resolution](fusion-occupancy-resolution.md) also
diagonalizes the physical part of the four-state nilpotent central seam.
Its leading coefficient has an exact rank, computable by a closed
22-state weighted automaton. This sharpens the earlier counting bound
and identifies its exponential scale. It does not prove a stronger
complete factoring bound or equality with the final contracted row rank.

## 1. Replace physical projectors by two matrix letters

Use the F4 auxiliary matrices Q,N from the
[nilpotent seam](nilpotent-seam-crossings.md), with delta=2 in F4.
The central gate, restricted to a vacuum pair or a triplet pair, is
respectively

```math
A_0(x)=I+x\bigl(Q-I+\delta(I+N)\bigr)+x^2\delta N,
\qquad A_2(x)=I+x(Q-I).                              \tag{1}
```

For an occupancy word s in {0,2}^n, let A_s be the ordered product
of these matrices in the central application order. Orthogonal
physical projectors give the exact identity

```math
\operatorname{rank}[x^a]A_{\rm central}(x)
 =\sum_{s\in\{0,2\}^n}d_{\#2(s)}
              \operatorname{rank}[x^a]A_s(x),         \tag{2}
```

where d_0=1 and d_k=Fibonacci(k-1) for k>=1. Rank here is on
the four-dimensional auxiliary space tensor the full 2n-strand
vacuum space. Extend to F16 to use the orthogonal fusion projectors;
matrix rank is unchanged by that field extension.

Equation (2) is useful beyond the leading degree: it replaces one
exponentially large matrix rank by weighted ranks of four-by-four
matrices. The leading degree admits a further finite-state reduction.

## 2. Four letters have degree at most seven

Only A0 has a quadratic term, delta N, and N^4=0. Thus every
four-letter product has degree at most seven. Put

```math
K_w=[x^7]A_{w_4}(x)A_{w_3}(x)A_{w_2}(x)A_{w_1}(x).
```

Eleven of the sixteen K_w vanish. The five remaining words have
zero or one triplet. Up to nonzero scalars their matrices are
the five generators recorded in the
[closure certificate](../evidence/nilpotent-leading-rank-automaton.json).
The all-vacuum generator is I+N^2 and has rank four; the four
one-triplet generators have ranks1,2,2,1.

At n=4m, the maximal allowed degree is7m. No four-letter block
can exceed seven, so the leading coefficient is exactly the ordered
product of the K_w. There are no omitted lower-degree contributions.
At n=4m+2, the last two letters must both be A0 to reach degree7m+4;
their leading product is delta^2 N^2. Thus the same block product
is followed on the left by N^2, up to a nonzero scalar.

Normalize a nonzero matrix by its first nonzero entry. Starting
from I, multiplication by the five generators reaches exactly22
projective matrices, including I, and is closed, allowing the zero
product as a discarded sink. This is checked on every transition,
not inferred from stabilization at a finite word length. Scalars
can be discarded because future zero tests and ranks are unchanged.

## 3. Weight the closed automaton by retained fusion

At each projective state keep two integer counters (v,c), counting
vacuum and charged multiplicities. A zero-triplet generator leaves
them unchanged. A one-triplet generator sends

```math
(v,c)\longmapsto(c,v+c),
```

the retained object2 fusion rule. Identical destination states add
their counters. Initialize (1,0) at I. The output is the sum of
v times the matrix rank at each state, or the rank after left
multiplication by N^2 for the two-letter tail.

Let R_m and S_m denote these two outputs. Exact generating functions are

```math
\sum_{m\ge0}R_m z^m
 =\frac{4-24z+55z^2-52z^3+22z^4}
 {(1-z)(1-3z+z^2)^2},\qquad
\sum_{m\ge0}S_m z^m=\frac{2-4z}{1-3z+z^2}.          \tag{3}
```

The first ranks R_m are4,4,15,53,178,576,1806. Also S_m=2v_m,
where v_m=Fibonacci(2m-1) for m>=1 and v_0=1. Consequently the
exact leading central auxiliary rank is

```math
\boxed{\operatorname{rank}C_{n,a_*}=
 \begin{cases}R_m,&n=4m,\\2v_m,&n=4m+2,\end{cases}
 \qquad a_*=2n-\lfloor n/4\rfloor.}                  \tag{4}
```

For completeness, (3) is certified as an all-length transfer identity.
The automaton has44 integer counters, transfer matrix T, initial
column u and output row o. For the first denominator, the checker
forms (T^5-7T^4+17T^3-17T^2+7T-I)u and verifies that oT^j
annihilates it for j=0,...,43. For the tail output it uses
(T^2-3T+I)u. Cayley--Hamilton then proves every later output relation.
The first coefficients give the numerators. This avoids interpreting
a recurrence fitted to initial values as a proof.

The double pole at z=phi^-2 has a nonzero numerator. Hence
R_m=Theta(m phi^(2m)), while S_m=Theta(phi^(2m)). With boundary
label count L=Theta(phi^n), both have exponent1/2 in L. The leading
central auxiliary bound therefore has this exact exponential scale.

## 4. What transfers to a whole crossing coefficient

The two outer boundary-degree halves contribute at most
h_(n,b)=min(b+1,n-b+1) convolutions. Equation (4) gives the refined
valid upper bound

```math
r(G_{a_*,b})\le h_{n,b}\operatorname{rank}C_{n,a_*}.  \tag{5}
```

This improves the earlier 4hW bound at its top degree. It does not
give a lower bound on the contracted row rank: at six sites all
leading boundary coefficients vanish although the central auxiliary
rank is two. Other coefficients and a complete favorable mixed
crossing sum remain open. The complete factoring score stays1/2.

The [rank formula](../tools/nilpotent_seam_top_rank.py) uses integer
recurrences. The [checker](../tools/check_nilpotent_seam_top_rank.py)
reconstructs every local generator, closes the automaton, checks the
44-counter identities, compares literal occupancy words, and checks
the refined bound against all complete small leading coefficients.

```sh
make verify-braid-leading-automaton
```

The counters have O(n) bits; transfer or recurrence evaluation takes
O(n) integer operations and O(n^2) elementary bit operations. No new
Lean formalization is used. The next constructive question is how to
couple such exact occupancy channels to the boundary ranks across
the full coefficient family, rather than improve only an edge.
