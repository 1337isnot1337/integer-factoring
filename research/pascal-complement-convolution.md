# Convolution bounds for arbitrary Pascal complements

## Result

The greedy Pascal complement of an arbitrary mode set satisfies a family
of sharp order inequalities. These inequalities use only addition and
ordering of the original modes. They hold over every field, including
characteristics other than two.

Let \(S=\{a_1<\cdots<a_n\}\) be a nonempty finite set of nonnegative
integers. Form the Pascal columns

~~~math
u_c=\bigl(\tbinom c{a_i}\bigr)_{i=1}^n
~~~

over a field \(K\), and let \(C=\{c_1<\cdots<c_n\}\) be the greedy
independent columns with \(c\notin S\). Write
\(\nu_1<\cdots<\nu_{2n}\) for the sorted union \(S\cup C\).

**Theorem.** The complement exists, and

~~~math
\boxed{\qquad
\nu_i\le a_i\quad(1\le i\le n),\qquad
\nu_{i+j}\le a_i+a_j+1\quad(1\le i,j\le n).
\qquad}                                               \tag{1}
~~~

In particular,

~~~math
\nu_{2i-1}\le2a_i,\qquad \nu_{2i}\le2a_i+1,\qquad
c_n\le2\max S+1.                                     \tag{2}
~~~

All the convolution bounds in (1) are equalities for the initial interval
\(S=\{0,\ldots,n-1\}\). In characteristic two, the
[Pascal complement theorem](schur-pascal-complement.md) identifies these
union orders with the exact normalized one-species infinity Smith orders.
Thus (1) is a uniform constraint on every such mode geometry.

Section 3 gives a stronger complement envelope that also retains the
known locations of \(S\). It yields a rigorous lower bound for the
common-height pressure. A general entropy floor for all mode sets is
not proved here. Finite numerical searches are not a substitute for that
remaining inequality. The derivation and checks used substantial OpenAI
assistant support.

## 1. Every translated consecutive Pascal block is invertible

For integers \(A,B\ge0\), consider the square matrix

~~~math
V_{a,j}=\binom{B+j}{a},\qquad 0\le a,j\le A.
~~~

Vandermonde's identity gives

~~~math
\binom{B+j}{a}=\sum_{h=0}^{a}\binom B{a-h}\binom jh.
                                                               \tag{3}
~~~

This factors \(V\) into a lower triangular and an upper triangular
integer matrix, each with diagonal entries one. Hence \(\det V=1\).
The matrix remains invertible over every field.

Consequently any \(r\) selected rows among \(0,\ldots,A\) have rank
\(r\) on the columns \(B,\ldots,B+A\). Removing \(k\) of those
columns reduces their rank by at most \(k\).

Taking \(A=\max S\) and \(B=A+1\) proves existence of the greedy
complement: this entire consecutive block lies outside \(S\) and has
row rank \(n\). The greedy scan has collected \(n\) columns by
\(2A+1\). No periodicity assumption or characteristic restriction is
needed.

## 2. Counting a prefix of the union

For an integer \(x\), let

~~~math
F(x)=|S\cap[0,x]|,\qquad
r(x)=\mathrm{rank}_K\{u_c:0\le c\le x,\ c\notin S\}.
~~~

The greedy scan selects exactly \(r(x)\) columns through \(x\).
Therefore

~~~math
|(S\cup C)\cap[0,x]|=F(x)+r(x).                       \tag{4}
~~~

Fix \(i,j\), and set \(x=a_i+a_j+1\). Apply Section 1 to rows
\(0,\ldots,a_i\) and columns \(a_j+1,\ldots,a_j+a_i+1\).
The first \(i\) selected rows, indexed by \(a_1,\ldots,a_i\), have
rank \(i\) on this block. Exactly \(F(x)-j\) of its columns belong
to \(S\). Deleting them leaves outside columns of rank at least

~~~math
i-(F(x)-j)=i+j-F(x).                                  \tag{5}
~~~

Their rank in the full set of \(n\) rows cannot be smaller. Thus
\(r(x)\ge i+j-F(x)\), and (4) gives at least \(i+j\) union elements
through \(x\). This is the second inequality in (1). The first follows
from \(S\subseteq S\cup C\). Setting \(j=i\) and using distinct
integer orders proves (2).

Equivalently, for every \(x\ge0\),

~~~math
|(S\cup C)\cap[0,x]|
\ge\max_{-1\le A\le x}\bigl(F(A)+F(x-A-1)\bigr),
\qquad F(-1)=0.                                      \tag{6}
~~~

One can prove (6) directly with the consecutive block
\(x-A,\ldots,x\). It is the counting-function version of the
convolution bounds. The inequalities do not determine the actual
complement in every characteristic; extra linear dependencies can make
the union orders smaller.

## 3. An explicit largest complement allowed by these inequalities

For every pair \(i,j\), put

~~~math
b_{ij}=a_i+a_j+1,\qquad
k_{ij}=\max\{0,i+j-F(b_{ij})\}.                       \tag{7}
~~~

At least \(k_{ij}\) complement columns must occur through \(b_{ij}\).
Define the deadline for the \(k\)-th complement element by

~~~math
d_k=\min\{b_{ij}:k_{ij}\ge k\},\qquad 1\le k\le n.
                                                               \tag{8}
~~~

Every set in (8) is nonempty: the pair \((n,n)\) has \(k_{nn}=n\).
The deadlines are nondecreasing, and the actual complement satisfies
\(c_k\le d_k\).

Construct \(\widehat C=\{\widehat c_1<\cdots<\widehat c_n\}\)
backward. Choose \(\widehat c_n\) as the largest nonnegative integer
outside \(S\) at most \(d_n\). Then choose \(\widehat c_k\) as the
largest nonnegative integer outside \(S\) satisfying

~~~math
\widehat c_k\le d_k,\qquad
\widehat c_k<\widehat c_{k+1}.                         \tag{9}
~~~

This procedure cannot fail, and

~~~math
\boxed{\qquad c_k\le\widehat c_k\quad(1\le k\le n).
\qquad}                                               \tag{10}
~~~

**Proof.** Induct backward. For \(k=n\), the actual \(c_n\) is an
available choice. If \(c_{k+1}\le\widehat c_{k+1}\), then
\(c_k<c_{k+1}\le\widehat c_{k+1}\) and \(c_k\le d_k\).
Thus \(c_k\) is an available choice at step \(k\), and the largest
choice is at least \(c_k\).

The envelope is computed from \(O(n^2)\) pair deadlines and elementary
integer comparisons. It uses no field arithmetic. Avoiding the finite
set \(S\) when applying (9) requires skipping at most \(n\) excluded
integers in total, even when the mode indices have large gaps.

The envelope need not be a Pascal basis. For example, over
\(\mathbb F_2\),

~~~math
S=\{0,1,2,4\},\qquad C=\{3,5,6,7\},\qquad
\widehat C=\{3,5,7,9\}.                               \tag{11}
~~~

For an initial interval, \(\widehat C=C=\{n,\ldots,2n-1\}\).
Thus both the exact theorem and this relaxation retain the interval
profile without loss.

## 4. A pressure bound with no polynomial elimination

In characteristic two, write \(U=S\cup C\) and
\(\widehat U=S\cup\widehat C\). For the common infinity height
\(h\) in the [Pascal complement note](schur-pascal-complement.md),
define

~~~math
\begin{aligned}
\mathcal B_{S,h}(s)
&=\prod_{a\in S}(1+s^{-(2a+1)})^2
  \prod_{\nu\in U}(1+s^{-(2\nu+1-2n+h)}),\\
\widehat{\mathcal B}_{S,h}(s)
&=\prod_{a\in S}(1+s^{-(2a+1)})^2
  \prod_{\nu\in\widehat U}(1+s^{-(2\nu+1-2n+h)}).
\end{aligned}                                         \tag{12}
~~~

For every real \(h\) and \(s>1\), (10) implies

~~~math
\mathcal B_{S,h}(s)\ge\widehat{\mathcal B}_{S,h}(s).
                                                               \tag{13}
~~~

The physical height product \(\mathcal Z_{S,h,R}(s)\) does not depend
on \(C\). Therefore a certified inequality
\(\widehat{\mathcal B}_{S,h}(s)\ge\mathcal Z_{S,h,R}(s)^R\)
rules out that pressure witness for the exact profile as well. A strict
inequality in the opposite direction for the envelope is only a reason
to compute the actual complement; it is not a crossing certificate.

The inequalities (1), (10), and (13) hold uniformly. This note does not
establish

~~~math
\widehat{\mathcal B}_{S,h}(s)
\ge\mathcal Z_{S,h,2\sqrt2}(s)^{2\sqrt2}
\quad\text{for every }S,h,s,                           \tag{14}
~~~

nor claim that initial intervals minimize the optimized pressure over
all mode sets. Equation (14) is a concrete remaining analytic question.
It would supply a floor for the stated common-height family. It would
still not cover independent mode heights or two occupied species.

## 5. Exact reproduction

```sh
python3 tools/check_pascal_complement_convolution.py \
  --output evidence/pascal-complement-convolution.json
```

The standard-library checker tests the integer Pascal factorization,
all nonempty subsets of a fixed small interval in characteristic two,
independent modular scans in characteristics three, five, and seven,
all interval sizes through 64, and mode indices with up to 256 bits.
It checks every inequality in (1) and (10) with integers. These checks
support the implementation; Sections 1--3 give the proofs for all sizes.
