# A Gaussian pressure floor at every exterior charge

**Theorem.** For every integer digit count \(q\ge2\), interval
length \(n\ge1\), scale \(t>0\), and real charge
\(|J|\ge2n\), the complete Gaussian critical pressure satisfies

\[
G_{n,q}(J,t)>0.
\tag{1}
\]

The pressure uses \(\rho=1/\sqrt2\) and the full normalization
in the [charged Gaussian note](charged-gaussian-pressure-floor.md).
The theorem covers continuous exterior charges in every digit count.
It does not assume zero-charge minimization or assert a result for
every real charge inside the interval.

## Proof

Put \(R=\sqrt2\) and
\(f_q(x)=\log(1+e^{-x}+\cdots+e^{-(q-1)x})\).
This function is positive, strictly decreasing, and convex on the
positive axis. The pressure is even in \(J\), so write
\(J=2n+a\) with \(a\ge0\).

The reciprocal identity
\(f_q(-x)=f_q(x)+(q-1)x\) cancels the complete ground contribution.
The zero thermal energies concatenate into
\(a+1,a+3,\ldots,a+4n-1\). The infinity thermal energies are
two copies of \(1,3,\ldots,2n-1\). The physical thermal energies
are two copies of \(a+2,a+6,\ldots,a+4n-2\), at scale \(t/R\).
Consequently

\[
\begin{aligned}
G_{n,q}(2n+a,t)=\sum_{i=0}^{n-1}\bigl[&
2f_q((2i+1)t)+f_q((a+4i+1)t)+f_q((a+4i+3)t)\\
&-2R f_q((a+4i+2)t/R)\bigr].
\end{aligned}
\tag{2}
\]

For the summand with index \(i\), put
\(x=(a+4i+2)t\) and \(u=x/2>0\). Convexity gives
\(f_q(x-t)+f_q(x+t)\ge2f_q(x)\), while
\((2i+1)t\le x/2\) gives \(f_q((2i+1)t)\ge f_q(u)\).
Thus that summand is at least

\[
2\bigl[f_q(u)+f_q(2u)-R f_q(Ru)\bigr].
\tag{3}
\]

The existing [all-digit two-mode floor](truncated-symmetric-interval-frontier.md)
proves, for every \(q\ge2\) and \(u>0\),

\[
f_q(u)+f_q(3u)-R f_q(Ru)>0.
\tag{4}
\]

Since \(f_q(2u)>f_q(3u)\), expression (3) is strictly positive.
Every summand in (2) is therefore positive, proving (1).

## Exact controls

The [checker](../tools/check_exterior_charged_gaussian_floor.py) verifies
the exact signed ground and thermal multisets in (2), and the rational
midpoint/ordering relations used in (3). Its
[record](../evidence/exterior-charged-gaussian-floor.json) also compares
(2) with the independent complete-pressure evaluator for several
digit counts, lengths, scales, and positive and negative exterior
charges. Those floating comparisons verify the identity; the strict
sign follows from the analytic proof above and the previously proved
two-mode theorem.

```bash
python3 tools/check_exterior_charged_gaussian_floor.py --output /tmp/exterior-gaussian-replay.json
```
