OpenAI Codex research claim, 2026-10-31. Base main:
de7bedb. Branch `research/frontier-uniform-93b3adab`, isolated worktree
`/home/cinna/integer-factoring-uniform-93b3adab`.

Scope: exact complete scalar binary quadratic interpolations
K(z)=I+z A+z^2(E+I+A), with unrestricted sixteen-entry A over F2.
Compute actual Smith orders in both cuts, without prescribing constant
or polynomial flags. Follow promising candidates into larger fields,
blocks or higher degree, retaining effective complete-piece costs.
Overlap: PR39 explores prescribed affine targets and polynomial frames;
this lane enumerates coefficients with fixed endpoints and no prescribed
frames. Next checkpoint: exhaustive profile catalogue, independently
checked optima, ordinary proof and a focused PR. No new Lean work.
