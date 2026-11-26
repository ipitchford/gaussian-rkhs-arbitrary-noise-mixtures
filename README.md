# Additive-noise rigidity, a sharp growth threshold and an exact nonlinear-mixture boundary for Gaussian kernel spaces

Version 0.1.0-candidate · DOI [10.5281/zenodo.24550034](https://doi.org/10.5281/zenodo.24550034) (all versions:
[10.5281/zenodo.24550033](https://doi.org/10.5281/zenodo.24550033)) · Creator: Anonymous · Text CC0-1.0, code MIT ·
Release page: https://evidencepress.org/releases/gaussian-rkhs-arbitrary-noise-mixtures/

**Status: unrefereed research candidate.** Written proofs with computational checks. Internal AI-assisted review only
(see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)); not peer review, independent reproduction or formal verification.

## What this is about

To learn how a random system evolves, kernel methods such as kernel extended dynamic mode decomposition (kernel EDMD)
average smooth "test functions" over one random step: Tf(x) = E[f(X') | X = x]. The analysis is cleanest when this
averaging keeps every test function inside the same space of smooth functions, most often the reproducing kernel Hilbert
space (RKHS) of a Gaussian kernel of width P. If it does not, the averaging operator is not even bounded on that space,
and error bounds that assume invariance do not apply.

An earlier Evidence Press release ([doi:10.5281/zenodo.24089772](https://doi.org/10.5281/zenodo.24089772)) found the
exact condition for linear systems and showed that "map plus noise" steps X' = F(x) + Z keep a Gaussian space only if F
is affine, but only for noise whose characteristic function is not too small. It left two questions: can that condition
be dropped, and can any transition whose randomness depends on the state keep the space? This paper answers both.

## What is shown

- **(A) No noise rescues nonlinear drift.** For every noise distribution, including noise whose frequency content is
  cut off entirely, boundedness between global Gaussian RKHSs forces F to be injective and affine, and the operator norm
  has an exact formula (Theorem 2.1).
- **(B) The threshold is sharp.** Composition into any Hilbert space of entire functions whose evaluation norms grow, in
  the imaginary direction, more slowly than every exponential of an exponential forces an affine map; for every faster
  growth rate there is a target where a nonlinear map x + c·sin x works (Theorem 3.1, Proposition 3.2, Corollary 3.3).
- **(C) State-dependent mixtures can pass, exactly.** A finite mixture of Gaussian branches with a common linear part,
  fixed covariances and state-dependent mixing probabilities is bounded exactly when each probability extends to an
  entire function within an explicit Gaussian growth bound (Theorems 5.1 and 5.5). With branch variance q: a bump gate
  exp(−αx²) works exactly when α ≤ 8q/(1+4q); a probit gate Φ(βx) exactly when β² ≤ 16q/(1+4q); a logistic gate never
  (Corollaries 6.1, 6.2). Probit gates can pass; logistic gates never can.
- **(D) Fock spaces.** The same Gram-matrix argument gives a several-variable no-cancellation theorem for sums of
  weighted composition operators on Fock spaces with a common invertible linear part (Corollary 5.3).

**Not new / limits.** The noise is formal at the operator level (it reduces to deterministic composition); the norm
formula, multiplier criterion and smoothing spectrum are classical. The one-variable Fock no-cancellation statement was
proved concurrently and more generally by Cho, Choe, Koo and Tien (arXiv:2611.26206) and by Lin, Lu, Yang and Zu
(arXiv:2611.19948). All results concern global RKHS norms only, not finite-data prediction or other norms. Open:
branches with different linear parts, state-dependent covariances, Lévy-type state-dependent kernels and general Markov
transitions.

## Contents

| Path | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf) | The paper (21 pages); source [manuscript.tex](manuscript.tex), which inputs [varying-covariances.tex](varying-covariances.tex) |
| [figures/phase-diagram.pdf](figures/phase-diagram.pdf), [.png](figures/phase-diagram.png) | Figure 1: the exact admissible (α, q) region for the bump gate |
| [code/verify.py](code/verify.py) | 97 checks: spectral constants, tube norm, scalar thresholds, Gram identities against quadrature, triangular-noise norms, warping identity, smoothing formulas |
| [code/check_extensions.py](code/check_extensions.py) | 63 checks in 50-digit arithmetic (one numpy Fock Gram test): a_*, warping growth bounds, probit bound and threshold, logistic pole, Fock Gram conjugation, threshold identities |
| [code/make_figures.py](code/make_figures.py) | Draws Figure 1 |
| [results/](results) | verification.json, extensions.json (approved outputs) |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md), [RESEARCH_GATES.json](RESEARCH_GATES.json) | Internal review summary; prior-art and limitation record |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment |
| [refs.bib](refs.bib), [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Bibliography, citation, licences, hashes |

## Replay (seconds)

From the package root (the scripts write to ./results and ./figures), with [uv](https://docs.astral.sh/uv/):

```sh
uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with sympy==1.14.0 python3 -I -B -O code/verify.py --output results/verification.json
uv run --python 3.12 --with numpy==2.3.5 --with sympy==1.14.0 --with mpmath==1.3.0 python3 -I -B -O code/check_extensions.py --output results/extensions.json
uv run --python 3.12 --with numpy==2.3.5 --with matplotlib==3.10.8 python3 -I -B code/make_figures.py
```

Expected: verify.py reports status PASS with 97 checks; check_extensions.py prints
`{"check_count": 63, "failed": [], "status": "PASS"}`. The scripts use no `assert`, so checks stay active under `-O`.
Comparison policy (GitHub CI on Linux): the interpreter-version field is excluded; check names, pass flags, summaries
and non-numeric fields must match exactly; numeric fields within absolute 1e-12 or relative 1e-9 (quadrature residuals
can differ at about 1e-15 across platforms), except scipy's own quadrature error estimates (`quad_error_estimate`,
diagnostics, not pass criteria) within relative 1e-4. On the release host both JSON files were byte-identical to the shipped ones.

## Cite

Anonymous (2026). *Additive-noise rigidity, a sharp growth threshold and an exact nonlinear-mixture boundary for
Gaussian kernel spaces* (Version 0.1.0-candidate). Evidence Press. https://doi.org/10.5281/zenodo.24550034
