# AI index: Additive-noise rigidity and the nonlinear-mixture boundary for Gaussian kernel spaces (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex) with
[varying-covariances.tex](varying-covariances.tex); Section 1 states the problem, what was known, Theorems A-C, what is
not new, scope and a roadmap; Section 7 lists the computational checks and open cases).
DOI 10.5281/zenodo.24550034 (concept 10.5281/zenodo.24550033). Follow-up to Evidence Press release
doi:10.5281/zenodo.24089772 (gaussian-rkhs-stochastic-invariance), answering its open projects (i) and (ii).
Status: unrefereed candidate; internal AI-assisted review, proof audits and computation only.

## Notation

k_P(x, y) = exp[−(x−y)ᵀP⁻¹(x−y)], Gaussian kernel with width matrix P ≻ 0; H_P(R^m) its RKHS on all of R^m (global
norms). T f(x) = E[f(X') | X = x]. Additive noise: X' = F(x) + Z, φ = characteristic function of Z. Mixture:
Σ_j p_j(x) N(Mx + b_j, Q_j). D_j = R⁻¹ − Mᵀ(P + 4Q_j)⁻¹M. Scalar example: P = R = 1, M = 1, branch variance q.
ev_z = evaluation at a complex point z = x + iy.

## Claims and where they are proved

| ID | Claim | Location | Status / dependence |
|---|---|---|---|
| A | For every law of Z, T: H_P(R^m) → H_R(R^d) is bounded iff F(x) = Ax + b, A of full column rank, plus a spectral domination condition; exact operator norm (eq. norm) | Theorem 2.1 (Theorem A); Remark 2.2 (reduction to composition); Lemmas 2.3, 2.4; Corollary 2.5 (essential norm) | written proof; noise is formal at operator level (Remark 2.2) |
| B | Stationary source with spectral measure dominating Lebesgue measure near zero, target a Hilbert space of entire functions with log‖ev_{x+iy}‖ of zero exponential type in y ⇒ F affine; for every type c > 0 a bandlimited source and stationary target admit a bounded composition with x + c'·sin x | Theorem 3.1, Proposition 3.2, Corollary 3.3 (Theorem B) | written proof; sharpness for Gaussian sources open |
| N | Noise with φ(ξ) = (1−|ξ|)_+: affine slopes allowed; ‖T_a‖² = 1/|a| iff |a| ≤ a_* ≈ 4.54349 (c_* ≈ 4.91081 at t ≈ 0.71533); slopes 1-4 give 1, 1/2, 1/3, 1/4 | Section 4 | written proof + 50-digit check |
| C | Finite mixture with common M, distinct (Q_j, b_j): bounded iff M full column rank and each p_j extends to an entire function with |p_j(x+iy)| ≤ C e^{2yᵀD_j y} | Theorem 5.1 (common covariance), Lemma 5.2 (multiplier criterion, classical form), Lemma 5.4, Theorem 5.5 (different covariances) (Theorem C) | written proof |
| F | Several-variable Fock no-cancellation: Σ_j W_{ψ_j, Az+b_j} on F²(C^d), A invertible, b_j distinct, bounded iff each summand bounded | Corollary 5.3 | written proof; one-variable case concurrent with Cho-Choe-Koo-Tien and Lin-Lu-Yang-Zu |
| G1 | Scalar bump gate e^{−αx²}: admissible iff α ≤ 8q/(1+4q), boundary included; equivalently q ≥ q_c = α/[4(2−α)]; q = 1/12, 1/4, 3/4 at α = 1/2, 1, 3/2 | Corollary 6.1, Figure 1 | written proof + exact checks |
| G2 | Probit gate Φ(βx): admissible iff β² ≤ 16q/(1+4q) (= 4D); logistic gate never (poles at iπ); kernel-width rule ℓ² ≥ σ²(σ²+4q)/(8q) | Corollary 6.2, Section 6 | written proof + 50-digit checks |
| S | Singular values of the scalar smoothing-and-gating operator | Proposition A.1 | classical calculation, checked |

## Evidence

| Check | Program | Expected output | Arithmetic |
|---|---|---|---|
| Spectral constants, tube norm, thresholds (G1), Gram identities (5.7) and unequal-covariance Gram against quadrature, triangular-noise norms (N), Bessel/warping identity (B), smoothing formulas (S) | [code/verify.py](code/verify.py) | [results/verification.json](results/verification.json): 97 checks, status PASS | sympy exact; numpy/scipy quadrature with per-check tolerances |
| c_*, t_*, a_* and slope norms (N); warping and dilated growth bounds (B); probit bound C₀ and threshold, divergence above it, logistic pole (G2); Fock Gram unitary conjugation in two variables (F); threshold identities (G1, G2) | [code/check_extensions.py](code/check_extensions.py) | [results/extensions.json](results/extensions.json): `{"check_count": 63, "failed": [], "status": "PASS"}` | mpmath 50 digits; one numpy Fock test (tolerance 1e-10, observed 3.9e-15) |
| Figure 1 | [code/make_figures.py](code/make_figures.py) | [figures/phase-diagram.png](figures/phase-diagram.png) | numpy, matplotlib |

Negative control: changing the probit threshold `mp.sqrt(4 * D)` to `mp.sqrt(2 * D)` in check_extensions.py makes six
checks fail and the program exit 1 under python -O ([REPLAY_RECEIPT.md](REPLAY_RECEIPT.md); also run in CI). The
finite calculations test constants and implementations; they do not certify the analytic theorems.

Numbers quoted on the release page: 8q/(1+4q); 16q/(1+4q); q = 1/4 at α = 1 (1/12 and 3/4 at α = 1/2 and 3/2);
a_* ≈ 4.54; 97 and 63 checks.

## Limits

The noise in A is formal at the operator level; the norm formula, multiplier criterion and smoothing spectrum are
classical. One-variable Fock no-cancellation is due concurrently and more generally to Cho, Choe, Koo and Tien
(arXiv:2611.26206) and Lin, Lu, Yang and Zu (arXiv:2611.19948). Results concern global RKHS norms only, not finite-data
prediction, other spaces or weaker norms. Open: branches with different linear parts, state-dependent covariances,
infinitely many state-dependent components, state-dependent Lévy kernels, general Markov transitions; sharpness of B
for Gaussian sources.

## Replay and records

Replay commands: [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art
and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Licences:
[LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes: [MANIFEST.sha256](MANIFEST.sha256).
