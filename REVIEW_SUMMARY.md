# Review summary (internal, AI-assisted)

This is a public summary of the internal checking of version 0.1.0-candidate, without scores. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (24-25 November 2026)

1. **A five-perspective Claude review panel** (editor, methodology, domain, perspective, devil's advocate) on the
   source draft. The methodology reviewer audited every original proof. No mathematical error was found; the panel
   asked for revisions of attribution, positioning and presentation, and suggested extensions.
2. **Clean replays** of both check programs from a fresh copy under `python -O`, with outputs identical to the shipped
   results, and a negative control.
3. **One revision round with extensions**: a context-first introduction; credit for concurrent and earlier work; the
   additive-noise growth hypothesis weakened to zero exponential type, with a nonlinear example for every positive type
   (sharpness); a several-variable Fock corollary; exact probit and logistic gate classification and a kernel-width
   rule; the exact triangular-noise threshold a_*.
4. **Verification**: a re-review of everything that changed and a proof audit of all new and rewritten proofs. No
   major or critical issue; the minor items were applied before release.

No Astra/Codex second-model review took place in this batch.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| Concurrent preprints (Cho-Choe-Koo-Tien; Lin-Lu-Yang-Zu) prove one-variable Fock no-cancellation | Both read and credited in the introduction and Section 5.3; the Gaussian-RKHS statement is positioned as a counterpart; the several-variable case with a common invertible linear part added (Corollary 5.3) |
| The noise is formal: T factors through deterministic composition | Remark 2.2; the introduction states that the new content is cross-space composition rigidity |
| Method credit to Ikeda-Ishikawa-Sawano | Credited before Lemma 2.3 and in the introduction |
| Missing citations and metadata | Tien-Khoi 2020, Hai-Huy-Tien, Sang-Zhao, Ishikawa added; bibliography copied from Crossref, DataCite and arXiv records |
| Presentation and process wording | Introduction rebuilt (problem, why it matters, what was known, what is new, what is not new, scope, roadmap, notation); section openers; author line "Anonymous"; AI-use statement |
| Triangular-noise threshold not sharp | Exact threshold a_* ≈ 4.54349, checked to 50 digits |
| Relevance for kernel EDMD | Probit/logistic gate classification and the kernel-width rule (Corollary 6.2) |
| Continuous-time remark overstated | Restricted to a constant jump rate |
| Branches with different linear parts | Attempted; branch features can coincide at complex states, so the uniform Gram bound fails; stated as open in Section 7 |

## What remains open

Mixtures with different linear parts, state-dependent covariances, infinitely many state-dependent components,
state-dependent Lévy kernels and general Markov transitions; sharpness of the zero-type condition for Gaussian sources.
