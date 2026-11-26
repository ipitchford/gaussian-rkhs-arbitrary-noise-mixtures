# Replay receipt (v0.1.0-candidate, 26 November 2026)

Release replay by the publisher from a fresh copy of the approved candidate (code, requirements.txt), macOS arm64,
working directory = package root, one CPU, a few seconds per program. Interpreters as printed inside the run (see
[ENVIRONMENT.txt](ENVIRONMENT.txt)): CPython 3.12.11 with the pinned packages.

| Program | Flags | Exit | Output |
|---|---|---|---|
| code/verify.py | -I -B -O | 0 | results/verification.json, byte-identical; status PASS, 97 of 97 checks |
| code/check_extensions.py | -I -B -O | 0 | results/extensions.json, byte-identical; `{"check_count": 63, "failed": [], "status": "PASS"}` |
| code/make_figures.py | -I -B | 0 | figures/phase-diagram.png byte-identical; phase-diagram.pdf differs only by its embedded creation data |
| Negative control | check_extensions.py with the probit threshold 4D changed to 2D (`mp.sqrt(D + D)`), -I -B -O | 1 | six probit checks fail (threshold formula and above-threshold growth at q = 0.25, 1, 3) |

stderr was empty for all passing runs. The outputs match the lead's replay of 25 November 2026. After the replay only
the title-page DOI line of manuscript.tex, the rebuilt PDF and the release records changed; code, results, figures and
requirements.txt are byte-identical to the approved candidate. GitHub CI on Linux reruns all three programs under the
comparison policy in ENVIRONMENT.txt, with the negative control.
