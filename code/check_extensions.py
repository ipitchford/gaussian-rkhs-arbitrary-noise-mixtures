#!/usr/bin/env python3
"""Finite checks for the constants and examples added in the revised manuscript.

Covers: the triangular-noise threshold (eq. astar), the growth bounds of the
warping target and its dilations (eq. typetwo, Corollary cor:sharp), the probit
growth bound and threshold (Corollary cor:gates), the unitary-conjugation form
of the Fock Gram matrix (Corollary cor:fock), and the scalar threshold
identities of the examples section.

These are finite calculations; they do not prove the theorems.  Every check
uses an explicit conditional (no assert), so `python -O` keeps them active.
The script is deterministic (fixed seed for the Fock Gram test).

Usage: python3 code/check_extensions.py --output results/extensions.json
"""

from __future__ import annotations

import argparse
import json
import platform
import sys

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 50
CHECKS: list[dict] = []


def check(name: str, passed: bool, **details) -> None:
    CHECKS.append({"name": name, "passed": bool(passed), **{k: str(v) for k, v in details.items()}})


def triangular() -> dict:
    f = lambda t: -2 * mp.log(1 - t) / t**2
    t_star = mp.findroot(lambda t: mp.diff(f, t), mp.mpf("0.7"))
    c_star = f(t_star)
    a_star = mp.sqrt(1 + 4 * c_star)
    check("c_star_value", abs(c_star - mp.mpf("4.91081")) < 5e-6, c_star=mp.nstr(c_star, 20))
    check("t_star_value", abs(t_star - mp.mpf("0.71533")) < 5e-6, t_star=mp.nstr(t_star, 20))
    check("a_star_value", abs(a_star - mp.mpf("4.54349")) < 5e-6, a_star=mp.nstr(a_star, 20))
    # c_star is a minimum: second derivative positive, and grid minimum agrees
    grid_min = min(f(mp.mpf(k) / 2000) for k in range(1, 2000))
    check("c_star_is_grid_minimum", grid_min >= c_star - mp.mpf("1e-30") and grid_min - c_star < 1e-5,
          grid_min=mp.nstr(grid_min, 15))

    def sup_maximand(a):
        c = (a**2 - 1) / 4
        g = lambda t: (1 - t) ** 2 * mp.e ** (c * t**2)
        vals = [g(mp.mpf(k) / 4000) for k in range(0, 4001)]
        return max(vals)

    results = {}
    for a in [1, 2, 3, 4]:
        m = sup_maximand(mp.mpf(a))
        results[a] = m / a
        check(f"triangular_norm_slope_{a}", abs(m - 1) < mp.mpf("1e-40"), squared_norm=mp.nstr(m / a, 20))
    below = sup_maximand(a_star - mp.mpf("1e-6"))
    above = sup_maximand(a_star + mp.mpf("1e-3"))
    check("triangular_threshold_below", abs(below - 1) < mp.mpf("1e-30"), sup=mp.nstr(below, 20))
    check("triangular_threshold_above", above > 1 + mp.mpf("1e-4"), sup=mp.nstr(above, 20))
    return {"c_star": mp.nstr(c_star, 25), "t_star": mp.nstr(t_star, 25), "a_star": mp.nstr(a_star, 25)}


def warping_growth() -> dict:
    kappa = 2 * mp.e - 1
    rows = []
    for a in [mp.mpf("0.1"), mp.mpf("0.25"), mp.mpf("0.49")]:
        for y in [mp.mpf(v) for v in ["0.01", "0.3", "1", "2.5", "4"]]:
            sq = (mp.sinh(2 * a * y) / (2 * a * y)) * (mp.e ** (mp.e ** (2 * y)) + mp.e ** (mp.e ** (-2 * y)) - 1) / kappa
            two_lambda = mp.log(sq)
            lower = mp.e ** (2 * y) - mp.log(kappa)
            upper = mp.e ** (2 * y) + 2 * a * y
            ok = lower <= two_lambda <= upper
            rows.append((float(a), float(y), ok))
            check(f"warping_growth_a={mp.nstr(a,3)}_y={mp.nstr(y,3)}", ok,
                  two_lambda=mp.nstr(two_lambda, 15), lower=mp.nstr(lower, 15), upper=mp.nstr(upper, 15))
            # direct spectral integral for the squared evaluation norm at iy
            if y <= 1:
                direct = mp.mpf(0)
                for n in range(-90, 91):  # tail beyond |n|=90 is below 1e-60 relative for y<=1
                    direct += mp.quad(lambda e: mp.e ** (-2 * e * y), [n - a, n + a]) / (2 * a * kappa * mp.factorial(abs(n)))
                check(f"warping_eval_norm_direct_a={mp.nstr(a,3)}_y={mp.nstr(y,3)}",
                      abs(direct - sq) / sq < mp.mpf("1e-25"), direct=mp.nstr(direct, 20), closed=mp.nstr(sq, 20))
    # dilation: Lambda_c(y) = Lambda_star(c y / 2) obeys the corollary bounds
    for c in [mp.mpf("0.2"), mp.mpf("1"), mp.mpf("3")]:
        a = mp.mpf("0.25")
        for y in [mp.mpf("0.5"), mp.mpf("2"), mp.mpf("5")]:
            yy = c * y / 2
            sq = (mp.sinh(2 * a * yy) / (2 * a * yy)) * (mp.e ** (mp.e ** (2 * yy)) + mp.e ** (mp.e ** (-2 * yy)) - 1) / kappa
            lam = mp.log(sq) / 2
            ok = (mp.e ** (c * y) - mp.log(kappa)) / 2 <= lam <= (mp.e ** (c * y) + a * c * y) / 2
            check(f"dilated_growth_c={mp.nstr(c,3)}_y={mp.nstr(y,3)}", ok, Lambda=mp.nstr(lam, 15))
    return {"grid_points": len(rows)}


def Phi(w):
    return (1 + mp.erf(w / mp.sqrt(2))) / 2


def probit() -> dict:
    C0 = mp.mpf(1) / 2 + (mp.sqrt(mp.e) + 1 + mp.sqrt(mp.pi / 2)) / mp.sqrt(2 * mp.pi)
    worst = mp.mpf(0)
    for u in [mp.mpf(v) / 4 for v in range(-40, 41)]:
        for v in [mp.mpf(k) / 4 for k in range(-40, 41)]:
            r = abs(Phi(mp.mpc(u, v))) / mp.e ** (v**2 / 2)
            worst = max(worst, r)
    check("probit_growth_bound_grid", worst <= C0, worst_ratio=mp.nstr(worst, 15), C0=mp.nstr(C0, 15))
    lower_ok = all(abs(Phi(mp.mpc(0, v))) >= mp.e ** ((abs(v) - 1) ** 2 / 2) / mp.sqrt(2 * mp.pi)
                   for v in [mp.mpf(k) / 2 for k in range(2, 41)])
    check("probit_lower_bound_imaginary_axis", lower_ok)
    # scalar threshold: weight ratio sup_y |Phi(i beta y)| e^{-2 D y^2}
    out = {}
    for q in [mp.mpf("0.25"), mp.mpf("1"), mp.mpf("3")]:
        D = 4 * q / (1 + 4 * q)
        beta_c = mp.sqrt(4 * D)
        check(f"probit_threshold_formula_q={mp.nstr(q,3)}", abs(beta_c**2 - 16 * q / (1 + 4 * q)) < mp.mpf("1e-45"))
        ys = [mp.mpf(k) / 2 for k in range(1, 61)]
        at = [abs(Phi(mp.mpc(0, beta_c * y))) * mp.e ** (-2 * D * y**2) for y in ys]
        above = [abs(Phi(mp.mpc(0, 1.05 * beta_c * y))) * mp.e ** (-2 * D * y**2) for y in ys]
        check(f"probit_at_threshold_bounded_q={mp.nstr(q,3)}", max(at) <= C0, max_ratio=mp.nstr(max(at), 10))
        check(f"probit_above_threshold_grows_q={mp.nstr(q,3)}", above[-1] > 1e6 and above[-1] > above[-11],
              last_ratio=mp.nstr(above[-1], 10))
        out[mp.nstr(q, 3)] = mp.nstr(beta_c, 20)
    # logistic pole
    pole = mp.mpc(0, mp.pi)
    check("logistic_pole_at_i_pi", abs(1 + mp.e ** (-pole)) < mp.mpf("1e-45"))
    return {"C0": mp.nstr(C0, 20), "probit_beta_c": out}


def fock_gram() -> dict:
    rng = np.random.default_rng(20261125)
    alpha = 0.7
    d, N = 2, 3
    L = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    b = rng.normal(size=(N, d)) + 1j * rng.normal(size=(N, d))

    def inner(z, w):  # <z, w> = sum z_k conj(w_k)
        return np.sum(z * np.conj(w))

    def gram(points):
        G = np.empty((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                wi, wj = points[i], points[j]
                G[i, j] = np.exp(alpha * inner(wi, wj) - alpha * (abs(inner(wi, wi)) + abs(inner(wj, wj))) / 2)
        return G

    G0 = gram(b)
    ev0 = np.sort(np.linalg.eigvalsh((G0 + G0.conj().T) / 2))
    check("fock_G0_positive_definite", ev0[0] > 0, min_eig=ev0[0])
    worst = 0.0
    for _ in range(25):
        z = rng.normal(size=d) * 0.8 + 1j * rng.normal(size=d) * 0.8
        pts = [L @ z + b[j] for j in range(N)]
        Gz = gram(pts)
        evz = np.sort(np.linalg.eigvalsh((Gz + Gz.conj().T) / 2))
        worst = max(worst, float(np.max(np.abs(evz - ev0))))
    check("fock_gram_unitarily_conjugate", worst < 1e-10, max_eigen_difference=worst)
    return {"min_eigenvalue_G0": float(ev0[0]), "max_eigen_difference": worst}


def scalar_identities() -> dict:
    q, s2 = sp.symbols("q sigma2", positive=True)
    D1 = 1 - 1 / (1 + 4 * q)
    check("D_scalar_identity", sp.simplify(D1 - 4 * q / (1 + 4 * q)) == 0)
    check("bump_threshold_is_2D", sp.simplify(2 * D1 - 8 * q / (1 + 4 * q)) == 0)
    check("probit_threshold_is_4D", sp.simplify(4 * D1 - 16 * q / (1 + 4 * q)) == 0)
    Dw = 1 / s2 - 1 / (s2 + 4 * q)
    check("width_D_identity", sp.simplify(Dw - 4 * q / (s2 * (s2 + 4 * q))) == 0)
    ell2 = s2 * (s2 + 4 * q) / (8 * q)
    check("width_gate_threshold", sp.simplify(1 / ell2 - 2 * Dw) == 0)
    a = sp.symbols("alpha", positive=True)
    qc = a / (4 * (2 - a))
    check("q_c_inverts_bump_threshold", sp.simplify(8 * qc / (1 + 4 * qc) - a) == 0)
    return {}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/extensions.json")
    args = parser.parse_args()
    sections = {}
    for name, fn in [("triangular", triangular), ("warping_growth", warping_growth), ("probit", probit),
                     ("fock_gram", fock_gram), ("scalar_identities", scalar_identities)]:
        try:
            sections[name] = fn()
        except Exception as error:  # record and fail
            check(f"section_execution_{name}", False, error=repr(error))
            sections[name] = {"execution_error": repr(error)}
    failed = [c["name"] for c in CHECKS if not c["passed"]]
    report = {
        "title": "Checks for the revised additions (threshold, sharpness, gates, Fock Gram)",
        "scope": "Finite checks of constants and examples; not proofs.",
        "versions": {"python": platform.python_version(), "mpmath": mp.__version__,
                     "numpy": np.__version__, "sympy": sp.__version__},
        "sections": sections,
        "checks": CHECKS,
        "summary": {"check_count": len(CHECKS), "failed": failed, "status": "PASS" if not failed else "FAIL"},
    }
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    print(json.dumps(report["summary"]))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
