#!/usr/bin/env python3
"""Independent finite checks for the Gaussian-kernel boundary manuscript.

This script checks constants, identities, and selected examples.  It does not
prove the classification or establish novelty.  All checks use explicit
conditionals, so Python's -O option does not disable them.

Usage:
    python code/verify.py
    python code/verify.py --output /absolute/path/verification.json
"""

from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import dblquad, quad
from scipy.optimize import minimize_scalar
from scipy.special import factorial, jv, ndtr
import sympy as sp


CHECKS: list[dict] = []


def check(name: str, passed: bool, **details) -> None:
    CHECKS.append({"name": name, "passed": bool(passed), **details})


def exact_check(name: str, left, right) -> None:
    residual = sp.simplify(left - right)
    check(name, residual == 0, residual=str(residual), expected=str(right))


def relative_error(actual, expected) -> float:
    return float(abs(actual - expected) / max(1.0, abs(expected)))


def symbolic_constants() -> dict:
    xi, x, y = sp.symbols("xi x y", real=True)
    spectral_rows = []
    tube_rows = []
    for p in [sp.Rational(1, 2), sp.Integer(1), sp.Integer(3)]:
        rho = sp.sqrt(p) / (2 * sp.sqrt(sp.pi)) * sp.exp(-p * xi**2 / 4)
        mass = sp.integrate(rho, (xi, -sp.oo, sp.oo))
        second = sp.integrate(xi**2 * rho, (xi, -sp.oo, sp.oo))
        exact_check(f"spectral_mass_P={p}", mass, 1)
        exact_check(f"spectral_second_moment_P={p}", second, 2 / p)
        spectral_rows.append({"P": str(p), "mass": str(mass), "second_moment": str(second)})
        for t in [sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2)]:
            beta = t / p
            fourier_amplitude = sp.exp(-xi**2 / (4 * beta)) / (2 * sp.sqrt(sp.pi * beta))
            spectral_norm = sp.integrate(fourier_amplitude**2 / rho, (xi, -sp.oo, sp.oo))
            x_integral = sp.integrate(sp.exp(-2 * beta * x**2), (x, -sp.oo, sp.oo))
            y_integral = sp.integrate(sp.exp((2 * beta - 4 / p) * y**2), (y, -sp.oo, sp.oo))
            tube_norm = sp.simplify(2 * x_integral * y_integral / (sp.pi * p))
            expected = 1 / sp.sqrt(beta * p * (2 - beta * p))
            exact_check(f"Gaussian_probe_spectral_P={p}_beta={beta}", spectral_norm, expected)
            exact_check(f"Gaussian_probe_tube_P={p}_beta={beta}", tube_norm, expected)
            tube_rows.append({"P": str(p), "beta": str(beta), "norm_squared": str(expected)})

    threshold_rows = []
    for alpha in [sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2)]:
        critical = alpha / (4 * (2 - alpha))
        critical_slack = sp.simplify(8 * critical / (1 + 4 * critical) - alpha)
        exact_check(f"mixture_critical_slack_alpha={alpha}", critical_slack, 0)
        below = sp.simplify(8 * (critical - sp.Rational(1, 100)) / (1 + 4 * (critical - sp.Rational(1, 100))) - alpha)
        above = sp.simplify(8 * (critical + sp.Rational(1, 100)) / (1 + 4 * (critical + sp.Rational(1, 100))) - alpha)
        check(f"mixture_threshold_sides_alpha={alpha}", below < 0 and above > 0,
              below_slack=str(below), above_slack=str(above))
        eta = sp.symbols("eta", real=True)
        delta = sp.symbols("delta", positive=True)
        q = critical + delta
        a = (2 - alpha) / (8 * alpha)
        dilation = 2 / (2 - alpha)
        original = -(xi - eta)**2 / (4 * alpha) + xi**2 / 8 - eta**2 / 8 - q * eta**2 / 2
        completed = -a * (xi - dilation * eta)**2 - delta * eta**2 / 2
        exact_check(f"modulation_complete_square_alpha={alpha}", original, completed)
        hs_original = 1 / sp.sqrt(alpha * (4 * q * (2 - alpha) - alpha))
        hs_reduced = 1 / (2 * sp.sqrt(alpha * (2 - alpha) * delta))
        exact_check(f"modulation_HS_constant_alpha={alpha}", hs_original, hs_reduced)
        threshold_rows.append({"alpha": str(alpha), "q_critical": str(critical)})

    alpha, q = sp.Rational(1, 2), sp.Rational(1, 4)
    qc = alpha / (4 * (2 - alpha))
    u = sp.sqrt(alpha * (2 - alpha) * (q - qc))
    lam = sp.sqrt(1 + u**2) - u
    exact_check("modulation_exact_lambda_alpha_half_q_quarter", lam, 1 / sp.sqrt(2))
    exact_check("modulation_exact_HS_alpha_half_q_quarter", 1 / (2 * u), sp.sqrt(2))
    return {
        "method": "SymPy exact rational Gaussian integration and symbolic simplification",
        "spectral_parameter_count": 3,
        "Gaussian_probe_parameter_count": 9,
        "threshold_parameter_count": 3,
        "spectral": spectral_rows,
        "Gaussian_probe_norms": tube_rows,
        "critical_thresholds": threshold_rows,
    }


def direct_gram(P: float, M: float, q_values: np.ndarray, shifts: np.ndarray, y: float):
    count = len(shifts)
    actual = np.empty((count, count), dtype=complex)
    predicted = np.empty_like(actual)
    max_quad_estimate = 0.0
    rho_coefficient = math.sqrt(P) / (2 * math.sqrt(math.pi))
    for i in range(count):
        for j in range(count):
            shift = shifts[i] - shifts[j]
            spectral_width = P + 2 * (q_values[i] + q_values[j])

            # Direct integration of h_i(z) conjugate(h_j(z)) rho_P.
            # The real-state phase cancels.  Combining real exponents avoids
            # overflowing individual features on the infinite integration range.
            def amplitude(t):
                source_log_density = -P * t * t / 4
                first_noise_log_character = -q_values[i] * t * t / 2
                second_noise_log_character = -q_values[j] * t * t / 2
                imaginary_state_factor = -2 * M * y * t
                return rho_coefficient * math.exp(source_log_density + first_noise_log_character
                                                 + second_noise_log_character + imaginary_state_factor)

            real, er = quad(lambda t: amplitude(t) * math.cos(shift * t), -np.inf, np.inf,
                            epsabs=2e-10, epsrel=2e-11, limit=200)
            imag, ei = quad(lambda t: amplitude(t) * math.sin(shift * t), -np.inf, np.inf,
                            epsabs=2e-10, epsrel=2e-11, limit=200)
            actual[i, j] = real + 1j * imag
            linear = -2 * M * y + 1j * shift
            predicted[i, j] = math.sqrt(P / spectral_width) * np.exp(linear**2 / spectral_width)
            max_quad_estimate = max(max_quad_estimate, er + ei)
    error = float(np.linalg.norm(actual - predicted, ord=2) / max(1.0, np.linalg.norm(predicted, ord=2)))
    diagonal = np.sqrt(np.real(np.diag(predicted)))
    normalised = predicted / diagonal[:, None] / diagonal[None, :]
    minimum_eigenvalue = float(np.linalg.eigvalsh(normalised)[0])
    return error, max_quad_estimate, minimum_eigenvalue


def gram_quadrature() -> dict:
    P, M = 1.0, 0.8
    shifts = np.array([-1.0, 0.0, 2.0])
    ys = [-1.3, 0.0, 0.6, 1.1]
    common_rows, varied_rows = [], []
    for q in [0.0, 0.2, 1.0]:
        for y in ys:
            error, estimate, smallest = direct_gram(P, M, np.full(3, q), shifts, y)
            check(f"common_covariance_Gram_q={q}_y={y}", error < 2e-9,
                  scaled_matrix_error=error, tolerance=2e-9)
            common_rows.append({"q": q, "y": y, "scaled_matrix_error": error,
                                "quad_error_estimate": estimate, "normalised_min_eigenvalue": smallest})
    varied_q = np.array([0.0, 0.25, 0.9])
    for y in [-2.0, 0.0, 0.75, 2.0]:
        error, estimate, smallest = direct_gram(P, M, varied_q, shifts, y)
        check(f"varied_covariance_Gram_y={y}", error < 2e-9 and smallest > 0,
              scaled_matrix_error=error, minimum_eigenvalue=smallest, tolerance=2e-9)
        varied_rows.append({"y": y, "scaled_matrix_error": error,
                            "quad_error_estimate": estimate, "normalised_min_eigenvalue": smallest})
    return {
        "method": "independent real and imaginary adaptive quadrature over the real line",
        "P": P, "M": M, "shifts": shifts.tolist(),
        "common_covariance_matrix_count": 12,
        "varied_covariance_matrix_count": 4,
        "complex_entry_integral_count": 144,
        "real_quad_call_count": 288,
        "varied_covariances": varied_q.tolist(),
        "common_covariance_results": common_rows,
        "varied_covariance_results": varied_rows,
        "limitation": "Positive sampled Gram eigenvalues do not prove a uniform lower bound over all imaginary states.",
    }


def triangular_noise() -> dict:
    rows = []
    for slope in [1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 20.0]:
        c = (slope * slope - 1) / 4

        def log_profile(t):
            return -np.inf if t >= 1 else 2 * np.log1p(-t) + c * t * t

        grid = np.linspace(0.0, 1.0, 20001)
        values = np.full_like(grid, -np.inf)
        values[:-1] = 2 * np.log1p(-grid[:-1]) + c * grid[:-1]**2
        index = int(np.argmax(values))
        left = grid[max(0, index - 2)]
        right = grid[min(len(grid) - 2, index + 2)]
        optimum = minimize_scalar(lambda t: -log_profile(t), bounds=(left, right), method="bounded",
                                  options={"xatol": 1e-14})
        candidates = [(0.0, 0.0), (float(optimum.x), float(-optimum.fun))]
        numeric_t, numeric_log_max = max(candidates, key=lambda item: item[1])

        # Independent stationary-point formula: c t (1-t)=1.  The smaller
        # root is a local minimum; compare the larger root with the endpoint.
        exact_candidates = [(0.0, 0.0)]
        if c >= 4:
            upper_root = (1 + math.sqrt(1 - 4 / c)) / 2
            exact_candidates.append((upper_root, float(log_profile(upper_root))))
        exact_t, exact_log_max = max(exact_candidates, key=lambda item: item[1])
        log_error = abs(numeric_log_max - exact_log_max)
        actual_squared = math.exp(numeric_log_max) / abs(slope)
        exact_squared = math.exp(exact_log_max) / abs(slope)
        check(f"triangular_noise_norm_slope={slope}", log_error < 2e-9,
              log_max_error=log_error, tolerance=2e-9)
        if slope <= 4:
            check(f"triangular_noise_rational_norm_slope={slope}", abs(actual_squared - 1 / slope) < 1e-12,
                  squared_norm=actual_squared, expected_squared_norm=1 / slope, tolerance=1e-12)
        rows.append({"slope": slope, "numeric_maximiser": numeric_t, "stationary_formula_maximiser": exact_t,
                     "log_max_error": log_error, "numeric_squared_norm": actual_squared,
                     "formula_squared_norm": exact_squared, "log_squared_norm": exact_log_max - math.log(abs(slope))})
    return {
        "characteristic_function": "max(1-abs(xi),0)",
        "P": 1, "R": 1,
        "parameter_count": 7,
        "grid_points_per_parameter": 20001,
        "method": "global dense scan plus bounded optimisation, compared with the stationary-point formula",
        "results": rows,
    }


def bessel_warping() -> dict:
    a, epsilon, cutoff = 0.4, 0.7, 12
    xi_nodes, xi_weights = np.polynomial.legendre.leggauss(96)
    xis = a * xi_nodes
    # Multiplication by a and the density 1/(2a) leave these weights /2.
    weights = xi_weights / 2
    orders = np.arange(-cutoff, cutoff + 1)
    factorial_weights = factorial(abs(orders), exact=False)
    theta_count = 2048
    thetas = 2 * math.pi * np.arange(theta_count) / theta_count
    oscillations = np.exp(1j * epsilon * xis[:, None] * np.sin(thetas)[None, :])
    fourier_factors = np.exp(-1j * orders[:, None] * thetas[None, :])
    # Direct periodic trapezoidal quadrature of the defining Fourier integral.
    coefficient_integrals = oscillations @ fourier_factors.T / theta_count
    library_coefficients = jv(orders[None, :], epsilon * xis[:, None])
    coefficient_error = float(np.max(abs(coefficient_integrals - library_coefficients)))
    check("Bessel_coefficients_against_Fourier_quadrature", coefficient_error < 2e-13,
          maximum_absolute_error=coefficient_error, tolerance=2e-13)

    normaliser = 2 * math.e - 1
    direct_weighted = normaliser * float(np.sum(weights[:, None] * factorial_weights[None, :] * abs(coefficient_integrals)**2))
    bessel_weighted = normaliser * float(np.sum(weights[:, None] * factorial_weights[None, :] * abs(library_coefficients)**2))
    norm_error = abs(direct_weighted - bessel_weighted)
    check("Bessel_weighted_norm_g_equal_one", norm_error < 2e-11,
          absolute_error=norm_error, tolerance=2e-11)

    t_max = epsilon * a
    r = t_max**2 / 4
    # Bound only the omitted Bessel-series tail, not floating-point/quadrature error.
    series_tail_bound = (normaliser * 2 * math.exp(t_max**2 / 2)
                         * r**(cutoff + 1) / math.factorial(cutoff + 1)
                         / (1 - r / (cutoff + 2)))
    pointwise_rows = []
    for x in [-10.0, -2.3, 0.0, 0.7, 4.0, 10.0]:
        reconstructed = np.sum(weights[:, None] * coefficient_integrals
                               * np.exp(1j * (xis[:, None] + orders[None, :]) * x))
        warped = x + epsilon * math.sin(x)
        expected = float(np.sinc(a * warped / math.pi))
        error = float(abs(reconstructed - expected))
        check(f"Bessel_warped_probe_x={x}", error < 2e-11,
              absolute_error=error, tolerance=2e-11)
        pointwise_rows.append({"x": x, "expected_sinc": expected, "real_reconstruction": float(reconstructed.real),
                               "imaginary_reconstruction": float(reconstructed.imag), "absolute_error": error})

    evaluation_rows = []
    for y in [0.0, 0.4, 0.8]:
        source_laplace = 1.0 if y == 0 else math.sinh(2 * a * y) / (2 * a * y)
        closed = source_laplace * (math.exp(math.exp(2 * y)) + math.exp(math.exp(-2 * y)) - 1) / normaliser
        terms = np.arange(-60, 61)
        numerical = source_laplace * float(np.sum(np.exp(-2 * terms * y) / factorial(abs(terms)))) / normaliser
        error = relative_error(numerical, closed)
        check(f"factorial_target_evaluation_y={y}", error < 2e-13,
              scaled_error=error, tolerance=2e-13)
        evaluation_rows.append({"y": y, "closed_squared_evaluation_norm": closed,
                                "spectral_sum_squared_evaluation_norm": numerical, "scaled_error": error})
    return {
        "a": a, "epsilon": epsilon, "input_feature": "g(xi)=1, with source norm squared 1",
        "source_quadrature_nodes": 96, "periodic_trapezoidal_nodes": theta_count,
        "coefficient_order_range": [-cutoff, cutoff], "coefficient_comparison_count": 96 * (2 * cutoff + 1),
        "maximum_coefficient_error": coefficient_error,
        "direct_Fourier_coefficient_weighted_norm_squared": direct_weighted,
        "Bessel_weighted_norm_squared": bessel_weighted,
        "weighted_norm_absolute_difference": norm_error,
        "analytic_omitted_Bessel_tail_upper_bound": series_tail_bound,
        "pointwise_reconstructions": pointwise_rows,
        "target_evaluation_checks": evaluation_rows,
        "limitation": "The analytic tail bound covers the Bessel-series truncation only. Quadrature and floating-point checks are numerical, not interval-certified.",
    }


def modulation_kernel(xi, eta, alpha: float, q: float):
    exponent = -(xi - eta)**2 / (4 * alpha) + xi**2 / 8 - eta**2 / 8 - q * eta**2 / 2
    return np.exp(exponent) / (2 * math.sqrt(math.pi * alpha))


def nystrom_singular_values(alpha: float, q: float, nodes: int, cutoff_sd: float):
    qc = alpha / (4 * (2 - alpha))
    delta = q - qc
    a = (2 - alpha) / (8 * alpha)
    dilation = 2 / (2 - alpha)
    eta_sd = math.sqrt(1 / (2 * delta))
    xi_sd = math.sqrt(dilation**2 / (2 * delta) + 1 / (4 * a))
    eta_limit, xi_limit = cutoff_sd * eta_sd, cutoff_sd * xi_sd
    roots, weights = np.polynomial.legendre.leggauss(nodes)
    xis, etas = xi_limit * roots, eta_limit * roots
    wx, we = xi_limit * weights, eta_limit * weights
    matrix = (np.sqrt(wx)[:, None]
              * modulation_kernel(xis[:, None], etas[None, :], alpha, q)
              * np.sqrt(we)[None, :])
    singular = np.linalg.svd(matrix, compute_uv=False)
    return singular[:6]**2, eta_limit, xi_limit


def gaussian_modulation() -> dict:
    rows = []
    for alpha, q in [(0.5, 0.25), (1.0, 0.75), (1.5, 1.25)]:
        qc = alpha / (4 * (2 - alpha))
        delta = q - qc
        u = math.sqrt(alpha * (2 - alpha) * delta)
        lam = math.sqrt(1 + u * u) - u
        hs_expected = 1 / (2 * u)
        a = (2 - alpha) / (8 * alpha)
        dilation = 2 / (2 - alpha)
        eta_sd = math.sqrt(1 / (2 * delta))
        xi_sd = math.sqrt(dilation**2 / (2 * delta) + 1 / (4 * a))
        eta_limit, xi_limit = 9 * eta_sd, 9 * xi_sd

        # Genuine two-dimensional integration of |K(xi,eta)|^2 on a rectangle.
        # After normalisation this integrand is a centred bivariate Gaussian.
        # The union bound on its two marginal tails bounds rectangle truncation.
        hs_numeric, quadrature_error = dblquad(
            lambda xi, eta: float(modulation_kernel(xi, eta, alpha, q)**2),
            -eta_limit, eta_limit, lambda eta: -xi_limit, lambda eta: xi_limit,
            epsabs=1e-9, epsrel=1e-9,
        )
        analytic_rectangle_tail = hs_expected * 4 * float(ndtr(-9.0))
        hs_error = relative_error(hs_numeric, hs_expected)
        check(f"modulation_HS_2D_alpha={alpha}_q={q}", hs_error < 2e-8,
              scaled_error=hs_error, tolerance=2e-8)

        coarse, coarse_eta, coarse_xi = nystrom_singular_values(alpha, q, 140, 8.0)
        fine, fine_eta, fine_xi = nystrom_singular_values(alpha, q, 200, 9.0)
        expected = lam**(2 * np.arange(6) + 1)
        eigen_errors = abs(fine - expected) / expected
        grid_differences = abs(fine - coarse) / expected
        for n in range(6):
            check(f"modulation_Nystrom_alpha={alpha}_q={q}_n={n}",
                  eigen_errors[n] < 2e-7 and grid_differences[n] < 2e-7,
                  relative_formula_error=float(eigen_errors[n]), relative_grid_difference=float(grid_differences[n]),
                  tolerance=2e-7)
        rows.append({
            "alpha": alpha, "q": q, "q_critical": qc, "delta": delta, "u": u, "lambda": lam,
            "formula_HS_norm_squared": hs_expected, "direct_2D_HS_norm_squared": hs_numeric,
            "HS_scaled_error": hs_error, "dblquad_error_estimate": quadrature_error,
            "analytic_rectangle_tail_upper_bound": analytic_rectangle_tail,
            "HS_rectangle": {"eta": [-eta_limit, eta_limit], "xi": [-xi_limit, xi_limit]},
            "squared_singular_values_formula": expected.tolist(),
            "squared_singular_values_coarse_Nystrom": coarse.tolist(),
            "squared_singular_values_fine_Nystrom": fine.tolist(),
            "relative_singular_value_errors": eigen_errors.tolist(),
            "relative_Nystrom_grid_differences": grid_differences.tolist(),
            "coarse_Nystrom": {"nodes_per_axis": 140, "eta_limit": coarse_eta, "xi_limit": coarse_xi},
            "fine_Nystrom": {"nodes_per_axis": 200, "eta_limit": fine_eta, "xi_limit": fine_xi},
        })
    return {
        "parameter_count": 3, "two_dimensional_integral_count": 3,
        "squared_singular_values_per_parameter": 6, "Nystrom_matrices": 6,
        "method": "2D adaptive quadrature and two Gauss-Legendre Nystrom discretisations of the original Fourier kernel",
        "results": rows,
        "limitation": "The finite discretisations check the first six squared singular values. They do not prove the complete Mehler expansion or the classification theorem.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).resolve().parents[1] / "results" / "verification.json"
    parser.add_argument("--output", type=Path, default=default, help="JSON output file")
    args = parser.parse_args()
    report = {
        "title": "Independent finite verification of Gaussian-kernel boundary constants and examples",
        "scope": "Exact symbolic checks and deterministic numerical examples; not a proof or a priority assessment.",
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "sympy": sp.__version__},
        "random_seed": None,
        "randomness": "None: all parameter choices, quadrature rules and grids are deterministic.",
        "sections": {},
    }
    for name, function in [
        ("symbolic_constants", symbolic_constants),
        ("Gaussian_Gram_quadrature", gram_quadrature),
        ("triangular_noise", triangular_noise),
        ("Bessel_warping", bessel_warping),
        ("Gaussian_modulation", gaussian_modulation),
    ]:
        try:
            report["sections"][name] = function()
        except Exception as error:
            check(f"section_execution_{name}", False, exception_type=type(error).__name__, exception=str(error))
            report["sections"][name] = {"execution_error": str(error)}
    failed = [item for item in CHECKS if not item["passed"]]
    report["checks"] = CHECKS
    report["summary"] = {"check_count": len(CHECKS), "passed_count": len(CHECKS) - len(failed),
                         "failed_count": len(failed), "status": "PASS" if not failed else "FAIL"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output.resolve()), **report["summary"]}, sort_keys=True))
    if failed:
        print(json.dumps({"failed_checks": failed}, indent=2, sort_keys=True))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
