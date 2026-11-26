#!/usr/bin/env python3
"""Plot the exact analytic boundary; no simulated or fitted observations."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / "figures"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.labelsize": 11, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    alpha = np.linspace(0, 1.94, 1400)
    q_critical = alpha / (4 * (2 - alpha))
    plum, moss, copper = "#674171", "#BCD19B", "#E3BEA0"
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    ax.fill_between(alpha, 0, q_critical, color=copper, alpha=0.75)
    ax.fill_between(alpha, q_critical, 1.5,
                    where=q_critical <= 1.5, color=moss, alpha=0.8)
    ax.plot(alpha, q_critical, color=plum, linewidth=2.8)
    for a, q, label in [(0.5, 1/12, r"$(1/2,\,1/12)$"),
                        (1, 1/4, r"$(1,\,1/4)$"),
                        (1.5, 3/4, r"$(3/2,\,3/4)$")]:
        ax.scatter([a], [q], s=38, color=plum, zorder=3)
        ax.annotate(label, (a, q), xytext=(-5, 14),
                    textcoords="offset points", ha="right", fontsize=10,
                    color="#3B2940")
    ax.text(0.53, 1.12, "Bounded operator", weight="bold", ha="center",
            color="#36442A", fontsize=12)
    ax.text(0.53, 0.99, "The boundary is included", ha="center",
            color="#36442A", fontsize=10)
    ax.text(1.43, 0.18, "Unbounded operator", weight="bold", ha="center",
            color="#704D37", fontsize=11)
    ax.text(1.23, 1.30, r"$q_c=\dfrac{\alpha}{4(2-\alpha)}$",
            color=plum, fontsize=15, ha="center")
    ax.set(xlim=(0, 1.94), ylim=(0, 1.5),
           xlabel=r"Weight concentration $\alpha$ in $w(x)=e^{-\alpha x^2}$",
           ylabel=r"Gaussian component variance $q$")
    ax.set_xticks([0, 0.5, 1, 1.5])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
    ax.set_title("Exact boundary for a nonlinear Markov transition",
                 loc="left", pad=13, fontsize=13, weight="bold")
    ax.grid(alpha=0.18, color="#84798A", linewidth=0.6)
    ax.set_axisbelow(True)
    for extension in ["pdf", "png"]:
        fig.savefig(output / f"phase-diagram.{extension}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
