"""Generate the M0 validation figure and print the numbers the note quotes.

Integrates the reference body (a declared placeholder, not Approach A's real
geometry) from rest for 5 s with RK4 at h = 1e-3 s, compares w(t) with the
closed form, and saves the figure under docs/notes/figures/.

Run from anywhere: .venv/bin/python scripts/plot_m0.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # file output only; no display needed
import matplotlib.pyplot as plt
import numpy as np

from rovsim.sim.integrate import integrate
from rovsim.sim.vertical import (
    VerticalBody,
    analytic_heave,
    derivative,
    terminal_velocity,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURE_PATH = REPO_ROOT / "docs" / "notes" / "figures" / "m0-terminal-velocity.png"

T_END = 5.0  # s
STEP = 1.0e-3  # s


def main() -> None:
    body = VerticalBody(mass=10.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05)

    times, states = integrate(derivative(body), np.array([0.0, 0.0]), (0.0, T_END), STEP)
    w_sim = states[:, 1]
    _, w_exact = analytic_heave(body, times)
    w_terminal = terminal_velocity(body)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(times, w_exact, color="tab:blue", linewidth=3.0, alpha=0.45,
            label="closed form: $w_t \\tanh(a_0 t / w_t)$")
    ax.plot(times, w_sim, color="tab:red", linewidth=1.2, linestyle=":",
            label=f"RK4 simulation, h = {STEP:g} s")
    ax.axhline(w_terminal, color="black", linestyle="--", linewidth=1.0,
               label=f"hand-calculated terminal velocity, {w_terminal:.6f} m/s")
    ax.set_xlabel("time t (s)")
    ax.set_ylabel("heave velocity w (m/s, positive down)")
    ax.set_title("Reference body released from rest: approach to terminal velocity")
    ax.set_xlim(0.0, T_END)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    plt.close(fig)

    w_final = float(w_sim[-1])
    relative_error = abs(w_final - w_terminal) / abs(w_terminal)
    print(f"Predicted terminal velocity (closed form): {w_terminal:.10f} m/s")
    print(f"Simulated w at t = {T_END:g} s:               {w_final:.10f} m/s")
    print(f"Relative error:                            {relative_error:.2e}")
    print(f"Figure saved to {FIGURE_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
