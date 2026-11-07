"""Generate the M0 validation figure and print the numbers the note quotes.

Integrates the reference body (a declared placeholder, not Approach A's real
geometry) from rest for 5 s with RK4 at h = 1e-3 s, compares w(t) with the
closed form, and saves the figure under docs/notes/figures/.

Run from anywhere: .venv/bin/python scripts/plot_m0.py
"""

import math
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
    net_static_force,
    terminal_velocity,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
FIGURE_PATH = REPO_ROOT / "docs" / "notes" / "figures" / "m0-terminal-velocity.png"

T_END = 5.0  # s
STEP = 1.0e-3  # s
STUDY_STEPS = (0.1, 0.05, 0.025)  # s, for the step-size study

# Terminal velocity worked by hand, to six figures: sqrt(19.6133 / 25.0).
# Derivation in docs/notes/m0-vertical-motion.md, section 2. Deliberately a
# constant rather than recomputed, so it stays an independent check.
HAND_CALCULATED_TERMINAL_VELOCITY = 0.885738  # m/s


def max_w_error(body: VerticalBody, h: float) -> float:
    """Max |w_sim - w_exact| over the run for step size h."""
    times, states = integrate(derivative(body), np.array([0.0, 0.0]), (0.0, T_END), h)
    _, w_exact = analytic_heave(body, times)
    return float(np.max(np.abs(states[:, 1] - w_exact)))


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
    # Distance from the asymptote w_t: the closed form alone sits this far away
    # at t = 5 s, so this is the physical transient, not integrator error.
    residual = abs(w_final - w_terminal) / abs(w_terminal)
    residual_predicted = 2.0 * math.exp(-2.0 * T_END / (body.mass * w_terminal / net_static_force(body)))
    integrator_error_final = abs(w_final - float(w_exact[-1])) / abs(w_terminal)
    rel_hand = abs(w_final - HAND_CALCULATED_TERMINAL_VELOCITY) / HAND_CALCULATED_TERMINAL_VELOCITY
    force = net_static_force(body)
    k = 0.5 * body.fluid_density * body.drag_coefficient * body.reference_area
    tau = body.mass * w_terminal / force
    z_exact, _ = analytic_heave(body, times)
    dev_z = float(np.max(np.abs(states[:, 0] - z_exact)))
    dev_w = float(np.max(np.abs(w_sim - w_exact)))
    errors = [max_w_error(body, h) for h in STUDY_STEPS]

    print(f"Net static force F:                        {force:.4f} N downward")
    print(f"Drag group k = 0.5*rho*Cd*A:               {k:.1f} N s^2/m^2")
    print(f"Hand-calculated terminal velocity:         {HAND_CALCULATED_TERMINAL_VELOCITY} m/s")
    print(f"Predicted terminal velocity (closed form): {w_terminal:.10f} m/s")
    print(f"Simulated w at t = {T_END:g} s:               {w_final:.10f} m/s")
    print(f"Distance from terminal velocity at t = {T_END:g} s:  {residual:.2e} (transient residual, relative)")
    print(f"  predicted analytically, 2*exp(-2t/tau):  {residual_predicted:.2e}")
    print(f"Integrator error at t = {T_END:g} s vs closed form: {integrator_error_final:.1e} (relative to w_t)")
    print(f"Relative error vs hand calculation:        {rel_hand:.2e}")
    print(f"Time constant tau = m*w_t/F:               {tau:.6f} s ({T_END / tau:.1f} tau in {T_END:g} s)")
    print(f"Depth at t = {T_END:g} s:                      {states[-1, 0]:.6f} m")
    print(f"Max trajectory deviation from closed form: {dev_z:.1e} m in z, {dev_w:.1e} m/s in w")
    print("Step-size study (max |w error| over the run):")
    for i, (h, e) in enumerate(zip(STUDY_STEPS, errors)):
        ratio = f"  ratio {errors[i - 1] / e:.1f}" if i else ""
        print(f"  h = {h:<6g} s: {e:.2e}{ratio}")
    print(f"Figure saved to {FIGURE_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
