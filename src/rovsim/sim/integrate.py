"""Fixed-step fourth-order Runge-Kutta (RK4) integrator.

This module knows nothing about physics. It takes a derivative function
f(t, y) = dy/dt and steps the state y forward in time. It is unit-agnostic:
whatever units f and y use, the results come back in the same units.
"""

import math
from collections.abc import Callable

import numpy as np

# A derivative function: given time t and state y, return dy/dt (same shape as y).
Derivative = Callable[[float, np.ndarray], np.ndarray]


def rk4_step(f: Derivative, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Advance the state y from time t to t + h with one classical RK4 step."""
    # Four slope estimates across the step:
    k1 = f(t, y)                          # k1 = f(t,       y)
    k2 = f(t + h / 2, y + h / 2 * k1)     # k2 = f(t + h/2, y + h/2 * k1)
    k3 = f(t + h / 2, y + h / 2 * k2)     # k3 = f(t + h/2, y + h/2 * k2)
    k4 = f(t + h, y + h * k3)             # k4 = f(t + h,   y + h * k3)

    # Weighted average of the slopes: y(t + h) = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(
    f: Derivative,
    y0: np.ndarray,
    t_span: tuple[float, float],
    h: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate dy/dt = f(t, y) from t_span[0] to t_span[1] with step size h.

    Returns (times, states). times has shape (n,) and states has shape
    (n, len(y0)); row i of states is the state at times[i]. The first sample
    is the initial condition. The last step is shortened if needed so the
    final sample lands exactly on t_span[1].

    Raises ValueError if h is not a positive finite number, if t_end is before
    t_start, or if the state becomes NaN or infinite (the solution has diverged).

    Note: the divergence check looks at the state after each full step, not at
    the intermediate slopes k1..k4. Overflow warnings inside a step are
    silenced, so an infinity absorbed mid-step (e.g. 1/inf) could in principle
    give a finite final state that the check does not catch.
    """
    if not (math.isfinite(h) and h > 0):
        raise ValueError(f"step size h must be positive and finite, got {h}")

    t_start, t_end = t_span
    if t_end < t_start:
        raise ValueError(f"t_end must not be before t_start, got t_span = {t_span}")
    y = np.asarray(y0, dtype=float)

    # Number of steps; the small tolerance stops float noise (e.g. 10.000000001)
    # from adding a spurious extra step when h divides the span exactly.
    n_steps = math.ceil((t_end - t_start) / h - 1e-9)

    times = [t_start]
    states = [y]
    for i in range(n_steps):
        t = t_start + i * h
        # The last step ends exactly on t_end, so it may be shorter than h.
        t_next = t_end if i == n_steps - 1 else t_start + (i + 1) * h
        # Overflow warnings are silenced here because the check below reports
        # divergence as a clear error instead.
        with np.errstate(over="ignore", invalid="ignore"):
            y = rk4_step(f, t, y, t_next - t)
        if not np.isfinite(y).all():
            raise ValueError(
                f"state is not finite after step {i + 1} (t = {t_next}); "
                "the solution has diverged"
            )
        times.append(t_next)
        states.append(y)

    return np.array(times), np.vstack(states)
