"""One-degree-of-freedom vertical (heave) model of a body in water.

Convention: z is positive DOWN (SNAME body-fixed convention), so the heave
velocity w = dz/dt is positive when the body sinks and negative when it rises.

Fluid: still, fresh water (1000 kg/m^3 by default). The target site is a
Florida sinkhole, not the sea, so the default is NOT the 1025 kg/m^3 of
seawater.

Force balance (Newton's second law along z):

    m * w_dot = (m - rho*V) * g  -  0.5 * rho * Cd * A * w * |w|

- (m - rho*V) * g is weight minus buoyancy: positive (downward) when the body
  is heavier than the water it displaces.
- 0.5 * rho * Cd * A * w * |w| is quadratic drag. It is written w*|w|, not
  w**2, so that it always opposes the motion (a rising body, w < 0, is pushed
  back down).

Added mass is deliberately absent until M2. Omitting it makes the transient
faster than reality but leaves the terminal velocity exact.

The reference body used by the tests (10 kg, 8 L, Cd = 1.0, A = 0.05 m^2) is a
declared placeholder chosen so the arithmetic can be done by hand. It is NOT
Approach A's real numbers.

This module is pure: no file, network or console I/O, and no global state.
"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

G: float = 9.80665  # standard gravity, m/s^2


@dataclass(frozen=True)
class VerticalBody:
    """Physical parameters of a body moving vertically through water."""

    mass: float  # body mass, kg (must be > 0)
    volume: float  # displaced volume, m^3 (must be >= 0)
    drag_coefficient: float  # quadratic drag coefficient Cd, dimensionless (must be >= 0)
    reference_area: float  # plan-view reference area A, m^2 (must be >= 0)
    fluid_density: float = 1000.0  # fluid density rho, kg/m^3 (must be > 0); fresh water

    def __post_init__(self) -> None:
        # Strictly positive: zero mass divides by zero in heave_acceleration,
        # zero fluid density is not a fluid.
        if self.mass <= 0.0:
            raise ValueError(f"mass must be > 0, got {self.mass}")
        if self.fluid_density <= 0.0:
            raise ValueError(f"fluid_density must be > 0, got {self.fluid_density}")
        # Zero is allowed (legitimate limiting cases, e.g. zero drag), negative is not.
        if self.volume < 0.0:
            raise ValueError(f"volume must be >= 0, got {self.volume}")
        if self.reference_area < 0.0:
            raise ValueError(f"reference_area must be >= 0, got {self.reference_area}")
        if self.drag_coefficient < 0.0:
            raise ValueError(
                f"drag_coefficient must be >= 0, got {self.drag_coefficient}"
            )


def net_static_force(body: VerticalBody) -> float:
    """Weight minus buoyancy in newtons, positive downward: (m - rho*V) * g."""
    return (body.mass - body.fluid_density * body.volume) * G


def drag_force(body: VerticalBody, w: float) -> float:
    """Quadratic drag in newtons for heave velocity w (m/s); opposes motion."""
    return (
        -0.5
        * body.fluid_density
        * body.drag_coefficient
        * body.reference_area
        * w
        * abs(w)
    )


def heave_acceleration(body: VerticalBody, w: float) -> float:
    """Heave acceleration w_dot in m/s^2 at heave velocity w (m/s)."""
    return (net_static_force(body) + drag_force(body, w)) / body.mass


def derivative(body: VerticalBody) -> Callable[[float, np.ndarray], np.ndarray]:
    """Return f(t, y) for the integrator, mapping state y = [z, w] to [w, w_dot].

    The model is autonomous, so the time argument t is unused.
    """

    def f(t: float, y: np.ndarray) -> np.ndarray:
        w = y[1]
        return np.array([w, heave_acceleration(body, w)])

    return f


def terminal_velocity(body: VerticalBody) -> float:
    """Terminal heave velocity in m/s, positive downward (z positive down).

    At terminal velocity the acceleration is zero, so the net static force F
    balances drag: F = k * w_t * |w_t| with k = 0.5 * rho * Cd * A.
    Hence w_t = sign(F) * sqrt(|F| / k).

    Degenerate cases: no drag (Cd * A == 0) never stops accelerating, so the
    result is +inf or -inf by the sign of F; F == 0 (neutral buoyancy) gives 0.0.
    """
    force = net_static_force(body)
    if force == 0.0:
        return 0.0
    # Same constant as in drag_force; keep the two in step if the drag model changes.
    k = 0.5 * body.fluid_density * body.drag_coefficient * body.reference_area
    if k == 0.0:
        return float(np.copysign(np.inf, force))
    return float(np.copysign(np.sqrt(abs(force) / k), force))


def analytic_heave(
    body: VerticalBody, t: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Exact (z, w) at times t for release from rest at z = 0.

    With a0 = F / m and terminal velocity w_t, the equation
    m * w_dot = F - k * w * |w| integrates to
        w(t) = w_t * tanh(a0 * t / w_t)
        z(t) = (w_t**2 / a0) * ln(cosh(a0 * t / w_t)),
    and dz/dt = (w_t**2 / a0) * (a0 / w_t) * tanh(...) = w(t), so the pair closes.

    Degenerate cases: F == 0 returns zero arrays (the body never moves); with
    no drag (w_t infinite) nothing balances F, so the acceleration never
    changes and the limit of the formulas above is constant acceleration:
    w = a0*t and z = 0.5*a0*t**2 (a literal tanh/cosh evaluation would give
    inf * 0 = NaN).
    """
    t = np.asarray(t, dtype=float)
    force = net_static_force(body)
    if force == 0.0:
        return np.zeros_like(t), np.zeros_like(t)
    a0 = force / body.mass
    w_t = terminal_velocity(body)
    if np.isinf(w_t):
        return 0.5 * a0 * t**2, a0 * t
    x = a0 * t / w_t
    # ln(cosh x) = logaddexp(x, -x) - ln 2, which does not overflow for large |x|.
    return (w_t**2 / a0) * (np.logaddexp(x, -x) - np.log(2.0)), w_t * np.tanh(x)
