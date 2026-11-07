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
