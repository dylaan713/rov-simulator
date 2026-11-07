import numpy as np
import pytest
from rovsim.sim.integrate import integrate
from rovsim.sim.vertical import (
    VerticalBody,
    analytic_heave,
    derivative,
    net_static_force,
    terminal_velocity,
)

REFERENCE = VerticalBody(
    mass=10.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
)
HAND_CALCULATED_TERMINAL_VELOCITY = 0.885738  # m/s down; see plan Task 4


def test_terminal_velocity_matches_hand_calculation():
    assert terminal_velocity(REFERENCE) == pytest.approx(
        HAND_CALCULATED_TERMINAL_VELOCITY, rel=1e-5
    )


def test_simulated_body_reaches_hand_calculated_terminal_velocity():
    # THE M0 VALIDATION. Released from rest, integrated for ~11 time constants.
    f = derivative(REFERENCE)
    _, states = integrate(f, np.array([0.0, 0.0]), (0.0, 5.0), 1e-3)
    assert states[-1, 1] == pytest.approx(
        HAND_CALCULATED_TERMINAL_VELOCITY, rel=1e-3
    )


def test_simulated_trajectory_matches_closed_form_everywhere():
    f = derivative(REFERENCE)
    times, states = integrate(f, np.array([0.0, 0.0]), (0.0, 5.0), 1e-3)
    z_exact, w_exact = analytic_heave(REFERENCE, times)
    assert np.max(np.abs(states[:, 1] - w_exact)) < 1e-9
    assert np.max(np.abs(states[:, 0] - z_exact)) < 1e-9


def test_step_size_controls_accuracy():
    # Why step size matters: the error must fall steeply as h shrinks.
    f = derivative(REFERENCE)

    def error(h):
        times, states = integrate(f, np.array([0.0, 0.0]), (0.0, 1.0), h)
        _, w_exact = analytic_heave(REFERENCE, times)
        return np.max(np.abs(states[:, 1] - w_exact))

    errors = [error(h) for h in (0.1, 0.05, 0.025)]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 1e-6
    # RK4 is fourth order: halving h cuts the error by 2**4 = 16. Measured
    # ratios are about 17 and 18; requiring > 14 leaves ~20% margin below that
    # and excludes an order-3 method (ratio near 8) decisively, not marginally.
    assert errors[0] / errors[1] > 14.0
    assert errors[1] / errors[2] > 14.0


def test_buoyant_body_rises_to_negative_terminal_velocity():
    light = VerticalBody(
        mass=6.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    w_t = terminal_velocity(light)
    assert w_t < 0.0
    # Same |F| and k as REFERENCE, so the hand value applies with a minus sign.
    assert w_t == pytest.approx(-HAND_CALCULATED_TERMINAL_VELOCITY, rel=1e-5)
    _, states = integrate(derivative(light), np.array([0.0, 0.0]), (0.0, 5.0), 1e-3)
    assert states[-1, 1] == pytest.approx(w_t, rel=1e-3)
    assert states[-1, 0] < 0.0  # shallower than it started


def test_body_released_above_terminal_speed_decelerates_to_it():
    f = derivative(REFERENCE)
    fast = 3.0 * HAND_CALCULATED_TERMINAL_VELOCITY
    _, states = integrate(f, np.array([0.0, fast]), (0.0, 5.0), 1e-3)
    steps = np.diff(states[:, 1])
    assert steps[0] < 0.0
    assert np.all(steps <= 0.0)  # monotone decay, never overshoots below w_t
    assert states[-1, 1] == pytest.approx(
        HAND_CALCULATED_TERMINAL_VELOCITY, rel=1e-3
    )


def test_zero_drag_gives_infinite_terminal_velocity():
    # degenerate case: zero drag
    dragless = VerticalBody(
        mass=10.0, volume=8.0e-3, drag_coefficient=0.0, reference_area=0.05
    )
    assert terminal_velocity(dragless) == float("inf")


def test_neutrally_buoyant_closed_form_is_zero_not_nan():
    # degenerate case: neutral buoyancy
    neutral = VerticalBody(
        mass=8.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    assert terminal_velocity(neutral) == 0.0
    z, w = analytic_heave(neutral, np.linspace(0.0, 5.0, 11))
    assert np.all(np.isfinite(z)) and np.all(np.isfinite(w))
    assert np.allclose(z, 0.0) and np.allclose(w, 0.0)


def test_zero_drag_gives_constant_acceleration_motion():
    # With no drag there is nothing to balance the net static force, so the
    # body never reaches a terminal velocity: w = a0*t and z = a0*t**2/2.
    dragless = VerticalBody(
        mass=10.0, volume=8.0e-3, drag_coefficient=0.0, reference_area=0.05
    )
    a0 = net_static_force(dragless) / dragless.mass  # 1.96133 m/s**2
    t = np.array([0.0, 1.0, 2.0, 3.0])
    assert a0 == pytest.approx(1.96133, rel=1e-5)
    z, w = analytic_heave(dragless, t)
    assert np.allclose(w, a0 * t, rtol=0.0, atol=1e-12)
    assert np.allclose(z, 0.5 * a0 * t**2, rtol=0.0, atol=1e-12)


def test_closed_form_stays_finite_at_large_time():
    # ln(cosh x) overflows for x > ~710 (t > ~320 s here); the closed form must not.
    z, w = analytic_heave(REFERENCE, np.array([1000.0]))
    assert np.all(np.isfinite(z)) and np.all(np.isfinite(w))
    assert w[0] == pytest.approx(HAND_CALCULATED_TERMINAL_VELOCITY, rel=1e-5)
