# tests/test_vertical_validation.py
import numpy as np
import pytest
from rovsim.sim.integrate import integrate
from rovsim.sim.vertical import (
    VerticalBody,
    analytic_heave,
    derivative,
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
    assert np.max(np.abs(states[:, 1] - w_exact)) < 1e-6
    assert np.max(np.abs(states[:, 0] - z_exact)) < 1e-6


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


def test_buoyant_body_rises_to_negative_terminal_velocity():
    light = VerticalBody(
        mass=6.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    w_t = terminal_velocity(light)
    assert w_t < 0.0
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
    # Review Focus 5.
    dragless = VerticalBody(
        mass=10.0, volume=8.0e-3, drag_coefficient=0.0, reference_area=0.05
    )
    assert terminal_velocity(dragless) == float("inf")


def test_neutrally_buoyant_closed_form_is_zero_not_nan():
    # Review Focus 5.
    neutral = VerticalBody(
        mass=8.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    assert terminal_velocity(neutral) == 0.0
    z, w = analytic_heave(neutral, np.linspace(0.0, 5.0, 11))
    assert np.all(np.isfinite(z)) and np.all(np.isfinite(w))
    assert np.allclose(z, 0.0) and np.allclose(w, 0.0)
