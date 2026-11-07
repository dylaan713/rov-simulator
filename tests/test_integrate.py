import numpy as np
import pytest
from rovsim.sim.integrate import integrate, rk4_step


def test_single_step_matches_exponential_decay():
    # dy/dt = -y, y(0) = 1 -> y(h) = exp(-h). RK4 local error is O(h^5):
    # exactly h^5/120 = 8.3e-8 here, so 2e-7 leaves margin without going slack.
    h = 0.1
    y = rk4_step(lambda t, y: -y, 0.0, np.array([1.0]), h)
    assert y[0] == pytest.approx(np.exp(-h), abs=2e-7)


def test_fourth_order_convergence():
    # Halving h must cut the global error by roughly 2**4 = 16.
    def error(h):
        _, states = integrate(lambda t, y: -y, np.array([1.0]), (0.0, 1.0), h)
        return abs(states[-1, 0] - np.exp(-1.0))

    assert error(0.1) / error(0.05) > 10.0


def test_trajectory_shape_and_endpoints():
    times, states = integrate(lambda t, y: -y, np.array([1.0]), (0.0, 1.0), 0.1)
    assert times.shape == (11,)
    assert states.shape == (11, 1)
    assert times[0] == 0.0
    assert times[-1] == 1.0


def test_final_step_is_shortened_to_land_on_t_end():
    # h does not divide the span evenly: the last step is shortened to end on t_end.
    times, _ = integrate(lambda t, y: -y, np.array([1.0]), (0.0, 1.0), 0.3)
    assert times[-1] == 1.0
    assert np.all(np.diff(times) > 0)


@pytest.mark.parametrize("bad_h", [0.0, -0.1, float("nan"), float("inf")])
def test_non_positive_or_non_finite_step_is_rejected(bad_h):
    # A zero, negative, NaN or infinite step size must be refused before stepping.
    with pytest.raises(ValueError):
        integrate(lambda t, y: -y, np.array([1.0]), (0.0, 1.0), bad_h)


def test_diverging_solution_raises_instead_of_returning_nan():
    # dy/dt = y**2 from y(0) = 1 blows up at t = 1; this must be an error, not NaN.
    with pytest.raises(ValueError, match=r"not finite after step \d+ \(t = "):
        integrate(lambda t, y: y**2, np.array([1.0]), (0.0, 2.0), 0.01)


def test_time_argument_is_passed_at_the_right_stage_times():
    # dy/dt = t**3, y(0) = 0 -> y(1) = 1/4. RK4 weights are Simpson's rule, so a
    # cubic in t is integrated exactly; a wrong stage time (e.g. k2 using t
    # instead of t + h/2) would give a different answer.
    _, states = integrate(lambda t, y: np.array([t**3]), np.array([0.0]), (0.0, 1.0), 0.25)
    assert states[-1, 0] == pytest.approx(0.25, abs=1e-12)


def test_two_component_state_harmonic_oscillator():
    # y = [position, velocity], y' = [velocity, -position] -> [cos t, -sin t].
    times, states = integrate(
        lambda t, y: np.array([y[1], -y[0]]), np.array([1.0, 0.0]), (0.0, 2.0), 1e-3
    )
    assert states.shape == (times.shape[0], 2)
    assert states[-1, 0] == pytest.approx(np.cos(2.0), abs=1e-9)
    assert states[-1, 1] == pytest.approx(-np.sin(2.0), abs=1e-9)


def test_end_time_before_start_time_is_rejected():
    with pytest.raises(ValueError, match="t_end"):
        integrate(lambda t, y: -y, np.array([1.0]), (1.0, 0.0), 0.1)
