# tests/test_integrate.py
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
    assert times[-1] == pytest.approx(1.0, abs=1e-12)


def test_final_step_is_shortened_to_land_on_t_end():
    # Review Focus 3: h does not divide the span evenly.
    times, _ = integrate(lambda t, y: -y, np.array([1.0]), (0.0, 1.0), 0.3)
    assert times[-1] == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.diff(times) > 0)


@pytest.mark.parametrize("bad_h", [0.0, -0.1, float("nan"), float("inf")])
def test_non_positive_or_non_finite_step_is_rejected(bad_h):
    # Review Focus 1.
    with pytest.raises(ValueError):
        integrate(lambda t, y: -y, np.array([1.0]), (0.0, 1.0), bad_h)


def test_diverging_solution_raises_instead_of_returning_nan():
    # Review Focus 2: dy/dt = y**2 from y(0) = 1 blows up at t = 1.
    with pytest.raises(ValueError, match="finite"):
        integrate(lambda t, y: y**2, np.array([1.0]), (0.0, 2.0), 0.01)
