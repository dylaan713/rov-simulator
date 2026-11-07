# tests/test_vertical.py
import numpy as np
import pytest
from rovsim.sim.vertical import (
    G,
    VerticalBody,
    derivative,
    drag_force,
    heave_acceleration,
    net_static_force,
)

REFERENCE = VerticalBody(
    mass=10.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
)


def test_heavy_body_has_net_downward_static_force():
    # (10 kg - 1000 kg/m3 * 0.008 m3) * 9.80665 = 2 kg * 9.80665 = 19.6133 N down
    assert net_static_force(REFERENCE) == pytest.approx(19.6133, abs=1e-4)


def test_buoyant_body_has_net_upward_static_force():
    light = VerticalBody(
        mass=6.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    assert net_static_force(light) == pytest.approx(-2.0 * G, abs=1e-9)


def test_neutrally_buoyant_body_is_in_equilibrium_at_rest():
    neutral = VerticalBody(
        mass=8.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    assert net_static_force(neutral) == pytest.approx(0.0, abs=1e-12)
    assert heave_acceleration(neutral, 0.0) == pytest.approx(0.0, abs=1e-12)


def test_drag_opposes_motion_in_both_directions():
    # k = 0.5 * rho * Cd * A = 0.5 * 1000 * 1.0 * 0.05 = 25.0 N s^2/m^2
    # At w = 2 m/s: 25.0 * 2**2 = 100 N, opposing the motion.
    assert drag_force(REFERENCE, 2.0) == pytest.approx(-100.0, abs=1e-9)
    assert drag_force(REFERENCE, -2.0) == pytest.approx(100.0, abs=1e-9)
    assert drag_force(REFERENCE, 0.0) == 0.0


def test_drag_is_quadratic_in_speed():
    assert drag_force(REFERENCE, 4.0) == pytest.approx(
        4.0 * drag_force(REFERENCE, 2.0), rel=1e-12
    )


def test_derivative_maps_state_to_velocity_and_acceleration():
    f = derivative(REFERENCE)
    dstate = f(0.0, np.array([3.0, 2.0]))
    assert dstate[0] == pytest.approx(2.0)
    assert dstate[1] == pytest.approx(heave_acceleration(REFERENCE, 2.0))


def test_zero_volume_area_and_drag_are_accepted():
    # Zero is a legitimate limiting case for these three: a body that displaces
    # no water, presents no area, or has no drag (zero drag gives an infinite
    # terminal velocity, which Task 4 tests). Only negative values are unphysical.
    body = VerticalBody(
        mass=10.0, volume=0.0, drag_coefficient=0.0, reference_area=0.0
    )
    assert body.volume == 0.0
    assert body.drag_coefficient == 0.0
    assert body.reference_area == 0.0


def test_heave_acceleration_at_rest_is_static_force_over_mass():
    # No drag at w = 0: 19.6133 N / 10 kg = 1.96133 m/s^2, downward
    assert heave_acceleration(REFERENCE, 0.0) == pytest.approx(1.96133, abs=1e-5)


def test_heave_acceleration_decelerates_when_drag_exceeds_net_weight():
    # At w = 2 m/s drag is 100 N up, net weight is 19.6133 N down:
    # (19.6133 - 100) / 10 = -8.03867 m/s^2, so the body is decelerating.
    assert heave_acceleration(REFERENCE, 2.0) == pytest.approx(-8.03867, abs=1e-5)


@pytest.mark.parametrize(
    "field,value",
    [
        ("mass", 0.0),
        ("mass", -1.0),
        ("volume", -1e-3),
        ("reference_area", -0.1),
        ("fluid_density", 0.0),
        ("drag_coefficient", -0.5),
    ],
)
def test_unphysical_parameters_are_rejected(field, value):
    # Parameters that are negative, or zero where zero is unphysical
    # (mass, fluid_density), are rejected.
    kwargs = dict(
        mass=10.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    kwargs[field] = value
    with pytest.raises(ValueError):
        VerticalBody(**kwargs)
