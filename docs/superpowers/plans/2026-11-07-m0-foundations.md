# M0 Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the repository and test harness, and simulate a body moving in one axis under weight, buoyancy and quadratic drag, validated against a hand-calculated terminal velocity and a closed-form transient.

**Architecture:** A fixed-step RK4 integrator that knows nothing about physics (`sim/integrate.py`, pure functions over NumPy state vectors) driving a 1-DOF heave model that knows nothing about integration (`sim/vertical.py`, pure functions over an immutable parameter object). The same two-sided purity is what makes the M2 6-DOF model a drop-in replacement for the derivative function, and it keeps `sim/` free of I/O, threads and clocks as the design requires.

**Tech Stack:** Python 3.12 (Homebrew, `/opt/homebrew/bin/python3.12`), NumPy, pytest, Matplotlib (plots only, never imported by `sim/`). SciPy is in the dependency list for later milestones but M0 does not use it — the integrator is hand-written because writing it is the point.

**Spec:** `docs/superpowers/specs/2026-11-05-rov-simulator-design.md`

## Global Constraints

- **Sign convention: SNAME body-fixed, z positive down.** Positive `z` is deeper, positive `w` is sinking, positive force is downward. Every sign in every test and docstring follows this.
- **Units are SI throughout, stated in every parameter's docstring.** m, kg, s, N, m/s, m², m³, kg/m³.
- `g = 9.80665 m/s²`, a single module-level constant `G` in `sim/vertical.py`. No other literal gravity anywhere.
- `rho = 1000.0 kg/m³` **fresh water** by default — the target site is a freshwater sinkhole. Overridable per body.
- **`sim/` is pure.** No file, network or console I/O; no threads; no reads of the wall clock; no global mutable state. Time is always a parameter.
- **Tests before implementation**, and each physics test asserts a number derived independently of the code under test (closed form or hand calculation), never a value captured from a previous run.
- **The owner must be able to read and explain every line.** Prefer an explicit loop over a clever one-liner; name variables after the physical quantity, not the symbol.
- Python >= 3.12. Type hints on every public function. No class where a function will do.

## Review Focus

Five conditions the design implies, that a reasonable person will hit, and that the happy-path tests would not catch. Each has a test pinned to the task that owns the code.

1. **Step size `h` <= 0 or non-finite** — silently returns the initial state or divides by zero; must raise `ValueError`. *(Task 2)*
2. **A diverging solution** — NaN/inf propagates and the run is reported as a success; must raise `ValueError` naming the step where the state stopped being finite. *(Task 2)*
3. **`t_span` not an integer multiple of `h`** — trajectory silently ends early or overshoots; the final step must be shortened so the last sample lands exactly on `t_end`. *(Task 2)*
4. **Non-positive mass, volume, area or fluid density** — produces a plausible-looking but meaningless trajectory; must raise `ValueError` at construction. *(Task 3)*
5. **Degenerate closed forms** — zero drag (`Cd` or `A` = 0) makes terminal velocity infinite, and a neutrally buoyant body makes the `tanh` solution divide by zero; both must return defined values (`inf`, and exact zeros) rather than NaN. *(Task 4)*

---

### Task 1: Repository, environment and test harness

**Files:**
- Create: `.gitignore`, `pyproject.toml`, `README.md`
- Create: `src/rovsim/__init__.py`, `src/rovsim/sim/__init__.py`
- Create: `tests/test_harness.py`

**Interfaces:**
- Consumes: nothing.
- Produces: an importable package `rovsim` with subpackage `rovsim.sim`; a `.venv` at the repo root; `pytest` runnable from the repo root with `src/` on the path.

**Requires explicit approval before running:** `git init` (first git write) and the dependency install.

- [ ] **Step 1: Initialise the repository**

```bash
cd /path/to/rov-simulator
git init -b main
```

- [ ] **Step 2: Write `.gitignore`**

Ignore: `.venv/`, `__pycache__/`, `*.pyc`, `.pytest_cache/`, `.DS_Store`, `*.egg-info/`, `build/`, `dist/`.
Do **not** ignore `docs/notes/figures/` — the validation plots are deliverables and are committed.

- [ ] **Step 3: Write `pyproject.toml`**

A single `[project]` table: name `rovsim`, version `0.1.0`, `requires-python = ">=3.12"`, dependencies `numpy`, `scipy`, `matplotlib`, and `[project.optional-dependencies] dev = ["pytest"]`. Build backend: setuptools with `[tool.setuptools] package-dir = {"" = "src"}`. A `[tool.pytest.ini_options]` table with `testpaths = ["tests"]` and `addopts = "-v"`.

- [ ] **Step 4: Create the virtual environment and install**

```bash
/opt/homebrew/bin/python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e ".[dev]"
```

Editable install is what puts `src/rovsim` on the path — no `sys.path` manipulation in tests, ever.

- [ ] **Step 5: Write the harness test**

```python
# tests/test_harness.py
def test_package_imports():
    import rovsim.sim
    assert rovsim.sim is not None


def test_numpy_available():
    import numpy as np
    assert int(np.__version__.split(".")[0]) >= 2
```

- [ ] **Step 6: Run the harness test**

Run: `.venv/bin/pytest`
Expected: 2 passed.

- [ ] **Step 7: Write `README.md`**

Title, one paragraph of what this is and why the simulator exists (compress the design's "Why this, and why the simulator"), a link to the design doc, how to set up and run the tests, and an empty **Predicted vs measured** table with columns `Quantity | Predicted (closed form) | Simulated | Error | Note`. Task 5 fills the first row.

- [ ] **Step 8: Commit**

```bash
git add .gitignore pyproject.toml README.md src tests docs HANDOFF.md
git commit -m "chore: initialise repository, venv and pytest harness"
```

---

### Task 2: Fixed-step RK4 integrator

**Files:**
- Create: `src/rovsim/sim/integrate.py`
- Create: `tests/test_integrate.py`

**Interfaces:**
- Consumes: nothing from Task 1 but the package path.
- Produces:
  - `rk4_step(f: Callable[[float, np.ndarray], np.ndarray], t: float, y: np.ndarray, h: float) -> np.ndarray`
  - `integrate(f: Callable[[float, np.ndarray], np.ndarray], y0: np.ndarray, t_span: tuple[float, float], h: float) -> tuple[np.ndarray, np.ndarray]` returning `(times, states)` where `times` has shape `(n,)`, `states` has shape `(n, len(y0))`, `times[0] == t_span[0]` and `times[-1] == t_span[1]` exactly.

- [ ] **Step 1: Write the failing tests**

```python
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `.venv/bin/pytest tests/test_integrate.py`
Expected: collection error — `ModuleNotFoundError: No module named 'rovsim.sim.integrate'`.

- [ ] **Step 3: Implement `rk4_step` and `integrate` in `src/rovsim/sim/integrate.py`**

`rk4_step` is the classical four-stage formula. Write it with the four slopes on separate named lines (`k1`..`k4`) and a comment giving the formula, because this function is the one the owner should be able to recite.

`integrate` validates `h` first, then steps from `t_span[0]`, shortening the last step so the final time equals `t_span[1]`, and raises `ValueError` naming the step index and time if `np.isfinite(y).all()` is false after any step. Accumulate into a Python list of arrays and `np.vstack` once at the end; do not grow an array in the loop.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `.venv/bin/pytest tests/test_integrate.py`
Expected: all pass (8 cases including the parametrised ones).

- [ ] **Step 5: Commit**

```bash
git add src/rovsim/sim/integrate.py tests/test_integrate.py
git commit -m "feat(sim): add fixed-step RK4 integrator with divergence guard"
```

---

### Task 3: Vertical body model — weight, buoyancy, quadratic drag

**Files:**
- Create: `src/rovsim/sim/vertical.py`
- Create: `tests/test_vertical.py`

**Interfaces:**
- Consumes: nothing (the model is independent of the integrator).
- Produces:
  - `G: float = 9.80665`
  - `VerticalBody` — frozen dataclass with fields `mass: float`, `volume: float`, `drag_coefficient: float`, `reference_area: float`, `fluid_density: float = 1000.0`, validating in `__post_init__`.
  - `net_static_force(body: VerticalBody) -> float` — weight minus buoyancy, positive downward: `(mass - fluid_density * volume) * G`.
  - `drag_force(body: VerticalBody, w: float) -> float` — `-0.5 * fluid_density * drag_coefficient * reference_area * w * abs(w)`, which opposes motion in both directions.
  - `heave_acceleration(body: VerticalBody, w: float) -> float` — `(net_static_force + drag_force) / mass`.
  - `derivative(body: VerticalBody) -> Callable[[float, np.ndarray], np.ndarray]` — returns the function `integrate` wants, mapping state `[z, w]` to `[w, heave_acceleration(body, w)]`.

The reference body used by every test in Tasks 3–5, chosen so the arithmetic can be done by hand:

| Field | Value | Note |
| --- | --- | --- |
| `mass` | 10.0 kg | |
| `volume` | 8.0e-3 m³ | 8 L displaced, so 2 kg net heavy in fresh water |
| `drag_coefficient` | 1.0 | blunt body, order-of-magnitude placeholder |
| `reference_area` | 0.05 m² | plan-view area |
| `fluid_density` | 1000.0 kg/m³ | fresh water |

These are **not** Approach A's real numbers and the docstring must say so.

**Deviation from the design, recorded deliberately.** The design's traceability table says M0 "produces the drag model the hardware track needs for `pc1`". This plan does not, for three reasons: no CAD or hull geometry exists to derive a coefficient from; `pc1` needs a design speed that is still TBD (it was to be set at the PDR, which was missed); and `pc1` is the owner's own calc under the mentor-only agreement that the software design explicitly leaves in force for the hardware track. M0 therefore validates the *model and the integrator* against a closed form using a declared placeholder body. The Approach-A coefficient derivation moves to M2, where the hull geometry it needs will exist. This goes in `docs/decisions.md` in Task 5.

- [ ] **Step 1: Write the failing tests**

```python
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
    # Review Focus 4.
    kwargs = dict(
        mass=10.0, volume=8.0e-3, drag_coefficient=1.0, reference_area=0.05
    )
    kwargs[field] = value
    with pytest.raises(ValueError):
        VerticalBody(**kwargs)
```

Note the deliberate asymmetry: `mass` and `fluid_density` must be strictly positive, while `volume`, `reference_area` and `drag_coefficient` may be zero (a zero-drag or zero-volume body is a legitimate limiting case to test against) but never negative. Encode exactly that.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `.venv/bin/pytest tests/test_vertical.py`
Expected: collection error — no module `rovsim.sim.vertical`.

- [ ] **Step 3: Implement `sim/vertical.py`**

Module docstring states the convention (z positive down), the fluid assumption and the force balance
`(m) * w_dot = (m - rho*V)*g - 0.5*rho*Cd*A*w*|w|`, plus one line saying added mass is deliberately absent until M2 and that its omission makes the transient faster than reality but leaves terminal velocity exact.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `.venv/bin/pytest tests/test_vertical.py`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/rovsim/sim/vertical.py tests/test_vertical.py
git commit -m "feat(sim): add 1-DOF heave model with weight, buoyancy and quadratic drag"
```

---

### Task 4: Closed-form validation — terminal velocity and transient

This is M0's gate: the simulator either matches an independently derived solution or M0 is not done.

**Files:**
- Modify: `src/rovsim/sim/vertical.py` (append the two closed-form functions)
- Create: `tests/test_vertical_validation.py`

**Interfaces:**
- Consumes: `VerticalBody`, `derivative`, `net_static_force` (Task 3); `integrate` (Task 2).
- Produces:
  - `terminal_velocity(body: VerticalBody) -> float` — signed, positive downward:
    `sign(F) * sqrt(2*|F| / (rho*Cd*A))` with `F = net_static_force(body)`; returns `+-inf` when `Cd*A == 0` and `0.0` when `F == 0`.
  - `analytic_heave(body: VerticalBody, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]` — exact `(z, w)` for release from rest at `z = 0`:
    `w(t) = w_t * tanh(a0 * t / w_t)` and `z(t) = (w_t**2 / a0) * log(cosh(a0 * t / w_t))`, where `a0 = F / m` is the initial acceleration and `w_t = terminal_velocity(body)`. Returns zero arrays when `w_t == 0`.

The hand calculation this task exists to confirm, done once, by hand, with the reference body:

```
F_net = (10.0 - 1000.0 * 0.008) * 9.80665      =  19.6133 N  (down)
k     = 0.5 * 1000.0 * 1.0 * 0.05              =  25.0  N s^2/m^2
w_t   = sqrt(F_net / k) = sqrt(0.784532)       =  0.885738 m/s
tau   = m * w_t / F_net = 10 * 0.885738/19.6133 = 0.451601 s (so 5 s is ~11 tau)
```

- [ ] **Step 1: Write the failing validation tests**

```python
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
```

- [ ] **Step 2: Run them to verify they fail**

Run: `.venv/bin/pytest tests/test_vertical_validation.py`
Expected: `ImportError: cannot import name 'terminal_velocity'`.

- [ ] **Step 3: Append `terminal_velocity` and `analytic_heave` to `src/rovsim/sim/vertical.py`**

Use `np.logaddexp`-free plain `np.log(np.cosh(x))` only if it does not overflow at the largest `x` the tests reach (`a0*t/w_t` is about 11 at t = 5 s, so `cosh` is ~3e4 — safe). Guard the degenerate branches before any division. Both docstrings must contain the derivation in two lines, because the note in Task 5 quotes them.

> **Superseded — the first sentence only.** The instruction to use plain `np.log(np.cosh(x))` was wrong, and its reasoning was the error: it scopes safety to the range the tests happen to reach rather than to the function's domain, and `cosh` overflows above argument 710 (about t = 320 s for the reference body). The implemented replacement is `logaddexp(x, -x) - ln 2`; see `docs/decisions.md` and entry 3 of `docs/ai-assistance-log.md`. The step above is left exactly as written, as a historical record. The rest of the step still stands.

- [ ] **Step 4: Run the whole suite**

Run: `.venv/bin/pytest`
Expected: every test passes. Record the simulated terminal velocity printed by Task 5's script for the README table.

- [ ] **Step 5: Commit**

```bash
git add src/rovsim/sim/vertical.py tests/test_vertical_validation.py
git commit -m "feat(sim): validate heave model against closed-form terminal velocity"
```

---

### Task 5: M0 deliverables — technical note, plot, decision log, AI-assistance log

The design makes documentation part of each milestone's definition of done. M0 is not done without this task.

**Files:**
- Create: `scripts/plot_m0.py`
- Create: `docs/notes/m0-vertical-motion.md`
- Create: `docs/decisions.md`
- Create: `docs/ai-assistance-log.md`
- Create: `docs/notes/figures/m0-terminal-velocity.png` (generated, committed)
- Modify: `README.md` (fill the predicted-vs-measured row)
- Modify: `HANDOFF.md` (state update)

**Interfaces:**
- Consumes: everything from Tasks 2–4. `scripts/` may import Matplotlib; `src/rovsim/sim/` still may not.
- Produces: the figure and the numbers quoted in the note and the README.

- [ ] **Step 1: Write `scripts/plot_m0.py`**

A `main()` that integrates the reference body for 5 s at h = 1e-3, plots simulated `w(t)` against the closed form plus a dashed horizontal line at the hand-calculated terminal velocity, labels axes with units, saves to `docs/notes/figures/m0-terminal-velocity.png` at 150 dpi, and prints the three numbers the note and README quote: predicted terminal velocity, simulated final velocity, and relative error. It takes no arguments and writes paths relative to the repo root resolved from `__file__`.

- [ ] **Step 2: Run it**

Run: `.venv/bin/python scripts/plot_m0.py`
Expected: the figure exists and the printed relative error is below 1e-3.

- [ ] **Step 3: Write `docs/notes/m0-vertical-motion.md`**

Sections, in this order: what is being modelled (the force balance, with the convention); the derivation of terminal velocity and of the `tanh` transient; what is deliberately *not* modelled yet (added mass, tether, thrusters) and which direction each omission biases the answer; the validation method; the figure; the result table; and what the next milestone needs from this one. Quote the real printed numbers from Step 2, not the plan's.

- [ ] **Step 4: Start `docs/decisions.md` and `docs/ai-assistance-log.md`**

`decisions.md`: a table `Date | Decision | Why | Rejected alternative`. Seed it with the decisions this plan made — Python 3.12 over system 3.9, hand-written RK4 over `scipy.integrate.solve_ivp`, fresh water over seawater, added mass deferred to M2, `src/` layout, the non-negative-vs-strictly-positive parameter rule, and the `pc1` deviation recorded in Task 3.

`ai-assistance-log.md`: a table `Date | What the assistant produced | What was wrong or missing | How it was caught`. Seed it with M0's actual record, written honestly — including anything the owner had to correct during execution. An empty or flattering log makes this deliverable worthless.

- [ ] **Step 5: Update `README.md` and `HANDOFF.md`**

README: fill the terminal-velocity row of the predicted-vs-measured table and link the note.
HANDOFF: `Done:` M0 complete with the validation result; `Next step:` M1 — thruster model, 3-thruster configuration matrix and pilot mixing.

- [ ] **Step 6: Run the full suite one last time and commit**

Run: `.venv/bin/pytest`
Expected: all pass.

```bash
git add scripts docs README.md HANDOFF.md
git commit -m "docs: add M0 technical note, validation plot, decision and AI-assistance logs"
```

---

## Open items, not part of this plan

- **Public GitHub repository.** The design assumes public from M0 and notes it is reversible in only one direction. `gh` is authenticated as `dylaan713`. No remote is created by this plan; it needs a separate yes.
- **Dependency install** (`numpy`, `scipy`, `matplotlib`, `pytest` into `.venv`) needs explicit approval before Task 1 Step 4 runs.
- **Hours per week** is still unestablished, so this plan carries no dates.
