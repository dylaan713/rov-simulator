# M1 — Thrusters and Allocation — Design

Date: 2026-11-19
Status: Draft for approval
Parent: `docs/superpowers/specs/2026-11-05-rov-simulator-design.md` (milestone M1)
Requirements served: R-3 (forward, backward, turn, rise, descend)

## Intent

Give the simulated vehicle thrusters, decide how pilot commands become thruster
commands, and let it drive. M1 is the first point where the simulator stops
being a falling object and becomes a vehicle that responds to a pilot.

The software is built by AI coding agents under the owner's direction
(`docs/decisions.md`, 2026-11-19), so this spec is written to be built from without
going back to the owner. Every value below that describes the real vehicle is a
**declared placeholder**: thruster ratings, positions and the hull's drag belong
to the hardware track, which is parked and stays the owner's.

## Goals

- A thruster model with **deadband, saturation and asymmetric reverse thrust**.
- A thruster configuration matrix built from thruster positions and directions,
  not typed in by hand.
- Allocation: a requested force and moment becomes three thruster commands, and
  a request the vehicle physically cannot produce is **rejected, not
  approximated**.
- Pilot mixing for R-3: forward/back, turn, rise/descend in, thruster commands
  out, with saturation that keeps the pilot's intended turn.
- A 4-DOF body (surge, heave, yaw, plus position on a map) so the vehicle can be
  driven and plotted, validated against closed forms.
- Technical note, plot, decision-log entries and AI-assistance-log entries, as
  M0 did.
- **The owner can explain every concept in M1** — deadband, saturation, the
  configuration matrix, underactuation, allocation, why the horizontal pair is
  scaled together — without needing to follow the code or redo the arithmetic.

## Non-goals

- 6-DOF, added mass, coupling between axes, restoring moments. All M2.
- Thruster spin-up lag (a first-order delay between command and thrust).
  Deferred to M3, where it matters for the control loop and the failsafe.
- Gamepad input, networking, the vehicle control loop. M3 and M4.
- Real thruster curves or geometry. Placeholders until the owner supplies them.

## Anti-goals

- **A lateral command silently turned into something else.** The vehicle has no
  sway thruster; asking for sideways force must fail loudly.
- **Saturation that changes what the pilot asked for.** Clipping each thruster
  independently at full forward-and-turn quietly kills the turn. That is the
  common bug this milestone must not have.
- **Placeholder numbers that look like the real vehicle.** Every placeholder is
  labelled as one in its docstring and in the note, as M0's reference body was.
- **The simulator becoming the project.** M1 stays small; the ground station
  (M4) is what satisfies most requirements, and the project deadline is fixed.

## Constraints

- Everything in the parent spec and in `docs/decisions.md` stands: hand-written
  RK4 (`sim/integrate.py`, reused unchanged), fresh water 1000 kg/m^3, SNAME
  body-fixed axes, SI units, `pytest -W error`, `sim/` pure (no I/O, threads or
  clock), Python 3.12, tests before implementation, every physics test asserting
  a number derived independently of the code.
- `sim/vertical.py` is not modified. M0's validation stays exactly as published.
- The three-layer boundary: thruster physics is `sim/`; allocation and pilot
  mixing are vehicle control software, `vehicle/`. `vehicle/` must not import
  `sim/`.

## Approaches considered

| | Approach | Trade-off |
|---|---|---|
| A | Thruster model and allocation only, no body | Smallest, and every check is algebra. But nothing moves, so there is nothing to watch and no closed-form driving test; the allocation is never exercised by dynamics until M2 |
| **B — chosen** | Thruster model, allocation, and a simple uncoupled 4-DOF body | Moderate size. The vehicle drives, and a steady turn has an exact closed-form radius to validate against. The body is not wasted at M2: a 6-DOF model with coupling switched off must reproduce it, so it becomes M2's test oracle |
| C | Skip ahead to full 6-DOF in M1 | Merges M1 and M2 into one long stretch with no checkpoint, and the allocation bugs would be tangled up with coupling bugs. Rejected |

B was the assistant's call; the owner handed over the how.

## Design

### Axes and signs (SNAME)

Body axes: x forward, y starboard, z down. Forces X, Y, Z; moments K (roll),
M (pitch), N (yaw). Positive N turns the nose to starboard (a right turn, seen
from above). Moment from a force F at position r is r × F, so a forward force
X_i at lateral offset y_i gives N = −y_i · X_i: the **port** thruster pushing
forward turns the vehicle **right**. World frame is NED: north, east, down, with
heading psi measured from north toward east.

### Thruster model — `sim/thruster.py`

Command u is dimensionless in [−1, 1] (the ESC command). Steady thrust:

    |u| <= d            ->  T = 0                       (deadband)
    s = (|u| - d)/(1 - d)                                (rescaled, 0..1)
    u > d               ->  T = +T_fwd_max * s^2
    u < -d              ->  T = -T_rev_max * s^2

Quadratic because propeller thrust goes as rpm squared and ESC command is
roughly proportional to rpm. |u| > 1 is saturated to 1 (a command the ESC
clips), **not** rejected. On the vehicle side, out-of-range commands are a bug
the vehicle software must never produce, and a test asserts it never does (see
*Allocation* for where that guarantee lives).

Placeholder thruster, identical for all three:
`T_fwd_max = 10.0 N`, `T_rev_max = 8.0 N`, `d = 0.05`.

The inverse (thrust wanted -> command) belongs to the vehicle side and lives in
`vehicle/`, with its own copy of the thruster parameters. Inverting through the
deadband is **deadband compensation**: a small non-zero thrust request maps to a
command just past the deadband, not into it. Zero maps to exactly zero.

### Geometry and the configuration matrix — `rovsim/geometry.py`

A `Thruster` has a position (m, body frame) and a unit direction. The 6×3
configuration matrix B has one column per thruster: `[direction; position ×
direction]`. Wrench = B · thrusts.

`geometry.py` sits outside both `sim/` and `vehicle/` because both need the type.
They hold separate instances: the simulator's is the truth, the vehicle's is
what the control software believes. They are equal for now; M6 is where they
start to differ.

Placeholder layout (`l = 0.15 m`):

| Thruster | Position (x, y, z) | Direction |
|---|---|---|
| port | (0, −0.15, 0) | (1, 0, 0) |
| starboard | (0, +0.15, 0) | (1, 0, 0) |
| vertical | (0, 0, 0) | (0, 0, 1) — positive thrust pushes **down** |

So B's rows are: X = [1, 1, 0], Y = [0, 0, 0], Z = [0, 0, 1],
K = [0, 0, 0], M = [0, 0, 0], N = [0.15, −0.15, 0].

### Allocation — `vehicle/allocation.py`

The vehicle controls exactly three axes: X, Z, N. The 3×3 block of B for those
rows is inverted directly (it is square and, for the placeholder, invertible:
T_port = X/2 + N/(2l), T_stbd = X/2 − N/(2l), T_vert = Z). Requests:

- A requested wrench with any non-zero Y, K or M raises `ValueError` naming the
  axis. Exact zero, not a tolerance: the pilot path never produces these, so a
  non-zero value is always a caller's bug.
- If the layout makes the controlled block singular, construction raises.
- Allocation also returns the **parasitic wrench**: what B produces on Y, K, M
  from the allocated thrusts. Zero for the placeholder; non-zero once a real
  layout has offsets, and reported rather than hidden.

Allocation is split into three steps, all in `vehicle/allocation.py`, because the
pilot path (M1) and the M3 controller both call them and must get the same
guarantee:

1. `allocate(wrench)` — the linear inverse above. **No magnitude check**: it may
   return thrusts beyond a thruster's limit. It answers "what would it take".
2. `limit(thrusts)` — the saturation policy: the horizontal pair scaled by one
   shared factor (see *Pilot mixing*), the vertical thruster clipped on its own.
   Its output is always within every thruster's directional limit.
3. `to_commands(thrusts)` — the inverse thruster model with deadband
   compensation. A thrust beyond its limit by more than 1e-9 relative raises
   `ValueError` (it means step 2 was skipped); within that tolerance it is
   clamped, so floating-point scaling landing a hair over 10 N cannot produce
   u = 1.0000000001.

The "never leaves [−1, 1]" guarantee is a property of `limit` followed by
`to_commands`, and is tested there, independent of the pilot.

*Assistant's call (builder check, 2026-11-19).*

### Pilot mixing and saturation — `vehicle/pilot.py`

Inputs: `surge`, `yaw`, `heave`, each in [−1, 1]; outside that range raises
(the "plausible answer to a nonsensical question" rule). There is no sway input,
so R-3's "no strafe" is enforced by the signature as well as by allocation.

Scaling to a wrench, using the vehicle's believed thruster limits:

    X = surge * 2*T_fwd_max   (surge >= 0)   or   surge * 2*T_rev_max  (surge < 0)
    N = yaw   * 2*l*T_rev_max                     (a full turn runs one thruster in reverse)
    Z = heave * T_fwd_max     (heave >= 0)   or   heave * T_rev_max    (heave < 0)

Then `allocate`, `limit`, `to_commands`. If either horizontal thruster exceeds
its limit in its direction, `limit` **scales both horizontal thrusts by the same
factor** so the worse one sits exactly at its limit. That keeps the ratio of X
to N and gives up speed instead. In this model the steady turn radius is
R = sqrt(X·k_r / (N·k_u)), so a fixed X/N fixes the steady radius. It does not
fix the transient path, which still changes with the scale; the technical note
must say so. The vertical thruster is independent.

Worked example (a required test): full surge + full yaw gives X = 20 N,
N = 2.4 N·m, so T_port = 18 N and T_stbd = 2 N. Port exceeds 10 N by 1.8×, so
both scale by 1/1.8: T_port = 10.0 N, T_stbd = 1.111 N, and X/N is unchanged at
8.333 m⁻¹. Steady radius R = sqrt(20·2.4 / (2.4·20)) = 1.000000 m, with or
without the scaling. Clipping per thruster instead gives 10 N and 2 N, so
X = 12 N, N = 1.2 N·m and R = sqrt(12·2.4 / (1.2·20)) = sqrt(1.2) = 1.095445 m:
a different turn.

### The 4-DOF body — `sim/planar.py`

State `[north, east, depth, psi, u, w, r]`. Uncoupled, no added mass (M2), no
sway (v = 0 holds exactly in this model; the real vehicle skids in a turn and
M2 adds that):

    m   * du/dt = X - k_u * u*|u|
    m   * dw/dt = (m - rho*V)*g + Z - k_w * w*|w|
    I_z * dr/dt = N - k_r * r*|r|
    d(north)/dt = u cos psi,  d(east)/dt = u sin psi,  d(depth)/dt = w,  dpsi/dt = r

Thruster commands come in as a function of time, `command(t)`, so the derivative
genuinely depends on t. That exercises RK4's mid-step time arguments, which
AI-assistance-log entry 2 showed the M0 suite could not see.

Placeholder body: `m = 10 kg`, `V = 10.2 L` (0.2 kg positively buoyant, as a
real ROV is trimmed, so it floats up when idle and the vertical thruster can
hold depth), `k_u = 20 N·s²/m²`, `k_w = 25 N·s²/m²`, `I_z = 0.25 kg·m²`,
`k_r = 2.4 N·m·s²`. Gravity and fluid density from the same constants as M0.

## Validation

Each number below is hand-calculated in this spec and must be re-derived
independently in the test, never captured from a run.

1. **Thruster:** deadband gives exactly 0 at |u| <= 0.05; u = 1 gives 10 N,
   u = −1 gives −8 N; u = 1.5 gives 10 N; monotonic over [−1, 1]; the inverse
   round-trips across the range, and 7 N maps to u = 0.05 + 0.95·sqrt(0.7) =
   0.844827.
2. **Configuration matrix:** equals the hand-written B above.
3. **Surge produces no yaw:** pure surge allocates equal port and starboard
   thrust and N = 0 exactly. (The design's headline M1 test.)
4. **Lateral is rejected:** non-zero Y, K, M each raise. Zero Y, K, M accepted.
5. **R-3:** each of forward, backward, turn left, turn right, rise and descend
   produces thrust of the right sign on the right thrusters, and nothing else.
6. **Saturation:** the worked example above, to the stated figures, including
   steady radius 1.000000 m (driven in the body) against 1.095445 m for
   per-thruster clipping; `limit` + `to_commands` never leave [−1, 1] over a
   sweep of wrenches far beyond the limits; `to_commands` raises on an
   un-limited over-limit thrust.
7. **Steady surge:** X = 10 N gives u_t = sqrt(10/20) = 0.707107 m/s, and the
   transient matches M0's tanh closed form with surge parameters.
8. **Steady turn:** X = 10 N, N = 0.6 N·m gives u_t = 0.707107 m/s and
   r_t = sqrt(0.6/2.4) = 0.5 rad/s, so the track converges to a circle of radius
   u_t/r_t = 1.414214 m, turning clockwise seen from above.
9. **Heave:** with no thrust, rises at w_t = −sqrt(1.96133/25) = −0.280095 m/s.
   With Z = +10 N, sinks at w_t = sqrt(8.03867/25) = 0.567051 m/s. With zero
   thrust and M0's body parameters, heave matches `vertical.py` to round-off.
10. **Time-varying command (the reason `command(t)` exists).** Surge with zero
    drag (`k_u = 0`, a legal limiting case like M0's zero drag), both horizontal
    thrusters commanded u(t) = 0.05 + 0.95·(t/2) for t in [0, 2] s, so s = t/2
    and each gives 10·(t/2)² = 2.5 t² N. Then X = 5 t², du/dt = 0.5 t²,
    u(t) = t³/6 and north(t) = t⁴/24. At t = 2 s: u = 1.333333 m/s,
    north = 0.666667 m. RK4 integrates these exactly, so the test holds to
    round-off at any step size, and moving RK4's k2 or k4 time argument breaks
    it. Neither this test nor item 11 may cross the deadband or reverse
    direction mid-run.
11. **Fourth order:** step-size error ratios near 16, run on item 8's spin-up
    from rest at constant thrust (smooth: u and r stay non-negative, so u|u| is
    u²). Not on a stepped command or a deadband crossing, which break the ratio.
12. **Mutation check:** the coordinator breaks the new code deliberately (swap
    port and starboard, drop the shared saturation scale, flip the N sign, break
    the deadband edge, move RK4's k2 and k4 time arguments) and confirms a test fails for each. Survivors are logged.

## Deliverables

- `src/rovsim/geometry.py`, `src/rovsim/sim/thruster.py`,
  `src/rovsim/sim/planar.py`, `src/rovsim/vehicle/allocation.py`,
  `src/rovsim/vehicle/pilot.py`, with tests for each.
- `scripts/plot_m1.py`: a top-down track of a scripted drive (straight, turn,
  straight) and the steady-turn circle against its closed-form radius.
- `docs/notes/m1-thrusters-allocation.md`: the physics, derivations,
  validation and plot, with the placeholder disclaimer up front. It **opens with
  a "Concepts" section**: each idea above in plain language, one or two
  paragraphs, with the physical intuition and why it matters on a real ROV, at
  the level of explaining it to another engineer. Derivations follow for anyone who
  wants them.
- New rows in `docs/decisions.md`; new entries in `docs/ai-assistance-log.md`
  for every failure found, however small.

## Left to the builder

- Function and class names beyond the module names above.
- How `command(t)` is represented (a callable is enough).
- Test file layout.
- Plot styling.

## Open for the owner

None blocking. Real thruster ratings, the real layout and the real hull drag
replace the placeholders whenever the hardware track supplies them; the
modules take them as parameters, so no code changes.
