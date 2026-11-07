# ROV Software Stack and Vehicle Simulator — Design

Date: 2026-11-05
Status: Draft for approval
Related: `~/Projects/Personal/_parked/rov-build` (hardware track, parked)
Requirements: "ROV Scope" (Google Doc, 2026-09-23) — R-1..R-24, C-1..C-2

## Intent

Build the surface control software and vehicle control code that the ROV Scope
requirements already call for, and build a physics simulator that stands in for
the vehicle so that software can be written and tested before any hardware
exists.

### Why this, and why the simulator

The requirements never ask for a simulator. They ask for a ground station:
gamepad control over the tether, live video, a live sensor display, data
logging at 1 Hz, a battery warning, and a thruster cut-off within one second of
losing the link.

None of that can be built against a vehicle that does not exist. No parts have
been ordered; no CAD exists. The simulator is not a side project and not a
vanity piece — **it is the stand-in vehicle that makes the required software
buildable now.** When hardware does arrive, the simulated vehicle is swapped
for the real one and the ground station already works, already logs, and has
already had its failsafe tested a thousand times in situations too dangerous to
stage in a sinkhole.

That is also why this is worth showing an employer. Developing vehicle software
against a simulated plant, then swapping in hardware, is how real vehicle
programs are run. It is a stronger story than either piece alone.

### Secondary intent

Two things the owner stated directly:

1. **Learn to build software properly.** He is a mechanical engineering
   undergraduate with strong engineering fundamentals — statics, dynamics,
   ODEs, applied statistics, MATLAB, one C++ course — and limited software
   experience. He directs; the assistant explains as it goes. **Teaching is a
   requirement of this project, not a courtesy.**
2. **Find the real limits of AI-assisted development.** The domain was chosen
   because it has a ground truth that cannot be bluffed. A web app rewards
   plausible-looking code because "it renders" feels like success. A dynamics
   model either matches a closed-form solution and a hand calculation or it
   does not.

## Working agreement — changed 2026-11-05

The parked hardware repo instructs the assistant to act as reviewer and mentor
only, and never to write the owner's calcs, CAD or code. **For software work,
that is superseded by the owner's decision on 2026-11-05: co-development, with
the assistant teaching as it goes.** The mentor-only agreement still stands for
the hardware track — calcs, CAD and fabrication decisions remain his.

## Goals

- Deliver the software the requirements call for: R-3, R-9, R-10, R-14, R-15,
  R-16, R-20, and scope items 4 and 5.
- A 6-DOF simulator faithful enough to develop and tune that software against.
- The owner can read, explain and defend every part of the system.
- Documentation produced as the work happens, not retrofitted — see below.

## Non-goals

- No CAD, fabrication, procurement or mechanical analysis. That track is parked
  and stays parked unless the owner restarts it.
- Depth hold and heading hold are **stretch goals**, per the owner's own scope
  decision. They are not requirements and do not gate anything.
- No autonomy, mapping or SLAM.
- Not photorealistic. Visualisation exists to make behaviour legible.
- Not hard real-time. Python with soft real-time timing is accepted.

## Anti-goals

Outcomes that would make this a failure even if the code runs:

- **Code the owner cannot read or explain.** If he cannot defend a subsystem,
  that subsystem has negative value against the employer goal.
- **The simulator becoming the project** while the required ground station
  never gets built. The simulator is a means.
- **The assistant writing everything while the owner watches.** That fails both
  stated goals at once — it neither teaches nor reveals where AI assistance
  actually breaks down.

## Constraints

- Python 3, NumPy, SciPy. Chosen because NumPy maps closely onto the MATLAB the
  owner knows, and a tethered ROV at these rates does not need C++ speed.
- The vehicle it models is **Approach A**, already decided: 3 flooded brushless
  thrusters with PETG propellers and ducts, acrylic tube housing, Raspberry Pi
  4, Ethernet tether to a laptop, depth from an automotive 0-30 psi sender,
  heading from a compass, 10 m rated depth, 30 min runtime, ~$492 of parts.
- **Three thrusters means the vehicle is underactuated.** R-3 asks only for
  forward, backward, turn, rise and descend - there is deliberately no sway. The
  allocation layer must not silently command lateral force the vehicle cannot
  produce.
- C-2: built, tested **and documented** by 2027-01-31.
- Owner's available hours per week: not established. Milestones are sized as
  testable units, not calendar estimates.

## Architecture - the decided parts

Three layers with a hard boundary.

    sim/        physics. Pure functions over state. No I/O, no threads, no clock.
    vehicle/    control software. Sensors in, thruster commands out. Runs on the Pi.
    ground/     surface station. Gamepad, display, logging, video.

**The central decision: the vehicle software must not know whether it is in
simulation or in water.** It talks to one `VehicleIO` interface with two
implementations - `SimIO` backed by the simulator, `HardwareIO` backed by the
Pi. Everything else follows from that seam. It is what lets the ground station
be built now, what makes the failsafe testable without drowning anything, and
what makes hardware bring-up a swap rather than a rewrite.

This mirrors the structure the owner arrived at independently in a previous
project: a pure domain layer kept strictly away from the UI.

### Physics scope

- 6-DOF rigid body, body-fixed SNAME convention.
- Rigid-body mass matrix plus **added mass** - the water accelerated along with
  the hull. For a blunt ROV this can rival the vehicle's own mass and is the
  main reason underwater vehicles handle unlike anything else.
- Hydrodynamic damping, linear and quadratic.
- Restoring forces from weight and buoyancy acting through separated centres of
  gravity and buoyancy. This is what makes the vehicle self-right in roll and
  pitch, and it is why three thrusters is a viable configuration at all.
- Thruster model including deadband and saturation. Pretending thrusters are
  linear is the most common cause of sim-to-real mismatch.
- Tether drag, simplified. Catenary behaviour deferred.
- Fixed-step RK4 integration.
- **Sensor models with noise, bias, quantisation and latency.** Without these,
  software tuned in simulation will not survive contact with hardware, which
  defeats the purpose of building the simulator.

Coefficients are derived from analytic primitives - cylinder and ellipsoid
added mass, flat-plate drag - using Approach A's intended geometry. This is the
only path that makes the closed-form validation in M2 possible, and it is the
same calculation the hardware track needs for `pc1`.

### Control scope

- Allocation across 3 thrusters, handling the underactuated case explicitly.
- Direct pilot mixing for R-3: forward, backward, turn, rise, descend.
- Link-loss failsafe: all thrusters stopped within 1 second (R-15).
- Depth hold and heading hold if time allows. Stretch, not core.

## Requirements traceability

| Milestone | Requirements served |
| --- | --- |
| M0 Foundations | none directly; teaches integration, and produces the drag model the hardware track needs for `pc1` (deferred to M2 — see docs/decisions.md) |
| M1 Thrusters and allocation | R-3 |
| M2 Six degrees of freedom | R-1, R-17 (analytically) |
| M3 Vehicle loop, sensors, failsafe | R-8, R-15, R-19 behaviour |
| M4 Ground station | R-9, R-10, R-14, R-16, R-20, scope 4 and 5 |
| M5 Video path | R-11, R-12, R-13 |
| M6 Hardware swap | gated on the hardware track restarting |

## Milestones

Each ends in something testable that teaches one specific idea.

- **M0 - Foundations.** Repository, test harness, a body moving in one axis
  under weight, buoyancy and drag. *Teaches:* numerical integration, terminal
  velocity, why step size matters. *Validation:* hand-calculated terminal
  velocity matches the simulator within a stated tolerance.
- **M1 - Thrusters and allocation.** Thruster model, 3-thruster configuration
  matrix, pilot mixing. *Teaches:* actuation, underactuation, why this vehicle
  cannot strafe. *Validation:* commanded surge produces no yaw; a lateral
  command is rejected rather than silently approximated; saturation is correct.
- **M2 - Six degrees of freedom.** Full rigid body with added mass, damping and
  restoring forces. *Validation:* energy conservation with damping disabled;
  decay rates match closed-form predictions; the vehicle self-rights from an
  initial roll.
- **M3 - Vehicle loop and failsafe.** The control loop, sensor models, telemetry
  out, and the 1-second link-loss cut-off. *Validation:* failsafe fires under
  injected link loss, dropped packets and delayed packets.
- **M4 - Ground station.** Gamepad input, live telemetry display, 1 Hz logging
  to a common format, tether-length entry, battery warning, mission replay.
  **This is the milestone that satisfies most of the actual requirements.**
- **M5 - Video path.** Capture, stream, record, and time-sync to the sensor log.
  Developed against a synthetic source until a camera exists.
- **M6 - Hardware swap.** `HardwareIO`, Pi bring-up, bench tests. Gated on the
  hardware track restarting.

## Documentation - a first-class deliverable

C-2 requires the project be documented, and the owner has asked specifically
that this produce something tangible to show. Documentation is therefore part
of each milestone's definition of done, not a phase at the end.

- Each milestone produces a short technical note: the physics being modelled,
  the derivation, the validation method, the plot, and the result. These
  assemble into the design report C-2 requires.
- A **decision log** recording why each design choice was made and what was
  rejected.
- An **AI-assistance log**: where the assistant was wrong, what it missed, what
  had to be corrected. This is the direct evidence for the owner's stated goal
  of finding the limits of AI-assisted development, and no one else's portfolio
  will have one.
- Public GitHub repository, so the work is linkable from an application. The
  repo README carries the predicted-versus-measured results as they arrive.

## Method

- Tests before implementation. Here tests are not hygiene - they are how a
  physical claim is proven. A test asserting a released body rises at the
  hand-calculated terminal velocity teaches more than any explanation.
- Each milestone is reviewed with the owner before the next begins.
- The assistant explains the physics and the reasoning, not only the code. Code
  the owner cannot read is a defect, not a style preference.

## Left to the builder

- Test framework (pytest assumed).
- Visualisation library for the ground station.
- Telemetry transport and message schema.
- Repository layout within the three-layer boundary.
- Whether the simulator runs in-process during early milestones or across the
  real transport from the start.

## Open questions

- Should the GitHub repository be public from M0? Assumed yes, for the
  employer goal. Reversible only in one direction.
- Hours per week, to set a realistic pace.
