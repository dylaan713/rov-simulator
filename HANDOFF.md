# Handoff
Last updated 2026-11-07.

## State Update
- **Done:** M0 complete (2026-11-07, branch `m0-foundations`). Hand-written RK4 and 1-DOF vertical model; simulated terminal velocity 0.8857381099 m/s vs closed form 0.8857381103 m/s, relative error 4.8e-10; step-size study confirms 4th order. 39 tests pass under `pytest -W error`. Note: `docs/notes/m0-vertical-motion.md`; logs: `docs/decisions.md`, `docs/ai-assistance-log.md`.
- **Next step:** M1 — thruster model, 3-thruster configuration matrix and pilot mixing.

## Context a new session needs
- Hardware track is **parked** at `~/Projects/Personal/_parked/rov-build` (no progress since 2026-09-24; PDR missed 2026-10-18). This project is software only.
- **Working agreement for software: co-development, assistant teaches as it goes** (Dylan's decision, 2026-11-05). This supersedes the mentor-only rule in the parked repo's CLAUDE.md, which still applies to hardware calcs/CAD.
- Public GitHub repo: approved; the push happens after M0 (not yet done).
- Still open: hours per week (unestablished, so no dates are planned).
- M0 uses a placeholder body, not Approach A's real geometry; real drag coefficients wait for M2 (see `docs/decisions.md`).
