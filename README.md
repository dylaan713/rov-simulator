# ROV Simulator

A 6-DOF underwater ROV simulator, plus the vehicle and ground-station software developed against it.

The project requirements ask for a ground station: gamepad control over the tether, live video, a live sensor display, 1 Hz data logging, a battery warning, and thruster cut-off within one second of losing the link. None of that can be built against a vehicle that does not exist yet. The simulator is the stand-in vehicle that makes the required software buildable now; when hardware arrives, the simulated vehicle is swapped for the real one.

Full design: [docs/superpowers/specs/2026-11-05-rov-simulator-design.md](docs/superpowers/specs/2026-11-05-rov-simulator-design.md)

## Setup and tests

Requires Python 3.12 or newer.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e ".[dev]"
.venv/bin/pytest
```

## Predicted vs measured

| Quantity | Predicted (closed form) | Simulated | Error | Note |
|---|---|---|---|---|
