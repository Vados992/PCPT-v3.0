# Computational reproducibility

## Environment and deterministic workflow

- Python >=3.10; dependencies pinned by bounded ranges in `pyproject.toml` (a production archive should additionally freeze a fully resolved lockfile).
- `python -m pytest -q` exercises the null, exact result, thermal screen, field QA and holdout estimates.
- `pcpt model --scenario configs/baseline.json` reports synthetic analytical and quadrature results in JSON.
- `pcpt benchmark` executes the legacy SciPy finite-volume pressure diffusion benchmark and generates the graphs and CSVs in `figures/` and `results/`.
- Every field CSV is hashed SHA-256 before analysis. Audit JSON contains the source identity and dataset quality checks.

## Expected synthetic benchmark

For `beta=2e-5 Pa^-1`, `A=10 kPa`, and scenario in `configs/baseline.json`, the exact cycle-mean flux should be `~5.0125e-9 m/s`, the linear `beta=0` model should yield zero cycle-mean flow, and the thermal Peclet estimate should be `~0.1053`.

## Known limitations

- The finite volume solver exercises the porous-flow **pressure** equation; it is not a nonlinear multi-phase / temperature solver.
- `field.py` uses a reduced-order proxy and a simple approximate autocorrelation adjustment; it can only screen candidate signals. Its output is **never** an acceptance certificate.
- No provided field CSV contains measurements; templates are empty by design.
- GitHub Actions may be unable to run due to GitHub account billing/spending-limit restrictions. Local `pytest` results do not imply hosted CI passed.
