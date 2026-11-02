# PCPT v3.0 — Mathematical & Experimental Validation

**Theory author / system architect:** Vadym Tsinderhoz (2026)  
**Status:** Computational benchmark available; **H-04/GRM and H-05 thermal effect not empirically validated**.  
**Copyright:** © 2026 Vadym Tsinderhoz. All rights reserved. See [LICENSE](LICENSE).

A reproducible research system for the **Planetary Coupled Pulse Theory (PCPT)**, focused on the proposed **Geological Rectification Mechanism (GRM)**. It implements a documented, open-boundary nonlinear Darcy-flow candidate and a zero-mean linear control, a finite-volume mass-conservation check, an illustrative thermal-advection screen, an optional browser-based research dashboard, and a strictly non-confirmatory field-data screening pipeline.

> **Scientific boundary:** The calculations are synthetic demonstrations. This repository does not prove a global planetary pulse, natural geological rectification, climate-significant petroleum cooling, or an extraction-driven loss of planetary cooling. The chosen pressure boundary is **assumed**, not inferred directly from lunar forcing. A non-zero analytical mean is not field evidence.

## Install and run

Requires Python **3.10+**. Locally:

```bash
git clone https://github.com/Vados992/PCPT-v3.0.git
cd PCPT-v3.0
python -m venv .venv
# macOS/Linux: source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,ui]'
python -m pytest -q
pcpt model --scenario configs/baseline.json --out outputs/model_summary.json
pcpt benchmark  # generates figures/ and results/; computationally heavier
streamlit run app.py
```

**Minimal CLI only:** `pip install -e .` (UI is optional). For Docker, `docker compose up --build` runs the local Streamlit dashboard at `http://localhost:8501`. `make check` runs tests and compiles Python; `make benchmark` executes the time-dependent finite-volume reference model.

The full PCPT v3.0 scientific manuscript and four precomputed synthetic figures are included under `assets/`. The local and GitHub versions of these assets have been uploaded as identical binary content.

## System layers

| Layer | Path | Status |
| --- | --- | --- |
| Baseline analytic nonlinear Darcy model and thermal screen | `pcpt/science.py` | Runnable; synthetic |
| Time-dependent conservative finite-volume check | `grm_benchmark.py` | Runnable; synthetic |
| CSV ingestion, quality gates, source hash and controlled holdout | `pcpt/field.py` | Runnable; screening, **not causality** |
| Command-line interface | `pcpt/cli.py` | Runnable |
| Local dashboard | `app.py` | Runnable with `ui` extras |
| Tests and CI definition | `tests/`, `.github/workflows/ci.yml` | Local tests; GitHub Actions account/billing may block runs |
| Scientific assumptions and falsification gates | `docs/` | Protocol / preregistration requirements |

## Reproduce the known benchmark

For parameters `k0=1e-13 m²`, `mu=1e-3 Pa s`, `L=10 m`, `beta=2e-5 Pa^-1`, and harmonic left pressure `A=10000 Pa` with an open right boundary, the exact cycle mean is

`q_bar = k0 * [I0(beta*A)-1] / (mu*L*beta)`.

The expected **synthetic** mean is approximately `5.0125e-9 m/s`. For `beta=0`, mean flux is zero. This does **not** imply that real tidal stress produces the assumed pressure boundary or permeability asymmetry. A **sealed** reservoir has no time-averaged boundary throughflow under periodic steady conditions.

## Field data (no fake observations)

Input schema: `data/templates/observations.csv`. Upload a **single, calibrated site/sensor series** with strictly increasing UTC timestamps, known uncertainty, independent negative-control flux and explicit operation flags. Measurement provenance is hashed (SHA-256). Runs with missing controls, recorded interventions, insufficient cycles, or inadequate coverage fail closed.

```bash
pcpt field YOUR_OBSERVATIONS.csv --period-s 44712 --out outputs/field_screen.json
```

The software compares a linear proxy versus an exponential-mobility proxy with a chronological held-out block, computes full-cycle means with an autocorrelation-aware approximate interval, and returns **`SCREENING_ONLY__INDEPENDENT_REPLICATION_REQUIRED`** even when a signal is statistically detectable. No live sensor integration or pre-registered external study has been performed.

## Falsification gates

1. **G0 — Energy/pressure scale**: demonstrate the proposed forcing can exceed instrumental resolution and an appropriate noise model.
2. **G1 — Known physics**: recover established tidal/ocean and poroelastic responses on an independent historical dataset.
3. **G2 — Pilot**: co-located pressure, displacement and temperature, with exact operational metadata.
4. **G3 — H-04**: demonstrate directional cycle-mean transport, conserved mass, suitable negative controls and blind replication.
5. **G4 — H-05**: independently quantify additional heat transport relative to a conduction-only null.
6. **G5 — Regional/global scaling**: require energy closure and conservative geographic bounds; terminate negligible effects.

Full experiment protocol: [docs/EXPERIMENT_PROTOCOL.md](docs/EXPERIMENT_PROTOCOL.md). Scientific claims: [docs/CLAIMS_AND_LIMITS.md](docs/CLAIMS_AND_LIMITS.md).

## Security and use

The dashboard binds to localhost by default; do **not** expose it publicly or upload confidential borehole/operational records without a suitable data governance plan. No credentials are required for the synthetic benchmark. All rights reserved: public readability is **not** permission for downstream reuse.
