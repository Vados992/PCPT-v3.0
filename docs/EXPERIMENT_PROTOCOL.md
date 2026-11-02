# PCPT v3.0 preregistered verification protocol (draft)

**Hypotheses:** H-04 (non-zero cycle-mean directed fluid transport) and H-05 (additional heat transport beyond the conduction-only null). The protocol is a design: no new real-world data have been collected.

## Phase A — Pre-registration

- Choose and document site, geometry, lithology, effective permeability, fluid viscosity, stress, pressure gradient, porosity, storage, temperature, and sensor calibration certificates **before** opening the confirmation window.
- Declare cycle periods, sampling cadence, source of ephemerides or pressure forcing, lag windows, exclusion criteria, operational intervention records, test statistics, primary outcome, and a practical minimum detectable effect in SI units.
- Freeze code, parameter estimation, numerical model, statistical thresholds and versioned analysis plan; hash inputs and configuration.
- Include a simultaneously measured **negative control** (sensor/instrument/site selected based on expected insensitivity) and a model of weather and industrial operations. A negative control alone does not prove absence of confounding.

## Phase B — Laboratory apparatus

1. Horizontal column with well-characterized porous medium and **calibrated** pressure/flow/temperature instrumentation.
2. Compare (a) linear-material calibration, (b) candidate stress-sensitive medium, and (c) sealed/no-throughflow control; include positive controls with imposed DC flow.
3. Prescribe sinusoidal pressure perturbation with measured boundary forcing, independent closed-loop clock, logged pumps/valves and controlled thermal boundary conditions.
4. Confirm **full-period integrated** net mass flux by direct mass/volume measurement or conservative tracers. Distinguish displacement, phase lag, and dispersion from directed transport.
5. Repeat with forcing amplitude, sign, period and sample orientation changed. Reversal and amplitude dependence should be predicted before the experiment.
6. Blind replicate on an independently operated apparatus. A result not exceeding combined sensor plus model uncertainty fails H-04.

## Phase C — Field observatory

- Minimum instrumentation: pore pressure, in-situ temperature, fluid flux/tracer, barometric pressure, geodetic strain/GNSS, well logs and operations; magnetic/oceanic inputs only when applicable.
- Prefer paired reference / production / depleted zones, **not** unpaired post-hoc comparisons.
- At least 12 complete candidate periods for the included screening example. Actual study duration must be justified by seasonal, operational and harmonic confounding.
- Hold out final 40% of complete cycles with preregistered parameters; handle cycle autocorrelation (screen provides a conservative AR(1) approximation, not a final inferential guarantee).
- Evaluate whether the nonlinear candidate model improves genuinely independent prediction against suitable physically informed nulls (linear Darcy, known pumping/pressure effects, weather, and tide/wellbore transfer).

## Gate decisions

**H-04 advancement requires**: observed directional mean transport, non-zero within a prespecified uncertainty budget, actual fluid-mass closure, confounder controls, physically credible energy source and blind independent replication. **H-05 advancement requires**: measured extra heat power with heat-storage/advection/conduction accounting, replicated comparison against thermal null. **G5 requires**: area-weighted power integration with conservative upper and lower bounds and a demonstrated link to environmental energy budget. Otherwise issue an upper bound or reject the stronger claim.

## Important methodological limitation

The program's CSV screening is a **diagnostic aid**, not a qualified scientific hypothesis test. Pre-registration, physical instrumentation, robust autocorrelation-aware modeling, intervention analysis, multiple comparison control and independent review remain external responsibilities.
