"""Conservative field-data analysis with explicit provenance and held-out evaluation.

No observational detection in this package is claimed as a natural validation.
Only a pre-registered, independently replicated experiment can do that.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import t as student_t

REQUIRED = ('timestamp_utc', 'site_id', 'sensor_id', 'pressure_Pa',
            'pressure_gradient_Pa_m', 'darcy_flux_m_s', 'control_flux_m_s',
            'flux_uncertainty_m_s', 'operation_flag', 'source_id')


class DataQualityError(ValueError):
    """Rejected measurement dataset (data-quality failure)."""


def load_observations(path: str | Path) -> tuple[pd.DataFrame, dict]:
    source = Path(path)
    raw = source.read_bytes()
    frame = pd.read_csv(source, dtype={'site_id': str, 'sensor_id': str,
                                       'source_id': str, 'operation_flag': str})
    absent = sorted(set(REQUIRED).difference(frame.columns))
    if absent:
        raise DataQualityError(f'Missing required columns: {absent}')
    if not len(frame):
        raise DataQualityError('Empty observations')
    # A single site/sensor is required to avoid silently averaging independent instruments.
    for column in ('site_id', 'sensor_id', 'source_id'):
        if frame[column].isna().any() or frame[column].nunique() != 1:
            raise DataQualityError(f'{column} requires exactly one non-null identifier')
    time = pd.to_datetime(frame['timestamp_utc'], utc=True, errors='coerce')
    if time.isna().any() or time.duplicated().any() or not time.is_monotonic_increasing:
        raise DataQualityError('Times must be valid, unique, UTC, and strictly ordered')
    frame['time'] = time
    for key in ('pressure_Pa', 'pressure_gradient_Pa_m', 'darcy_flux_m_s',
                'control_flux_m_s', 'flux_uncertainty_m_s'):
        frame[key] = pd.to_numeric(frame[key], errors='coerce')
        if not np.isfinite(frame[key].to_numpy(dtype=float)).all():
            raise DataQualityError(f'Non-finite or missing {key}')
    if (frame.flux_uncertainty_m_s <= 0).any():
        raise DataQualityError('flux_uncertainty_m_s must be >0')
    # Explicit operational flag permits only known unperturbed readings.
    flags = frame.operation_flag.fillna('').astype(str).str.strip().str.lower()
    if not flags.isin(['0', '1']).all():
        raise DataQualityError('operation_flag must be 0/1; unknown operations not silently imputed')
    if (flags == '1').any():
        raise DataQualityError('Dataset contains operational interventions; segment them BEFORE confirmatory evaluation')
    elapsed = frame.time.diff().dropna().dt.total_seconds().to_numpy()
    if len(elapsed) == 0 or (elapsed <= 0).any():
        raise DataQualityError('At least two valid chronological observations are required')
    return frame, {'sha256': hashlib.sha256(raw).hexdigest(), 'path': source.name,
                   'n_rows': len(frame), 'site_id': str(frame.site_id.iloc[0]),
                   'sensor_id': str(frame.sensor_id.iloc[0]),
                   'source_id': str(frame.source_id.iloc[0]),
                   'interval_median_s': float(np.median(elapsed)),
                   'synthetic_indicator': 'SYNTHETIC' in str(frame.source_id.iloc[0]).upper()}


def _predict(theta, gradient, pressure):
    # Restrict to a physically simplified Darcy-proxy form, NOT a full poroelastic inversion.
    scale, beta = theta
    return scale * gradient * np.exp(np.clip(beta * pressure, -8, 8))


def evaluate(path: str | Path, *, period_s: float = 44712.0,
             holdout_fraction: float = 0.4, min_cycles: int = 12) -> dict:
    if not (period_s > 0 and 0.2 <= holdout_fraction <= 0.6 and min_cycles >= 6):
        raise ValueError('Invalid fixed protocol: positive period, holdout 0.2..0.6, min_cycles>=6')
    df, provenance = load_observations(path)
    t = (df.time - df.time.iloc[0]).dt.total_seconds().to_numpy()
    cycle = np.floor(t / period_s).astype(int)
    n_cycles = int(cycle.max())  # Last partial cycle discarded
    if n_cycles < min_cycles:
        raise DataQualityError(f'Insufficient complete cycles: {n_cycles} < {min_cycles}')
    # Cover each cycle evenly and exclude large data gaps.
    cadence = provenance['interval_median_s']
    if cadence > period_s / 8:
        raise DataQualityError('Sampling too sparse: fewer than 8 samples per cycle')
    d = df.loc[cycle < n_cycles].copy()
    d['cycle'] = cycle[cycle < n_cycles]
    count = d.groupby('cycle').size()
    expected = period_s / cadence
    if len(count) != n_cycles or (count < 0.8 * expected).any():
        raise DataQualityError('Missing cycles or less than 80% cycle sampling coverage')
    n_train_cycles = int(n_cycles * (1 - holdout_fraction))
    if n_train_cycles < 4 or n_cycles - n_train_cycles < 4:
        raise DataQualityError('Insufficient full cycles in train and holdout')
    train = d[d.cycle < n_train_cycles]
    test = d[d.cycle >= n_train_cycles]
    gtr = train.pressure_gradient_Pa_m.to_numpy(float)
    ptr = train.pressure_Pa.to_numpy(float)
    ytr = train.darcy_flux_m_s.to_numpy(float)
    gt = test.pressure_gradient_Pa_m.to_numpy(float)
    pt = test.pressure_Pa.to_numpy(float)
    yt = test.darcy_flux_m_s.to_numpy(float)
    # Predictors are fit exclusively on train. Nonlinear parameter bounded to avoid overflow.
    def residual(theta):
        return (_predict(theta, gtr, ptr) - ytr) / train.flux_uncertainty_m_s.to_numpy(float)
    # pressure gradient and pressure must vary sufficiently for parameter identification.
    if np.std(gtr) <= 1e-10 or np.std(ptr) <= 1e-6:
        raise DataQualityError('Insufficient forcing variability for identification')
    linear = float(np.dot(gtr, ytr) / max(np.dot(gtr, gtr), 1e-200))
    lin_pred = linear * gt
    try:
        fit = least_squares(residual, x0=[max(linear, 1e-13), 0],
                            bounds=([0., -1e-3], [1e-3, 1e-3]),
                            x_scale=[1e-10, 1e-5], max_nfev=2000)
    except ValueError as exc:
        raise DataQualityError(f'Model fit failed: {exc}') from exc
    if not fit.success:
        raise DataQualityError('Nonlinear model fit did not converge')
    nonlin_pred = _predict(fit.x, gt, pt)
    rmse_linear = float(np.sqrt(np.mean((yt - lin_pred)**2)))
    rmse_grm = float(np.sqrt(np.mean((yt - nonlin_pred)**2)))
    # Cycle-level, not sample-level inference to avoid inflated n from autocorrelation.
    cycle_means = test.groupby('cycle')[['darcy_flux_m_s', 'control_flux_m_s']].mean()
    x = cycle_means.darcy_flux_m_s.to_numpy(float)
    controls = cycle_means.control_flux_m_s.to_numpy(float)
    if len(x) < 4 or np.std(x, ddof=1) == 0:
        raise DataQualityError('Insufficient between-cycle variation for honest uncertainty estimate')
    def robust_ci(v):
        n = len(v); sd = np.std(v, ddof=1)
        lag = float(np.corrcoef(v[1:], v[:-1])[0,1]) if n > 3 else 0.
        lag = float(np.clip(np.nan_to_num(lag), -0.9, 0.9))
        # Conservative effective n for positive serial correlation, never increase n for negative rho.
        effective_n = max(2.01, n * (1-max(0.,lag))/(1+max(0.,lag)))
        radius = float(student_t.ppf(0.975, max(effective_n-1, 1)) * sd / np.sqrt(effective_n))
        return {'mean_m_s': float(np.mean(v)), 'ci95_m_s': [float(np.mean(v)-radius), float(np.mean(v)+radius)],
                'cycle_lag1_rho': lag, 'effective_n_approx': effective_n}
    target = robust_ci(x); ctrl = robust_ci(controls) if np.std(controls, ddof=1)>0 else None
    observed = target['ci95_m_s'][0] > 0 or target['ci95_m_s'][1] < 0
    control_ok = bool(ctrl is not None and ctrl['ci95_m_s'][0] <= 0 <= ctrl['ci95_m_s'][1])
    fit_gain = float((rmse_linear - rmse_grm) / max(rmse_linear, 1e-30))
    return {'research_status': 'SYNTHETIC_INPUT' if provenance['synthetic_indicator'] else 'OBSERVATIONAL_ASSOCIATION_NOT_CAUSAL_PROOF',
            'provenance': provenance, 'protocol': {'period_s': period_s, 'n_complete_cycles': n_cycles,
                'train_cycles': n_train_cycles, 'holdout_cycles': n_cycles-n_train_cycles},
            'linear_control': {'holdout_rmse_m_s': rmse_linear, 'train_coefficient': linear},
            'candidate_grm': {'holdout_rmse_m_s': rmse_grm, 'train_scale': float(fit.x[0]),
                              'train_beta_Pa_inv': float(fit.x[1]), 'relative_rmse_gain': fit_gain},
            'cycle_mean_holdout': target, 'negative_control_mean_holdout': ctrl,
            'screening_flags': {'nonzero_cycle_mean_ci_excludes_zero': bool(observed),
                'negative_control_ci_includes_zero': control_ok,
                'incremental_holdout_gain_at_least_5pct': bool(fit_gain>=0.05)},
            'decision': 'SCREENING_ONLY__INDEPENDENT_REPLICATION_REQUIRED',
            'limitations': ['Pre-specified sampling and holdout required', 'Operations excluded, not modelled',
                'Control site comparability must be documented', 'Proxy fit is not mechanistic attribution',
                'Independent field replication and calibrated sensors remain mandatory']}


def write_result(result: dict, path: str | Path):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
