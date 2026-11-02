"""Independent analytical checks of the GRM pressure-boundary benchmark.

Model assumptions: horizontal 1D saturated single phase, imposed left pressure,
constant right pressure, k(p)=k0 exp(beta*p) in a bounded test regime.
Never interpret these calculations as geological observations.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import math

import numpy as np
from scipy.special import i0


@dataclass(frozen=True)
class Scenario:
    length_m: float = 10.0
    k0_m2: float = 1e-13
    viscosity_Pa_s: float = 1e-3
    storage_Pa_inv: float = 2e-8
    beta_Pa_inv: float = 2e-5
    pressure_amplitude_Pa: float = 1e4
    period_s: float = 12.42 * 3600
    thermal_conductivity_W_m_K: float = 2.0
    volumetric_heat_capacity_J_m3_K: float = 4.2e6
    temperature_difference_K: float = 0.3

    def __post_init__(self):
        positive = ['length_m', 'k0_m2', 'viscosity_Pa_s', 'storage_Pa_inv',
                    'period_s', 'thermal_conductivity_W_m_K', 'volumetric_heat_capacity_J_m3_K']
        if any(not math.isfinite(getattr(self, n)) or getattr(self, n) <= 0 for n in positive):
            raise ValueError('Positive finite values required for length, permeability, viscosity, storage, period, thermal properties')
        if not math.isfinite(self.beta_Pa_inv) or not math.isfinite(self.pressure_amplitude_Pa):
            raise ValueError('Non-finite input parameter')
        if self.pressure_amplitude_Pa < 0 or abs(self.beta_Pa_inv * self.pressure_amplitude_Pa) > 8:
            raise ValueError('Amplitude must be >=0 and |beta * amplitude| <= 8 in this screening model')

    @classmethod
    def read(cls, source: str | Path):
        obj = json.loads(Path(source).read_text(encoding='utf-8'))
        allowed = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in obj.items() if k in allowed})


def mean_flux(s: Scenario) -> float:
    """Exact <q> = (k0/(mu L beta)) [I0(beta A)-1]."""
    beta, amplitude = s.beta_Pa_inv, s.pressure_amplitude_Pa
    if beta == 0 or amplitude == 0:
        return 0.0
    return float(s.k0_m2 / (s.viscosity_Pa_s * s.length_m * beta) * (i0(beta * amplitude) - 1))


def pressure(t_s, s: Scenario):
    return s.pressure_amplitude_Pa * np.sin(2 * np.pi * np.asarray(t_s) / s.period_s)


def flux(t_s, s: Scenario):
    p = pressure(t_s, s)
    if s.beta_Pa_inv == 0:
        return s.k0_m2 / (s.viscosity_Pa_s * s.length_m) * p
    return (s.k0_m2 / (s.viscosity_Pa_s * s.length_m) *
            np.expm1(s.beta_Pa_inv * p) / s.beta_Pa_inv)


def energy_screen(s: Scenario, q: float) -> dict:
    """Steady advection-diffusion at fixed boundary temperatures; NOT coupled transients."""
    pe = s.volumetric_heat_capacity_J_m3_K * q * s.length_m / s.thermal_conductivity_W_m_K
    baseline = -s.thermal_conductivity_W_m_K * s.temperature_difference_K / s.length_m
    ratio = pe / math.expm1(pe) if abs(pe) > 1e-8 else 1 - pe / 2 + pe * pe / 12
    advective_result = baseline * ratio
    return {'peclet': pe, 'conduction_only_W_m2': baseline,
            'heat_flux_at_inlet_W_m2': advective_result,
            'increment_vs_conduction_W_m2': advective_result - baseline,
            'note': 'boundary-referenced synthetic steady thermal screen, not measured heat flow'}


def calculate(s: Scenario) -> dict:
    value = mean_flux(s)
    t = np.linspace(0, s.period_s, 4096, endpoint=False)
    quadrature = float(np.mean(flux(t, s)))
    return {'research_status': 'SYNTHETIC_MODEL_ONLY', 'parameters': asdict(s),
            'analytic_mean_flux_m_s': value, 'numerical_mean_flux_m_s': quadrature,
            'quadrature_abs_error_m_s': abs(value - quadrature),
            'linear_null_mean_flux_m_s': 0.0,
            'thermal_screen': energy_screen(s, value),
            'limitations': ['Not empirical validation', 'Uncalibrated constitutive relation',
                            'Externally imposed oscillatory pressure and open boundaries']}
