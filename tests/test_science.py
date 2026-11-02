from dataclasses import replace
import numpy as np
import pytest
from pcpt.science import Scenario, calculate, mean_flux, flux, energy_screen


def test_exact_mean_and_quadrature():
    s=Scenario(); r=calculate(s)
    assert r['research_status']=='SYNTHETIC_MODEL_ONLY'
    assert mean_flux(s)==pytest.approx(5.0125e-9,rel=1e-4)
    assert r['quadrature_abs_error_m_s']<1e-15


def test_linear_and_zero_amplitude():
    s=Scenario(beta_Pa_inv=0);t=np.arange(512)*s.period_s/512
    assert mean_flux(s)==0
    assert abs(np.mean(flux(t,s)))<1e-18
    assert mean_flux(Scenario(pressure_amplitude_Pa=0))==0


def test_sign_flip_and_thermal_screen():
    s=Scenario(); assert mean_flux(s)==pytest.approx(-mean_flux(replace(s,beta_Pa_inv=-s.beta_Pa_inv)),abs=1e-19)
    h=energy_screen(s, mean_flux(s))
    assert h['peclet']==pytest.approx(.1052628,rel=1e-4)
    assert h['increment_vs_conduction_W_m2']>0


def test_invalid_parameter_rejected():
    with pytest.raises(ValueError): Scenario(length_m=0)
    with pytest.raises(ValueError): Scenario(pressure_amplitude_Pa=-1)
    with pytest.raises(ValueError): Scenario(beta_Pa_inv=100)
