"""Optional local research dashboard: streamlit run app.py."""
from pathlib import Path
import json
import tempfile

import numpy as np
import pandas as pd
import streamlit as st

from pcpt.science import Scenario, calculate, flux, pressure
from pcpt.field import evaluate, DataQualityError

st.set_page_config(page_title='PCPT v3.0 Research Lab', page_icon='🌍', layout='wide')
st.title('PCPT v3.0 | Mathematical & Experimental Validation')
st.warning('Research prototype. Synthetic model outputs are NOT empirical confirmation of GRM, thermal or climate claims.')
mode = st.radio('Workspace', ['Synthetic GRM benchmark', 'Field CSV screening'], horizontal=True)
if mode == 'Synthetic GRM benchmark':
    c1, c2, c3 = st.columns(3)
    with c1: amplitude = st.number_input('Pressure amplitude [Pa]', min_value=0., max_value=100000., value=10000., step=500.)
    with c2: beta = st.number_input('Permeability sensitivity [1/Pa]', min_value=-1e-4, max_value=1e-4, value=2e-5, format='%.7f')
    with c3: length = st.number_input('Column length [m]', min_value=.1, max_value=1000., value=10.)
    scenario = Scenario(pressure_amplitude_Pa=amplitude, beta_Pa_inv=beta, length_m=length)
    result = calculate(scenario)
    st.json({k: v for k,v in result.items() if k != 'parameters'})
    t = np.linspace(0, scenario.period_s, 300)
    st.line_chart(pd.DataFrame({'Linear null [m/s]': flux(t, Scenario(beta_Pa_inv=0, pressure_amplitude_Pa=amplitude, length_m=length)),
                                'GRM candidate [m/s]': flux(t, scenario)}, index=t/scenario.period_s), x_label='Cycle fraction')
else:
    st.markdown('Upload your own calibrated, time-synchronised measurements. Missing fields, interventions, bad timestamps, and insufficient cycles are rejected. Screening never produces "theory proven" status.')
    st.code(','.join(['timestamp_utc','site_id','sensor_id','pressure_Pa','pressure_gradient_Pa_m', 'darcy_flux_m_s', 'control_flux_m_s', 'flux_uncertainty_m_s', 'operation_flag', 'source_id']))
    upload = st.file_uploader('Measurements CSV', type=['csv'])
    period = st.number_input('Pre-registered forcing period [seconds]', min_value=60., value=44712.)
    if upload is not None and st.button('Run conservative screening'):
        with tempfile.TemporaryDirectory() as root:
            p = Path(root)/'observations.csv'
            p.write_bytes(upload.getvalue())
            try:
                output = evaluate(p, period_s=period)
                st.json(output)
                st.download_button('Download audit JSON', json.dumps(output, indent=2), 'pcpt_screen.json')
            except (DataQualityError, ValueError, OSError) as exc:
                st.error(f'Input rejected: {exc}')
