from datetime import datetime, timezone, timedelta
import csv
import numpy as np
import pytest
from pcpt.field import load_observations, evaluate, DataQualityError, REQUIRED


def make_field(tmp_path, cycles=16, samples=12, with_ops=False):
    output=tmp_path/'field.csv';period=1200.;dt=period/samples
    rng=np.random.default_rng(12)
    with output.open('w',newline='') as f:
        w=csv.writer(f); w.writerow(REQUIRED)
        for i in range(cycles*samples+1):
            pressure=2000*np.sin(2*np.pi*i/samples)
            gradient=pressure/10
            q=1e-9*gradient*np.exp(2e-4*pressure)+rng.normal(0,2e-8)
            control=rng.normal(0,2e-8)
            w.writerow([(datetime(2026,1,1,tzinfo=timezone.utc)+timedelta(seconds=i*dt)).isoformat(),
                        'SITE-1','SENSOR-1',pressure,gradient,q,control,3e-8,
                        '1' if with_ops and i==5 else '0','SYNTHETIC-test-fixture'])
    return output


def test_valid_synthetic_returns_screen_only(tmp_path):
    p=make_field(tmp_path)
    r=evaluate(p,period_s=1200)
    assert r['research_status']=='SYNTHETIC_INPUT'
    assert r['protocol']['n_complete_cycles']==16
    assert r['decision']=='SCREENING_ONLY__INDEPENDENT_REPLICATION_REQUIRED'
    assert len(r['provenance']['sha256'])==64


def test_operation_rejected(tmp_path):
    p=make_field(tmp_path,with_ops=True)
    with pytest.raises(DataQualityError,match='interventions'):
        evaluate(p, period_s=1200)


def test_bad_or_sparse_data_rejected(tmp_path):
    p=make_field(tmp_path,cycles=4)
    with pytest.raises(DataQualityError,match='Insufficient'):
        evaluate(p,period_s=1200)
    p.write_text('timestamp_utc,darcy_flux_m_s\n2026-01-01T00:00:00Z,1\n')
    with pytest.raises(DataQualityError,match='Missing required'):
        load_observations(p)
