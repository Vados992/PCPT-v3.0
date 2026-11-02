"""Synthetic analytical and numerical checks for PCPT v3.0."""
import unittest
import numpy as np
import grm_benchmark as b
class TestGRMBenchmark(unittest.TestCase):
    def test_linear_zero_mean(self):
        self.assertEqual(b.mean_exact(0.0),0.0)
        t=np.arange(512)*b.P['period_s']/512
        self.assertLess(abs(float(np.mean(b.flux_steady(t,0)))),1e-18)
    def test_sign_reversal(self):
        self.assertAlmostEqual(b.mean_exact(2e-5),-b.mean_exact(-2e-5),delta=1e-17)
    def test_quadrature(self):
        t=np.arange(512)*b.P['period_s']/512
        self.assertAlmostEqual(float(np.mean(b.flux_steady(t,2e-5))),b.mean_exact(2e-5),delta=1e-15)
    def test_thermal_energy(self):
        h=b.thermal_solution(b.mean_exact(2e-5))
        self.assertAlmostEqual(h['Pe'],0.10526279184902743,delta=1e-10)
        self.assertAlmostEqual(h['incremental_total_heat_W_m2'],0.003102492707038093,delta=1e-10)
    def test_fv_one_step_mass_conservation(self):
        p=np.zeros(32)
        new,q=b.nonlinear_timestep(p,350.0,350.0,b.P['beta_Pa_inv'])
        lhs=b.P['C_Pa_inv']*b.P['L_m']*float(np.mean(new-p))
        rhs=(q[0]-q[-1])*350.
        self.assertLess(abs(lhs-rhs),1e-9)
if __name__=='__main__':unittest.main()
