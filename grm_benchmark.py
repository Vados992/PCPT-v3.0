#!/usr/bin/env python3
"""PCPT v3.0: independently rerunnable synthetic GRM verification benchmark.

Research-status notice: all time series are SYNTHETIC model solutions, not field data.
MIT/BSD-free? Copyright (c) 2026 Vadym Tsinderhoz. All rights reserved.
Supplementary independently implemented numerical benchmark generated for PCPT v3.0.
"""
import csv
import json
import math
import platform
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import i0
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'; FIG=ROOT/'figures'
RESULTS.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)

P={
    'L_m':10.,'k0_m2':1e-13,'mu_Pa_s':1e-3,
    'C_Pa_inv':2e-8,'beta_Pa_inv':2e-5,
    'pressure_amplitude_Pa':10000.,'period_s':12.42*3600,
    'water_heat_capacity_vol_J_m3_K':4.2e6,
    'thermal_conductivity_W_m_K':2.,
    'temperature_difference_K':0.30,
    'reference_planetary_heat_flow_W_m2':0.092,
}


def primitive(p, beta, k0=P['k0_m2'], mu=P['mu_Pa_s']):
    """F'(p)=k(p)/mu. Kirchhoff transform, stable at beta=0."""
    if abs(beta)<1e-25:
        return (k0/mu)*np.asarray(p)
    return (k0/mu)*np.expm1(beta*np.asarray(p))/beta


def derivative(p, beta, k0=P['k0_m2'], mu=P['mu_Pa_s']):
    return (k0/mu)*np.exp(beta*np.asarray(p))


def forcing(t, amplitude=P['pressure_amplitude_Pa']):
    return amplitude*np.sin(2*np.pi*np.asarray(t)/P['period_s'])


def flux_steady(t,beta,amplitude=P['pressure_amplitude_Pa']):
    return primitive(forcing(t,amplitude),beta)/P['L_m']


def mean_exact(beta, amplitude=P['pressure_amplitude_Pa']):
    if abs(beta)<1e-25:return 0.0
    return P['k0_m2']/(P['mu_Pa_s']*P['L_m']*beta)*(i0(beta*amplitude)-1.)


def nonlinear_timestep(previous, t,dt,beta, iterations=12):
    """Backward Euler finite volume, Newton linearized Kirchhoff-flux residual.

    C dp/dt + dq/dx=0, q=-dF(p)/dx; boundary pressures Dirichlet.
    Integrated continuity is respected to Newton residual tolerance.
    """
    x=np.asarray(previous,float).copy(); n=len(x)
    dx=P['L_m']/n; c=P['C_Pa_inv']; left=float(forcing(t)); right=0.
    # 1D finite-volume end-cell face distances are dx/2.
    for _ in range(iterations):
        f=primitive(x,beta); fp=derivative(x,beta)
        q=np.empty(n+1)
        q[0]=2*(float(primitive(left,beta))-f[0])/dx
        q[1:-1]=(f[:-1]-f[1:])/dx
        q[-1]=2*(f[-1]-float(primitive(right,beta)))/dx
        res=c*(x-previous)/dt + (q[1:]-q[:-1])/dx
        if max(abs(res))<1e-12:break
        coef=np.full(n,2.); coef[[0,-1]]=3.
        diag=c/dt +coef*fp/dx**2
        upper=-fp[1:]/dx**2
        lower=-fp[:-1]/dx**2
        ab=np.zeros((3,n));ab[1]=diag;ab[0,1:]=upper;ab[2,:-1]=lower
        delta=solve_banded((1,1),ab,-res)
        x=x+delta
        if max(abs(delta))<1e-8: break
    else:
        raise RuntimeError(f'Newton did not converge at t={t}: {max(abs(res))}')
    f=primitive(x,beta)
    q=np.empty(n+1); q[0]=2*(float(primitive(left,beta))-f[0])/dx
    q[1:-1]=(f[:-1]-f[1:])/dx
    q[-1]=2*(f[-1]-float(primitive(right,beta)))/dx
    r=c*(x-previous)/dt +(q[1:]-q[:-1])/dx
    if max(abs(r))>1e-10:
        raise RuntimeError(f'Continuity residual too large: {max(abs(r))}')
    return x,q


def simulate(n=32,dt=120.,cycles=10,beta=P['beta_Pa_inv']):
    period=P['period_s'];steps_per_period=round(period/dt);dt=period/steps_per_period
    steps=steps_per_period*cycles
    previous=np.zeros(n)
    # keep last two cycles, avoid storing enormous intermediates
    series=[];mass_max=0.;
    for z in range(1,steps+1):
        t=z*dt
        x,q=nonlinear_timestep(previous,t,dt,beta)
        # Model mass balance per unit cross-section, integral(Cp dx) [m].
        mass_change=P['C_Pa_inv']*P['L_m']*np.mean(x-previous)
        boundary_integral=(q[0]-q[-1])*dt
        err=mass_change-boundary_integral
        mass_max=max(mass_max,abs(err))
        if z>(cycles-2)*steps_per_period:
            series.append([t,t/period,float(forcing(t)),q[0],q[-1],float(np.mean(x))])
        previous=x
    arr=np.asarray(series)
    last=arr[-steps_per_period:]
    # compare spatial average q over last cycle; for periodic steady state inlet~outlet.
    qin=float(np.mean(last[:,3]));qout=float(np.mean(last[:,4]));
    return {'rows':arr,'mean_inlet_m_s':qin,'mean_outlet_m_s':qout,
            'boundary_mean_difference_m_s':qin-qout,
            'max_step_mass_residual_m':mass_max,
            'steps_per_cycle':steps_per_period,'dt_s':dt,
            'final_mean_p_Pa':float(np.mean(previous))}


def thermal_solution(q):
    L=P['L_m'];rho_cp=P['water_heat_capacity_vol_J_m3_K'];lam=P['thermal_conductivity_W_m_K'];
    delta=P['temperature_difference_K'];pe=rho_cp*q*L/lam
    xi=np.linspace(0.,1.,150)
    if abs(pe)<1e-9:
        profile=delta*xi;grad0=delta/L
    else:
        profile=delta*np.expm1(pe*xi)/np.expm1(pe)
        grad0=delta/L*pe/np.expm1(pe)
    # heat flux referenced to inlet T(0), *not* naive q*rho_cp*dT
    total_heat_left=-lam*grad0
    conduction_only=-lam*delta/L
    advective_scale=rho_cp*q*delta
    return {'Pe':pe,'thermal_x_m':xi,'thermal_deltaT_K':profile,
            'reference_conduction_W_m2':conduction_only,
            'total_heat_left_W_m2':total_heat_left,
            'incremental_total_heat_W_m2':total_heat_left-conduction_only,
            'rhoCp_q_dT_scale_W_m2':advective_scale}


def main():
    beta=P['beta_Pa_inv']; period=P['period_s']; q_exact=float(mean_exact(beta))
    table=[]
    for N in [8,16,32,64,128,256,512,1024]:
        t=np.arange(N)*period/N
        sampled=float(np.mean(flux_steady(t,beta)))
        table.append({'samples_per_cycle':N,'quadrature_qmean_m_s':sampled,
                      'analytic_qmean_m_s':q_exact,'abs_error_m_s':abs(sampled-q_exact)})
    for beta_test in [0.,1e-5,2e-5,-2e-5]:
        assert abs(np.mean(flux_steady(np.arange(4096)*period/4096,beta_test))-
                   mean_exact(beta_test)) < 5e-16
    assert abs(mean_exact(0.))<1e-22
    assert abs(mean_exact(beta)+mean_exact(-beta))<1e-21
    assert abs(mean_exact(beta,0.))<1e-21
    sim=simulate();
    thermal=thermal_solution(q_exact)
    # A reproducible baseline finite-volume simulation for beta = 0 (linear null).
    linear=simulate(cycles=8,beta=0.)
    # Sensitivity, beta(A) two-point alpha, expected ~A^2 in weak nonlinearity
    sensitivity=[]
    for b in [-3e-5,-2e-5,-1e-5,0,1e-5,2e-5,3e-5]:
        for a in [2500.,5000.,10000.,15000.]:
            val=mean_exact(b,a)
            sensitivity.append({'beta_per_Pa':b,'amplitude_Pa':a,'qmean_m_s':val,
                                'per_m2_liters_per_day':val*86400*1000})

    with (RESULTS/'quadrature_convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    with (RESULTS/'sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(sensitivity[0]));w.writeheader();w.writerows(sensitivity)
    with (RESULTS/'transient_last_two_cycles.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['t_s','period_count','left_pressure_Pa','qin_m_s','qout_m_s','mean_pressure_Pa']);w.writerows(sim['rows'])
    vals={'research_status':'SYNTHETIC / MODEL RESULTS ONLY; NO NEW FIELD DATA',
          'scenario':P,'analytic_mean_q_m_s':q_exact,
          'analytic_mean_q_L_m2_day':q_exact*86400*1000,
          'analytic_mean_q_m_year':q_exact*365.25*86400,
          'FV_n_cells':32,'FV_cycles':10,'FV_dt_s':sim['dt_s'],
          'FV_mean_inlet_m_s':sim['mean_inlet_m_s'],
          'FV_mean_outlet_m_s':sim['mean_outlet_m_s'],
          'FV_rel_to_analytic':sim['mean_outlet_m_s']/q_exact-1,
          'FV_cycle_mass_error_bound_m_s':sim['boundary_mean_difference_m_s'],
          'FV_max_step_mass_closure_abs_m':sim['max_step_mass_residual_m'],
          'linear_FV_mean_outlet_m_s':linear['mean_outlet_m_s'],
          'linear_FV_mean_inlet_m_s':linear['mean_inlet_m_s'],
          'linear_FV_max_mass_error_m':linear['max_step_mass_residual_m'],
          'heat_Pe':thermal['Pe'],
          'heat_conduction_only_W_m2':thermal['reference_conduction_W_m2'],
          'heat_with_flow_inlet_W_m2':thermal['total_heat_left_W_m2'],
          'heat_increment_W_m2':thermal['incremental_total_heat_W_m2'],
          'heat_advective_scale_W_m2':thermal['rhoCp_q_dT_scale_W_m2'],
          'runtime':{'python':platform.python_version(),'numpy':np.__version__,
                     'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},
          'validation_conclusion':'Synthetic verification does not validate hypothesis in nature.'}
    (RESULTS/'summary.json').write_text(json.dumps(vals,indent=2),encoding='utf8')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.dpi':130,'axes.spines.top':False,'axes.spines.right':False})
    # 1 cycle q vs p for steady analytic case
    phase=np.linspace(0,2*np.pi,900);T=phase/(2*np.pi)*period
    fig,ax=plt.subplots(figsize=(8.0,3.7));
    ax.plot(phase/(2*np.pi),flux_steady(T,0)*1e9,label='linear null: beta = 0',lw=2)
    ax.plot(phase/(2*np.pi),flux_steady(T,beta)*1e9,label='nonlinear test: beta = 2e-5 1/Pa',lw=2)
    ax.axhline(q_exact*1e9,color='grey',ls='--',lw=1.5,label='nonlinear cycle mean')
    ax.set_xlabel('Fraction of imposed cycle');ax.set_ylabel('Darcy flux [nm/s]')
    ax.set_title('Synthetic boundary-driven GRM benchmark: oscillation versus non-zero mean')
    ax.legend(loc='upper right',frameon=False,fontsize=8);ax.grid(alpha=.16);fig.tight_layout()
    fig.savefig(FIG/'01_flux_waveform.png',dpi=175);plt.close(fig)
    # transient q In q Out
    arr=sim['rows'];fig,ax=plt.subplots(figsize=(8.0,3.7))
    ax.plot(arr[:,1],arr[:,3]*1e9,label='inlet flux',lw=1.5)
    ax.plot(arr[:,1],arr[:,4]*1e9,label='outlet flux',lw=1.5,alpha=.9)
    ax.axhline(q_exact*1e9,color='grey',ls='--',label='analytic cycle mean')
    ax.set_xlabel('Elapsed cycles');ax.set_ylabel('Darcy flux [nm/s]')
    ax.set_title('Finite-volume solution: diffusive delay and the conserved cycle-mean flux')
    ax.legend(loc='upper right',frameon=False,fontsize=8);ax.grid(alpha=.16);fig.tight_layout()
    fig.savefig(FIG/'02_transient_fv.png',dpi=175);plt.close(fig)
    # beta sensitivity curves
    fig,ax=plt.subplots(figsize=(8.0,3.7))
    bgrid=np.linspace(-3e-5,3e-5,300)
    for a in [2500.,5000.,10000.,15000.]:
        ax.plot(bgrid*1e5,np.array([mean_exact(b,a) for b in bgrid])*1e9,label=f'A={a/1000:g} kPa',lw=1.8)
    ax.axhline(0,color='grey',lw=1)
    ax.set_xlabel('Permeability-pressure sensitivity beta [1e-5 /Pa]')
    ax.set_ylabel('Cycle-mean Darcy flux [nm/s]')
    ax.set_title('Synthetic sensitivity: sign, forcing amplitude and nonlinear strength')
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.16);fig.tight_layout()
    fig.savefig(FIG/'03_sensitivity.png',dpi=175);plt.close(fig)
    # thermal profile
    fig,ax=plt.subplots(figsize=(8.0,3.7))
    ax.plot(thermal['thermal_x_m'],thermal['thermal_deltaT_K'],lw=2,label='advection + conduction')
    ax.plot(thermal['thermal_x_m'],P['temperature_difference_K']*thermal['thermal_x_m']/P['L_m'],ls='--',lw=2,label='conduction only')
    ax.set_xlabel('Distance x [m]');ax.set_ylabel('Temperature above left boundary [K]')
    ax.set_title(f'Illustrative steady temperature profile, Pe={thermal["Pe"]:.3f}')
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.16);fig.tight_layout()
    fig.savefig(FIG/'04_thermal_profile.png',dpi=175);plt.close(fig)
    print(json.dumps(vals,indent=2))

if __name__=='__main__':main()
